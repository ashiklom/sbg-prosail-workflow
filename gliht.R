library(raster)
library(sf)
library(fs)
library(here)
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)

dir_walk(here("R"), source)

# Duke SPP1 coordinates
lat <- 35.97
lon <- -79.10

# Solar conditions
doy <- 297
# TODO: Assuming time is 10:30 AM, but is this right?
local_time_hr <- 10.5
cosz <- cos_solar_zenith_angle(doy, lat, lon, local_time_hr)
zenith_angle_deg <- acos(cosz) * 180 / pi

indir <- here("data", "g-liht", "Duke_AM_SPP1_24Oct2013")

hdrfiles <- dir_ls(indir, glob = "*.hdr")
reflfile <- dir_ls(indir, glob = "*refl*.hdr")

stopifnot(length(reflfile) == 1)
reflfile_nohdr <- str_remove(reflfile, "\\.hdr$")

refl <- brick(reflfile_nohdr)

bands <- gliht_band2num(names(refl))
bands_int <- round(bands)

rgb <- c(700, 510, 460)
irgb <- which_closest(rgb, bands)

samplefile <- path(indir, "mysamples.csv")
if (!file_exists(samplefile)) {
  mask_zero <- refl[[irgb[1]]] > 0
  refl_rgb <- mask(refl[[irgb]], mask_zero, maskvalue = FALSE)

  plotRGB(refl_rgb, stretch = "hist")

  tidy_samples <- function(samp_df, ...) {
    .dots <- rlang::list2(...)
    samp_df %>%
      as_tibble() %>%
      pivot_longer(-cell, names_to = "wavelength", values_to = "reflectance") %>%
      mutate(
        wavelength = gliht_band2num(wavelength),
        reflectance = reflectance / 10000,
        !!!.dots
      )
  }

  # Select some samples interactively to inspect
  darksamples <- tidy_samples(click(refl, cell = TRUE, show = FALSE),
                              category = "dark")
  lightsamples <- tidy_samples(click(refl, cell = TRUE, show = FALSE),
                               category = "light")
  redsamples <- tidy_samples(click(refl, cell = TRUE, show = FALSE),
                             category = "red")
  all_tidy <- bind_rows(darksamples, lightsamples, redsamples)
  write_csv(all_tidy, samplefile)
} else {
  all_tidy <- read_csv(samplefile)
}

cells <- unique(all_tidy$cell)

prosail2gliht <- function(params) {
  # Params:
  # 1-6 -- N, Cab, Car, Cbrown, Cw, Cm,
  # 7-8 -- LIDFa, LIDFb,
  # 9-11 -- LAI, q, psoil
  sail_args <- sail_defaults
  sail_args["tts"] <- zenith_angle_deg
  # Observer zenith is 0. To a first approximation...
  sail_args["tto"] <- 0
  # Sun-sensor azimuth angle is 0. Probably not, but OK for now...
  sail_args["psi"] <- 0
  # Leaf traits, LIDF
  sail_args["N"] <- params[1]
  sail_args["Cab"] <- params[2]
  sail_args["Car"] <- params[3]
  sail_args["Cw"] <- params[4]
  sail_args["Cm"] <- params[5]
  sail_args["LAI"] <- params[6]
  sail_args["psoil"] <- params[7]

  sail_all <- PEcAnRTM::pro4sail(sail_args)
  # Clear sky conditions: Assume diffuse contribution is small, so use BDR (4)
  # TODO: More sophisticated treatment?
  # TODO: Grab closest bands. In reality, should use FWHM...
  sail <- as.numeric(sail_all[, 4][[bands_int]])
  sail
}

fit_prosail_cell <- function(cell, dat = all_tidy) {
  datsub <- dplyr::filter(dat, cell == !!cell)
  gliht <- datsub %>%
    dplyr::arrange(wavelength) %>%
    dplyr::pull(reflectance)

  optfun <- function(params) {
    sail <- tryCatch(
      prosail2gliht(params),
      error = function(e) {
        message("Hit the following error: ",
                conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(sail)) return(1e9)
    if (any(!is.finite(sail))) return(1e9)
    sum((prosail2gliht(params) - gliht) ^ 2)
  }

  start <- sail_defaults[c("N", "Cab", "Car", "Cw", "Cm",
                           "LAI", "psoil")]

  ## Test that parameters work
  xx <- start
  test1 <- optfun(xx)
  stopifnot(is.numeric(test1), test1 < 1e7)
  xx["LAI"] <- 5
  test2 <- optfun(xx)
  stopifnot(is.numeric(test2), test2 < 1e7)

  # With base R optimization
  fit <- optim(
    start, optfun,
    method = "L-BFGS-B",
    lower = c(1, 0, 0, 0, 0, 0, 0),
    upper = c(10, 150, 50, 0.1, 0.1, 10, 1)
  )

  # With differential evolution algorithm
  ## fit <- DEoptim::DEoptim(
  ##   optfun,
  ##   lower = c(1, 0, 0, 0, 0, 0, 0),
  ##   upper = c(10, 150, 50, 0.1, 0.1, 10, 1),
  ##   control = list(itermax = 1000)
  ## )

  fit
}

# Fit PRO4SAIL to each cell
fits <- list()
for (cell in cells) {
  message("Cell ", cell)
  result <- tryCatch(
    fit_prosail_cell(cell),
    error = function(e) {
      message("Hit error: ", conditionMessage(e))
      return(list(par = NULL))
    }
  )
  fits[[as.character(cell)]] <- result$par
}

fit_params <- fits %>%
  map_dfr(enframe, name = "parameter", .id = "cell") %>%
  mutate(cell = as.numeric(cell))

ggplot(fit_params) +
  aes(x = "", y = value) +
  geom_violin() +
  geom_jitter() +
  facet_wrap(vars(parameter), scales = "free_y")

sim_fits <- map(fits, prosail2gliht) %>%
  map(as_tibble) %>%
  map(function(x) mutate(x, wavelength = bands)) %>%
  bind_rows(.id = "cell") %>%
  mutate(cell = as.numeric(cell))

fit_obs <- sim_fits %>%
  left_join(all_tidy, c("cell", "wavelength"))

ggplot(fit_obs) +
  aes(x = wavelength) +
  geom_line(aes(y = value, color = "PRO4SAIL BDR")) +
  geom_line(aes(y = reflectance, color = "G-LiHT")) +
  facet_wrap(vars(cell))

ggsave("figures/gliht-sail-duke-AM.png",
       width = 13.5, height = 8, dpi = 300)

fit_obs %>%
  mutate(model_bias = value - reflectance) %>%
  ggplot() +
  aes(x = wavelength, y = model_bias, group = cell) +
  geom_line(color = "gray30") +
  geom_hline(yintercept = 0, color = "red")

ggsave("figures/gliht-sail-duke-AM-bias.png",
       width = 7, height = 4, dpi = 300)
