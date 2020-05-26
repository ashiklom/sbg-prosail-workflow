library(data.table)
library(fs)
library(ggplot2)
library(ncdf4)

plsrdir <- path("data", "plsr-models", "PLSR_Coefficients 2")
plsrfile <- path(plsrdir, "FFT_AVIRIS_Mean_Canopy_N_PLSR_Coefficients.csv")
stopifnot(file_exists(plsrfile))

dat <- fread(plsrfile)
dat[, wavelength := as.numeric(gsub("Wave_", "", V1))]

dat_sub <- dat[!is.na(wavelength)][wavelength >= 400]
intercept <- dat[V1 == "(Intercept)", Perc_N]
waves <- dat_sub[["wavelength"]]
iwaves <- waves - 399
coefs <- dat_sub[["Perc_N"]]

# Look at the coefficients
ggplot(dat[!is.na(wavelength)]) +
  aes(x = wavelength, y = Perc_N) +
  geom_line()

# Read CLM spectra
clmfile <- path("data", "outputs", "clm-monthly.nc")
stopifnot(file_exists(clmfile))

nc <- nc_open(clmfile)

bhr <- ncvar_get(nc, "bhr", start = c(1, 1, 7, 1), count = c(-1, -1, 1, -1))
lat <- ncvar_get(nc, "lat")
lon <- ncvar_get(nc, "lon")
bhr_sub <- bhr[,,iwaves]

# Remove flat spectra (NA)
naspec <- apply(bhr_sub, c(1, 2), function(x) all(x == 0))
bhr_sub[naspec] <- NA

do_nitrogen <- function(x) {
  intercept + sum(10000 * x * coefs)
}

nitrogen <- apply(bhr_sub, c(1, 2), do_nitrogen)
nitrogen[nitrogen < 0] <- 0

nitrogen_df <- expand.grid(lon = lon, lat = lat)
nitrogen_df$nitrogen <- c(nitrogen)
setDT(nitrogen_df)

ggplot(nitrogen_df) +
  aes(x = lon, y = lat, fill = nitrogen) +
  geom_tile() +
  scale_fill_viridis_c()

ggsave("~/Pictures/clm-fake-n-map.png")
