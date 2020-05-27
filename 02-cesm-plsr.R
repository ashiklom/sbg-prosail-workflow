library(data.table)
library(fs)
library(ggplot2)
library(ncdf4)
library(magrittr)
library(sf)

plsrdir <- path("data", "plsr-models", "PLSR_Coefficients 2")
plsrfile <- path(plsrdir, "FFT_AVIRIS_Mean_Canopy_N_PLSR_Coefficients.csv")
stopifnot(file_exists(plsrfile))

dat <- fread(plsrfile)
dat[, wavelength := as.numeric(gsub("Wave_", "", V1))]

# PROSAIL is only 400-2500nm. This removes coefficients for wavelengths < 400
# nm. Those coefficients are all zero anyway, so it shouldn't affect the results.
dat_sub <- dat[!is.na(wavelength)][wavelength >= 400]
intercept <- dat[V1 == "(Intercept)", Perc_N]
waves <- dat_sub[["wavelength"]]
# PROSAIL is 400-2500 nm. This converts wavelengths into integer indices (i.e.
# 400 becomes 1, 401 becomes 2) for easier subsetting.
iwaves <- waves - 399
coefs <- dat_sub[["Perc_N"]]

# Plot the coefficients
ggplot(dat[!is.na(wavelength)]) +
  aes(x = wavelength, y = Perc_N) +
  geom_line()

# Read CLM spectra
clmfile <- path("data", "outputs", "clm-monthly-withprospect.nc")
stopifnot(file_exists(clmfile))

nc <- nc_open(clmfile)

# NOTE: Just using hemispherical-directional reflectance here for now. In
# reality, should probably use an average of hemispherical-directional and
# bi-directional reflectance weighted by direct/diffuse incident solar
# radiation, or something similar.
refl <- ncvar_get(nc, "hdr", start = c(1, 1, 7, 1), count = c(-1, -1, 1, -1))
lat <- ncvar_get(nc, "lat")
lon <- ncvar_get(nc, "lon")
refl_sub <- refl[,,iwaves]

# Set flat spectra (all reflectance = 0) to NA
naspec <- apply(refl_sub, c(1, 2), function(x) all(x == 0))
refl_sub[naspec] <- NA

do_nitrogen <- function(x) {
  # Model assumes that reflectance is in AVIRIS units, which are reflectance x
  # 10000 (so they can be distributed as integers).
  intercept + sum(10000 * x * coefs)
}

nitrogen <- apply(refl_sub, c(1, 2), do_nitrogen)
nitrogen[nitrogen < 0] <- 0

nitrogen_df <- expand.grid(lon = lon, lat = lat)
nitrogen_df$nitrogen <- c(nitrogen)
setDT(nitrogen_df)

# Land polygons from Natural Earth. This should be downloaded by the
# `01-cesm-prosail.R`.
landmask_file <- path("data", "landmask", "ne_110m_land.shp")
landmask_sf <- read_sf(landmask_file)

ggplot() +
  geom_sf(data = landmask_sf) +
  geom_tile(aes(x = lon, y = lat, fill = nitrogen),
            data = nitrogen_df) +
  scale_fill_viridis_c(na.value = NA)

ggsave("figures/clm-plsr-nitrogen-map.png", width = 10, height = 7,
       dpi = 300)
