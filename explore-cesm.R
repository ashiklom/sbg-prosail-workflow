library(fs)
library(here)
library(ncdf4)
library(raster)
library(magrittr, include.only = "%>%")
library(sf)

# Load helper functions
dir_walk(here("R"), source)

stopifnot(
  requireNamespace("rnaturalearth", quietly = TRUE),
  requireNamespace("rnaturalearthdata", quietly = TRUE),
  requireNamespace("rgdal", quietly = TRUE),
  requireNamespace("PEcAnRTM", quietly = TRUE)
)

test_outdir <- dir_create(here("data", "cesm", "monthly_test_run"))
fname <- "CLM5_1584188701_f45_f45_mg37.clm2.h0.1857-02-01-00000.nc"
fpath <- path(test_outdir, fname)

# Download CLM
if (!file_exists(fpath)) {
  fosf_id <- "27a8v"
  download.file(path("https://osf.io", "download", fosf_id), fpath)
}

clm_input_path <- path(test_outdir, "surfdata_4x5_hist_78pfts_CMIP6_simyr2000_c190214.nc")
# Download CLM
if (!file_exists(clm_input_path)) {
  fosf_id <- "ftdzh"
  download.file(path("https://osf.io", "download", fosf_id), clm_input_path)
}

innc <- ncdf4::nc_open(clm_input_path)

wgs84 <- crs("+init=epsg:4326")

landmask_file <- path("data", "landmask", "ne_110m_land.shp")
if (!file_exists(landmask_file)) {
  rnaturalearth::ne_download(
    "small", "land", "physical",
    returnclass = "sf", destdir = path_dir(landmask_file),
    load = FALSE
  )
}

landmask <- read_sf(landmask_file) %>%
  st_combine() %>%
  as_Spatial()

elai_b <- brick(fpath, varname = "ELAI", crs = wgs84) %>%
  # Convert 0-360 to -180-180
  rotate()

elai_v <- extract(elai_b, landmask, cellnumbers = TRUE)[[1]]
notna <- apply(is.na(elai_v), 1, sum) == 0
xy <- xyFromCell(elai_b, elai_v[notna, "cell"])
elai <- elai_v[notna, -1]
stopifnot(nrow(xy) == nrow(elai))

ncell <- nrow(elai)
nmonth <- ncol(elai)
nn <- ncell * nmonth
wl <- 400:2500

nwl <- length(wl)

# Set PROSPECT parameters
# TODO: Cab, Car should vary seasonally
N <- matrix(1 + rlnorm(nn, 0.3, 0.2), ncell, nmonth)
Cab <- matrix(rtnorm(nn, 40, 10), ncell, nmonth)
Car <- matrix(rlnorm(nn, 1, 1), ncell, nmonth)
Cbrown <- matrix(rbeta(nn, 0.9, 1), ncell, nmonth)
Cw <- matrix(rlnorm(nn, log(0.01), 1), ncell, nmonth)
Cm <- matrix(rlnorm(nn, log(0.01), 1), ncell, nmonth)
lidfa <- matrix(runif(nn, -0.4, -0.3), ncell, nmonth)
lidfb <- matrix(runif(nn, -0.2, -0.1), ncell, nmonth)

# Now, run PROSPECT at each value
result <- array(numeric(), c(ncell, nmonth, nwl, 4))
pb <- progress::progress_bar$new(total = nn)

for (i in seq_len(ncell)) {
  for (j in seq_len(nmonth)) {
    result[i, j, ,] <- PEcAnRTM::pro4sail(c(
      N = N[i, j], Cab = Cab[i, j], Car = Car[i, j], Cbrown = Cbrown[i, j],
      Cw = Cw[i, j], Cm = Cm[i, j],
      LIDFa = lidfa[i, j], LIDFb = lidfb[i, j], TypeLIDF = 1,
      LAI = elai[i, j],
      q = 0.01, tts = 30, tto = 10,
      psi = 0, psoil = 0.7
    ))
    pb$tick()
  }
}

# Create output NetCDF file
times <- colnames(elai) %>%
  gsub("^X", "", .) %>%
  strptime("%Y.%m.%d", tz = "UTC") %>%
  as.POSIXct() %>%
  as.numeric()

londim <- ncdim_def("lon", "degrees_east", unique(xy[, "x"]))
latdim <- ncdim_def("lat", "degrees_north", unique(xy[, "y"]))
timedim <- ncdim_def("time", posixct_cf_unit, times)
wavedim <- ncdim_def("wavelength", "nanometers", wl)

fillval <- -999

bhr <- ncvar_def("bhr", "0-1", list(londim, latdim, timedim, wavedim), fillval,
                 "bi-hemispherical reflectance")
hdr <- ncvar_def("hdr", "0-1", list(londim, latdim, timedim, wavedim), fillval,
                 "hemispherical-directional reflectance")
dhr <- ncvar_def("dhr", "0-1", list(londim, latdim, timedim, wavedim), fillval,
                 "directional-hemispherical reflectance")
bdr <- ncvar_def("bdr", "0-1", list(londim, latdim, timedim, wavedim), fillval,
                 "bi-directional reflectance")

outdir <- dir_create(here("data", "outputs"))
outnc <- nc_create(
  path(outdir, "clm-monthly.nc"),
  list(bhr, hdr, dhr, bdr),
  force_v4 = TRUE
)

# Write results to file
# TODO: This step is really slow. Is there a way to optimize it?
message("Writing output.")
pb <- progress::progress_bar$new(total = ncell)
for (icell in seq_len(ncell)) {
  ilon <- which(londim$vals == xy[icell, "x"])
  ilat <- which(latdim$vals == xy[icell, "y"])
  ncvar_put(outnc, bhr, result[icell,,,1],
            start = c(ilon, ilat, 1, 1), count = c(1, 1, -1, -1))
  ncvar_put(outnc, hdr, result[icell,,,2],
            start = c(ilon, ilat, 1, 1), count = c(1, 1, -1, -1))
  ncvar_put(outnc, dhr, result[icell,,,3],
            start = c(ilon, ilat, 1, 1), count = c(1, 1, -1, -1))
  ncvar_put(outnc, bdr, result[icell,,,4],
            start = c(ilon, ilat, 1, 1), count = c(1, 1, -1, -1))
  pb$tick()
}

nc_close(outnc)

saveRDS(
  list(xy = xy,
       N = N, Cab = Cab, Car = Car, Cbrown, Cw = Cw, Cm = Cm,
       lidfa = lidfa, lidfb = lidfb, lai = elai),
  path(outdir, "clm-monthly-true-values.rds")
)

# Test the output file
if (interactive()) {
  nc2 <- nc_open(path(outdir, "clm-monthly.nc"))
  zz <- ncvar_get(nc2, "bhr")
  nc_close()

  col <- colorRampPalette(c("blue3", "green3", "red3"))(12)
  matplot(400:2500, t(zz[1,8,,]), type = 'l', col = col, lty = 1,
          xlab = "Wavelength (nm)", ylab = "BHR")
  legend("topright", as.character(1:12), lty = 1, col = col,
         ncol = 2, title = "Month")

}
