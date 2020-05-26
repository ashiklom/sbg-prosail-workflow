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
if (!file_exists(clm_input_path)) {
  fosf_id <- "ftdzh"
  download.file(path("https://osf.io", "download", fosf_id), clm_input_path)
}

clm_param_path <- path(test_outdir, "clm5_params.c171117.nc")
if (!file_exists(clm_param_path)) {
  fosf_id <- "t84u9"
  download.file(path("https://osf.io", "download", fosf_id), clm_param_path)
}

ncparam <- nc_open(clm_param_path)
nat_pfts <- trimws(ncvar_get(ncparam, "pftname")[1:15])

innc <- nc_open(clm_input_path)
nat_pft <- ncvar_get(innc, "PCT_NAT_PFT")

pft_prospect <- read.csv("data/clm-pft-prosail-param.csv")
pft_prospect_mat <- rbind(rep(0, 6), as.matrix(pft_prospect[,-1]))

prospect_grid <- array(numeric(), c(dim(nat_pft)[1:2], ncol(pft_prospect_mat)))
for (i in seq_len(nrow(nat_pft))) {
  for (j in seq_len(ncol(nat_pft))) {
    prospect_grid[i, j, ] <- (nat_pft[i, j, ] / 100) %*% pft_prospect_mat
  }
}

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

# Now, run PROSPECT at each value
result <- array(numeric(), c(ncell, nmonth, nwl, 4))
pb <- progress::progress_bar$new(total = nn)
for (i in seq_len(ncell)) {
  for (j in seq_len(nmonth)) {
    ilat <- which(unique(xy[, "y"]) == xy[i, "y"])
    ilon <- which(unique(xy[, "x"]) == xy[i, "x"])
    prospect_params <- prospect_grid[ilon, ilat, ]
    if (all(prospect_params == 0)) {
      result[i, j, , ] <- NA
    } else {
      result[i, j, ,] <- PEcAnRTM::pro4sail(c(
        # PROSPECT params: N Cab Car Cbrown Cw Cm
        prospect_params,
        # SAIL params
        LIDFa = -0.35, LIDFb = -0.15, TypeLIDF = 1,
        LAI = elai[i, j],
        q = 0.01, tts = 30, tto = 10,
        psi = 0, psoil = 0.7
      ))
    }
    pb$tick()
  }
}

# Test the output
## iii <- which(!is.na(result[, 7, 5, 2]))
## matplot(t(result[sample(iii, 100),7,,2]), type = "l")

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

prospect_N <- ncvar_def("N", "1-Inf", list(londim, latdim), fillval,
                        "effective number of leaf layers")
prospect_Cab <- ncvar_def("Cab", "ug cm-2", list(londim, latdim), fillval,
                          "chlorophyll content")
prospect_Car <- ncvar_def("Car", "ug cm-2", list(londim, latdim), fillval,
                          "carotenoid content")
prospect_Cbrown <- ncvar_def("Cbrown", "0-1", list(londim, latdim), fillval,
                             "brown matter content")
prospect_Cw <- ncvar_def("Cw", "g cm-2", list(londim, latdim), fillval,
                         "leaf water content")
prospect_Cm <- ncvar_def("Cm", "ug cm-2", list(londim, latdim), fillval,
                         "dry matter content")
elai_nc <- ncvar_def("ELAI", "unitless", list(londim, latdim, timedim), fillval,
                     "Effective leaf area index")
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
  path(outdir, "clm-monthly-withprospect.nc"),
  list(prospect_N, prospect_Cab, prospect_Car, prospect_Cbrown,
       prospect_Cw, prospect_Cm, elai_nc,
       bhr, hdr, dhr, bdr),
  force_v4 = TRUE
)

# Write results to file
# TODO: This step is really slow. Is there a way to optimize it?
message("Writing output.")
pb <- progress::progress_bar$new(total = ncell)
for (icell in seq_len(ncell)) {
  ilon <- which(londim$vals == xy[icell, "x"])
  ilat <- which(latdim$vals == xy[icell, "y"])
  # PROSPECT params
  ncvar_put(outnc, prospect_N, prospect_grid[ilon, ilat, 1],
            start = c(ilon, ilat), count = c(1, 1))
  ncvar_put(outnc, prospect_Cab, prospect_grid[ilon, ilat, 2],
            start = c(ilon, ilat), count = c(1, 1))
  ncvar_put(outnc, prospect_Car, prospect_grid[ilon, ilat, 3],
            start = c(ilon, ilat), count = c(1, 1))
  ncvar_put(outnc, prospect_Cbrown, prospect_grid[ilon, ilat, 4],
            start = c(ilon, ilat), count = c(1, 1))
  ncvar_put(outnc, prospect_Cw, prospect_grid[ilon, ilat, 5],
            start = c(ilon, ilat), count = c(1, 1))
  ncvar_put(outnc, prospect_Cm, prospect_grid[ilon, ilat, 6],
            start = c(ilon, ilat), count = c(1, 1))
  # LAI from CLM
  ncvar_put(outnc, elai_nc, elai[icell,],
            start = c(ilon, ilat, 1), count = c(1, 1, -1))
  # SAIL outputs
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

# Test the output file
if (interactive()) {
  nc2 <- nc_open(path(outdir, "clm-monthly-withprospect.nc"))
  zz <- ncvar_get(nc2, "hdr")
  nc_close()

  col <- colorRampPalette(c("blue3", "green3", "red3"))(12)
  matplot(400:2500, t(zz[1,8,,]), type = 'l', col = col, lty = 1,
          xlab = "Wavelength (nm)", ylab = "BHR")
  legend("topright", as.character(1:12), lty = 1, col = col,
         ncol = 2, title = "Month")

}
