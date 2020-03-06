# Example PRO4SAIL simulation workflow with simulated LAI
# Author: Alexey Shiklomanov
library(raster)

# Set RNG seed for reproducibility...
set.seed(8675309)

# Multivariate normal distribution
rmvn <- function(n, mu = 0, Sigma = diag(length(mu))) {
  p <- length(mu)
  stopifnot(all(dim(Sigma) == p))
  if (!all(dim(Sigma) == p)) stop("Dimension mismatch")
  D <- chol(Sigma)
  matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p))
}

simgrid <- expand.grid(1:50, 1:50)
n <- nrow(simgrid)

distance <- as.matrix(dist(simgrid))
phi <- 0.05

# Use lognormal distribution here
# Mean LAI is exp(0.5) ~~ 1.65
X <- exp(rmvn(1, rep(0.5, n), 0.2 * exp(-phi * distance)))
## hist(X, xlab = "LAI")

Xraster <- rasterFromXYZ(cbind(simgrid - 0.5, t(X)))
## plot(Xraster)

lai <- Xraster

# Now, try using this LAI field to run PRO4SAIL
wl <- 400:2500
# Last dimension:
# 1 - Bi-hemispherical
# 2 - Hemispherical-directional
# 3 - Directional-hemispherical
# 4 - Bi-directional
result <- array(numeric(), c(NROW(lai), NCOL(lai), length(wl), 4))
pb <- progress::progress_bar$new(total = length(lai))
for (i in 1:NROW(lai)) {
  for (j in 1:NCOL(lai)) {
    result[i, j, ,] <- PEcAnRTM::pro4sail(c(
      N = 1.4, Cab = 40, Car = 8, Cbrown = 0,
      Cw = 0.01, Cm = 0.01,
      LIDFa = -0.35, LIDFb = -0.15, TypeLIDF = 1,
      LAI = lai[i,j],
      q = 0.01, tts = 30, tto = 10,
      psi = 0, psoil = 0.7
    ))
    pb$tick()
  }
}

# Create the binary file
# BSQ: row x column x band, which matches the array format
outdir <- file.path("data", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
writeBin(c(result[,,,4]), file.path(outdir, "prosail-bdr.bsq"))
writeBin(c(result[,,,3]), file.path(outdir, "prosail-hdr.bsq"))

# Now, write the header files
hdr_bdr <- c(
  "ENVI",
  "description = {PROSAIL-simulated reflectance}",
  glue::glue("lines = {NCOL(result)}"),
  glue::glue("samples = {NROW(result)}"),
  glue::glue("bands = {length(wl)}"),
  "header offset = 0",
  "file type = ENVI Standard",
  "data type = 5",
  "interleave = bsq",
  "sensor type = Unknown",
  "byte order = 0",
  "wavelength units = Nanometers",
  glue::glue("wavelength = { <<paste(wl, collapse = ', ')>> }",
             .open = "<<", .close = ">>")
)
writeLines(hdr_bdr, file.path(outdir, "prosail-bdr.bsq.hdr"))
writeLines(hdr_bdr, file.path(outdir, "prosail-hdr.bsq.hdr"))

## result_nir <- result[,, 380, 4]
## result_r <- result[,, 270, 4]
## result_ndvi <- (result_nir - result_r) / (result_nir + result_r)
## result_ndvi_raster <- raster(result_ndvi, xmn = 1, xmx = 50, ymn = 1, ymx = 50)
## par(mfrow = c(1, 3))
## plot(result_ndvi_raster, main = "NDVI")
## plot(Xraster, main = "True LAI")
## plot(getValues(result_ndvi_raster), getValues(Xraster),
##      xlab = "NDVI", ylab = "LAI")
