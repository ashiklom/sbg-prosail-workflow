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

lai <- rasterFromXYZ(cbind(simgrid - 0.5, t(X)))

# Random truncated normal
rtnorm <- function(n, mu, sigma) {
  x <- rnorm(n, mu, sigma)
  while (any(x < 0)) {
    x[x < 0] <- rnorm(sum(x < 0), mu, sigma)
  }
  x
}

# PROSPECT parameter distributions
drawN <- function(n = 1) 1 + rlnorm(n, 0.3, 0.2)
drawCab <- function(n = 1) rtnorm(n, 40, 10)
drawCar <- function(n = 1) rlnorm(n, 1, 1)
drawCbrown <- function(n = 1) rbeta(n, 0.9, 1)
drawCw <- function(n = 1) rlnorm(n, log(0.01), 1)
drawCm <- function(n = 1) rlnorm(n, log(0.01), 1)

# Draw random values for PROSPECT parameters
ni <- NROW(lai)
nj <- NCOL(lai)
nn <- ni * nj
N <- matrix(drawN(nn), ni, nj)
Cab <- matrix(drawCab(nn), ni, nj)
Car <- matrix(drawCar(nn), ni, nj)
Cbrown <- matrix(drawCbrown(nn), ni, nj)
Cw <- matrix(drawCw(nn), ni, nj)
Cm <- matrix(drawCm(nn), ni, nj)

lidfa <- matrix(runif(nn, -0.4, -0.3), ni, nj)
lidfb <- matrix(runif(nn, -0.2, -0.1), ni, nj)

# Now, try using this LAI field to run PRO4SAIL
wl <- 400:2500
# Last dimension:
# 1 - Bi-hemispherical
# 2 - Hemispherical-directional
# 3 - Directional-hemispherical
# 4 - Bi-directional
result <- array(numeric(), c(ni, nj, length(wl), 4))
pb <- progress::progress_bar$new(total = length(lai))
for (i in seq_len(ni)) {
  for (j in seq_len(nj)) {
    result[i, j, ,] <- PEcAnRTM::pro4sail(c(
      ## N = 1.4, Cab = 40, Car = 8, Cbrown = 0,
      ## Cw = 0.01, Cm = 0.01,
      N = N[i, j], Cab = Cab[i, j], Car = Car[i, j], Cbrown = Cbrown[i, j],
      Cw = Cw[i, j], Cm = Cm[i, j],
      LIDFa = lidfa[i, j], LIDFb = lidfb[i, j], TypeLIDF = 1,
      LAI = lai[i,j],
      q = 0.01, tts = 30, tto = 10,
      psi = 0, psoil = 0.7
    ))
    pb$tick()
  }
}

# Create the binary file

outdir <- file.path("data", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Swap the rows and columns of the PROSAIL array. Not sure why this is
# necessary, but maybe it's something to do with column- vs. row-major order of
# memory storage? Regardless, this makes the Python image appear correct.
result_t <- array(numeric(), c(nj, ni, length(wl), 4))
for (b in seq_len(4)) {
  for (a in seq_along(wl)) {
    result_t[,, a, b] <- t(result[,, a, b])
  }
}
# BSQ: row x column x band, which matches the array format
writeBin(c(result_t[,,,4]), file.path(outdir, "prosail-bdr.bsq"))
writeBin(c(result_t[,,,3]), file.path(outdir, "prosail-hdr.bsq"))

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

# ...and save the true values for validation, etc.
true_values <- list(
  lai = as.matrix(lai), N = N, Cab = Cab, Car = Car, Cbrown = Cbrown,
  Cw = Cw, Cm = Cm, lidfa = lidfa, lidfb = lidfb
)
saveRDS(true_values, file.path(outdir, "true-values.rds"))

result_nir <- result[,, 380, 4]
result_r <- result[,, 270, 4]
result_ndvi <- (result_nir - result_r) / (result_nir + result_r)
result_ndvi_raster <- raster(result_ndvi, xmn = 1, xmx = 50, ymn = 1, ymx = 50)

dir.create("figures", showWarnings = FALSE)
pdf(file.path("figures", "random-lai.pdf"), width = 14, height = 7)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))
plot(lai, main = "True LAI")
plot(result_ndvi_raster, main = "NDVI")
plot(getValues(result_ndvi_raster), getValues(lai),
     xlab = "NDVI (780 nm / 670 nm)", ylab = "LAI")
dev.off()
