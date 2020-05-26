library(data.table)
library(fs)
library(ggplot2)
library(PEcAnRTM)

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

do_nitrogen <- function(x) intercept + sum(10000 * x * coefs)

prosail_range <- function(mmin, mmax, var) {
  out <- array(numeric(), c(2101, 4, nn))
}

# LAI variation
nn <- 50
lai_seq <- seq(0, 5, length.out = nn)
prosail_lai <- array(numeric(), c(2101, 4, nn))
for (i in seq_len(nn)) {
  lai <- lai_seq[i]
  prosail_lai[,,i] <- pro4sail(c(1.4, 40, 8, 0, 0.01, 0.01,
                                 -0.35, -0.15, 1,
                                 lai, 0.01, 30, 10, 0, 0.7))
}

chl_seq <- seq(10, 100, length.out = nn)
prosail_chl <- array(numeric(), c(2101, 4, nn))
for (i in seq_len(nn)) {
  chl <- chl_seq[i]
  prosail_chl[,,i] <- pro4sail(c(1.4, chl, 8, 0, 0.01, 0.01,
                                 -0.35, -0.15, 1,
                                 3, 0.01, 30, 10, 0, 0.7))
}

plai_sub <- prosail_lai[iwaves,2,]
lai_nitrogen <- apply(plai_sub, 2, do_nitrogen)
pchl_sub <- prosail_chl[iwaves,2,]
chl_nitrogen <- apply(ps_sub, 2, do_nitrogen)

png("~/Pictures/nitrogen-sensitivity.png")
par(mfrow = c(1, 2))
plot(lai_seq, lai_nitrogen, xlab = "LAI")
plot(chl_seq, chl_nitrogen, xlab = "Chl")
dev.off()

png("~/Pictures/spec-variability.png", width = 800)
par(mfrow = c(1, 2))
matplot(prosail_lai[,2,], type = "l")
matplot(prosail_chl[,2,], type = "l")
dev.off()

matplot(prosail_out[,,15], type = "l")
legend("topright", c("bhr", "hdr", "dhr", "bdr"),
       col = 1:4, lty = 1:4)
