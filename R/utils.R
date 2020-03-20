library(magrittr, include.only = "%>%")

posixct_cf_unit <- "seconds since 1970-01-01 00:00:00 -00:00"

which_closest <- function(x, from) {
  vapply(x, function(i) which.min(abs(i - from)), numeric(1))
}

sail_defaults <- c(
  N = 1.5, Cab = 40, Car = 8, Cbrown = 0, Cw = 0.01, Cm = 0.009,
  LIDFa = -0.35, LIDFb = -0.15, TypeLIDF = 1,
  LAI = 3, q = 0.01,
  tts = 30, tto = 10, psi = 0,
  psoil = 1
)
