# Random truncated normal
rtnorm <- function(n, mu, sigma) {
  x <- rnorm(n, mu, sigma)
  while (any(x < 0)) {
    x[x < 0] <- rnorm(sum(x < 0), mu, sigma)
  }
  x
}
