gliht_band2num <- function(bandname) {
  out <- bandname %>%
    stringr::str_remove("^X") %>%
    stringr::str_remove("\\.Nanometers$") %>%
    as.numeric()
}
