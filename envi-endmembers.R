library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

fname2 <- "./data/hypertrace-endmembers/multi_endmember_roi_w_bad_bands.txt"
file_raw <- readLines(fname2)
header_raw <- grep("^;", file_raw, value = TRUE)
roi_names_raw <- grep("ROI name:", header_raw, value = TRUE)
roi_type <- str_match(roi_names_raw, ": (.*)")[,2]
col_names <- str_subset(header_raw, "^; +ID") %>%
  str_remove("^; +") %>%
  str_replace_all(" {2,}", ",") %>%
  str_split(",") %>%
  .[[1]]

tab2 <- read.table(fname2, comment.char = ";")
colnames(tab2) <- col_names

dat_wide <- as_tibble(tab2) %>%
  mutate(class = roi_type, i = row_number()) %>%
  select(i, class, starts_with("B"))

dat <- pivot_longer(
  dat_wide,
  starts_with("B"),
  names_to = "band",
  values_to = "value"
)

envi_file <- "./data/hypertrace-endmembers/ang20180815t203458_corr_v2r2_img.hdr"
envi_raw <- readLines(envi_file)

replace_header_tag <- function(lines, tag, value) {
  out <- lines
  newline <- glue::glue("{tag} = {value}")
  pattern <- glue::glue("^ *{tag} *=")
  out[str_detect(out, pattern)] <- newline
  out
}

waves <- str_subset(envi_raw, "^wavelength +=") %>%
  str_match("\\{(.*)\\}") %>%
  .[,2] %>%
  str_split(",") %>%
  .[[1]] %>%
  as.numeric()

bands_waves <- tibble(band = paste0("B", seq_along(waves)), wave = waves)

dat2 <- left_join(dat, bands_waves, "band")

plt <- dat2 %>%
  mutate(value = if_else(wave > 1350 & wave < 1450, NA_real_, value),
         value = if_else(wave > 1800 & wave < 1920, NA_real_, value)) %>%
  ggplot() +
  aes(x = wave, y = value, color = class, group = i) +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  theme_bw()

ggsave(
  "figures/envi-endmembers.png", plt,
  width = 7, height = 7, dpi = 300
)

envi_endmember <- function(m) {
  endmember <- unique(m$class)
  stopifnot(length(endmember) == 1)

  # Convert data to a matrix
  m2 <- m %>%
    select(starts_with("B")) %>%
    as.matrix()

  marr <- array(m2, c(1, nrow(m2), ncol(m2)))
  dimnames(marr) <- list(NULL, NULL, paste0("Wave_", waves))

  mr <- raster::brick(marr)
  raster::metadata(mr) <- list(Wavelength = waves)

  outdir <- file.path("data", "outputs", "envi-endmembers")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(outdir, paste0(endmember, ".envi"))
  raster::writeRaster(mr, outfile, overwrite = TRUE)

  # Add other metadata to header
  outfile_hdr <- file.path(outdir, paste0(endmember, ".hdr"))
  md_fields <- paste("wavelength", "correction factors", "fwhm", "bbl",
                     "smoothing factors", "data ignore value",
                     sep = "|")
  other_metadata <- str_subset(
    envi_raw,
    glue::glue("^({md_fields}) *=")
  )
  lapply(other_metadata, write, file = outfile_hdr, append = TRUE)
  outfile
}

dat_wide %>%
  group_split(class) %>%
  lapply(envi_endmember)
