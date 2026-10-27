suppressPackageStartupMessages({
  library(admiral)
  library(dplyr)
  library(tidyr)
  library(haven)
  library(lubridate)
  library(stringr)
})

read_sdtm <- function(domain) {
  read_xpt(file.path("data-raw", paste0(domain, ".xpt"))) |> convert_blanks_to_na()
}

save_adam <- function(df, name) saveRDS(df, file.path("adam", paste0(name, ".rds")))

# SAS rounds half away from zero, R rounds half to even
round_sas <- function(x, digits = 0) {
  m <- 10^digits
  sign(x) * floor(abs(x) * m + 0.5 + 1e-9) / m
}
