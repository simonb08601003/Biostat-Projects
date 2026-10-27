# ADQSADAS: ADAS-Cog(11) total by visit, windowed, with LOCF
source("R/00_setup.R")

adsl <- readRDS("adam/adsl.rds")
qs <- read_sdtm("qs") |> filter(QSCAT == "ALZHEIMER'S DISEASE ASSESSMENT SCALE")

# The study codes 14 items; ADAS-Cog(11) drops items 03, 09 and 10.
# 1-3 missing items: prorate to 70 using the max score of the items observed.
item_max <- tibble(
  QSTESTCD = sprintf("ACITM%02d", c(1, 2, 4:8, 11:14)),
  max = c(10, 5, 5, 5, 5, 8, 12, 5, 5, 5, 5)
)
totals <- qs |>
  inner_join(item_max, by = "QSTESTCD") |>
  summarise(
    n_obs = sum(!is.na(QSSTRESN)),
    total = sum(QSSTRESN, na.rm = TRUE),
    max_obs = sum(max[!is.na(QSSTRESN)]),
    .by = c(USUBJID, VISITNUM, VISIT, QSDTC)
  ) |>
  mutate(AVAL = case_when(n_obs == 11 ~ total,
                          n_obs >= 8 ~ total * 70 / max_obs))

windows <- tribble(
  ~AVISIT,    ~AVISITN, ~AWLO, ~AWHI, ~AWTARGET,
  "Baseline", 0,        -Inf,  1,     1,
  "Week 8",   8,        2,     84,    56,
  "Week 16",  16,       85,    140,   112,
  "Week 24",  24,       141,   Inf,   168
)

adqs <- totals |>
  left_join(qs |> filter(QSTESTCD == "ACTOT") |> select(USUBJID, VISITNUM, QSSEQ),
            by = c("USUBJID", "VISITNUM")) |>
  mutate(PARAMCD = "ACTOT", PARAM = "Adas-Cog(11) Subscore", PARAMN = 1) |>
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs(USUBJID),
    new_vars = exprs(TRTSDT, TRTP = TRT01P, TRTPN = TRT01PN, SITEGR1, AGE,
                     AGEGR1, SEX, RACE, ITTFL, EFFFL, COMP24FL)
  ) |>
  derive_vars_dt(new_vars_prefix = "A", dtc = QSDTC) |>
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) |>
  filter(!is.na(AVAL)) |>
  cross_join(windows) |>
  filter(ADY >= AWLO, ADY <= AWHI) |>
  mutate(AWTDIFF = abs(ADY - AWTARGET)) |>
  # closest to target day; ties go to the later visit
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, AVISITN),
    order = exprs(AWTDIFF, desc(ADY)),
    new_var = ANL01FL,
    mode = "first"
  ) |>
  select(-TRTSDT)

anl <- adqs |>
  filter(ANL01FL == "Y") |>
  mutate(ABLFL = if_else(AVISITN == 0, "Y", NA_character_)) |>
  derive_var_base(by_vars = exprs(USUBJID, PARAMCD), source_var = AVAL, new_var = BASE)

# LOCF into every missed post-baseline window, for all ITT subjects
expected <- adsl |>
  filter(ITTFL == "Y") |>
  select(USUBJID) |>
  cross_join(windows |> filter(AVISITN > 0) |> select(AVISIT, AVISITN)) |>
  mutate(PARAMCD = "ACTOT")

anl <- anl |>
  mutate(ADT = as.character(ADT)) |>
  derive_locf_records(
    dataset_ref = expected,
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(AVISITN),
    keep_vars = exprs(PARAM, PARAMN, BASE, TRTP, TRTPN, SITEGR1, AGE, AGEGR1,
                      SEX, RACE, ITTFL, EFFFL, COMP24FL)
  ) |>
  mutate(
    ADT = as.Date(ADT),
    ANL01FL = if_else(is.na(DTYPE), ANL01FL, "Y"),
    CHG = if_else(AVISITN > 0, AVAL - BASE, NA_real_),
    PCHG = if_else(AVISITN > 0 & BASE != 0, 100 * CHG / BASE, NA_real_)
  )

adqsadas <- bind_rows(anl, filter(adqs, is.na(ANL01FL))) |>
  arrange(USUBJID, AVISITN, ADY)

save_adam(adqsadas, "adqsadas")
message("ADQSADAS: ", nrow(adqsadas), " records")
