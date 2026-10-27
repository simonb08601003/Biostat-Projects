# ADSL
source("R/00_setup.R")

dm <- read_sdtm("dm")
ex <- read_sdtm("ex")
ds <- read_sdtm("ds")
sv <- read_sdtm("sv")
qs <- read_sdtm("qs")
vs <- read_sdtm("vs")
sc <- read_sdtm("sc")

adsl <- dm |>
  filter(ARMCD != "Scrnfail") |>
  mutate(
    TRT01P = ARM,
    TRT01A = ACTARM,
    TRT01PN = case_match(ARM, "Placebo" ~ 0, "Xanomeline Low Dose" ~ 54,
                         "Xanomeline High Dose" ~ 81),
    TRT01AN = TRT01PN
  ) |>
  derive_vars_merged(
    dataset_add = derive_vars_dt(ex, new_vars_prefix = "EXST", dtc = EXSTDTC),
    by_vars = exprs(USUBJID),
    new_vars = exprs(TRTSDT = EXSTDT),
    order = exprs(EXSTDT, EXSEQ),
    mode = "first"
  ) |>
  derive_vars_merged(
    dataset_add = derive_vars_dt(ex, new_vars_prefix = "EXEN", dtc = EXENDTC),
    by_vars = exprs(USUBJID),
    new_vars = exprs(TRTEDT = EXENDT),
    order = exprs(EXENDT, EXSEQ),
    mode = "last"
  ) |>
  derive_vars_dt(new_vars_prefix = "RFEN", dtc = RFENDTC) |>
  # last patch has no end date -> use end of study
  mutate(TRTEDT = if_else(is.na(TRTEDT) & !is.na(TRTSDT), RFENDT, TRTEDT)) |>
  derive_var_trtdurd() |>
  mutate(
    AGEGR1 = case_when(AGE < 65 ~ "<65", AGE <= 80 ~ "65-80", TRUE ~ ">80"),
    AGEGR1N = case_when(AGE < 65 ~ 1, AGE <= 80 ~ 2, TRUE ~ 3),
    RACEN = case_match(RACE, "WHITE" ~ 1, "BLACK OR AFRICAN AMERICAN" ~ 2,
                       "AMERICAN INDIAN OR ALASKA NATIVE" ~ 6),
    ITTFL = "Y",
    SAFFL = if_else(!is.na(TRTSDT), "Y", "N")
  ) |>
  derive_vars_merged(
    dataset_add = filter(ds, DSCAT == "DISPOSITION EVENT"),
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCDECOD = DSDECOD)
  ) |>
  mutate(DISCONFL = if_else(DCDECOD != "COMPLETED", "Y", NA_character_))

# completers: attended week 8/16/24 (VISITNUM 8/10/12)
visit_flag <- function(adsl, visnum, flag) {
  derive_var_merged_exist_flag(
    adsl, dataset_add = sv, by_vars = exprs(USUBJID),
    new_var = !!sym(flag), condition = VISITNUM == !!visnum,
    false_value = "N", missing_value = "N"
  )
}
adsl <- adsl |>
  visit_flag(8, "COMP8FL") |>
  visit_flag(10, "COMP16FL") |>
  visit_flag(12, "COMP24FL")

# efficacy population: safety + at least one post-baseline ADAS-Cog and CIBIC+
has_post <- function(cat) {
  qs |> filter(QSCAT == cat, VISITNUM > 3, !is.na(QSSTRESN)) |> pull(USUBJID) |> unique()
}
adas_post <- has_post("ALZHEIMER'S DISEASE ASSESSMENT SCALE")
cibic_post <- has_post("CLINICIAN'S INTERVIEW-BASED IMPRESSION OF CHANGE (CIBIC+)")
adsl <- adsl |>
  mutate(EFFFL = if_else(SAFFL == "Y" & USUBJID %in% adas_post & USUBJID %in% cibic_post,
                         "Y", "N"))

# baseline characteristics
mmse <- qs |>
  filter(QSCAT == "MINI-MENTAL STATE") |>
  summarise(MMSETOT = sum(QSSTRESN, na.rm = TRUE), .by = USUBJID)

adsl <- adsl |>
  derive_vars_merged(
    dataset_add = filter(vs, VSTESTCD == "HEIGHT", VISITNUM == 1),
    by_vars = exprs(USUBJID), new_vars = exprs(HEIGHTBL = VSSTRESN)
  ) |>
  derive_vars_merged(
    dataset_add = filter(vs, VSTESTCD == "WEIGHT", VISITNUM == 3),
    by_vars = exprs(USUBJID), new_vars = exprs(WEIGHTBL = VSSTRESN)
  ) |>
  mutate(
    HEIGHTBL = round_sas(HEIGHTBL, 1),
    WEIGHTBL = round_sas(WEIGHTBL, 1),
    BMIBL = round_sas(WEIGHTBL / (HEIGHTBL / 100)^2, 1),
    BMIBLGR1 = case_when(BMIBL < 25 ~ "<25", BMIBL < 30 ~ "25-<30", BMIBL >= 30 ~ ">=30")
  ) |>
  derive_vars_merged(
    dataset_add = filter(sc, SCTESTCD == "EDLEVEL"),
    by_vars = exprs(USUBJID), new_vars = exprs(EDUCLVL = SCSTRESN)
  ) |>
  left_join(mmse, by = "USUBJID")

# pool sites with < 3 subjects in any arm
small_sites <- adsl |>
  count(SITEID, TRT01P) |>
  complete(SITEID, TRT01P, fill = list(n = 0)) |>
  filter(n < 3) |>
  pull(SITEID)
adsl <- mutate(adsl, SITEGR1 = if_else(SITEID %in% small_sites, "900", SITEID))

save_adam(adsl, "adsl")
message("ADSL: ", nrow(adsl), " subjects")
