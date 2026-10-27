# ADAE and ADTTE (time to first dermatologic event)
source("R/00_setup.R")

adsl <- readRDS("adam/adsl.rds")

adae <- read_sdtm("ae") |>
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs(USUBJID),
    new_vars = exprs(TRTSDT, TRTEDT, RFENDT, TRTA = TRT01A, TRTAN = TRT01AN,
                     SAFFL, AGE, SEX, RACE)
  ) |>
  # only a missing day is imputed; year-only dates stay missing
  derive_vars_dt(new_vars_prefix = "AST", dtc = AESTDTC,
                 highest_imputation = "D", date_imputation = "first") |>
  derive_vars_dt(new_vars_prefix = "AEN", dtc = AEENDTC) |>
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ASTDT, AENDT)) |>
  mutate(
    TRTEMFL = if_else(!is.na(ASTDT) & ASTDT >= TRTSDT, "Y", "N"),
    # dermatologic events of special interest, per the SAP
    CQ01NAM = if_else(
      str_detect(AEDECOD, "APPLICATION|DERMATITIS|ERYTHEMA|BLISTER") |
        (AEBODSYS == "SKIN AND SUBCUTANEOUS TISSUE DISORDERS" &
           !AEDECOD %in% c("COLD SWEAT", "HYPERHIDROSIS", "ALOPECIA")),
      "DERMATOLOGIC EVENTS", NA_character_
    )
  ) |>
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(by_vars = exprs(USUBJID), order = exprs(ASTDT, AESEQ),
                  new_var = AOCC01FL, mode = "first"),
    filter = TRTEMFL == "Y" & !is.na(CQ01NAM)
  )
save_adam(adae, "adae")

# event = first dermatologic TEAE, censored at end of study
derm_event <- event_source(
  dataset_name = "adae",
  filter = AOCC01FL == "Y",
  date = ASTDT,
  set_values_to = exprs(EVNTDESC = "Dematologic Event Occured",  # spelling as in the reference ADaM
                        SRCDOM = "ADAE", SRCVAR = "ASTDT", SRCSEQ = AESEQ)
)
eos <- censor_source(
  dataset_name = "adsl",
  date = RFENDT,
  set_values_to = exprs(EVNTDESC = "Study Completion Date", SRCDOM = "ADSL", SRCVAR = "RFENDT")
)

adtte <- derive_param_tte(
  dataset_adsl = filter(adsl, SAFFL == "Y"),
  start_date = TRTSDT,
  event_conditions = list(derm_event),
  censor_conditions = list(eos),
  source_datasets = list(adsl = adsl, adae = adae),
  set_values_to = exprs(PARAMCD = "TTDE", PARAM = "Time to First Dermatologic Event")
) |>
  derive_vars_duration(new_var = AVAL, start_date = STARTDT, end_date = ADT) |>
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs(USUBJID),
    new_vars = exprs(TRTA = TRT01A, TRTAN = TRT01AN, AGE, AGEGR1, SEX, RACE, SITEGR1, SAFFL)
  )
save_adam(adtte, "adtte")

message("ADAE: ", nrow(adae), " records; ADTTE: ", nrow(adtte), " subjects, ",
        sum(adtte$CNSR == 0), " events")
