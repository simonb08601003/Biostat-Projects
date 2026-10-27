# QC: compare derived ADaM against the CDISC reference ADaM, variable by variable
source("R/00_setup.R")

same <- function(a, b) {
  if (inherits(a, "Date") || inherits(b, "Date")) { a <- as.numeric(a); b <- as.numeric(b) }
  if (is.character(a) || is.character(b)) {
    a <- na_if(as.character(a), ""); b <- na_if(as.character(b), "")
    eq <- a == b
  } else {
    eq <- abs(a - b) < 1e-6
  }
  (is.na(a) & is.na(b)) | (!is.na(eq) & eq)
}

compare <- function(dataset, ours, ref, by, vars) {
  j <- inner_join(ours, ref, by = by, suffix = c(".o", ".r"))
  tibble(
    dataset = dataset,
    variable = vars,
    records = nrow(j),
    matched = sapply(vars, \(v) sum(same(j[[paste0(v, ".o")]], j[[paste0(v, ".r")]])))
  )
}

ours <- function(x) readRDS(file.path("adam", paste0(x, ".rds")))
ref <- function(x) read_xpt(file.path("data-reference", paste0(x, ".xpt"))) |> convert_blanks_to_na()

qc <- bind_rows(
  compare("ADSL", ours("adsl"), rename(ref("adsl"), TRTDURD = TRTDUR), "USUBJID",
          c("TRT01P", "TRTSDT", "TRTEDT", "TRTDURD", "AGEGR1", "SAFFL", "ITTFL", "EFFFL",
            "COMP8FL", "COMP16FL", "COMP24FL", "DCDECOD", "HEIGHTBL", "WEIGHTBL", "BMIBL",
            "EDUCLVL", "MMSETOT", "SITEGR1")),
  compare("ADQSADAS", filter(ours("adqsadas"), ANL01FL == "Y"),
          filter(ref("adqsadas"), PARAMCD == "ACTOT", ANL01FL == "Y"),
          c("USUBJID", "AVISIT", "DTYPE"), c("AVAL", "BASE", "CHG")),
  compare("ADAE", ours("adae"), ref("adae"), c("USUBJID", "AESEQ"),
          c("ASTDT", "TRTEMFL", "CQ01NAM", "AOCC01FL")),
  compare("ADTTE", ours("adtte"), ref("adtte"), "USUBJID", c("STARTDT", "ADT", "CNSR", "AVAL"))
) |>
  mutate(match_rate = matched / records)

write.csv(qc, "outputs/qc_adam_vs_reference.csv", row.names = FALSE)
print(qc, n = Inf)
cat("\nAll variables match:", all(qc$matched == qc$records), "\n")
