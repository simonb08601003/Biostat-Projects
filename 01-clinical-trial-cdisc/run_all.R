# SDTM -> ADaM -> QC -> analyses -> tables
for (s in c("01_adsl", "02_adqsadas", "03_adae_adtte", "07_qc",
            "04_efficacy", "05_survival", "06_tlf")) {
  message("== ", s)
  source(file.path("R", paste0(s, ".R")), local = new.env())
}
