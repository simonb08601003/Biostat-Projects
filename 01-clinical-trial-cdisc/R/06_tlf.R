# Demographics and TEAE tables (rtables / tern)
source("R/00_setup.R")
suppressPackageStartupMessages({ library(rtables); library(tern) })

arms <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
adsl <- readRDS("adam/adsl.rds") |>
  mutate(
    TRT01P = factor(TRT01P, levels = arms),
    TRT01A = factor(TRT01A, levels = arms),
    AGEGR1 = factor(AGEGR1, levels = c("<65", "65-80", ">80")),
    SEX = factor(SEX, levels = c("F", "M"), labels = c("Female", "Male")),
    RACE = factor(str_to_title(RACE)),
    BMIBLGR1 = factor(BMIBLGR1, levels = c("<25", "25-<30", ">=30"))
  )

lyt_dm <- basic_table(title = "Table 14-2.01 Demographic and Baseline Characteristics",
                      subtitles = "Population: Intent-to-treat") |>
  split_cols_by("TRT01P") |>
  add_overall_col("Total") |>
  add_colcounts() |>
  analyze_vars(
    vars = c("AGE", "AGEGR1", "SEX", "RACE", "BMIBL", "BMIBLGR1", "EDUCLVL", "MMSETOT"),
    var_labels = c("Age (years)", "Age group", "Sex", "Race", "Baseline BMI (kg/m2)",
                   "BMI group", "Years of education", "MMSE total"),
    .stats = c("n", "mean_sd", "median", "range", "count_fraction")
  )
tbl_dm <- build_table(lyt_dm, adsl |> filter(ITTFL == "Y"))

adae <- readRDS("adam/adae.rds") |>
  filter(SAFFL == "Y", TRTEMFL == "Y") |>
  mutate(TRTA = factor(TRTA, levels = arms))
attr(adae$AEDECOD, "label") <- "Preferred term"
adsl_saf <- adsl |> filter(SAFFL == "Y") |> mutate(TRTA = TRT01A)

lyt_ae <- basic_table(title = "Table 14-5.02 Incidence of Treatment-Emergent Adverse Events by SOC and PT",
                      subtitles = "Population: Safety") |>
  split_cols_by("TRTA") |>
  add_colcounts() |>
  analyze_num_patients(
    vars = "USUBJID",
    .stats = c("unique", "nonunique"),
    .labels = c(unique = "Subjects with at least one TEAE", nonunique = "Total number of TEAEs")
  ) |>
  split_rows_by("AEBODSYS", child_labels = "visible", nested = FALSE,
                split_fun = drop_split_levels, label_pos = "topleft",
                split_label = "System organ class") |>
  summarize_num_patients(var = "USUBJID", .stats = "unique",
                         .labels = c(unique = "Subjects with at least one event")) |>
  count_occurrences(vars = "AEDECOD", .indent_mods = -1L) |>
  append_varlabels(adae, "AEDECOD", indent = 1L)
tbl_ae <- build_table(lyt_ae, adae, alt_counts_df = adsl_saf) |>
  prune_table() |>
  sort_at_path(path = "AEBODSYS", scorefun = cont_n_allcols) |>
  sort_at_path(path = c("AEBODSYS", "*", "AEDECOD"), scorefun = score_occurrences)

writeLines(toString(tbl_dm), "outputs/t_14_2_01_demographics.txt")
writeLines(toString(tbl_ae), "outputs/t_14_5_02_teae_soc_pt.txt")
saveRDS(list(dm = tbl_dm, ae = tbl_ae), "outputs/tables.rds")
cat(toString(tbl_dm))
cat(head(strsplit(toString(tbl_ae), "\n")[[1]], 30), sep = "\n")
