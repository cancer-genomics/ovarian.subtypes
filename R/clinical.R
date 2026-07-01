clean_colnames <- function(x) {
  rename <- dplyr::rename
  nms <- colnames(x)
  nms2 <- nms %>%
    tolower(.) %>%
    str_replace_all(., "(\\.)\\1+", "_") %>%
    str_replace_all(., "\\.", "_") %>%
    str_replace_all(., "_$", "")
  colnames(x) <- nms2
  x <- x %>%
    rename(
      years = years_from_diagnosis,
      overall_survival = overall_survival_status_0_alive_1_dead,
      histology = histological_tumor_type
    ) %>%
    mutate(histology = str_replace_all(histology, "ovarian", "Ovarian"))
  x
}

.tumor_clinical_data <- function(obj) {
  if (nrow(obj) == 1) {
    return(obj)
  }
  ntypes <- length(unique(obj$tumor_type))
  if (ntypes > 1) {
    obj$discordant_tumor_type <- TRUE
  }
  obj2 <- filter(obj, tumor.normal == "tumor") %>%
    distinct()
  if (nrow(obj2) == 1) {
    return(obj2)
  }
  obj2
}

clinical_data <- function(x) {
  x2 <- select(
    x,
    lab_id,
    subject_id2,
    tumor.normal,
    tumor_type,
    sex,
    age_at_diagnosis_or_surgery,
    stage_at_first_diagnosis,
    smoker,
    discordant_tumor_type
  )
  x3 <- ungroup(x2) %>%
    group_by(subject_id2) %>%
    nest()
  x3$data <- x3$data %>%
    map(.tumor_clinical_data)
  x4 <- unnest(x3, "data")
  x4
}

clean_sdata <- function(sdat) {
  molecular <- sdat[, 17:26]
  sdat2 <- sdat[, 1:16]
  sdat3 <- sdat2[, -c(3, 4, 6, 11, 16)]
  varnames <- c(
    "lab_id", "contact",
    "age_dx",
    "age_surgery",
    "sex",
    "tumor_type",
    "stage",
    "pfs",
    "is_alive",
    "os",
    "days_dx"
  )
  colnames(sdat3) <- varnames
  age.sx <- sdat3$age_surgery
  age.sx[age.sx == ""] <- NA
  age.sx <- sapply(strsplit(sdat3$age_surgery, "/"), "[", 1)
  age.sx <- as.integer(age.sx)
  sdat3$age_surgery <- age.sx
  sdat3
}

clindata_description <- function(varnames) {
  clindata.descr <- tibble(
    varname = varnames,
    description = c(
      "Unique lab identifier for the sample",
      "Contact PI providing sample",
      "Age at diagnosis (years)",
      "Age at surgery (years)",
      "Sex",
      "Histological tumor type",
      "FIGO stage (1988)",
      "Progression-free survival from diagnosis (days)",
      "Is alive (TRUE, FALSE)",
      "Overall survival from diagnosis (days)",
      "Days from diagnosis"
    )
  )
  clindata.descr
}

join_clinical_data <- function(clinical.data, sdat3) {
  cdat <- clinical.data %>%
    left_join(sdat3, join_by(lab_id, sex, tumor_type)) %>%
    mutate(age = ifelse(age_at_diagnosis_or_surgery == "NA", NA,
      as.integer(age_at_diagnosis_or_surgery)
    )) %>%
    select(-age_at_diagnosis_or_surgery) %>%
    mutate(stage = ifelse(is.na(stage), stage_at_first_diagnosis, stage)) %>%
    select(-stage_at_first_diagnosis) %>%
    mutate(age_surgery = ifelse(is.na(age_surgery), age, age_surgery)) %>%
    select(-age) %>%
    mutate(
      is_alive = ifelse(is_alive == 0, TRUE, FALSE),
      pfs = ifelse(pfs == "-", NA, as.numeric(pfs))
    )
  cdat
}

overall_survival <- function(cdat) {
  os <- cdat$os
  months <- os[grepl("m", os)] %>%
    strsplit("m") %>%
    sapply("[", 1) %>%
    as.numeric()
  days <- months * 30
  os[grepl("m", os)] <- days
  date.range <- os[grepl("/", os)] %>%
    str_split("-") %>%
    unlist() %>%
    matrix(3, 2, byrow = TRUE) %>%
    set_colnames(c("start", "end")) %>%
    as_tibble() %>%
    mutate(
      start = ymd(start),
      end = ymd(end)
    ) %>%
    mutate(days = end - start) %>%
    mutate(days = as.numeric(days))
  os[grepl("/", os)] <- date.range$days
  os <- as.numeric(os)
  return(os)
}

format_cancer_stage <- function(cdat) {
  stage <- cdat$stage
  stage[grepl("n/a", stage)] <- NA
  stage <- ifelse(stage == "NA", NA, stage) %>%
    ifelse(. == "", NA, .) %>%
    str_replace_all("\u2161", "II") %>%
    str_replace_all("\u2160", "I") %>%
    str_replace_all("1", "I") %>%
    str_replace_all("2", "II") %>%
    str_replace_all("3", "III") %>%
    str_replace_all("4", "IV") %>%
    toupper()
  stage
}

check_lab_ids <- function(sdat3, manifest2) {
  stopifnot(all(sdat3$lab_id %in% manifest2$lab_id))
  test <- select(
    manifest2, lab_id, stage_at_first_diagnosis,
    age_at_diagnosis_or_surgery, sex, tumor_type
  ) %>%
    mutate(
      age.y = age_at_diagnosis_or_surgery,
      age.y = ifelse(age.y == "NA", NA, age.y),
      age.y = as.integer(age.y)
    ) %>%
    select(-age_at_diagnosis_or_surgery)
  check <- left_join(sdat3, test, by = "lab_id") %>%
    select(
      lab_id, age_dx, age.y, age_surgery,
      sex.x, sex.y,
      tumor_type.x, tumor_type.y, stage, stage_at_first_diagnosis
    )
  TRUE
}

select_clinical_columns <- function(cdat) {
  cdat2 <- cdat %>%
    select(
      subject_id, lab_id, sex,
      tumor.normal, tumor_type,
      stage, smoker, age_dx,
      age_surgery, days_dx,
      is_alive, pfs, os, contact,
      discordant_tumor_type,
      ethn_race,
      country_sample_collection, year_of_dx, year_of_surgery
    ) %>%
    mutate(sex = ifelse(is.na(sex) & grepl("ovarian", tumor_type),
      "Female", sex
    )) %>%
    ungroup()
  cdat2
}

clean_clinical_data <- function(sdat, manifest2, country_dx_tx) {
  rename <- dplyr::rename
  clinical.data <- manifest2 %>%
    mutate(discordant_tumor_type = FALSE) %>%
    clinical_data() %>%
    rename(subject_id = subject_id2)

  sdat3 <- clean_sdata(sdat)
  clindata.descr <- clindata_description(colnames(sdat3))
  stopifnot(check_lab_ids(sdat3, manifest2))
  clinical.data$sex[clinical.data$lab_id == "CGOV104T_Rep"] <- "Female"
  stopifnot(all(sdat3$lab_id %in% clinical.data$lab_id))
  cdat <- join_clinical_data(clinical.data, sdat3)
  cdat <- join_country_dx_tx(cdat, country_dx_tx)
  cdat$os <- overall_survival(cdat)
  cdat$stage <- format_cancer_stage(cdat)
  cdat[cdat == "NA"] <- NA
  cdat2 <- select_clinical_columns(cdat)
  cdat2
}

read_sdata <- function(sfile) {
  sdat <- readRDS(sfile) %>%
    as_tibble()
  sdat
}

read_countrycoll <- function(cfile) {
  countrycoll <- read_csv(cfile)
  countrycoll <- setNames(
    countrycoll,
    c("subject_id", "ethn_race", "country_sample_collection")
  )
  countrycoll
}

read_dx_tx_dates <- function(dtfile) {
  dx_tx_dates <- read_csv(dtfile)
  dt <- dx_tx_dates %>%
    mutate(
      `Date of Diagnosis` = as.Date(`Date of Diagnosis`, format = "%m/%d/%y"),
      `Debulking surgery-Date` = as.Date(`Debulking surgery-Date`, format = "%m/%d/%y"),
      year_of_dx = year(`Date of Diagnosis`),
      year_of_surgery = year(`Debulking surgery-Date`)
    )
  dt2 <- dt %>%
    mutate(CGID = str_extract(CGID, "[A-Z]{1,10}[0-9]{1,5}"))
  dt3 <- dt2 %>%
    group_by(CGID) %>%
    summarize(
      year_of_dx = pick_notna(year_of_dx),
      year_of_surgery = pick_notna(year_of_surgery)
    )
  select(dt3, subject_id = CGID, year_of_dx, year_of_surgery)
}

join_dx_tx_dates <- function(countrycoll, dx_tx_dates) {
  country_dx_tx <- full_join(countrycoll, dx_tx_dates, by = "subject_id")
  country_dx_tx
}

join_country_dx_tx <- function(cdat, country_dx_tx) {
  cdat2 <- cdat %>%
    left_join(country_dx_tx, by = "subject_id")
  cdat2
}
