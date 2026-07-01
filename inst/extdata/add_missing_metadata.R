library(tidyverse)
library(tidyxl)

# Read the excel file with newly added missing data
df <- xlsx_cells("Added metadata 20250614.xlsx")

# Find the cells that are marked in red
formats <- xlsx_formats("Added metadata 20250614.xlsx")
red_cells_format <- which(formats$local$font$color$rgb == "FFFF0000")
# Note that the cell contents are all characters
red_cells <- subset(df, local_format_id %in% red_cells_format,
                    c("row", "col", "character"))

# Find the IDs that correspond to these red cells
ids_flt <- subset(df, row %in% red_cells$row & col == 1)
ids_sel <- ids_flt[, c("row", "character")]
ids <- ids_sel$character
names(ids) <- ids_sel$row

# Split into sex, ethnicity, and year of surgery (this may need padding to be
# accepted into previously added data)
categories <- split(red_cells, red_cells$col)
categories_ids <- lapply(categories, function(x) {
  bind_cols(x, tibble(id = ids[as.character(x$row)]))
})

# Save country (local_format_id 7) in center_info.csv
if (!is.null(categories_ids$`7`)) {
  center_info <- read_csv("center_info.csv")
  missing_center_info <- categories_ids$`7`
  for (i in seq_len(nrow(missing_center_info))) {
    id <- missing_center_info$id[i]
    country_sample_collection <- missing_center_info$character[i]
    print(id)
    print("Previous entry in center_info.csv")
    print(center_info %>% filter(`Subject ID` == id) %>% pull(`Country of Sample Collection`))
    print("Updated entry")
    print(country_sample_collection)
    print("Replace entry")
    center_info[center_info$`Subject ID` == id, "Country of Sample Collection"] <- country_sample_collection
  }

  write_csv(center_info, "center_info.csv", na = "")
}

# Save ethnicity (local_format_id 6) in center_info.csv
if (!is.null(categories_ids$`6`)) {
  center_info <- read_csv("center_info.csv")
  missing_center_info <- categories_ids$`6`
  for (i in seq_len(nrow(missing_center_info))) {
    id <- missing_center_info$id[i]
    ethn_race <- missing_center_info$character[i]
    print(id)
    print("Previous entry in center_info.csv")
    print(center_info %>% filter(`Subject ID` == id) %>% pull(`Ethnicity/Race`))
    print("Updated entry")
    print(ethn_race)
    print("Replace entry")
    center_info[center_info$`Subject ID` == id, "Ethnicity/Race"] <- ethn_race
  }

  write_csv(center_info, "center_info.csv", na = "")
}

# Save year_of_surgery (local_format_id 9) in diagnosis_surgery_dates.csv
if (!is.null(categories_ids$`9`)) {
  diagnosis_surgery_dates <- read_csv("diagnosis_surgery_dates.csv")
  missing_dx_tx_dates <- categories_ids$`9`
  for (i in seq_len(nrow(missing_dx_tx_dates))) {
    id <- missing_dx_tx_dates$id[i]
    tx_year <- missing_dx_tx_dates$character[i]
    # Create a dummy date to keep the same format
    tx_date <- paste(1, 1, substr(tx_year, 3, 4), sep = "/")
    print(id)
    print("Previous entry in diagnosis_surgery_dates.csv")
    print(diagnosis_surgery_dates %>% filter(CGID == id) %>% pull(`Debulking surgery-Date`))
    print("Updated entry")
    print(tx_date)
    print("Replace entry")
    diagnosis_surgery_dates[diagnosis_surgery_dates$CGID == id, "Debulking surgery-Date"] <- tx_date
  }

  write_csv(diagnosis_surgery_dates, "diagnosis_surgery_dates.csv", na = "")
}

# Save sex (local format_id 5) in manifest.rds
if (!is.null(categories_ids$`5`)) {
  manifest <- readRDS("manifest.rds")
  missing_sex <- categories_ids$`5`
  for (i in seq_len(nrow(missing_sex))) {
    id <- missing_sex$id[i]
    sex <- missing_sex$character[i]
    print(id)
    print("Previous entry in manifest.rds")
    print(manifest %>% filter(subject_id == id) %>% pull(sex))
    print("Updated entry")
    print(sex)
    print("Replace entry")
    manifest[manifest$subject_id == id, "sex"] <- sex
  }

  saveRDS(manifest, "manifest.rds")
}
