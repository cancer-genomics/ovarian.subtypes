# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.
# Set target options:
tar_option_set(packages = c("tidyverse", "here", "lubridate", "magrittr", "fs"), # Packages that your targets need for their tasks.
               format="rds")
## Run the R scripts in the R/ folder with your custom functions:
lapply(list.files("R", full.names=TRUE), source)
##tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.
# Replace the target list below with your own:
list(
    tar_target(file, here("inst", "extdata", "manifest.rds")),
    tar_target(sfile, here("inst", "extdata", "sdat.rds")),
    tar_target(manifest0, readRDS(file)),
    tar_target(sdat, read_sdata(sfile)),
    tar_target(facets.file, file.path("..", "output", "facets",
                                      "merge-facets-tables.R",
                                      "all-segments.txt")),
    tar_target(facets, read_facets(facets.file)),
    tar_target(facets.dir, file.path("..", "output", "facets-trellis/jhpce_directories")),
    tar_target(directory.listing, dir_ls(facets.dir, type="directory")),
    tar_target(manifest2, subject_id2(manifest0)),
    tar_target(cdat2, clean_clinical_data(sdat, manifest2)),
    tar_target(manifest.list, clean_manifest(manifest2, cdat2)),
    tar_target(manifest8, join_with_facets(manifest.list, facets, directory.listing)),
    tar_target(manifest, filter_discordant_tumors(manifest8)),
    tar_target(clinical, filter(cdat2, lab_id %in% manifest$lab_id))
)
