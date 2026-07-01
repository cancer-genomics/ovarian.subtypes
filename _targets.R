# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.
# Set target options:
tar_option_set(packages = c("tidyverse", "here",
                            "lubridate", "magrittr",
                            "fs",
                            "SummarizedExperiment"), # Packages that your targets need for their tasks.
               format="rds")
## Run the R scripts in the R/ folder with your custom functions:
lapply(list.files("R", full.names=TRUE), source)
## Columns stripped from published .rda objects before saving (Phase 2, 2026-06-23)
PHI_COLS <- c("pgdx_id", "bamfile", "bam_local", "size", "genotype_id", "facet_id")
##tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.
# Replace the target list below with your own:
list(
    tar_target(file, here("inst", "extdata", "manifest.rds"), format="file"),
    tar_target(sfile, here("inst", "extdata", "sdat.rds"), format="file"),
    tar_target(cfile, here("inst", "extdata", "center_info.csv"), format="file"),
    tar_target(dtfile, here("inst", "extdata", "diagnosis_surgery_dates.csv"), format="file"),
    tar_target(manifest0, readRDS(file)),
    tar_target(sdat, read_sdata(sfile)),
    tar_target(countrycoll, read_countrycoll(cfile)),
    tar_target(dx_tx_dates, read_dx_tx_dates(dtfile)),
    tar_target(country_dx_tx, join_dx_tx_dates(countrycoll, dx_tx_dates)),
    tar_target(facets.file, file.path("..", "output", "facets",
                                      "merge-facets-tables.R",
                                      "all-segments.txt")),
    tar_target(facets, read_facets(facets.file)),
    tar_target(facets.dir, file.path("..", "output", "facets-trellis/jhpce_directories")),
    tar_target(directory.listing, dir_ls(facets.dir, type="directory")),
    tar_target(manifest2, subject_id2(manifest0)),
    tar_target(cdat2, clean_clinical_data(sdat, manifest2, country_dx_tx)),
    tar_target(manifest.list, clean_manifest(manifest2, cdat2)),
    tar_target(manifest8, join_with_facets(manifest.list, facets, directory.listing)),
    tar_target(manifest9, filter_discordant_tumors(manifest8)),
    tar_target(clinical, filter(cdat2, lab_id %in% manifest9$lab_id)),
    tar_target(discordant, filter(manifest8, discordant_tumor_type)),
    tar_target(facets.file1, file.path("..", "output", "facets",
                                       "merge-facets-tables.R",
                                       "summary-stats.txt")),
    tar_target(facets.file2, file.path("..", "output", "facets-trellis", "summary-stats.txt")),
    tar_target(facets.purity, read_facets2(facets.file1, facets.file2)),
    tar_target(manifest, add_facets_purity(manifest9, facets.purity)),
    tar_target(save.manifest,
               save_object(select(manifest, -any_of(PHI_COLS)), "manifest"),
               format="file"),
    tar_target(save.clinical, save_object(clinical, "clinical"),
               format="file"),
    tar_target(save.discordant,
               save_object(select(discordant, -any_of(PHI_COLS)), "discordant"),
               format="file"),
    tar_target(endo.manifest, subset_endo(manifest)),
    tar_target(pathway.file, file.path("..", "output",
                                       "gene_pathway.rmd",
                                       "gene_pathway.csv")),
    tar_target(pathways, read_pathways(pathway.file, "endometrioid")),
    tar_target(idat.file, file.path("..", "output",
                                    "05-data_integration.rmd",
                                    "integrated_data.rds")),
    tar_target(idat, read_idat(idat.file, endo.manifest, pathways)),
    tar_target(idat.endometrioid, order_endo(idat, endo.manifest)),
    tar_target(save.idat.endo, save_object(idat.endometrioid, "idat.endometrioid"),
               format="file"),
    tar_target(gi.manifest, subset_gi(manifest)),
    tar_target(gi.pathways, read_pathways(pathway.file, "mucinous")),
    tar_target(idat.gi, read_idat(idat.file, gi.manifest, gi.pathways)),
    tar_target(save.idat.gi, save_object(idat.gi, "idat.gi")),
    tar_target(muc.manifest, subset_mucinous(manifest)),
    tar_target(muc.pathways.file, file.path("..", "output", "05-data_integration.rmd",
                                            "mucinous_factor_levels.rds")),
    tar_target(muc.levels, mucinous_pathways(muc.pathways.file)),
    tar_target(idat.muc, read_idat(idat.file, muc.manifest, gi.pathways)),
    tar_target(idat.mucinous, order_mucinous(idat.muc, muc.manifest, muc.levels$pathway)),
    tar_target(save.idat.muc, save_object(idat.mucinous, "idat.mucinous")),
    ##
    ## Methylation analyses
    ##
    tar_target(meth.file, file.path("..", "output", "methylation.Rmd", "se.rds")),
    tar_target(methylation, read_methylation_se(meth.file, manifest, discordant)),
    tar_target(save.meth,
               save_object(select(methylation, -any_of(PHI_COLS)), "methylation")),
    ## summarized experiment
    tar_target(tcga.file, file.path("..", "extdata", "se_lab_tcga.rds")),
    tar_target(metadata.file, file.path("..", "extdata", "combmetadata.rds")),
    tar_target(meth.se, read_methylation_data(metadata.file, tcga.file)),
    tar_target(methylation_se_pre, check_against_manifest(meth.se, manifest, discordant)),
    tar_target(match.file, file.path("inst", "extdata", "match_table.csv"),
               format="file"),
    tar_target(match_table, read.csv(match.file)),
    ##trace(update_tcga_barcodes, browser)
    tar_target(methylation_se_with_signet, update_tcga_barcodes(methylation_se_pre, match_table)),
    tar_target(signet.file, file.path("inst", "extdata", "stomach_muc_signet.csv"),
               format="file"),
    tar_target(signet.cases, get_signetring(signet.file, match_table)),
    tar_target(methylation_se, drop_signet(methylation_se_with_signet,
                                           signet.file,
                                           match_table)),
    ## Drop signet ring cases
    tar_target(save.meth.se, {
               se <- methylation_se
               colData(se) <- colData(se)[, setdiff(colnames(colData(se)), "basename")]
               save_object(se, "methylation_se")
               })
)
