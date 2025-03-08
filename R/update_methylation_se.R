# This code is used to update the methylation_se object saved in data to assign
# the actual barcode IDs to the TCGA samples instead of internal IDs
# To see how the map table was generated, see:
# match_barcodes_by_beta_corr.R in
# /dcs07/scharpf/data/skoul/Projects/ovarian_subtypes_tcga_methyl

###############################################################################
# This code is run in R directory
###############################################################################
# Do not rerun the code without the correct source data
# Last commit 4319329436b8aea5a0183205097850662c108085
if (FALSE) {
  # Load the methylation_se object and save it with different name for
  # "easy" comparison
  library(SummarizedExperiment)
  load("../data/methylation_se.rda")
  save(methylation_se, file = "../data/methylation_se_v0.rda")

  # Load the match_table.csv downloaded from the cluster
  match_table <- read.csv("../data/match_table.csv")
  barcodes <- match_table$Barcode.ID
  sum(is.na(barcodes)) # Only one is missing - check!
  # [1] 1
  # TCGA failed to update this sample
  # See https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/ for
  # more information
  barcodes[is.na(barcodes)] <- "TCGA-D5-6930-01A-11D-1926-05"
  names(barcodes) <- c(1:164)
  tissue_source <- match_table$Project.ID
  for (i in seq_along(tissue_source)) {
    if (is.na(tissue_source[i])) {
      # This is the NA observation
      tissue_source[i] <- "Colorectal mucinous"
    } else if (tissue_source[i] == "TCGA-COAD") {
      tissue_source[i] <- "Colorectal mucinous"
    } else if (tissue_source[i] == "TCGA-PAAD") {
      tissue_source[i] <- "Pancreas mucinous"
    } else if (tissue_source[i] == "TCGA-STAD") {
      tissue_source[i] <- "Stomach mucinous"
    } else if (tissue_source[i] == "TCGA-UCEC") {
      tissue_source[i] <- "Uterine endometrioid"
    } else {
      stop("Error! Wrong tissue source.")
      print(tissue_source[i])
    }
  }
  tissue_type <- match_table$Sample.Type
  tissue_type <- gsub("Primary ", "", tissue_type)
  tissue_type <- gsub("Solid Tissue ", "", tissue_type)
  tissue_type[is.na(tissue_type)] <- "Tumor"
  tissue_type_short <- sapply(tissue_type, function(x) substr(x, 1, 1))

  # Update the colData and colnames for the TCGA samples
  # Columns to update
  # lab_id, diagnosis, tumor, t.n
  # COAD = Colorectal mucinous
  # PAAD = Pancreas mucinous
  # STAD = Stomach mucinous
  # UCEC = Uterine endometrioid
  allcols <- colnames(methylation_se)
  levels(colData(methylation_se)$diagnosis) <- c("Uterine endometrioid",
                                                 levels(colData(methylation_se)$diagnosis))
  newlevels <- levels(colData(methylation_se)$diagnosis)[-match("Uterine endometrial", levels(colData(methylation_se)$diagnosis))]
  colData(methylation_se)$diagnosis[match(c(1:164), allcols)] <- tissue_source
  colData(methylation_se)$diagnosis[colData(methylation_se)$diagnosis == "Uterine endometrial"] <- "Uterine endometrioid"
  colData(methylation_se)$diagnosis <- factor(colData(methylation_se)$diagnosis, levels = newlevels)
  colData(methylation_se)$tumor[match(c(1:164), allcols)] <- tissue_type
  colData(methylation_se)$t.n[match(c(1:164), allcols)] <- tissue_type_short
  colData(methylation_se)$lab_id[match(c(1:164), allcols)] <- barcodes
  colnames(methylation_se)[match(c(1:164), allcols)] <- barcodes
  save(methylation_se, file = "../data/methylation_se.rda")
}
