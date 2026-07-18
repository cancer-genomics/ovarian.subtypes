## Step 2 of the methylation provenance plan: regenerate se_lab_tcga.rds from
## its grandparent inputs (bValsselect.rds, combmetadata.rds, se.rds).
## Ported from archive methylation_summarized_experiment.Rmd.
build_se_lab_tcga <- function(bValsselect_file, combmetadata_file, se_jhu_file,
                               manifest) {
  dlevels <- c("Uterine endometrial", "Ovarian endometrioid",
               "Ovarian mucinous", "Colorectal mucinous",
               "Pancreas mucinous", "Stomach mucinous")
  manifest2 <- manifest %>%
    cancer_names() %>%
    dplyr::select(-tumor_type) %>%
    dplyr::rename(tumor_type = tumor, tumor = tumor.normal) %>%
    dplyr::mutate(study = "JHU", diagnosis = factor(tumor_type, dlevels),
                  tumor = Hmisc::capitalize(tumor))

  lab.and.tcga <- readRDS(bValsselect_file) %>% t()
  metadata <- readRDS(combmetadata_file) %>% tibble::as_tibble()

  df <- tibble::tibble(
    lab_id    = metadata$Sample_Name,
    diagnosis = metadata$Diagnosis,
    tissue    = metadata$Tissue,
    type      = metadata$T.N,
    study     = metadata$batch
  ) %>%
    dplyr::mutate(
      study = ifelse(study %in% 1:2, "JHU", "TCGA"),
      tumor = ifelse(type == "T", "Tumor", "Normal"),
      diagnosis = factor(diagnosis, dlevels)
    ) %>%
    dplyr::select(lab_id, diagnosis, tumor, study)

  se.lab.tcga <- SummarizedExperiment::SummarizedExperiment(
    assays  = S4Vectors::SimpleList(beta = lab.and.tcga),
    colData = df
  )
  colnames(se.lab.tcga) <- df$lab_id

  ## Step 2 passes a file path (frozen se.rds); Step 4 passes the freshly
  ## regenerated `orse` object in memory. Accept either.
  se.jhu <- if (is.character(se_jhu_file)) readRDS(se_jhu_file) else se_jhu_file
  rowindex <- rownames(se.jhu) %in% rownames(se.lab.tcga)
  colindex <- which(!colnames(se.jhu) %in% colnames(se.lab.tcga))
  se.jhu2 <- se.jhu[rowindex, colindex]
  SummarizedExperiment::assays(se.jhu2) <- S4Vectors::SimpleList(
    beta = SummarizedExperiment::assay(se.jhu2) / 1000
  )
  se.jhu3 <- se.jhu2[rownames(se.lab.tcga), ]
  df.addl.jhu <- dplyr::filter(manifest2, lab_id %in% colnames(se.jhu3)) %>%
    dplyr::select(lab_id, diagnosis, tumor, study) %>%
    dplyr::arrange(match(lab_id, colnames(se.jhu3)))
  SummarizedExperiment::colData(se.jhu3) <- S4Vectors::DataFrame(df.addl.jhu,
                                                                   row.names = colnames(se.jhu3))

  cbind(se.lab.tcga, se.jhu3)
}

## ---------------------------------------------------------------------------
## Step 4 of the methylation provenance plan: regenerate the JHU `orse`
## (output/methylation.Rmd/se.rds, the grandparent input to build_se_lab_tcga)
## from the two batches' bVals/mVals matrices + sample sheets + frozen manifest.
##
## Ported verbatim from code/methylation.Rmd (the manual Step-4 producer). The
## SE-building primitives (createseobject/createoriginalSE and the merge/clean
## helpers) already live in R/assay_functions.R; these functions are the thin
## orchestration the Rmd performed inline. CLUSTER-ONLY: every input below is
## on JHPCE, not in the local checkout.
## ---------------------------------------------------------------------------

## Integer-encode a beta/M matrix, preserving dimnames. Verbatim from the
## methylation.Rmd `list_assays` chunk (defined inline there). beta uses
## scale=1000, M scale=100, the binary assays scale=1 -- build_se_lab_tcga()
## later divides the beta assay back by 1000, so this scaling is load-bearing.
tointeger <- function(x, scale = 1) {
  dns <- dimnames(x)
  if (scale == 1) {
    xi <- matrix(as.integer(x), nrow(x), ncol(x))
  } else {
    xx <- round(x * scale, 0)
    xi <- matrix(as.integer(xx), nrow(x), ncol(x))
  }
  dimnames(xi) <- dns
  xi
}

## Build the batch-2 `targets` sample table (the methylation.Rmd `create_targets`
## chunk, which its header flags "does not work"). Ported faithfully -- the
## historical failure mode is unknown and can only be diagnosed against the live
## batch-2 sample sheet on the cluster. VERIFY ON CLUSTER: (a) read.metharray.sheet
## finds exactly one CSV under baseDir and parses it; (b) the merge with
## methdat082620.csv preserves the expected batch-2 samples; (c) the resulting
## column names/order line up with the batch-1 `targets1` so mergetargets()'s
## rbind() succeeds (mergetargets reorders targets1 columns to match `targets`).
build_batch2_targets <- function(baseDir, methdat_file) {
  targets <- minfi::read.metharray.sheet(baseDir)
  sampdat <- utils::read.csv(methdat_file)
  sampdat <- sampdat[, c(1, 3:7)]
  colnames(sampdat)[1] <- "Sample_Name"
  targets <- merge(targets, sampdat, by = "Sample_Name")
  targets$sampletype <- paste0(targets$Diagnosis, "_", targets$T.N)
  tibble::as_tibble(targets)
}

## Regenerate `orse` exactly as code/methylation.Rmd's list_assays + results
## chunks do. Returns the object that was saved to se.rds (assays trimmed to
## beta+M, integer-encoded). The intermediate tumor-only `se` the Rmd also
## builds is dead code for se.rds (only fed the eval=FALSE proptest) and is
## omitted. NOTE cutoff.beta = 0.3 -- the Rmd overrides the 0.2 default.
build_orse <- function(bVals081820_file, mVals081820_file,
                       bVals_file, mVals_file, targets1_file,
                       methmanifest_file, baseDir, methdat_file) {
  ## Step 4 passes file paths (frozen matrices); Step 5 passes the freshly
  ## regenerated matrices/targets in memory. Accept either.
  rp <- function(x) if (is.character(x)) readRDS(x) else x
  bVals    <- rp(bVals081820_file)   # batch 2 beta
  mVals    <- rp(mVals081820_file)   # batch 2 M
  ann850k  <- minfi::getAnnotation(
    IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  mVals1   <- rp(mVals_file)         # batch 1 M
  bVals1   <- rp(bVals_file)         # batch 1 beta
  targets1 <- rp(targets1_file)      # batch 1 targets (raw read.metharray.sheet)
  man1     <- rp(methmanifest_file)  # frozen; no code producer
  targets  <- build_batch2_targets(baseDir, methdat_file)

  assays <- list(bVals, mVals, ann850k, baseDir,
                 mVals1, bVals1, targets1, targets, man1)
  names(assays) <- c("bVals", "mVals", "ann850k", "baseDir",
                     "mVals1", "bVals1", "targets1", "targets", "man1")

  assays2 <- createseobject(assays, cutoff.beta = 0.3, cutoff.M = -1.5)
  orse <- createoriginalSE(assays2)

  SummarizedExperiment::assays(orse)[[3]] <- tointeger(SummarizedExperiment::assays(orse)[[3]])
  SummarizedExperiment::assays(orse)[[4]] <- tointeger(SummarizedExperiment::assays(orse)[[4]])
  SummarizedExperiment::assays(orse)[[1]] <- tointeger(SummarizedExperiment::assays(orse)[[1]], scale = 1000)
  SummarizedExperiment::assays(orse)[[2]] <- tointeger(SummarizedExperiment::assays(orse)[[2]], scale = 100)

  SummarizedExperiment::assays(orse) <- SummarizedExperiment::assays(orse)[1:2]
  orse
}

## ---------------------------------------------------------------------------
## Step 5 of the methylation provenance plan: regenerate the bVals/mVals
## matrices from the JHU IDATs. Ported from the recovered batch scripts
## (code/methylation/{Methylation_data_script.R, differential_methylation_061119.R,
## preprocess_meth_0717.R}). Both batches share ONE filter chain (below); they
## differ only at I/O (batch 1 uses read.metharray.exp(force=TRUE); batch 2 does
## not). No set.seed is needed -- the matrices are deterministic given a fixed
## minfi/annotation version. CLUSTER-ONLY. minfi is NOT bit-reproducible across
## versions, so Step-5 matrix comparisons use TOLERANCE (see compare_matrix()),
## and pinned minfi/annotation versions must be recorded before running.
## ---------------------------------------------------------------------------

## Shared cold path: funnorm -> detP(<0.01 in ALL samples) -> drop chrX/chrY ->
## dropLociWithSnps(defaults) -> drop cross-reactive probes -> getBeta/getM.
## Filter ORDER is load-bearing for probe membership. `xreactive_file` is the
## 450K 48639-non-specific-probes-Illumina450k.csv from the
## methylationArrayAnalysis package extdata (applied to EPIC data in the
## original, faithfully reproduced here).
funnorm_filter_matrices <- function(RGSet, xreactive_file) {
  ann850k <- minfi::getAnnotation(
    IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  detP <- minfi::detectionP(RGSet)
  mSetSF <- minfi::preprocessFunnorm(RGSet)                    # default args
  detP <- detP[match(minfi::featureNames(mSetSF), rownames(detP)), ]
  keep <- rowSums(detP < 0.01) == ncol(mSetSF)                 # p<0.01 in all samples
  mSetSqFlt <- mSetSF[keep, ]
  sexprobes <- ann850k$Name[ann850k$chr %in% c("chrX", "chrY")]
  mSetSqFlt <- mSetSqFlt[!(minfi::featureNames(mSetSqFlt) %in% sexprobes), ]
  mSetSqFlt <- minfi::dropLociWithSnps(mSetSqFlt)              # default args
  xReactiveProbes <- utils::read.csv(xreactive_file)
  mSetSqFlt <- mSetSqFlt[!(minfi::featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID), ]
  list(bVals = minfi::getBeta(mSetSqFlt), mVals = minfi::getM(mSetSqFlt))
}

## Batch 1 (endomuc/methylation IDATs): read.metharray.exp(force=TRUE). Returns
## the raw sheet `targets` (=targets.rds, consumed by build_orse) alongside the
## matrices. WATCH-OUT: createseobject()'s clean.targets() applies POSITIONAL
## fixes (rows 3/12/31/38), so the regenerated sheet row order must match the
## frozen targets.rds -- verify on cluster.
build_batch1_matrices <- function(baseDir, xreactive_file) {
  targets <- minfi::read.metharray.sheet(baseDir)
  RGSet <- minfi::read.metharray.exp(targets = targets, force = TRUE)
  ## The frozen batch-1 bVals.rds/mVals.rds carry CLEANED CG-lab-ID colnames
  ## (e.g. CGCRC330N, with positional fixes CGCRC330Tl->CGCRC330T1 etc.).
  ## createseobject() consumes the batch-1 MATRICES already cleaned -- it only
  ## re-cleans targets1 -- so mergebVals()/mergemVals() align columns by the
  ## cleaned Sample_Name. Name the matrix columns with clean.targets() here to
  ## reproduce that invariant; otherwise the raw N_CGID_slide_array names fail
  ## the match() column selection in mergebVals() ("undefined columns selected").
  ## The RAW sheet is still returned as `targets` (createseobject cleans it).
  ## clean.targets() applies POSITIONAL fixes (rows 3/12/31/38), so this relies
  ## on read.metharray.sheet() preserving the canonical sheet row order.
  minfi::sampleNames(RGSet) <- clean.targets(targets)$Sample_Name
  mats <- funnorm_filter_matrices(RGSet, xreactive_file)
  c(list(targets = targets), mats)
}

## Batch 2 (meth_081720 IDATs): read.metharray.exp() with NO force. Batch-2
## targets are built separately by build_batch2_targets() inside build_orse, so
## only the matrices are returned here.
build_batch2_matrices <- function(baseDir, xreactive_file) {
  targets <- minfi::read.metharray.sheet(baseDir)
  RGSet <- minfi::read.metharray.exp(targets = targets)
  minfi::sampleNames(RGSet) <- targets$Sample_Name
  funnorm_filter_matrices(RGSet, xreactive_file)
}

read_methylation_se <- function(file, manifest, discordant) {
  rename <- dplyr::rename
  se <- readRDS(file)
  meth <- tibble(lab_id = colnames(se)) %>%
    mutate(platform = "methylation") %>%
    filter(!lab_id %in% discordant$lab_id)
  meth2 <- dplyr::select(manifest, -platform) %>%
    left_join(meth, by = "lab_id") %>%
    filter(!is.na(platform))
  any(!meth$lab_id %in% meth2$lab_id)
  meth.fuzzymatch <- filter(meth, !lab_id %in% meth2$lab_id) %>%
    rename(alt_id = lab_id) %>%
    mutate(lab_id = c(
      "CGCRC330N_1",
      "CGCRC330T_1",
      "CGCRC330T1_1",
      "CGOV177T_2",
      "CGOV179T_Rpt",
      "CGOV186T_2",
      "CGOV188N_2",
      "CGOV188T_2",
      "CGST1N_1",
      "CGST1T_2",
      "CGST2T_2"
    )) %>%
    dplyr::select(-alt_id)
  meth3 <- dplyr::select(manifest, -platform) %>%
    left_join(meth.fuzzymatch, by = "lab_id") %>%
    filter(!is.na(platform))
  meth4 <- bind_rows(meth2, meth3)
  methylation <- meth4
  methylation
}

read_methylation_data <- function(file, tcga.file) {
  rename <- dplyr::rename
  se <- readRDS(tcga.file)
  metadata <- readRDS(file) %>%
    as_tibble() %>%
    filter(grepl("^C", Sample_Name)) %>%
    rename(lab_id = Sample_Name)
  metadata2 <- colData(se) %>%
    as_tibble()
  md <- left_join(
    metadata2, metadata,
    join_by(lab_id)
  ) %>%
    dplyr::select(-c(Diagnosis, sampletype, Sample)) %>%
    set_colnames(tolower(colnames(.))) %>%
    mutate(t.n = substr(tumor, 1, 1))
  md2 <- as(md, "DataFrame")
  colData(se) <- md2
  colnames(se) <- se$lab_id
  se
}

#' @importFrom BiocGenerics cbind
check_against_manifest <- function(se, manifest, discordant) {
  is_jhu <- se$study == "JHU"
  jhu <- se[, is_jhu]
  in_manifest <- colnames(jhu) %in% manifest$lab_id
  notin_manifest <- colnames(jhu)[!in_manifest]
  notin_manifest2 <- notin_manifest[!notin_manifest %in% discordant$lab_id]
  dat <- tibble(
    to_map = notin_manifest2,
    lab_id = NA
  )
  for (i in seq_len(nrow(dat))) {
    id <- dat$to_map[i]
    if (id == "CGCRC330T") {
      dat$lab_id[i] <- "CGCRC330T_1"
      next()
    }
    ix <- grep(id, manifest$lab_id)
    if (length(ix) == 1) {
      dat$lab_id[i] <- manifest$lab_id[ix]
      next()
    }
    stop()
  }
  ix <- match(dat$to_map, colnames(jhu))
  colnames(jhu)[ix] <- dat$lab_id
  jhu$lab_id <- colnames(jhu)
  jhu2 <- jhu[, !colnames(jhu) %in% discordant$lab_id]
  stopifnot(all(colnames(jhu2) %in% manifest$lab_id))

  tcga <- se[, se$study == "TCGA"]
  methylation_se <- cbind(jhu2, tcga)
  methylation_se
}

update_tcga_barcodes <- function(meth.se, match_table) {
  barcodes <- match_table$Barcode.ID
  sum(is.na(barcodes))
  barcodes[is.na(barcodes)] <- "TCGA-D5-6930-01A-11D-1926-05"
  names(barcodes) <- c(1:164)
  tissue_source <- match_table$Project.ID
  for (i in seq_along(tissue_source)) {
    if (is.na(tissue_source[i])) {
      tissue_source[i] <- "Colorectal mucinous"
    } else if (tissue_source[i] == "TCGA-COAD") {
      tissue_source[i] <- "Colorectal mucinous"
    } else if (tissue_source[i] == "TCGA-PAAD") {
      tissue_source[i] <- "Pancreatic mucinous"
    } else if (tissue_source[i] == "TCGA-STAD") {
      tissue_source[i] <- "Stomach mucinous"
    } else if (tissue_source[i] == "TCGA-UCEC") {
      tissue_source[i] <- "Uterine endometrial"
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

  allcols <- colnames(meth.se)
  dx_levels <- levels(colData(meth.se)$diagnosis)
  new_dx_levels <- dx_levels
  new_dx_levels[match("Pancreas mucinous", new_dx_levels)] <- "Pancreatic mucinous"
  colData(meth.se)$diagnosis <- as.character(colData(meth.se)$diagnosis)
  is_panc_muc <- colData(meth.se)$diagnosis == "Pancreas mucinous"
  colData(meth.se)$diagnosis[is_panc_muc] <- "Pancreatic mucinous"
  colData(meth.se)$diagnosis[match(c(1:164), allcols)] <- tissue_source
  colData(meth.se)$diagnosis <- factor(colData(meth.se)$diagnosis,
    levels = new_dx_levels
  )
  colData(meth.se)$tumor[match(c(1:164), allcols)] <- tissue_type
  colData(meth.se)$t.n[match(c(1:164), allcols)] <- tissue_type_short
  colData(meth.se)$lab_id[match(c(1:164), allcols)] <- barcodes
  colnames(meth.se)[match(c(1:164), allcols)] <- barcodes
  return(meth.se)
}

get_signetring <- function(signet.file, mt) {
  signet.cases <- read_csv(signet.file, show_col_types = FALSE) %>%
    filter(mucinous == 0) %>%
    mutate(sample.id = paste0(sample.id, "A")) %>%
    dplyr::rename(Sample.ID = sample.id) %>%
    left_join(dplyr::select(mt, Sample.ID, Barcode.ID), by = "Sample.ID")
  signet.cases
}

drop_signet <- function(meth.se, signet.file, match_table) {
  mt <- as_tibble(match_table)
  signet.cases <- get_signetring(signet.file, mt)
  meth.se2 <- meth.se[, !colnames(meth.se) %in% signet.cases$Barcode.ID]
  return(meth.se2)
}

#' @export
pairedMeth <- function(methprop, manifest) {
  manifest2 <- manifest %>%
    dplyr::select(
      subject_id, lab_id, tumor_type,
      tumor.normal
    ) %>%
    distinct()
  tumors <- filter(manifest, tumor.normal == "tumor")
  tumortypes <- tumors %>%
    dplyr::select(subject_id, lab_id, tumor_type) %>%
    ungroup() %>%
    distinct()
  tt <- dplyr::select(tumortypes, -lab_id) %>%
    distinct()
  meth2 <- methprop %>%
    dplyr::select(Sample_Name, propmeth) %>%
    rename(lab_id = Sample_Name) %>%
    left_join(
      manifest2,
      join_by(lab_id)
    ) %>%
    dplyr::select(-tumor_type) %>%
    left_join(tt, by = "subject_id") %>%
    group_by(tumor_type) %>%
    nest()
  meth.matrix.list <- meth2$data %>%
    map(tumor_normal_matrix)
  nr <- sapply(meth.matrix.list, length)
  meth.matrix.list2 <- meth.matrix.list[nr > 0]
  meth.matrix.list2
}

tumor_normal_matrix <- function(x) {
  x.nested <- x %>%
    group_by(subject_id) %>%
    nest()
  nr <- map_int(x.nested$data, nrow)
  if (length(nr) < 4) {
    return(NULL)
  }
  x.nested2 <- x.nested[nr == 2, ]
  x.nested2$data %>%
    map_dfr(function(x) x) %>%
    pull(propmeth) %>%
    matrix(nc = 2, byrow = TRUE)
}

#' @export
project.cancer <- function(se.jhu, pc.tcga, ld.tcga, cancertype, use_pcs = 1:5) {
  se.jhu2 <- se.jhu[, se.jhu$diagnosis == cancertype]
  jhu.meth <- assays(se.jhu2)[[1]] %>%
    t()
  jhu.pcs <- predict(pc.tcga, newdata = jhu.meth) %>%
    as_tibble() %>%
    dplyr::select(paste0("PC", use_pcs)) %>%
    mutate(dx = se.jhu2$diagnosis)
  jhu.class.predictions <- predict(ld.tcga, newdata = jhu.pcs)
  jhu.x <- jhu.class.predictions$x[, c("LD1", "LD2")] %>%
    as_tibble() %>%
    mutate(
      dx = as.character(se.jhu2$diagnosis),
      tumor = factor(se.jhu2$tumor, c("Normal", "Tumor"))
    ) %>%
    rename(
      Groups = dx,
      tumor.normal = tumor
    ) %>%
    mutate(lab = "JHU")
  jhu.x
}

#' @export
project.cancer.prob <- function(se.jhu, pc.tcga, ld.tcga, cancertype, use_pcs = 1:5) {
  se.jhu2 <- se.jhu[, se.jhu$diagnosis == cancertype]
  jhu.meth <- assays(se.jhu2)[[1]] %>%
    t()
  jhu.pcs <- predict(pc.tcga, newdata = jhu.meth) %>%
    as_tibble() %>%
    dplyr::select(paste0("PC", use_pcs)) %>%
    mutate(dx = se.jhu2$diagnosis)
  jhu.class.predictions <- predict(ld.tcga, newdata = jhu.pcs)
  jhu.x <- jhu.class.predictions$posterior %>%
    as_tibble() %>%
    mutate(
      dx = as.character(se.jhu2$diagnosis),
      tumor = factor(se.jhu2$tumor, c("Normal", "Tumor"))
    ) %>%
    rename(
      Groups = dx,
      tumor.normal = tumor
    ) %>%
    mutate(lab = "JHU")
  jhu.x
}

#' Wrapper for principal component analysis of SummarizedExperiment object
#'
#' @export
mypca <- function(se, scale = FALSE, center = TRUE, rk) {
  x <- t(assays(se)[[1]])
  prcomp(x, scale = scale, center = center, rank. = rk)
}

#' Prepare TCGA and JHU methylation subsets for LDA training
#'
#' Reads signet ring exclusion IDs, splits methylation_se by study, removes
#' signet ring cases from TCGA, and optionally excludes a diagnosis and
#' restricts JHU to tumor samples only.
#'
#' @return Named list with elements \code{se.tcga} and \code{se.jhu}.
#' @export
filter_lda_samples <- function(methylation_se, signet_ring_file,
                               exclude_tcga_diagnosis = "Pancreatic mucinous",
                               jhu_tumor_only = TRUE) {
  pull_id <- function(x) {
    stringr::str_extract(x, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[A-Z0-9]{2}")
  }
  signet_ring <- readr::read_csv(signet_ring_file, show_col_types = FALSE) %>%
    dplyr::filter(mucinous == 0) %>%
    dplyr::pull(sample.id)

  se.jhu <- methylation_se[, methylation_se$study == "JHU"]
  if (jhu_tumor_only) se.jhu <- se.jhu[, se.jhu$tumor == "Tumor"]

  se.tcga <- methylation_se[, methylation_se$study == "TCGA"]
  se.tcga <- se.tcga[, !(pull_id(se.tcga$lab_id) %in% signet_ring)]
  if (!is.null(exclude_tcga_diagnosis)) {
    se.tcga <- se.tcga[, !se.tcga$diagnosis %in% exclude_tcga_diagnosis]
  }

  list(se.tcga = se.tcga, se.jhu = se.jhu)
}

#' Fit PCA + LDA model on a TCGA methylation SummarizedExperiment
#'
#' Runs PCA on the assay matrix, builds a feature tibble from the top
#' \code{n_pcs} principal components, and fits a linear discriminant model
#' predicting cancer diagnosis.
#'
#' @return Named list: \code{pc} (prcomp object), \code{ld} (lda object),
#'   \code{features} (tibble of PC scores + dx label used to train the LDA).
#' @export
fit_lda_model <- function(se, n_pcs = 5L) {
  pc <- mypca(se, rk = n_pcs)
  features <- pc$x[, seq_len(n_pcs)] %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dx = se$diagnosis) %>%
    dplyr::mutate(dx = droplevels(dx))
  ld <- MASS::lda(dx ~ ., features)
  list(pc = pc, ld = ld, features = features)
}

#' Project TCGA samples into LDA space and compute ellipses
#'
#' @return Named list: \code{obs} (tibble of LD scores, Groups, tumor.normal,
#'   lab_id, lab="TCGA"), \code{ell} (ellipse polygons), \code{axes} (axis
#'   label strings with variance explained).
#' @export
project_tcga_to_lda <- function(ld, features, se) {
  obs <- my.ggord.lda(ld, features$dx) %>%
    dplyr::mutate(tumor.normal = se$tumor, lab_id = colnames(se))
  list(obs = obs, ell = my.ellipse(obs), axes = axis.labels(ld))
}

#' Project JHU cancer types into TCGA LDA space
#'
#' For each cancer type in \code{cancer_types}, projects JHU samples and
#' appends \code{tcga_obs} as background. Returns a combined tibble suitable
#' for faceted visualisation, with a final TCGA background group
#' (\code{Groups = "TCGA"}). The \code{Groups} column is returned as a factor
#' with levels ordered for visualisation.
#'
#' @return Tibble with columns LD1, LD2, Groups (factor), tumor.normal, lab,
#'   lab_id.
#' @export
project_jhu_to_lda <- function(se.jhu, pc, ld, tcga_obs, cancer_types) {
  group_levels <- c("TCGA", "Colorectal mucinous", "Ovarian endometrioid",
                    "Uterine endometrial", "Ovarian mucinous",
                    "Pancreatic mucinous", "Stomach mucinous")
  jhu_proj <- lapply(cancer_types, function(ct) {
    project.cancer(se.jhu, pc, ld, ct) %>%
      dplyr::mutate(lab_id = colnames(se.jhu)[se.jhu$diagnosis == ct]) %>%
      dplyr::bind_rows(tcga_obs)
  })
  jhu_all <- dplyr::bind_rows(jhu_proj)
  tcga_background <- dplyr::mutate(tcga_obs, Groups = "TCGA")
  dplyr::bind_rows(jhu_all, tcga_background) %>%
    dplyr::mutate(Groups = factor(Groups, levels = group_levels))
}

#' Compute LDA posterior probabilities for JHU cancer types
#'
#' For each cancer type in \code{cancer_types}, calls
#' \code{project.cancer.prob} and attaches \code{lab_id} from the sample
#' names.
#'
#' @return Tibble with posterior probability columns, Groups, tumor.normal,
#'   lab, lab_id.
#' @export
project_jhu_posteriors <- function(se.jhu, pc, ld, cancer_types) {
  lapply(cancer_types, function(ct) {
    project.cancer.prob(se.jhu, pc, ld, ct) %>%
      dplyr::mutate(lab_id = colnames(se.jhu)[se.jhu$diagnosis == ct])
  }) %>%
    dplyr::bind_rows()
}
