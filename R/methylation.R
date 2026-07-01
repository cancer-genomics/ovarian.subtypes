read_methylation_se <- function(file, manifest, discordant) {
  rename <- dplyr::rename
  se <- readRDS(file)
  meth <- tibble(lab_id = colnames(se)) %>%
    mutate(platform = "methylation") %>%
    filter(!lab_id %in% discordant$lab_id)
  meth2 <- select(manifest, -platform) %>%
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
    select(-alt_id)
  meth3 <- select(manifest, -platform) %>%
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
    select(-c(Diagnosis, sampletype, Sample)) %>%
    set_colnames(tolower(colnames(.))) %>%
    mutate(t.n = substr(tumor, 1, 1))
  md2 <- as(md, "DataFrame")
  colData(se) <- md2
  colnames(se) <- se$lab_id
  se
}

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
    left_join(select(mt, Sample.ID, Barcode.ID), by = "Sample.ID")
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
    select(
      subject_id, lab_id, tumor_type,
      tumor.normal
    ) %>%
    distinct()
  tumors <- filter(manifest, tumor.normal == "tumor")
  tumortypes <- tumors %>%
    select(subject_id, lab_id, tumor_type) %>%
    ungroup() %>%
    distinct()
  tt <- select(tumortypes, -lab_id) %>%
    distinct()
  meth2 <- methprop %>%
    select(Sample_Name, propmeth) %>%
    rename(lab_id = Sample_Name) %>%
    left_join(
      manifest2,
      join_by(lab_id)
    ) %>%
    select(-tumor_type) %>%
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
    select(paste0("PC", use_pcs)) %>%
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
    select(paste0("PC", use_pcs)) %>%
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
