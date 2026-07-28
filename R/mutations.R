jags_data <- function(gdat) {
  gdat <- gdat %>%
    readRDS() %>%
    as_tibble() %>%
    set_colnames(gsub(" ", "_", colnames(.))) %>%
    set_colnames(gsub("\\.", "_", colnames(.))) %>%
    set_colnames(tolower(colnames(.)))
  ngenes <- length(levels(gdat$gene_symbol))
  nid <- length(levels(gdat$internal_id))
  ntumors <- gdat %>%
    group_by(tumor_type) %>%
    summarize(n = length(unique(internal_id)))
  gdat2 <- gdat %>%
    mutate(
      mutation = ifelse(!is.na(mutation), 1L, 0L),
      methylation = ifelse(!is.na(methylation), 1L, 0L),
      fusion = ifelse(!is.na(fusion), 1L, 0L),
      copynumber = ifelse(!is.na(copynumber), 1L, 0L)
    ) %>%
    mutate(
      total = mutation + methylation + fusion + copynumber,
      is_altered = ifelse(total > 0, 1L, 0L)
    ) %>%
    dplyr::select(gene_symbol, internal_id, tumor_type, is_altered) %>%
    mutate(
      tumor_type = ifelse(tumor_type == "Ovarian endometrioid",
        "ovarian",
        "uterine"
      ),
      internal_id = as.character(internal_id),
      gene_symbol = as.character(gene_symbol)
    ) %>%
    group_by(internal_id, gene_symbol) %>%
    summarize(
      tumor_type = unique(tumor_type),
      number_altered = sum(is_altered),
      is_altered = as.integer(number_altered > 0),
      .groups = "drop"
    )
  tumortype <- group_by(gdat2, internal_id) %>%
    summarize(
      tumor_type = unique(tumor_type),
      .groups = "drop"
    )
  gdat3 <- gdat2 %>%
    dplyr::select(gene_symbol, internal_id, is_altered) %>%
    spread(gene_symbol, is_altered) %>%
    left_join(tumortype, by = "internal_id")
  X <- dplyr::select(ungroup(gdat3), -internal_id) %>%
    dplyr::select(-tumor_type) %>%
    as.matrix()
  Y <- ifelse(gdat3$tumor_type == "ovarian", 1, 0)
  is_hypermutator <- X[, "hypermutator"] == 1
  X <- X[, -match("hypermutator", colnames(X))]
  Y <- Y[!is_hypermutator]
  X <- X[!is_hypermutator, ]
  X <- X[, colSums(X) >= 3]
  X <- cbind(1, X)
  colnames(X)[1] <- "intercept"
  list(X = X, Y = Y)
}

contingency_table <- function(x, Ns) {
  if (sum(x$number_altered) == 0) {
    return(NULL)
  }
  m <- x %>%
    mutate(tumor_type = as.character(tumor_type)) %>%
    group_by(tumor_type, is_altered) %>%
    summarize(
      n = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = "is_altered", values_from = "n") %>%
    set_colnames(c("tumor_type", "wt", "mt")) %>%
    mutate(n = mt + wt)
  m2 <- left_join(Ns, m, by = c("tumor_type", "n")) %>%
    dplyr::select(tumor_type, mt, n)
  m2[is.na(m2)] <- 0
  m <- dplyr::select(m2, -tumor_type) %>%
    as.matrix()
  rownames(m) <- m2$tumor_type
  m
}

stan_inputs <- function(ct) {
  y <- ct[, 1]
  n <- ct[, 2]
  J <- length(y)
  X <- diag(4)
  X <- X[, -1]
  list(y = y, n = n, J = J, K = 6, X = X)
}

#' Create list of y, n, J, and x for stan
#'
#' @export
inputs_endo <- function(ct) {
  y <- ct[1:2, 1]
  n <- ct[1:2, 2]
  J <- length(y)
  x <- c(0, 1)
  list(y = y, n = n, J = J, x = x)
}

#' Create list of y, n, J, and x for stan
#'
#' @export
inputs_mucinous <- function(ct) {
  y <- ct[3:4, 1]
  n <- ct[3:4, 2]
  J <- length(y)
  x <- c(0, 1)
  list(y = y, n = n, J = J, x = x)
}

sampling2 <- function(data, model, params, ...) {
  sampling(model,
    data = data,
    iter = params$iter,
    thin = params$thin,
    chains = params$chains,
    warmup = params$warmup,
    control = params$control, ...
  )
}

slice_params <- function(x) {
  nms <- rownames(x)
  x %>%
    as_tibble() %>%
    mutate(parameter = nms) %>%
    filter(grepl("^beta", parameter) | grepl("^theta", parameter)) %>%
    dplyr::select(
      parameter, mean, se_mean, sd, `2.5%`,
      `5%`, `50%`, `95%`, `97.5%`, `n_eff`, Rhat
    )
}

#' Create contingency table for comparing differences in mutation rates
#'
#' @export
complete_table <- function(x, Ns, tumor_order) {
  x2 <- full_join(Ns, x, by = c("tumor_type", "n"))
  x2[is.na(x2)] <- 0L
  x3 <- dplyr::select(x2, tumor_type, mt, n)
  x4 <- left_join(tumor_order, x3, by = "tumor_type")
  x5 <- as.matrix(x4[, 2:3])
  rownames(x5) <- x4$tumor_type
  x5
}

#' Stan output
#'
#' @export
stan_output <- function(data.list, model,
                        params,
                        summaryfun,
                        probs = c(
                          0.025, 0.05, 0.1,
                          0.5, 0.9, 0.95,
                          0.975
                        )) {
  tmp <- data.list %>%
    map(sampling2, model, params) %>%
    map(summaryfun, probs = probs) %>%
    map(1) %>%
    map(slice_params)
  return(tmp)
}

#' Label Pancreas, Stomach, and Colorectal mucinous cancers as GI mucinous
#'
#' @export
collapse_gi <- function(dat) {
  dat2 <- dat %>%
    mutate(tumor_type = as.character(tumor_type)) %>%
    mutate(tumor_type = case_when(
      tumor_type == "Pancreas mucinous" ~ "GI mucinous",
      tumor_type == "Stomach mucinous" ~ "GI mucinous",
      tumor_type == "Colorectal mucinous" ~ "GI mucinous",
      TRUE ~ tumor_type
    ))
  dat2
}

#' Provides full names for cancer subtypes
#'
#' @export
cancer_names <- function(x) {
  x2 <- x %>%
    mutate(tumor = Hmisc::capitalize(tumor_type))
  x3 <- x2 %>%
    mutate(tumor = case_when(
      tumor == "Colorectal" ~ "Colorectal mucinous",
      tumor == "Pancreas" ~ "Pancreas mucinous",
      tumor == "Stomach" ~ "Stomach mucinous",
      TRUE ~ tumor
    ))
  x3
}

#' Provide mucinous cancers
#' @export
muc <- function() c("Colorectal mucinous", "Ovarian mucinous", "Pancreaas mucinous", "Stomach mucinous")

#' Provide endometrioid/endometrial cancers
#' @export
endo <- function() c("Ovarian endometrioid", "Uterine endometrioid")

#' x contains lab_id
#' manifest contains lab_id and lab_id2
#' replace lab_id with lab_id2 when not equal
#' @export
swap_lab_id <- function(x, y) {
  y <- dplyr::select(y, lab_id, lab_id2)
  isf <- is.factor(x$lab_id)
  if (isf) {
    levs <- tibble(lab_id = levels(x$lab_id)) %>%
      inner_join(y, by = "lab_id")
    levs2 <- levs$lab_id2
  }
  x.y <- inner_join(x, y, by = "lab_id") %>%
    mutate(lab_id = ifelse(lab_id == lab_id2, lab_id, lab_id2)) %>%
    dplyr::select(-lab_id2)
  if (isf) {
    x.y$lab_id <- factor(x.y$lab_id, levs2)
  }
  return(x.y)
}

order_samples <- function(x, gene.levels) {
  x2 <- dplyr::select(x, lab_id, gene_symbol) %>%
    mutate(
      hypermutator = gene_symbol == "hypermutator",
      gene1.alt = gene_symbol == gene.levels[1],
      gene2.alt = gene_symbol == gene.levels[2],
      gene3.alt = gene_symbol == gene.levels[3],
      gene4.alt = gene_symbol == gene.levels[4],
      gene5.alt = gene_symbol == gene.levels[5]
    ) %>%
    group_by(lab_id) %>%
    summarize(
      hypermut = any(hypermutator),
      gene1 = any(gene1.alt),
      gene2 = any(gene2.alt),
      gene3 = any(gene3.alt),
      gene4 = any(gene4.alt),
      gene5 = any(gene5.alt)
    )
  x3 <- x2 %>%
    arrange(
      hypermut,
      -gene1,
      -gene2,
      -gene3,
      -gene4,
      -gene5
    )
  x3
}

read_pathways <- function(pathway.file, tumor.type) {
  rename <- dplyr::rename
  if (is.character(pathway.file)) {
    pathways <- read_csv(pathway.file, show_col_types = FALSE) %>%
      rename(gene_symbol = gene.symbol)
  } else {
    pathways <- rename(pathway.file, gene_symbol = gene.symbol)
  }
  pathways <- pathways %>%
    filter(tumor_type == tumor.type) %>%
    dplyr::select(-tumor_type) %>%
    mutate(
      pathway = str_replace_all(pathway, "TGFBR pathway", "TGFBR"),
      pathway = str_replace_all(pathway, "BRCA", "DNA repair")
    )
  if (tumor.type == "mucinous") {
    ix  <- which(pathways$gene_symbol == "JAK1"        & pathways$pathway == "Cell cycle")
    ix2 <- which(pathways$gene_symbol == "MED1-STAT5B" & pathways$pathway == "Other")
    pathways2 <- pathways[-c(ix, ix2), ]
  } else {
    pathways2 <- pathways
  }
  pathways2
}

read_integrated_data <- function(manifest, pathways) {
  tumortypes <- dplyr::select(manifest, lab_id, tumor_type) %>%
    ungroup() %>%
    distinct()
  idat <- here(
    "output", "01-data_integration.rmd",
    "integrated_data.rds"
  ) %>%
    readRDS() %>%
    mutate(gene_symbol = gene) %>%
    filter(lab_id %in% manifest$lab_id) %>%
    left_join(pathways, by = "gene_symbol") %>%
    left_join(tumortypes, by = "lab_id") %>%
    distinct()
  idat
}

read_idat <- function(idat.file, manifest, pathways) {
  rename <- dplyr::rename
  tumortypes <- dplyr::select(manifest, lab_id, tumor_type) %>%
    ungroup() %>%
    distinct()
  idat <- idat.file %>%
    readRDS() %>%
    filter(lab_id %in% manifest$lab_id) %>%
    left_join(pathways, by = c("gene" = "gene_symbol")) %>%
    left_join(tumortypes, by = "lab_id") %>%
    distinct() %>%
    filter(!is.na(pathway)) %>%
    mutate(alteration = ifelse(type == "mutation",
      "mutation", alteration
    )) %>%
    cancer_names() %>%
    dplyr::select(-tumor_type) %>%
    rename(tumor_type = tumor)

  idat
}

gene_list <- function(idat, pathway.levels2, tumor.levels) {
  gene.list <- idat %>%
    mutate(
      pathway = factor(pathway, pathway.levels2),
      tumor_type = factor(tumor_type, tumor.levels)
    ) %>%
    group_by(gene, pathway, tumor_type) %>%
    summarize(
      n = length(unique(lab_id)),
      .groups = "drop"
    ) %>%
    arrange(pathway, tumor_type, n) %>%
    group_by(pathway) %>%
    nest()
  gene.list
}

remove_duplicated_genes <- function(gene.list) {
  gl <- filter(gene.list, pathway != "Hypermutator") %>%
    pull(data) %>%
    map(function(x) {
      filter(x, !duplicated(gene)) %>%
        arrange(n)
    })
  gl
}

gene_levels <- function(gene.list) {
  gene.levels <- unnest(gene.list, "data") %>%
    pull(gene) %>%
    unique()
}

endo_order <- function(idat) {
  genes.for.sample.order <- idat %>%
    filter(
      pathway == "PI3K",
      tumor_type == "Ovarian endometrioid"
    ) %>%
    group_by(gene) %>%
    summarize(
      n = length(unique(lab_id)),
      .groups = "drop"
    ) %>%
    arrange(-n)
  genes.for.sample.order
}

muc_order <- function(idat) {
  genes.for.sample.order <- idat %>%
    filter(
      pathway == "Ras and TK receptors",
      tumor_type == "Ovarian mucinous"
    ) %>%
    group_by(gene) %>%
    summarize(
      n = length(unique(lab_id)),
      .groups = "drop"
    ) %>%
    arrange(-n)
}

endo_id_levels <- function(idat, genes.for.sample.order) {
  ovarian.order <- filter(idat, tumor_type == "Ovarian endometrioid") %>%
    mutate(gene_symbol = gene) %>%
    order_samples(gene.levels = genes.for.sample.order$gene)
  uterine.order <- filter(idat, tumor_type == "Uterine endometrioid") %>%
    mutate(gene_symbol = gene) %>%
    order_samples(gene.levels = genes.for.sample.order$gene)
  id.levels <- c(ovarian.order$lab_id, uterine.order$lab_id)
}

muc_id_levels <- function(idat, genes.for.sample.order) {
  ovarian.order <- filter(idat, tumor_type == "Ovarian mucinous") %>%
    mutate(gene_symbol = gene) %>%
    order_samples(gene.levels = genes.for.sample.order$gene)
  crc.order <- filter(idat, tumor_type == "Colorectal mucinous") %>%
    mutate(gene_symbol = gene) %>%
    order_samples(gene.levels = genes.for.sample.order$gene)
  id.levels <- c(ovarian.order$lab_id, crc.order$lab_id)
}

order_idat_endo <- function(idat, tumor.levels, pathway.levels, sample.order) {
  gene.list <- gene_list(idat, pathway.levels, tumor.levels)
  gl <- remove_duplicated_genes(gene.list)
  gene.list$data[gene.list$pathway != "Hypermutator"] <- gl
  gene.levels <- gene_levels(gene.list)
  id.levels <- endo_id_levels(idat, sample.order)
  plevels <- pathway.levels
  gene.levels2 <- gene.levels[gene.levels != "hypermutator"]
  idat2 <- idat %>%
    filter(gene != "hypermutator") %>%
    mutate(
      lab_id = factor(lab_id, id.levels),
      gene = factor(gene, gene.levels2),
      pathway = factor(pathway, plevels),
      tumor_type = factor(tumor_type, tumor.levels)
    )
  idat2
}

order_idat_mucinous <- function(idat, tumor.levels, pathway.levels, sample.order) {
  gene.list <- gene_list(idat, pathway.levels, tumor.levels)
  gl <- remove_duplicated_genes(gene.list)
  gene.list$data[gene.list$pathway != "Hypermutator"] <- gl
  gene.levels <- gene_levels(gene.list)
  id.levels <- muc_id_levels(idat, sample.order)
  plevels <- pathway.levels
  gene.levels2 <- gene.levels[gene.levels != "hypermutator"]
  idat2 <- idat %>%
    filter(gene != "hypermutator") %>%
    mutate(
      lab_id = factor(lab_id, id.levels),
      gene = factor(gene, gene.levels2),
      pathway = factor(pathway, plevels),
      tumor_type = factor(tumor_type, tumor.levels)
    )
  idat2
}

order_endo <- function(idat, manifest) {
  tumor.levels <- c("Ovarian endometrioid", "Uterine endometrioid")
  pathway.levels <- pathway_levels()
  pathway.levels2 <- c("Hypermutator", pathway.levels)
  sample.order <- endo_order(idat)
  idat2 <- order_idat_endo(idat, tumor.levels, pathway.levels2, sample.order)
  idat3 <- remove_duplicate_samples(idat2, manifest)
  idat3$pathway <- droplevels(idat3$pathway)
  idat3
}

exclude_cgcrc254_only_genes <- function(idat) {
  genes.to.drop <- idat %>%
    group_by(gene) %>%
    summarize(ids = paste(unique(lab_id), collapse = ",")) %>%
    filter(ids == "CGCRC254T")
  idat2 <- filter(idat, !gene %in% genes.to.drop$gene)
  idat2
}

order_mucinous <- function(idat, manifest, pathway.levels) {
  tumor.levels <- c("Ovarian mucinous", "Colorectal mucinous")
  sample.order <- muc_order(idat)
  pathway.levels2 <- c("Hypermutator", pathway.levels)
  idat2 <- order_idat_mucinous(idat, tumor.levels, pathway.levels2, sample.order)
  idat3 <- remove_duplicate_samples(idat2, manifest)
  idat4 <- exclude_cgcrc254_only_genes(idat3)
  idat4$pathway <- droplevels(idat4$pathway)
  idat4$gene <- droplevels(idat4$gene)
  idat4
}

order_gi <- function(idat, manifest) {
  tumor.levels <- c("Colorectal mucinous", "Stomach mucinous", "Pancreas mucinous")
  idat2 <- order_idat(idat, tumor.levels)
  idat2
}

remove_duplicate_samples <- function(idat2, manifest) {
  dup.samples <- idat2 %>%
    dplyr::select(lab_id) %>%
    distinct() %>%
    left_join(dplyr::select(manifest, subject_id, lab_id), by = "lab_id") %>%
    group_by(subject_id) %>%
    nest()
  nr <- map_dbl(dup.samples$data, nrow)
  dup.samples2 <- dup.samples[nr > 1, ]
  drop.samples <- dup.samples$data %>% map_dfr(function(x) x[-1, ])
  if (any("CGOV141T_1" %in% drop.samples$lab_id)) {
    drop.samples$lab_id[match("CGOV141T_1", drop.samples$lab_id)] <- "CGOV141T"
  }
  idat2 <- filter(idat2, !lab_id %in% drop.samples$lab_id)
  id.levels <- levels(idat2$lab_id)
  id.levels <- id.levels[!id.levels %in% drop.samples$lab_id]
  idat2$lab_id <- factor(idat2$lab_id, id.levels)
  idat2
}

## ── 05-data_integration functions ────────────────────────────────────────────

#' Expand multi-gene entries to one row per gene
#'
#' @param x Tibble with columns \code{lab_id} and \code{gene_symbol}.
#' @param sep Separator used in multi-gene \code{gene_symbol} strings.
#' @export
expand_genes <- function(x, sep = ", ") {
  genes  <- strsplit(x$gene_symbol, sep)
  ngenes <- purrr::map_int(genes, length)
  tibble::tibble(lab_id      = rep(x$lab_id, ngenes),
                 gene_symbol = unlist(genes))
}

#' Read and combine WES and WGS mutation CSVs into a single table
#' Read the consolidated mutations file (extdata/mutations.tsv)
#'
#' @param path Path to mutations.tsv (the canonical consolidated mutations file).
#' @export
read_mutations <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

#' Compute mutation spectra summary from a combined mutation table
#'
#' Hypermutators are not filtered here — they appear as samples with no
#' substitution rows in the input and thus produce empty columns in the
#' spectra panel, matching the original published CRC analysis.
#'
#' @param mt  Mutation table from \code{\link{read_mutations}}.
#' @export
build_mutation_spectra <- function(mt) {
  ## mutation format: chrN_start-end_REF_ALT; type column is lowercase
  subs1 <- dplyr::filter(mt, tolower(type) == "substitution",
                         grepl("[ATCG]_[ATCG]", mutation))
  subs2 <- dplyr::filter(mt, tolower(type) == "substitution",
                         grepl("[ATCG]/[ATCG]", mutation)) %>%
    dplyr::mutate(mutation = stringr::str_replace_all(mutation, "/", "_"))
  subs <- dplyr::bind_rows(subs1, subs2)
  complement <- c("A", "T", "C", "G") %>% stats::setNames(c("T", "A", "G", "C"))
  ref_alt  <- stringr::str_match(subs$mutation, "([ATCG])_([ATCG])$")
  ref_base <- ref_alt[, 2]
  alt_base <- ref_alt[, 3]
  mutspectra <- subs %>%
    dplyr::mutate(
      base1 = ifelse(ref_base %in% c("C", "T"), ref_base, complement[ref_base]),
      base2 = ifelse(ref_base %in% c("C", "T"), alt_base, complement[alt_base]),
      mutation = paste0(base1, ">", base2)
    ) %>%
    dplyr::select(lab_id, gene, mutation)
  nsubs <- mutspectra %>%
    dplyr::group_by(lab_id) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop")
  mutspectra %>%
    dplyr::left_join(nsubs, by = "lab_id") %>%
    dplyr::group_by(lab_id, mutation) %>%
    dplyr::summarize(number_substitutions = unique(n), n = unique(n),
                     n.type = dplyr::n(), percent = n.type / n,
                     .groups = "drop")
}

#' Read WES copy-number table (table_S5.tsv) and expand to one row per gene
#'
#' @param s5_path Path to table_S5.tsv.
#' @export
read_wes_cnv <- function(s5_path) {
  s5 <- readr::read_tsv(s5_path, show_col_types = FALSE) %>%
    clean_colnames3()
  s52 <- s5 %>%
    dplyr::filter(is_focal_cnv, cnv_type %in% c("LOSS", "GAIN") | is_loh) %>%
    dplyr::select(lab_id, ploidy, cnv_type, is_loh, total_copy_number,
                  clinically_relevant_genes, biologically_relevant_genes,
                  chromosome, start, end) %>%
    dplyr::rename(biol_gene = biologically_relevant_genes) %>%
    dplyr::mutate(
      biol_gene = ifelse(biol_gene == "-", NA_character_, biol_gene),
      biol_gene = ifelse(nchar(biol_gene) == 0L, NA_character_, biol_gene),
      type  = ifelse(cnv_type == "GAIN", "amplification", "deletion"),
      start = as.integer(start)
    ) %>%
    dplyr::rename(chrom = chromosome)
  genelist <- purrr::map(s52$biol_gene, ~ strsplit(.x, ",")[[1]])
  ngenes   <- purrr::map_int(genelist, length)
  tibble::tibble(lab_id      = rep(s52$lab_id, ngenes),
                 gene_symbol = unlist(genelist),
                 chrom       = rep(s52$chrom, ngenes),
                 start       = rep(s52$start, ngenes),
                 type        = rep(s52$type,  ngenes))
}

#' Read WGS deletion table (table_s8.rmd/table_S8.csv)
#'
#' Returns the raw cleaned-column deletion table; rows are one entry per sample.
#' @param del_path Path to table_S8.csv.
#' @export
read_wgs_deletions <- function(del_path) {
  readr::read_csv(del_path, show_col_types = FALSE) %>%
    clean_colnames3()
}

#' Read and process WGS amplicon table (table_s6.rmd/table_s6.tsv)
#'
#' @param amp_path Path to table_s6.tsv.
#' @export
read_wgs_amplicons <- function(amp_path) {
  readr::read_tsv(amp_path, show_col_types = FALSE) %>%
    clean_colnames3() %>%
    dplyr::mutate(gene_symbol = cancer_connection) %>%
    dplyr::filter(!is.na(gene_symbol)) %>%
    expand_genes() %>%
    dplyr::mutate(type = "copynumber", alteration = "amplification") %>%
    dplyr::distinct() %>%
    dplyr::rename(gene = gene_symbol)
}

#' Compute marginal alteration frequencies per sample
#'
#' @param mt       Mutation table from \code{\link{read_mutation_table}}.
#' @param cnv_wes  WES CNV table from \code{\link{read_wes_cnv}}.
#' @param del      Deletion table from \code{\link{read_wgs_deletions}}.
#' @param amp      Amplicon table from \code{\link{read_wgs_amplicons}}.
#' @export
build_marginal_frequencies <- function(mt, cnv_wes, del, amp) {
  mut_marginal <- mt %>%
    dplyr::group_by(lab_id) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(alteration = "mutation", type = "Mutation")
  wes_marginal <- cnv_wes %>%
    dplyr::rename(alteration = type) %>%
    dplyr::mutate(type = "Copy number") %>%
    dplyr::select(lab_id, gene_symbol, type, alteration) %>%
    dplyr::group_by(lab_id, alteration) %>%
    dplyr::summarize(n = dplyr::n(), alteration = unique(alteration),
                     type = unique(type), .groups = "drop")
  del_marginal <- del %>%
    dplyr::select(lab_id) %>%
    dplyr::group_by(lab_id) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(alteration = "deletion", type = "Copy number")
  amp_marginal <- dplyr::select(amp, lab_id) %>%
    dplyr::group_by(lab_id) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(alteration = "amplification", type = "Copy number")
  dplyr::bind_rows(mut_marginal, wes_marginal,
                   dplyr::bind_rows(del_marginal, amp_marginal))
}

## ── mutational-signature computation (BASE-02) ──────────────────────────────
##
## Ports `mut.to.sigs.input()` + `whichSignatures()` (package `deconstructSigs`,
## COSMIC v2 legacy 30-signature reference) forward from the archived 2019
## script into named functions. See `provenance/signature_matrices_provenance.md`
## for the full derivation this is based on, and `cards/BASE-02-port-mutational-
## signatures.md`'s Notes/log for the schema investigation and the environment
## blocker that has prevented end-to-end verification of this port so far.
##
## `extdata/mutations.tsv`'s coordinates are hg18 throughout (see
## `code/mutations.rmd`'s column-description comment) -- no liftover step is
## needed, unlike the archived script's mixed-build source table. Use
## `BSgenome.Hsapiens.UCSC.hg18`, never a different build: trinucleotide
## context is derived from the genome sequence at each call's coordinates, so
## the wrong build silently produces wrong signature weights.

#' Parse chr/pos/ref/alt from a mutations.tsv coordinate string
#'
#' `mutations.tsv`'s `mutation` column encodes each call in one of two
#' formats (both documented in `code/mutations.rmd`'s column-description
#' comment): \code{chrN_start-end_REF_ALT} (the majority; \code{start == end}
#' for substitutions) or \code{chrN:pos_REF/ALT} (a minority of WGS Strelka
#' calls). A handful of rows (all from CGCRC254T) match neither and are
#' returned as \code{NA} rows for the caller to drop.
#'
#' @param mutation Character vector, the `mutation` column of
#'   \code{\link{read_mutations}}'s output.
#' @export
parse_mutation_coords <- function(mutation) {
  fmt_a <- stringr::str_match(
    mutation, "^chr([0-9XYM]+)_([0-9]+)-[0-9]+_([ACGT])_([ACGT])$"
  )
  fmt_b <- stringr::str_match(
    mutation, "^chr([0-9XYM]+):([0-9]+)_([ACGT])/([ACGT])$"
  )
  chr_num <- dplyr::coalesce(fmt_a[, 2], fmt_b[, 2])
  tibble::tibble(
    chr = ifelse(is.na(chr_num), NA_character_, paste0("chr", chr_num)),
    pos = as.integer(dplyr::coalesce(fmt_a[, 3], fmt_b[, 3])),
    ref = dplyr::coalesce(fmt_a[, 4], fmt_b[, 4]),
    alt = dplyr::coalesce(fmt_a[, 5], fmt_b[, 5])
  )
}

#' Build a deconstructSigs-ready mutation table for one tumor arm
#'
#' Filters \code{mutations} to single-base substitutions (deconstructSigs'
#' trinucleotide-context method is defined for SNVs only; Deletion/Insertion/
#' Indel/Complex-indel rows are excluded here, matching the archived script's
#' use of the WES/WGS "Mutation.Position" columns which were substitution-only)
#' for samples in \code{lab_ids}, parses genomic coordinates, and drops
#' unparseable rows. Returns the `Sample`/chr/pos/ref/alt columns
#' \code{deconstructSigs::mut.to.sigs.input()} expects.
#'
#' @param mutations Mutation table from \code{\link{read_mutations}}.
#' @param lab_ids   Character vector of lab IDs defining the tumor arm.
#' @export
signature_input_table <- function(mutations, lab_ids) {
  subs <- dplyr::filter(
    mutations, tolower(type) == "substitution", lab_id %in% lab_ids
  )
  coords <- parse_mutation_coords(subs$mutation)
  dplyr::bind_cols(lab_id = subs$lab_id, coords) %>%
    dplyr::filter(!is.na(chr), !is.na(pos), !is.na(ref), !is.na(alt))
}

#' Fit COSMIC v2 mutational signatures per sample
#'
#' Runs \code{deconstructSigs::mut.to.sigs.input()} (hg18 reference) followed
#' by \code{deconstructSigs::whichSignatures(signatures.ref = signatures.cosmic,
#' tri.counts.method = "exome")} for every sample in \code{sig_input}, then
#' assembles the per-sample \code{$weights} into a signatures-by-samples
#' matrix with all-zero signature rows dropped -- matching the shape of the
#' committed \code{extdata/{endosigs,mucsigs}.rds} (signatures.cosmic is
#' bundled reference data shipped with `deconstructSigs`; no external
#' download needed, but the fitted weights depend on the installed package
#' version, hence pinning it in `renv.lock`).
#'
#' @param sig_input Output of \code{\link{signature_input_table}}.
#' @export
fit_mutational_signatures <- function(sig_input) {
  if (!requireNamespace("deconstructSigs", quietly = TRUE)) {
    stop(
      "Package 'deconstructSigs' is required but not installed. See ",
      "provenance/signature_matrices_provenance.md and cards/BASE-02-port-",
      "mutational-signatures.md for why this is a real, required dependency ",
      "(not optional dev tooling) and the environment blocker encountered ",
      "when installing it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg18", quietly = TRUE)) {
    stop(
      "Package 'BSgenome.Hsapiens.UCSC.hg18' is required but not installed. ",
      "extdata/mutations.tsv uses hg18 coordinates throughout (see ",
      "code/mutations.rmd's column-description comment) -- do not substitute ",
      "a different genome build; that would silently produce wrong ",
      "trinucleotide contexts and wrong signature weights.",
      call. = FALSE
    )
  }
  bsg <- getExportedValue(
    "BSgenome.Hsapiens.UCSC.hg18", "BSgenome.Hsapiens.UCSC.hg18"
  )
  sigs_input <- deconstructSigs::mut.to.sigs.input(
    mut.ref   = as.data.frame(sig_input),
    sample.id = "lab_id",
    chr       = "chr",
    pos       = "pos",
    ref       = "ref",
    alt       = "alt",
    bsg       = bsg
  )
  signatures_cosmic <- get(
    "signatures.cosmic", envir = asNamespace("deconstructSigs")
  )
  output <- lapply(rownames(sigs_input), function(id) {
    deconstructSigs::whichSignatures(
      tumor.ref         = sigs_input,
      signatures.ref    = signatures_cosmic,
      sample.id         = id,
      contexts.needed   = TRUE,
      tri.counts.method = "exome"
    )
  })
  names(output) <- rownames(sigs_input)
  weights <- lapply(output, function(x) x$weights)
  sig_matrix <- t(as.matrix(do.call(rbind, weights)))
  sig_matrix[rowSums(sig_matrix) > 0, , drop = FALSE]
}

#' Compute the endometrioid-arm mutational-signature matrix
#'
#' Reproduces \code{extdata/endosigs.rds}: filters \code{mutations} to lab IDs
#' with \code{manifest$tumor_type \%in\% c("ovarian endometrioid", "uterine
#' endometrioid")}, fits COSMIC v2 signatures per sample via
#' \code{\link{fit_mutational_signatures}}. See
#' \code{provenance/signature_matrices_provenance.md} for the full derivation
#' and known caveats (e.g. CGOV163T's very low mutation count).
#'
#' @param mutations Mutation table from \code{\link{read_mutations}}.
#' @param manifest  Manifest tibble (\code{data/manifest.rda}) with
#'   \code{lab_id}/\code{tumor_type} columns.
#' @export
endo_signature_matrix <- function(mutations, manifest) {
  lab_ids <- manifest$lab_id[
    manifest$tumor_type %in% c("ovarian endometrioid", "uterine endometrioid")
  ]
  fit_mutational_signatures(signature_input_table(mutations, lab_ids))
}

#' Compute the mucinous-arm mutational-signature matrix
#'
#' Reproduces \code{extdata/mucsigs.rds}: filters \code{mutations} to lab IDs
#' with \code{manifest$tumor_type \%in\% c("ovarian mucinous", "colorectal",
#' "pancreas", "stomach")} -- the four tumor types combined in the committed
#' matrix's columns (CGOV/CGCRC/CGPA/CGST prefixes) -- and fits COSMIC v2
#' signatures per sample via \code{\link{fit_mutational_signatures}}.
#'
#' Note: this uses the raw \code{manifest$tumor_type} values directly rather
#' than the existing \code{\link{muc}} helper, because \code{muc()} contains
#' a pre-existing typo ("Pancreaas mucinous") and expects the
#' \code{\link{cancer_names}}-capitalized vocabulary, neither of which matches
#' \code{manifest$tumor_type}'s raw lower-case values. See this card's
#' Notes/log for the full investigation.
#'
#' @param mutations Mutation table from \code{\link{read_mutations}}.
#' @param manifest  Manifest tibble (\code{data/manifest.rda}) with
#'   \code{lab_id}/\code{tumor_type} columns.
#' @export
muc_signature_matrix <- function(mutations, manifest) {
  lab_ids <- manifest$lab_id[
    manifest$tumor_type %in% c("ovarian mucinous", "colorectal", "pancreas", "stomach")
  ]
  fit_mutational_signatures(signature_input_table(mutations, lab_ids))
}

## ── Raw mutation-call consolidation (MAN-01) ────────────────────────────────
##
## Ports the former `code/mutations.rmd`'s raw-source consolidation forward
## into named, composable functions. `code/mutations.rmd` now just calls
## these and writes `output/mutations.rds` / `extdata/mutations.tsv` for
## anyone who wants the files on disk; the manuscript pipeline's
## `mutations_tbl` target computes the same table in memory via
## `build_mutations_tbl()` below, from the raw source files (PGDx/Strelka
## Excel reports, per-sample Strelka rerun TSVs, and the private manifest),
## rather than reading the frozen `extdata/mutations.tsv`.

#' WGS tumor samples excluded from the Strelka rerun calls
#'
#' These 12 lab IDs are present in both \code{extdata/strelka_reruns/} and
#' the current manifest, but were absent from the manifest when the original
#' published analysis (the former \code{02-04-mutations.rmd}) was run.
#' Including them would add ~1,548 rows to the mutations table and change
#' the manuscript's published mutation counts. This is a deliberate exclusion
#' to preserve reproducibility of the published analysis -- **not a bug** --
#' and a candidate for inclusion in a future manuscript revision.
#'
#' @export
EXCLUDED_STRELKA_RERUN_SAMPLES <- c(
  "CGOV463T", "CGOV467T", "CGOV469T", "CGOV470T", "CGOV471T", "CGOV474T",
  "CGOV477T", "CGOV478T", "CGOV480T", "CGOV484T", "CGOV485T", "CGOV487T"
)

#' Patch the CGPA367T rows in the PGDx Excel report
#'
#' \code{pgdx-compiled.xlsx} has shifted Type/Consequence/Context/MAF columns
#' for the report rows belonging to lab ID CGPA367T. Correct values were
#' verified against the archived \code{pgdx-compiled-patched.xlsx} (no longer
#' on disk). The affected rows are identified via a \code{manifest$lab_id}
#' lookup rather than a hardcoded vendor ID, so no vendor IDs appear in this
#' file (card \code{REL-01}).
#'
#' @param pgdx PGDx report tibble (must have \code{Prefix} and \code{Gene}
#'   columns).
#' @param manifest Private manifest tibble with \code{pgdx_id}/\code{lab_id}
#'   columns (\code{ovarian.subtypes/inst/extdata/manifest.rds}).
#' @export
patch_cgpa367t <- function(pgdx, manifest) {
  bad_id <- manifest$pgdx_id[manifest$lab_id == "CGPA367T"]
  stopifnot(length(bad_id) == 1L, !is.na(bad_id))
  patch <- tibble::tibble(
    Gene        = c("FRMD4B", "GGN",    "GNAQ",   "GRIK1"),
    Type        = "Substitution",
    Consequence = c("Nonsynonymous coding", "Nonsynonymous coding",
                    "Nonsynonymous coding", "Nonsense"),
    Context     = c("AAGTCNTTGCT", "GGGCCNGGTCG", "GAGTGNGTCCA", "AGGTTNATGTG"),
    MAF         = c(0.277778, 0.138889, 0.192308, 0.153846)
  )
  bad_rows <- grepl(bad_id, pgdx$Prefix, fixed = TRUE)
  stopifnot(sum(bad_rows) == nrow(patch))
  pgdx[bad_rows, c("Type", "Consequence", "Context", "MAF")] <-
    patch[match(pgdx$Gene[bad_rows], patch$Gene),
          c("Type", "Consequence", "Context", "MAF")]
  pgdx
}

#' Read and combine the PGDx and original-Strelka mutation-call Excel reports
#'
#' @param report_dir Directory containing \code{pgdx-compiled.xlsx} and
#'   \code{strelka-compiled.xlsx} (\code{extdata/mutation_reports}).
#' @param manifest Private manifest tibble (needed for the CGPA367T patch,
#'   see \code{\link{patch_cgpa367t}}).
#' @export
read_mutation_reports <- function(report_dir, manifest) {
  pgdx <- readxl::read_excel(file.path(report_dir, "pgdx-compiled.xlsx"), sheet = 1) %>%
    dplyr::select(-c("Original Nuc Change", "Original Coordinates")) %>%
    dplyr::mutate(caller = "PGDx")
  pgdx <- patch_cgpa367t(pgdx, manifest)

  strelka <- readxl::read_excel(file.path(report_dir, "strelka-compiled.xlsx"), sheet = 1) %>%
    dplyr::select(c("Prefix", "PGDx-like Change ID",
                    "Type", "Annotation", "Context",
                    "Gene Name",
                    "Feature ID (CCDS)", "HGVS.c",
                    "Tumor Distinct Read Depth (tier1)",
                    "Normal Distinct Read Depth (tier1)")) %>%
    dplyr::mutate(
      MAF = `Tumor Distinct Read Depth (tier1)` / `Normal Distinct Read Depth (tier1)`,
      MAF = round(MAF, 3)
    ) %>%
    dplyr::select(-c(`Tumor Distinct Read Depth (tier1)`,
                     `Normal Distinct Read Depth (tier1)`)) %>%
    dplyr::rename(
      `hg18 Nuc Change` = `PGDx-like Change ID`,
      CCDS              = `Feature ID (CCDS)`,
      `AA Change`       = "HGVS.c",
      Gene              = `Gene Name`,
      Consequence       = Annotation
    ) %>%
    dplyr::mutate(caller = "Strelka") %>%
    dplyr::select(colnames(pgdx))

  dplyr::bind_rows(pgdx, strelka)
}

#' Root a PGDx report \code{Prefix} string down to its manifest-matchable form
#' @keywords internal
mutation_report_root_id <- function(x) {
  tmp <- strsplit(x, "_Cancer")
  nm  <- sapply(tmp, "[", 1)
  stringr::str_replace_all(nm, "WGS_Ex", "WGS")
}

#' Match PGDx-prefixed report rows to the manifest's \code{pgdx_id}
#'
#' Splits \code{combined} into PGDx-prefixed (\code{Prefix} starting with
#' \code{P}/\code{L}/\code{V}) and non-PGDx-prefixed rows. For the
#' PGDx-prefixed rows, resolves each \code{Prefix}'s manifest \code{pgdx_id}
#' via a chain of successive heuristics -- the report \code{Prefix} strings
#' do not match the manifest's \code{pgdx_id} verbatim.
#'
#' @param combined Combined PGDx + Strelka report tibble from
#'   \code{\link{read_mutation_reports}}.
#' @param manifest Private manifest tibble.
#' @return List with \code{pgdx} (matched, \code{is_pgdx = TRUE}, has
#'   \code{pgdx_id}) and \code{notpgdx} (\code{is_pgdx = FALSE}) tibbles.
#' @export
match_pgdx_ids <- function(combined, manifest) {
  combined.pgdx <- dplyr::filter(combined, grepl("^[PLV]", Prefix)) %>%
    dplyr::mutate(root = mutation_report_root_id(Prefix),
                  root = stringr::str_replace_all(root, "__", ""))

  x          <- dplyr::select(combined.pgdx, root, Prefix) %>% dplyr::distinct()
  notmatched <- dplyr::filter(x, !root %in% manifest$pgdx_id)
  matched    <- dplyr::filter(x,  root %in% manifest$pgdx_id)

  newid <- rep(NA, nrow(notmatched))
  for (i in seq_len(nrow(notmatched))) {
    tmp <- strsplit(notmatched$root[i], "_")[[1]][1]
    ix  <- grep(tmp, manifest$pgdx_id)
    if (length(ix) > 1) {
      ix <- ix[grep("[0-9]T", manifest$pgdx_id[ix])]
      if (length(ix) > 1) stop("Still more than one hit")
    }
    if (length(ix) == 0) {
      if (grepl("^LP", tmp)) {
        lp.ids <- manifest$pgdx_id[grep("^LP", manifest$pgdx_id)]
        tmp    <- paste0(notmatched$root[i], "_WGS")
        ix     <- match(tmp, lp.ids)
        if (length(ix) == 0) ix <- match(tmp, manifest$pgdx_id)
        if (is.na(ix)) {
          tmp <- stringr::str_replace(tmp, "LP600", "LP00")
          ix  <- match(tmp, manifest$pgdx_id)
        }
      }
    }
    if (length(ix) == 0) {
      tmp <- stringr::str_replace(notmatched$root[i], "2WGS", "2_WGS")
      ix  <- match(tmp, manifest$pgdx_id)
    }
    if (is.na(manifest$pgdx_id[ix])) {
      tmp <- stringr::str_replace(tmp, "^Victor_", "")
      tmp <- stringr::str_replace(tmp, "_hg19", "")
      ix  <- match(tmp, manifest$pgdx_id)
    }
    newid[i] <- manifest$pgdx_id[ix]
  }

  matched$root2    <- matched$root
  notmatched$root2 <- newid
  x2 <- dplyr::bind_rows(matched, notmatched) %>%
    dplyr::select(-root) %>%
    dplyr::rename(pgdx_id = root2)

  combined.pgdx2   <- dplyr::left_join(combined.pgdx, x2, by = "Prefix") %>%
    dplyr::mutate(is_pgdx = TRUE)
  combined.notpgdx <- dplyr::filter(combined, !grepl("^[PLV]", Prefix)) %>%
    dplyr::mutate(is_pgdx = FALSE)

  list(pgdx = combined.pgdx2, notpgdx = combined.notpgdx)
}

#' Tumor-sample BAM basename lookup table derived from the manifest
#' @keywords internal
tumor_bam_lookup <- function(manifest) {
  dplyr::select(manifest, subject_id, lab_id, bam_local, tumor.normal, platform) %>%
    dplyr::filter(tumor.normal == "tumor") %>%
    dplyr::mutate(bam = basename(bam_local)) %>%
    dplyr::select(-bam_local)
}

#' Join mutation-report rows to the manifest for subject/lab IDs
#'
#' Joins the PGDx-matched rows (via \code{pgdx_id}) and non-PGDx rows (via
#' \code{Prefix == lab_id}) to the manifest, then, for any report rows that
#' still fail to match (some non-PGDx rows carry a BAM-filename-derived
#' \code{Prefix} rather than a lab ID -- two of these also carry a trailing
#' run suffix absent from the BAM filename), resolves the match via the BAM
#' basename instead.
#'
#' @param matched_ids Output of \code{\link{match_pgdx_ids}} (list with
#'   \code{pgdx}/\code{notpgdx} elements).
#' @param manifest Private manifest tibble.
#' @export
join_mutation_manifest <- function(matched_ids, manifest) {
  manifest2 <- dplyr::select(manifest, subject_id, lab_id, pgdx_id,
                             tumor.normal, platform) %>%
    dplyr::filter(tumor.normal == "tumor")

  m <- tumor_bam_lookup(manifest)

  merged1       <- dplyr::left_join(matched_ids$pgdx,    manifest2, by = "pgdx_id")
  merged2       <- dplyr::left_join(matched_ids$notpgdx, manifest2, by = c("Prefix" = "lab_id"))
  merge.attempt <- dplyr::bind_rows(merged1, merged2)

  merged     <- dplyr::filter(merge.attempt, !is.na(subject_id))
  not.merged <- dplyr::filter(merge.attempt,  is.na(subject_id))

  todo <- dplyr::select(not.merged, Prefix, root) %>%
    dplyr::distinct() %>%
    dplyr::mutate(bam = NA)

  for (i in seq_len(nrow(todo))) {
    if (!is.na(todo$root[i])) {
      ## Two report prefixes carry a trailing run suffix ("_" or "_A") that
      ## is absent from the BAM filename; strip it so the grep below
      ## resolves. Equivalent to the former hardcoded lookup: of the 14
      ## roots reaching this loop, only those two end in "_" or "_A".
      todo$root[i] <- stringr::str_replace(todo$root[i], "_A?$", "")
      ix <- grep(todo$root[i], m$bam)
      if (length(ix) != 1) stop()
      todo$bam[i] <- m$bam[ix]
    }
  }
  todo <- dplyr::left_join(todo, m, by = "bam")

  not.merged2 <- dplyr::left_join(
    dplyr::select(not.merged, -c(subject_id, lab_id, tumor.normal, platform)),
    todo, by = c("Prefix", "root")
  ) %>%
    dplyr::select(-bam) %>%
    dplyr::select(colnames(merged))

  dplyr::bind_rows(merged, not.merged2) %>% dplyr::arrange(platform)
}

#' Fill in \code{lab_id} for report rows the manifest join left missing
#'
#' A handful of rows join to a subject/platform via BAM matching but still
#' have no \code{lab_id} (the manifest's \code{lab_id} column is \code{NA}
#' for that BAM); for those, the report's own \code{Prefix} value serves as
#' \code{lab_id} directly.
#'
#' @param merged_original Output of \code{\link{join_mutation_manifest}}.
#' @param manifest Private manifest tibble.
#' @export
fix_missing_lab_id <- function(merged_original, manifest) {
  missing.lab <- dplyr::filter(merged_original,  is.na(merged_original$lab_id))
  notmissing  <- dplyr::filter(merged_original, !is.na(merged_original$lab_id))

  m <- tumor_bam_lookup(manifest)

  missing.lab2 <- dplyr::select(missing.lab, -c(subject_id, lab_id, tumor.normal, platform)) %>%
    dplyr::mutate(lab_id = Prefix)
  missing.lab3 <- dplyr::left_join(missing.lab2, m, by = "lab_id") %>%
    dplyr::select(colnames(notmissing))

  dplyr::bind_rows(notmissing, missing.lab3)
}

#' Read per-sample Strelka rerun TSVs
#'
#' Files are individual per-sample Strelka rerun reports, rsynced from
#' \code{/dcs05/scharpf/data/skoul/Projects/ovarian_subtypes_strelka/strelka-pipeline/outDir/Reports/filtered/}
#' into \code{extdata/strelka_reruns/} (one
#' \code{<lab_id>.small_variants_plv0.1.12_filtered.txt} file per lab ID).
#'
#' @param rerun_dir Directory of per-sample Strelka rerun TSVs
#'   (\code{extdata/strelka_reruns}).
#' @param manifest Private manifest tibble (restricts to lab IDs present in
#'   the manifest).
#' @export
read_strelka_reruns <- function(rerun_dir, manifest) {
  reports <- list.files(rerun_dir, full.names = TRUE)

  read_one_report <- function(fname) {
    tmp        <- readr::read_tsv(fname, show_col_types = FALSE)
    tmp$lab_id <- strsplit(basename(fname), "\\.")[[1]][1]
    tmp
  }

  purrr::map_dfr(reports, read_one_report) %>%
    dplyr::select(c(lab_id,
                    "Variant (hg18)", "Type", "Annotation", "Context",
                    "Gene Name", "Feature ID (CCDS)", "HGVS.c",
                    "Tumor Distinct Read Depth (tier1)",
                    "Normal Distinct Read Depth (tier1)")) %>%
    dplyr::mutate(
      MAF = `Tumor Distinct Read Depth (tier1)` / `Normal Distinct Read Depth (tier1)`,
      MAF = round(MAF, 3)
    ) %>%
    dplyr::select(-c(`Tumor Distinct Read Depth (tier1)`,
                     `Normal Distinct Read Depth (tier1)`)) %>%
    dplyr::rename(
      `hg18 Nuc Change` = `Variant (hg18)`,
      CCDS              = `Feature ID (CCDS)`,
      `AA Change`       = "HGVS.c",
      Gene              = `Gene Name`,
      Consequence       = Annotation
    ) %>%
    dplyr::mutate(caller = "Strelka", platform = "WGS") %>%
    dplyr::filter(lab_id %in% manifest$lab_id)
}

#' Replace original Strelka calls with rerun calls for corrected samples
#'
#' Excludes \code{\link{EXCLUDED_STRELKA_RERUN_SAMPLES}} from the rerun calls
#' (12 WGS tumor samples absent from the manifest when the original
#' published analysis was run; see that constant's documentation), then, for
#' every remaining rerun sample, drops its original calls and substitutes the
#' rerun's.
#'
#' @param mut_original Output of \code{\link{fix_missing_lab_id}}.
#' @param mut_rerun Output of \code{\link{read_strelka_reruns}}.
#' @export
integrate_mutation_calls <- function(mut_original, mut_rerun) {
  mut_rerun <- dplyr::filter(mut_rerun, !lab_id %in% EXCLUDED_STRELKA_RERUN_SAMPLES)

  dplyr::select(mut_original, colnames(mut_rerun)) %>%
    dplyr::filter(!lab_id %in% EXCLUDED_STRELKA_RERUN_SAMPLES) %>%
    dplyr::filter(!lab_id %in% mut_rerun$lab_id) %>%
    dplyr::bind_rows(mut_rerun)
}

#' Consolidate PGDx + Strelka (original & rerun) mutation calls
#'
#' Top-level orchestrator that reproduces the former \code{code/mutations.rmd}'s
#' \code{mutations} object (the table saved to \code{output/mutations.rds})
#' from raw sources: the PGDx/Strelka Excel reports, the per-sample Strelka
#' rerun TSVs, and the private manifest. Composes, in order,
#' \code{\link{read_mutation_reports}}, \code{\link{match_pgdx_ids}},
#' \code{\link{join_mutation_manifest}}, \code{\link{fix_missing_lab_id}},
#' \code{\link{read_strelka_reruns}}, and \code{\link{integrate_mutation_calls}}
#' (which applies \code{\link{EXCLUDED_STRELKA_RERUN_SAMPLES}}).
#'
#' @param report_dir Directory with \code{pgdx-compiled.xlsx}/
#'   \code{strelka-compiled.xlsx} (\code{extdata/mutation_reports}).
#' @param rerun_dir Directory of per-sample Strelka rerun TSVs
#'   (\code{extdata/strelka_reruns}).
#' @param manifest Private manifest tibble
#'   (\code{ovarian.subtypes/inst/extdata/manifest.rds}, loaded directly --
#'   **not** \code{\link{get_manifest}}'s PHI-stripped \code{data(manifest)},
#'   which lacks the \code{pgdx_id}/\code{bam_local} columns this port needs).
#' @export
consolidate_mutation_calls <- function(report_dir, rerun_dir, manifest) {
  combined     <- read_mutation_reports(report_dir, manifest)
  matched_ids  <- match_pgdx_ids(combined, manifest)
  merged_orig  <- join_mutation_manifest(matched_ids, manifest)
  mut_original <- fix_missing_lab_id(merged_orig, manifest)
  mut_rerun    <- read_strelka_reruns(rerun_dir, manifest)
  integrate_mutation_calls(mut_original, mut_rerun)
}

#' Format the consolidated mutations table for the canonical TSV export
#'
#' Drops rows with no platform (unmatched/ambiguous report rows) and renames
#' to the canonical \code{extdata/mutations.tsv} column names. See the
#' provenance comment in \code{code/mutations.rmd} for full column
#' documentation.
#'
#' \code{write_tsv()}/\code{read_tsv()}'s round trip -- the path
#' \code{\link{read_mutations}} has always used to serve this table -- treats
#' readr's default missing-value markers (\code{""} and the literal string
#' \code{"NA"}) as true \code{NA} on read-back. \code{aa_change} (sourced from
#' the report's \code{HGVS.c} column) is the one column where the raw report
#' data currently contains that literal string rather than a genuine missing
#' value, so it is normalized here to match \code{read_mutations()}'s existing,
#' already-in-production behavior -- this in-memory path must return the exact
#' same table \code{read_mutations(extdata/mutations.tsv)} has always returned.
#'
#' @param mutations Output of \code{\link{consolidate_mutation_calls}}.
#' @export
format_mutations_export <- function(mutations) {
  dplyr::filter(mutations, !is.na(platform)) %>%
    dplyr::select(
      lab_id,
      mutation    = `hg18 Nuc Change`,
      gene        = Gene,
      ccds        = CCDS,
      aa_change   = `AA Change`,
      type        = Type,
      consequence = Consequence,
      context     = Context,
      maf         = MAF,
      caller,
      platform
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(.x, "NA")))
}

#' Consolidate and format mutation calls into the canonical export shape
#'
#' Composes \code{\link{consolidate_mutation_calls}} with
#' \code{\link{format_mutations_export}} to produce the table in the exact
#' column layout of \code{extdata/mutations.tsv} -- the shape
#' \code{\link{read_mutations}} returns when reading that file from disk.
#' Used by the manuscript pipeline's \code{mutations_tbl} target to compute
#' the table in memory rather than reading the frozen TSV file.
#'
#' @inheritParams consolidate_mutation_calls
#' @export
build_mutations_tbl <- function(report_dir, rerun_dir, manifest) {
  format_mutations_export(consolidate_mutation_calls(report_dir, rerun_dir, manifest))
}
