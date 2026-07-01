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
    select(gene_symbol, internal_id, tumor_type, is_altered) %>%
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
    select(gene_symbol, internal_id, is_altered) %>%
    spread(gene_symbol, is_altered) %>%
    left_join(tumortype, by = "internal_id")
  X <- select(ungroup(gdat3), -internal_id) %>%
    select(-tumor_type) %>%
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
    select(tumor_type, mt, n)
  m2[is.na(m2)] <- 0
  m <- select(m2, -tumor_type) %>%
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
    select(
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
  x3 <- select(x2, tumor_type, mt, n)
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
  y <- select(y, lab_id, lab_id2)
  isf <- is.factor(x$lab_id)
  if (isf) {
    levs <- tibble(lab_id = levels(x$lab_id)) %>%
      inner_join(y, by = "lab_id")
    levs2 <- levs$lab_id2
  }
  x.y <- inner_join(x, y, by = "lab_id") %>%
    mutate(lab_id = ifelse(lab_id == lab_id2, lab_id, lab_id2)) %>%
    select(-lab_id2)
  if (isf) {
    x.y$lab_id <- factor(x.y$lab_id, levs2)
  }
  return(x.y)
}

order_samples <- function(x, gene.levels) {
  x2 <- select(x, lab_id, gene_symbol) %>%
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
    select(-tumor_type) %>%
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
  tumortypes <- select(manifest, lab_id, tumor_type) %>%
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
  tumortypes <- select(manifest, lab_id, tumor_type) %>%
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
    select(-tumor_type) %>%
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
    select(lab_id) %>%
    distinct() %>%
    left_join(select(manifest, subject_id, lab_id), by = "lab_id") %>%
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
