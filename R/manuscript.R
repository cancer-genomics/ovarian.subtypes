use_listing <- function(x, sid, lid) {
  ix <- match(sid, x$subject_id)
  if (is.na(ix)) stop("Failed to match subject id")
  dat <- x$data[[ix]]
  index <- match(lid, dat$lab_id)
  if (is.na(index)) stop("Failed to match lab id")
  dat <- dat[index, ]
  x$data[[ix]] <- dat
  x
}

keep_first_if_multiple <- function(tumors) {
  tmp <- tumors %>%
    group_by(subject_id) %>%
    nest()
  tmp2 <- use_listing(tmp, "CGOV141", "CGOV141T_1")
  tmp3 <- use_listing(tmp2, "CGOV197", "CGOV197T_1")
  tmp4 <- use_listing(tmp3, "CGOV375", "CGOV375T")
  tmp5 <- use_listing(tmp4, "CGCRC330", "CGCRC330T_1")
  tmp5$data <- tmp5$data %>%
    map(function(x) x[1, ])
  tumors2 <- unnest(tmp5, "data")
  return(tumors2)
}

number_tumor_types <- function(mfest, hmut) {
  w.hmut <- mfest %>%
    group_by(tumor_type) %>%
    summarize(n = length(unique(subject_id)))
  hmut2 <- left_join(hmut, select(mfest, subject_id, lab_id),
    by = "lab_id"
  )
  wo.hmut <- mfest %>%
    filter(!subject_id %in% hmut2$subject_id) %>%
    group_by(tumor_type) %>%
    summarize(n = length(unique(subject_id)))
  num <- left_join(w.hmut, wo.hmut,
    by = "tumor_type",
    suffix = c("", ".no_hypmut")
  )
  num
}

list_cancers <- function() {
  cancers <- list(
    "ovarian endometrioid",
    "ovarian mucinous",
    "uterine endometrioid",
    c("stomach", "colorectal", "pancreas")
  )
  cancers
}

types_as_vector <- function(types) {
  types2 <- setNames(types$n, types$tumor_type)
  types3 <- types2[c(
    "ovarian endometrioid", "ovarian mucinous",
    "uterine endometrioid",
    "colorectal", "pancreas", "stomach"
  )]
  types3
}

num <- function(manifest, tumor.abbrev) {
  tumors <- filter(manifest, tumor.normal == "tumor")
  tmp <- tumors %>%
    group_by(subject_id) %>%
    nest()
  tmp$data <- tmp$data %>%
    map(function(x) x[1, ])
  tumors <- unnest(tmp, "data")
  result <- tumors %>%
    group_by(tumor_type) %>%
    summarize(n = length(unique(subject_id))) %>%
    mutate(abbrev = c(
      "crc",
      "oe",
      "om",
      "pa",
      "st",
      "ue"
    )) %>%
    filter(abbrev == tumor.abbrev) %>%
    pull(n)
  result
}

#' Collapse GI mucinous cancer counts into a single named vector entry
#'
#' Unlike collapse_gi() which operates on a data frame's tumor_type column,
#' this function takes a named numeric vector of tumor counts and sums
#' colorectal + pancreas + stomach into a single "GI" entry.
collapse_gi_counts <- function(types2) {
  gi.cancers <- c("colorectal", "pancreas", "stomach")
  ov.ut.cancers <- c(
    "ovarian endometrioid",
    "ovarian mucinous",
    "uterine endometrioid"
  )
  types4 <- c(types2[ov.ut.cancers], sum(types2[gi.cancers]))
  names(types4) <- c(ov.ut.cancers, "GI")
  n.types2 <- setNames(
    types4,
    c("OE", "OM", "UE", "GI")
  )
  n.types2
}

numbers_in_parenthesis <- function(types2) {
  gi.cancers <- c("colorectal", "pancreas", "stomach")
  ov.ut.cancers <- c(
    "ovarian endometrioid",
    "ovarian mucinous",
    "uterine endometrioid"
  )
  types4 <- c(types2[ov.ut.cancers], sum(types2[gi.cancers]))
  names(types4) <- c(ov.ut.cancers, "GI")
  n.types <- paste0("(n=", types4, ")")
  names(n.types) <- names(types4)
  n.types
}

#' @export
y.of.n <- function(y, manifest, tumor.abbrev) {
  n <- num(manifest, tumor.abbrev)
  paste(y, "of", n)
}

#' @export
pctfun <- function(y, manifest, tumor.abbrev) {
  n <- num(manifest, tumor.abbrev)
  round(y / n * 100, 0)
}

get_tumors <- function(manifest) {
  tumors <- filter(manifest, tumor.normal == "tumor") %>%
    keep_first_if_multiple()
  tumors
}

n_subjects <- function(tumors) {
  N <- tumors %>%
    pull(subject_id) %>%
    unique() %>%
    length()
  N
}

ov_muc_subjects <- function(tumors) {
  ovmuc_ids <- filter(
    tumors,
    tumor_type %in% "ovarian mucinous"
  ) %>%
    pull(subject_id)
  ovmuc_ids
}

endometrial_cancers <- function() {
  c(
    "ovarian endometrioid",
    "uterine endometrioid"
  )
}

mucinous_cancers <- function() {
  c(
    "colorectal",
    "ovarian mucinous",
    "pancreas",
    "stomach"
  )
}

mucinous_labels <- function() c("CRC", "OM", "PM", "SM")

#' @export
num_hypermut <- function(tumor.type, stab3, stab4,
                         manifest,
                         cutoff = 500) {
  tumors <- filter(manifest, tumor.normal == "tumor")
  tmp <- tumors %>%
    group_by(subject_id) %>%
    nest()
  tmp$data <- tmp$data %>%
    map(function(x) x[1, ])
  tumors <- unnest(tmp, "data")
  tumortypes <- tumors %>%
    select(subject_id, lab_id, tumor_type) %>%
    distinct()
  nhypermut.wes <- stab3 %>%
    group_by(lab_id, tumor_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(
      n > cutoff,
      tumor_type %in% endometrial_cancers()
    ) %>%
    arrange(-n)
  nhypermut.wgs <- stab4 %>%
    group_by(lab_id, tumor_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(
      n > cutoff,
      tumor_type %in% endometrial_cancers()
    ) %>%
    arrange(-n)
  nhypermut <- bind_rows(nhypermut.wes, nhypermut.wgs) %>%
    arrange(-n) %>%
    group_by(tumor_type) %>%
    summarize(n = length(unique(lab_id))) %>%
    mutate(abbrev = c("oe", "ue"))
  filter(nhypermut, abbrev == tumor.type) %>%
    pull(n)
}

freqfun <- function(stab, gene_) {
  num <- stab %>%
    group_by(tumor_type, lab_id) %>%
    summarize(
      has_driver = any(gene %in% gene_),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    group_by(tumor_type) %>%
    summarize(n = sum(has_driver))
  freq <- paste0(num$n[1], "|", num$n[2])
  freq
}

#' @export
mutations_overall <- function(stab3, stab4, endo = TRUE, cutoff = 500) {
  if (endo) {
    tumor.types <- endo_cancers()
    labels <- endo_labels()
  } else {
    tumor.types <- mucinous_cancers()
    labels <- mucinous_labels()
  }
  names(labels) <- tumor.types
  nmut.wes <- stab3 %>%
    group_by(lab_id, tumor_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(tumor_type %in% tumor.types) %>%
    arrange(-n)
  nmut.wgs <- stab4 %>%
    group_by(lab_id, tumor_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(tumor_type %in% tumor.types) %>%
    arrange(-n)
  mut.overall <- bind_rows(nmut.wes, nmut.wgs) %>%
    arrange(-n) %>%
    filter(n <= cutoff) %>%
    group_by(tumor_type) %>%
    summarize(
      med = median(n),
      .groups = "drop"
    ) %>%
    mutate(tumor_type = labels[.$tumor_type])
  mut.oall <- setNames(
    mut.overall$med,
    mut.overall$tumor_type
  )
  mut.oall
}

#' @export
get_drivers <- function() {
  drivers <- tibble(gene = c(
    "PTEN",
    "PIK3CA",
    "KRAS",
    "PIK3R1",
    "NF1",
    "ARID1A"
  ))
  drivers
}

get_stats <- function(hyp.ue, hyp.oe) {
  a <- b <- 1
  ue <- rbeta(10e3, a + hyp.ue[1], b + hyp.ue[2])
  oe <- rbeta(10e3, a + hyp.oe[1], b + hyp.oe[2])
  post <- oe / ue
  fold <- round(mean(post), 2)
  cred.int <- quantile(post, c(0.025, 0.975)) %>%
    round(2) %>%
    paste(collapse = "-")
  stats <- paste0(
    fold, "-fold decrease, 95\\% posterior credible interval (CI): ",
    cred.int
  )
  stats
}

#' @export
endo_cancers <- function() {
  c(
    "ovarian endometrioid",
    "uterine endometrioid"
  )
}
endo_labels <- function() c("OE", "UE")

#' @export
rbind_endo <- function(stab3, stab4) {
  e.cancers <- endo_cancers()
  stab <- bind_rows(stab3, stab4) %>%
    filter(tumor_type %in% e.cancers)
  stab
}

tumors_with_drivers <- function(stab, manifest, drivers, tumor.types) {
  drivers <- filter(stab, gene %in% drivers$gene) %>%
    group_by(tumor_type, lab_id) %>%
    summarize(
      has_driver = any(gene %in% drivers$gene),
      .groups = "drop"
    ) %>%
    ungroup()
  drivers2 <- select(manifest, lab_id, tumor_type) %>%
    filter(tumor_type %in% tumor.types) %>%
    left_join(drivers, by = c("lab_id", "tumor_type")) %>%
    mutate(has_driver = ifelse(is.na(has_driver), FALSE, has_driver)) %>%
    group_by(tumor_type) %>%
    summarize(
      have_driver = sum(has_driver),
      n = n()
    ) %>%
    mutate(proportion = have_driver / n)
  drivers2
}

#' @export
tumors_with_drivers2 <- function(...) {
  e.labels <- c("OE", "UE")
  tumors_with_drivers(...) %>%
    pull(have_driver) %>%
    setNames(e.labels)
}

#' @export
prop_of_tumors_with_driver <- function(stab, manifest, drivers,
                                       tumor.types,
                                       labels = c("OE", "UE")) {
  drivers2 <- tumors_with_drivers(stab, manifest, drivers, tumor.types)
  prop <- drivers2$proportion %>%
    "*"(100) %>%
    round(1) %>%
    setNames(labels)
  prop
}

#' @export
get_esr1 <- function(stab3) {
  esr1 <- stab3 %>%
    filter(
      gene == "ESR1",
      tumor_type == "uterine endometrioid"
    ) %>%
    group_by(subject_id, lab_id) %>%
    summarize(
      gene = unique(gene),
      protein = unique(amino_acid_change),
      .groups = "drop"
    )
  esr1
}

#' @export
number_esr1 <- function(esr1) {
  n.esr1 <- length(unique(esr1$subject_id)) %>%
    broman::spell_out()
  n.esr1
}

#' @export
number_esr1_hotspots <- function(esr1) {
  n.hotspot <- esr1 %>%
    filter(grepl("^537", protein)) %>%
    group_by(subject_id) %>%
    summarize(n = n()) %>%
    nrow() %>%
    broman::spell_out()
  n.hotspot
}

#' @export
mucinous_mutations <- function(stab3, stab4) {
  cutoff <- 500
  nmut.wes <- stab3 %>%
    group_by(lab_id, tumor_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(
      n > cutoff,
      tumor_type %in% mucinous_cancers()
    ) %>%
    arrange(-n)
  nmut.wgs <- stab4 %>%
    group_by(lab_id, tumor_type) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(
      n > cutoff,
      tumor_type %in% mucinous_cancers()
    ) %>%
    arrange(-n)
  nmut <- bind_rows(nmut.wes, nmut.wgs) %>%
    arrange(-n)
  nmut
}

#' @export
get_ccnd1 <- function(idat.mucinous) {
  ccnd1 <- idat.mucinous %>%
    filter(gene == "CCND1") %>%
    group_by(alteration, tumor_type) %>%
    summarize(
      n = length(unique(lab_id)),
      .groups = "drop"
    )
  ccnd1
}

#' @export
get_erbb2 <- function(idat.mucinous, tumors) {
  tumor.levels <- c("Ovarian mucinous", "Colorectal mucinous")
  tumortypes <- tumors %>%
    select(subject_id, lab_id, tumor_type) %>%
    ungroup() %>%
    distinct()
  erbb2 <- idat.mucinous %>%
    filter(gene == "ERBB2") %>%
    group_by(tumor_type) %>%
    summarize(
      n = length(unique(lab_id)),
      .groups = "drop"
    )
  erbb2
}

#' @export
erbb2_amp <- function(erbb2) {
  erbb2.n <- filter(erbb2, tumor_type == "Ovarian mucinous") %>%
    pull(n)
  erbb2.n.amp <- filter(
    idat.mucinous,
    gene == "ERBB2",
    tumor_type == "Ovarian mucinous",
    alteration == "amplification"
  ) %>%
    nrow()
  erbb2.n.amp
}

#' @export
get_cdkn2 <- function(idat.mucinous) {
  cdkn2 <- idat.mucinous %>%
    filter(grepl("CDKN2[AB]", gene)) %>%
    group_by(gene, alteration, tumor_type) %>%
    summarize(
      n = length(unique(lab_id)),
      .groups = "drop"
    )
  cdkn2
}

#' @export
get_crc_drivers <- function(idat.mucinous) {
  crc <- idat.mucinous %>%
    filter(
      tumor_type == "Colorectal mucinous",
      gene %in% c("CDKN2A", "ERBB2", "CCND1"),
      type == "copynumber"
    )
  crc
}

#' @export
number_crc_not_hypermutated <- function(idat.mucinous, hypermut) {
  total.crc <- idat.mucinous %>%
    filter(tumor_type == "Colorectal mucinous") %>%
    select(lab_id) %>%
    distinct() %>%
    filter(!lab_id %in% hypermut$lab_id) %>%
    nrow()
  total.crc
}

#' @export
number_cdkn2a_mutations <- function(cdkn2a) {
  nmut <- filter(
    cdkn2a, tumor_type == "Ovarian mucinous",
    alteration == "mutation"
  ) %>%
    pull(n)
  nmut
}

#' @export
unique_tumortypes <- function(erbb2, ccnd1, cdkn2a) {
  tumor.types <- c(
    erbb2$tumor_type,
    ccnd1$tumor_type,
    cdkn2a$tumor_type
  ) %>%
    unique()
  tumor.types
}

#' @export
crc_credible_interval <- function(crc, types2, nmut) {
  a <- 0.5
  b <- 0.5
  prev.crc <- rbeta(
    10e3, nrow(crc) + a,
    types2["colorectal"] - nrow(crc) + b
  )
  prev.ov <- rbeta(
    10e3, nmut + a,
    types2["ovarian mucinous"] - nmut + b
  )
  fold.lower <- sprintf("%.2f", round(mean(prev.crc / prev.ov), 2))
  ci <- round(quantile(prev.crc, c(0.025, 0.975)), 2)
  ci2 <- sprintf("%.2f", ci)
  ci3 <- paste(ci2, collapse = "-")
}

#' @export
unique_alterations <- function(idat, gene_symbol, alt) {
  x <- filter(idat, gene == gene_symbol) %>%
    filter(alteration %in% alt) %>%
    filter(tumor_type == "Ovarian mucinous") %>%
    group_by(lab_id) %>%
    summarize(mutation = unique(alteration))
  x
}

#' @export
ras_alterations <- function(idat.mucinous) {
  ras <- idat.mucinous %>%
    filter(
      pathway %in% c("Ras and TK receptors"),
      type %in% c("mutation", "copynumber"),
      tumor_type == "Ovarian mucinous"
    ) %>%
    group_by(lab_id) %>%
    summarize(
      alteration = paste(unique(alteration), collapse = ","),
      .groups = "drop"
    )
  ras
}

#' @export
ras_prev_ci <- function(ras, types2) {
  ratio <- paste0(nrow(ras), " of ", types2["ovarian mucinous"])
  pct <- scales::percent(nrow(ras) / types2["ovarian mucinous"])
  ras.txt <- paste0(ratio, " (", pct, ")")
  ras.prev <- rbeta(
    10e3, nrow(ras) + 0.5,
    types2["ovarian mucinous"] - nrow(ras) + 0.5
  )
  mn.prev <- round(mean(ras.prev), 2)
  ci.ras <- quantile(ras.prev, c(0.025, 0.975)) %>%
    round(2) %>%
    paste(collapse = "-")
  ci.ras
}

read_s2 <- function(s2file) {
  stab2 <- read_csv(s2file, show_col_types = FALSE) %>%
    clean_colnames3() %>%
    set_colnames(c(
      "subject_id", "lab_id", "tumor.normal",
      "platform", "reads", "mapq"
    ))
  stab2
}

subset_to_ovmuc <- function(stab2, ovmuc_ids) {
  ovmuc <- filter(stab2, subject_id %in% ovmuc_ids) %>%
    select(subject_id, lab_id, tumor.normal) %>%
    filter(tumor.normal == "tumor")
}

#' @export
read_bayes <- function(bfile) {
  bayes <- readRDS(bfile) %>%
    filter(gene != "hypermutator") %>%
    unite("uid", c(gene, tumors), remove = FALSE)
  bayes
}

#' @export
posterior_mean <- function(bayes, gene.id, comparison) {
  post.mean <- bayes %>%
    filter(gene == gene.id, tumors == comparison) %>%
    select(mean) %>%
    round(2) %>%
    pull(mean)
  post.mean
}

posterior_ci <- function(bayes, gene.id, comparison) {
  post.ci <- bayes %>%
    filter(gene == gene.id, tumors == comparison) %>%
    select(`5%`, `95%`) %>%
    round(2) %>%
    paste(collapse = "-")
  post.ci
}

posterior_ci2 <- function(bayes, gene.id, comparison) {
  post.ci <- bayes %>%
    filter(gene == gene.id, tumors == comparison) %>%
    select(`5%`, `95%`) %>%
    mutate(across(everything(), function(x) -round(x, 2))) %>%
    rev() %>%
    paste(collapse = "-")
  post.ci
}

#' @export
cdkn2a_ci <- function(bayes) {
  posterior_ci2(bayes, "CDKN2A", "Ovarian mucinous vs GI mucinous")
}

#' @export
tp53_ci <- function(bayes) posterior_ci2(bayes, "TP53", "Ovarian mucinous vs GI mucinous")
#' @export
kras_ci <- function(bayes) posterior_ci2(bayes, "KRAS", "Ovarian mucinous vs GI mucinous")
#' @export
erbb2_ci <- function(bayes) posterior_ci2(bayes, "ERBB2", "Ovarian mucinous vs GI mucinous")
apc_ci <- function(bayes) posterior_ci(bayes, "APC", "Ovarian mucinous vs GI mucinous")

#' @export
update_comparison_labels <- function(bayes, nwrap = 25) {
  ue.v.oe <- "Uterine endometrioid - Ovarian endometrioid" %>%
    str_wrap(nwrap)
  gi.v.om <- "GI mucinous - Ovarian mucinous" %>%
    str_wrap(nwrap)
  bayes %>%
    filter(gene != "hypermutator") %>%
    unite("uid", c(gene, tumors), remove = FALSE) %>%
    mutate(
      tumors = str_replace_all(tumors, "vs", "-"),
      tumors = case_when(
        tumors == "Ovarian endometrioid - Uterine endometrioid" ~ ue.v.oe,
        tumors == "Ovarian mucinous - GI mucinous" ~ gi.v.om,
        TRUE ~ tumors
      ),
      tumors = str_wrap(tumors, 25)
    ) %>%
    filter(tumors %in% c(ue.v.oe, gi.v.om))
  bayes
}

#' @export
number_signif_differences <- function(bayes) {
  stats2 <- bayes %>%
    mutate(
      lor = mean,
      low = `5%`,
      high = `95%`
    )
  nsignif <- stats2 %>%
    group_by(tumors) %>%
    summarize(
      nsignif = sum(`5%` > 0 | `95%` < 0),
      n = n()
    )
  nsignif
}

#' @export
number_signif_random <- function(nsignif, alpha = 0.1) {
  ndrgenes <- nsignif$n
  rndm <- colSums(matrix(
    runif(1000 * ndrgenes[1]),
    ndrgenes[1],
    1000
  ) < alpha)
  rndm
}

subset_comparisons <- function(bayes.out, comparisons) {
  bayes.out2 <- remove_brackets_parameter(bayes.out) %>%
    filter(tumors %in% comparisons)
  bayes.out2
}

remove_brackets_parameter <- function(bayes.out) {
  bayes.out2 <- bayes.out %>%
    mutate(
      tmp = str_replace_all(parameter, "\\[", ""),
      tmp = str_replace_all(tmp, "\\]", ""),
      parameter = tmp
    ) %>%
    select(-tmp)
  bayes.out2
}

summarize_group1_mutation_rate <- function(bayes.out2) {
  group1 <- filter(bayes.out2, parameter == "theta1") %>%
    select(
      gene, pathway, ct,
      "5%", "50%", "95%"
    ) %>%
    mutate(
      cancertype = map_chr(ct, function(x) rownames(x)[1]),
      empirical = map_dbl(ct, function(x) x[1, 1] / x[1, 2]),
      label = map_chr(ct, function(x) paste0(x[1, 1], "/", x[1, 2]))
    )
  group1
}

summarize_group2_mutation_rate <- function(bayes.out2) {
  group2 <- filter(bayes.out2, parameter == "theta2") %>%
    select(
      gene, pathway, ct,
      "5%", "50%", "95%"
    ) %>%
    mutate(
      cancertype = map_chr(ct, function(x) rownames(x)[2]),
      empirical = map_dbl(ct, function(x) x[2, 1] / x[2, 2]),
      label = map_chr(ct, function(x) paste0(x[2, 1], "/", x[2, 2]))
    )
  group2
}

num_subj_w_meth <- function(manifest) {
  select(
    methylation, subject_id, lab_id,
    tumor.normal, tumor_type
  ) %>%
    filter(tumor.normal == "tumor") %>%
    group_by(tumor_type) %>%
    summarize(n = length(unique(subject_id))) %>%
    mutate(is_gi = tumor_type %in% c(
      "colorectal", "pancreas",
      "stomach"
    ))
}

wilcoxTest <- function(x) {
  broom::tidy(wilcox.test(x[, 1], x[, 2], paired = TRUE))
}

#' @export
prop_methylated_anova <- function(mfile) {
  methprop <- readRDS(mfile)
  meth.matrix.list <- pairedMeth(methprop, manifest)
  wilcox.stats <- meth.matrix.list %>%
    map_dfr(wilcoxTest)
  minp <- round(min(wilcox.stats$p.value), 2)
  anovastat <- summary(aov(propmeth ~ Diagnosis, data = methprop))
  list(minp = minp, anovastat = anovastat)
}

#' @export
read_lda <- function(lda.file, manifest) {
  mfest <- select(manifest, lab_id, subject_id)
  probLDA <- read_csv(lda.file, show_col_types = FALSE) %>%
    mutate(tumor.normal = factor(tumor.normal, c("Normal", "Tumor"))) %>%
    left_join(mfest, by = "lab_id") %>%
    mutate(prob.mucinous = `Colorectal mucinous` + `Stomach mucinous`)
  probLDA
}

#' Process LDA posteriors tibble into the probLDA format
#'
#' Accepts the raw posterior probability tibble from
#' \code{project_jhu_posteriors}, factors \code{tumor.normal}, joins
#' \code{subject_id} from the manifest, and adds \code{prob.mucinous}.
#' Use this instead of \code{read_lda} when the posteriors are already
#' in memory (e.g., as a targets pipeline value).
#'
#' @export
process_lda_posteriors <- function(posteriors, manifest) {
  mfest <- select(manifest, lab_id, subject_id)
  posteriors %>%
    mutate(tumor.normal = factor(tumor.normal, c("Normal", "Tumor"))) %>%
    left_join(mfest, by = "lab_id") %>%
    mutate(prob.mucinous = `Colorectal mucinous` + `Stomach mucinous`)
}

#' @export
ovarian_muc_lda <- function(probLDA) {
  filter(probLDA, Groups == "Ovarian mucinous") %>%
    filter(tumor.normal == "Tumor")
}

#' @export
endometrial_like <- function(ov.muc) {
  endo.like <- filter(ov.muc, `Uterine endometrial` > 0.6) %>%
    select(`Uterine endometrial`, tumor.normal, lab_id, subject_id)
}

#' @export
mucinous_like <- function(ov.muc) {
  filter(ov.muc, prob.mucinous >= 0.6) %>%
    select(lab_id, subject_id)
}

#' @export
summarize_wnt_pi3k <- function(clinical, idat.endometrioid,
                               idat.mucinous) {
  idat <- bind_rows(idat.endometrioid, idat.mucinous)
  idat2 <- idat %>%
    group_by(lab_id) %>%
    summarize(
      wnt = ifelse(any(pathway == "WNT"), 1, 0),
      pi3k = ifelse(any(gene != "PPP2R1A" & pathway == "PI3K"), 1, 0)
    )
  cdat <- clinical %>%
    left_join(
      idat2,
      join_by(lab_id)
    ) %>%
    filter(!is.na(os)) %>%
    select(lab_id, tumor_type, os, sex, stage, is_alive, wnt, pi3k) %>%
    mutate(
      event = ifelse(is_alive, 0, 1),
      wnt = ifelse(is.na(wnt), 0, wnt),
      pi3k = ifelse(is.na(pi3k), 0, pi3k),
      os = os / 365
    )
  cdat
}

matched_tumor_normal_pairs <- function(manifest) {
  subjects <- filter(manifest, platform == "WGS") %>%
    filter(tumor.normal == "tumor") %>%
    pull(subject_id)
  nsubj <- length(subjects)
  matched <- filter(manifest, platform == "WGS") %>%
    group_by(subject_id) %>%
    nest()
  matched1 <- matched[map_int(matched$data, nrow) > 1, ]
  isPair <- function(x) all(c("normal", "tumor") %in% x$tumor.normal)
  if (!all(map_lgl(matched1$data, isPair))) stop("subjects with missing paired data")
  matched1
}

#' @export
number_tcga_cancers <- function(se.tcga) {
  nt <- colData(se.tcga) %>%
    as_tibble() %>%
    filter(diagnosis != "Pancreatic mucinous") %>%
    mutate(diagnosis = as.character(diagnosis)) %>%
    pull(diagnosis) %>%
    table()
  nt2 <- paste(as.integer(nt), tolower(names(nt))) %>%
    paste(collapse = ", ") %>%
    str_replace("46", "and 46")
  nt2
}

#' @export
purity_exclusion <- function(manifest, threshold = 0.2) {
  exclusion1 <- filter(manifest, is.na(purity), is_na_purity)
  exclusion2 <- filter(manifest, !is.na(purity), purity <= 0.2)
  excluded <- bind_rows(exclusion1, exclusion2)
  excluded
}

#' @export
cellularity <- function(s3file, s4file) {
  tabs3 <- read_tsv(s3file,
    show_col_types = FALSE
  )
  tabs4 <- read_tsv(s4file,
    show_col_types = FALSE
  )
  wes <- tabs3 %>%
    group_by(`Lab id`) %>%
    summarize(max.maf = max(MAF))
  wgs <- tabs4 %>%
    group_by(`Lab id`) %>%
    summarize(max.maf = max(MAF))
  dat <- bind_rows(wes, wgs)
  return(dat)
}

#' @export
get_idat_muc <- function(manifest) {
  idat.muc <- filter(idat.mucinous, lab_id %in% manifest$lab_id)
  return(idat.muc)
}

#' @export
get_idat_endo <- function(manifest) {
  idat.endo <- filter(idat.endometrioid, lab_id %in% manifest$lab_id)
  return(idat.endo)
}

#' @export
get_hypermut <- function() {
  return(hypermut)
}

get_methylation <- function(manifest) {
  keep <- colnames(methylation) %in% manifest$lab_id | methylation$study == "TCGA"
  meth <- methylation[, keep]
  return(meth)
}

#' @export
get_methylation_se <- function(manifest) {
  keep <- colnames(methylation_se) %in% manifest$lab_id | methylation_se$study == "TCGA"
  meth.se <- methylation_se[, keep]
  return(meth.se)
}

#' @export
get_clinical <- function(manifest) {
  clindat <- filter(clinical, lab_id %in% manifest$lab_id)
  return(clindat)
}

odds_ratio <- function(gene, x) {
  gene1 <- filter(x, gene_symbol == gene)
  or <- round(gene1$median_odds, 1)
  ci <- c(gene1$lower_odds, gene1$upper_odds) %>%
    round(1) %>%
    paste(collapse = "-")
  list(or = or, ci = ci)
}

#' @export
number_wes_tumors <- function(manifest) {
  filter(manifest, platform == "WES") %>%
    select(subject_id, lab_id, tumor.normal) %>%
    distinct() %>%
    filter(tumor.normal == "tumor") %>%
    nrow()
}

#' @export
number_wgs_tumors <- function(manifest) {
  filter(manifest, platform == "WGS") %>%
    select(subject_id, lab_id, tumor.normal) %>%
    distinct() %>%
    filter(tumor.normal == "tumor") %>%
    nrow()
}

low_cellularity <- function(mafs, min.maf = 0.2) {
  low.cell <- filter(mafs, max.maf < min.maf) %>%
    distinct()
  low.cell
}

#' @export
median_coverage <- function(dat) {
  cov <- dat %>%
    group_by(Platform) %>%
    summarize(
      median = median(`Expected coverage`),
      median = round(median, 0),
      median = paste0(median, "X")
    )
  cov
}

#' @export
tumor_types <- function(tumors) {
  tumortypes <- tumors %>%
    select(subject_id, lab_id, tumor_type, tumor.normal) %>%
    distinct()
  tumortypes
}

#' @export
read_stable <- function(filename, tumortypes, manifest) {
  dat <- read_tsv(filename, show_col_types = FALSE) %>%
    clean_colnames3() %>%
    left_join(tumortypes, by = "lab_id") %>%
    filter(lab_id %in% manifest$lab_id)
  dat
}

#' @export
number_survival <- function(clinical) {
  n.surv <- clinical %>%
    select(subject_id, lab_id, os) %>%
    filter(!is.na(os)) %>%
    nrow()
  n.surv
}

crc_subjects <- function(tumors) {
  filter(tumors, grepl("colorectal", tumor_type)) %>%
    pull(subject_id)
}

#' @export
summarize_subtypes <- function(mfest, hmut) {
  mfest.tumors <- get_tumors(mfest)
  types <- number_tumor_types(mfest.tumors, hmut)
  N <- sum(types$n)
  tumortypes <- tumor_types(mfest.tumors)
  types2 <- types_as_vector(types)
  ncrc <- types$n[types$tumor_type == "colorectal"]
  ovmuc_ids <- ov_muc_subjects(mfest.tumors)
  crc_ids <- crc_subjects(mfest.tumors)
  n.types <- numbers_in_parenthesis(types2)
  n.types2 <- collapse_gi_counts(types2)
  result <- list(
    tumors = mfest.tumors,
    types = types,
    tumortypes = tumortypes,
    N = N,
    types2 = types2,
    ncrc = ncrc,
    ovmuc_ids = ovmuc_ids,
    crc_ids = crc_ids,
    n.types = n.types,
    n.types2 = n.types2
  )
  return(result)
}

#' @export
purity_numbermut <- function(stab3, stab4, mfest) {
  mfest2 <- mfest %>%
    filter(tumor.normal == "tumor") %>%
    select(lab_id, purity)
  stab <- bind_rows(stab3, stab4) %>%
    group_by(lab_id) %>%
    summarize(nmut = n())
  stab2 <- left_join(mfest2, stab, by = "lab_id") %>%
    mutate(nmut = ifelse(is.na(nmut), 0, nmut))
  spearman <- cor(stab2$nmut, stab2$purity, use = "complete")
  spearman
}

fig5_manifest <- function(manifest) {
  manifest2 <- manifest %>%
    cancer_names() %>%
    select(-tumor_type) %>%
    rename(tumor_type = tumor) %>%
    filter(
      !lab_id %in% hypermut$lab_id,
      tumor.normal == "tumor"
    ) %>%
    collapse_gi()
  tmp <- manifest2 %>%
    group_by(subject_id) %>%
    nest()
  tmp$data <- tmp$data %>%
    map(function(x) x[1, ])
  manifest2 <- unnest(tmp, "data")
  manifest2
}

get_denominators <- function(mfest) {
  mfest %>%
    group_by(tumor_type) %>%
    summarize(n = length(unique(subject_id)))
}
