## Data integrity: PHI-free columns, expected dimensions, valid values.
## These run against the committed .rda files — a regression guard ensuring
## that rerunning the package pipeline never reintroduces stripped columns.

library(ovarian.subtypes)
library(SummarizedExperiment)

PHI_COLS <- c("pgdx_id", "bamfile", "bam_local", "size", "genotype_id", "facet_id")

# ── manifest ─────────────────────────────────────────────────────────────────

test_that("manifest has no PHI columns", {
  data(manifest, package = "ovarian.subtypes")
  expect_false(any(PHI_COLS %in% names(manifest)))
})

test_that("manifest has expected core columns", {
  data(manifest, package = "ovarian.subtypes")
  expect_true(all(c(
    "lab_id", "subject_id", "tumor.normal",
    "platform", "tumor_type", "purity"
  ) %in% names(manifest)))
})

test_that("all manifest lab_ids are CG identifiers", {
  data(manifest, package = "ovarian.subtypes")
  expect_true(all(grepl("^CG", manifest$lab_id)))
})

test_that("manifest purity values are in (0, 1] or NA", {
  data(manifest, package = "ovarian.subtypes")
  p <- manifest$purity[!is.na(manifest$purity)]
  expect_true(all(p > 0 & p <= 1))
})

test_that("manifest has expected tumor types", {
  data(manifest, package = "ovarian.subtypes")
  expected <- c(
    "ovarian endometrioid", "ovarian mucinous",
    "uterine endometrioid", "colorectal", "pancreas", "stomach"
  )
  tumors <- unique(manifest$tumor_type[manifest$tumor.normal == "tumor"])
  expect_true(all(expected %in% tumors))
})

test_that("manifest has both tumor and normal samples", {
  data(manifest, package = "ovarian.subtypes")
  expect_true("tumor" %in% manifest$tumor.normal)
  expect_true("normal" %in% manifest$tumor.normal)
})

test_that("manifest has at least 200 rows", {
  data(manifest, package = "ovarian.subtypes")
  expect_gte(nrow(manifest), 200)
})

# ── discordant ───────────────────────────────────────────────────────────────

test_that("discordant has no PHI columns", {
  data(discordant, package = "ovarian.subtypes")
  expect_false(any(PHI_COLS %in% names(discordant)))
})

test_that("discordant lab_ids are CG identifiers", {
  data(discordant, package = "ovarian.subtypes")
  expect_true(all(grepl("^CG", discordant$lab_id)))
})

test_that("discordant samples are not in manifest", {
  data(manifest, package = "ovarian.subtypes")
  data(discordant, package = "ovarian.subtypes")
  # discordant tumors should be absent from the main manifest
  disc_tumors <- discordant$lab_id[discordant$discordant_tumor_type]
  expect_false(any(disc_tumors %in% manifest$lab_id))
})

# ── methylation (data frame) ──────────────────────────────────────────────────

test_that("methylation has no PHI columns", {
  data(methylation, package = "ovarian.subtypes")
  expect_false(any(PHI_COLS %in% names(methylation)))
})

test_that("methylation lab_ids are CG identifiers", {
  data(methylation, package = "ovarian.subtypes")
  expect_true(all(grepl("^CG", methylation$lab_id)))
})

test_that("methylation lab_ids are a subset of manifest lab_ids", {
  data(manifest, package = "ovarian.subtypes")
  data(methylation, package = "ovarian.subtypes")
  expect_true(all(methylation$lab_id %in% manifest$lab_id))
})

# ── methylation_se (SummarizedExperiment) ─────────────────────────────────────

test_that("methylation_se colData has no basename column", {
  data(methylation_se, package = "ovarian.subtypes")
  expect_false("basename" %in% names(SummarizedExperiment::colData(methylation_se)))
})

test_that("methylation_se has expected colData columns", {
  data(methylation_se, package = "ovarian.subtypes")
  cols <- names(SummarizedExperiment::colData(methylation_se))
  expect_true(all(c("lab_id", "diagnosis", "study", "tumor") %in% cols))
})

test_that("methylation_se has JHU and TCGA samples", {
  data(methylation_se, package = "ovarian.subtypes")
  study <- SummarizedExperiment::colData(methylation_se)$study
  expect_true("JHU" %in% study)
  expect_true("TCGA" %in% study)
})

test_that("methylation_se JHU sample names are CG identifiers", {
  data(methylation_se, package = "ovarian.subtypes")
  study <- SummarizedExperiment::colData(methylation_se)$study
  jhu_ids <- colnames(methylation_se)[study == "JHU"]
  expect_true(all(grepl("^CG", jhu_ids)))
})

# ── idat objects ─────────────────────────────────────────────────────────────

test_that("idat.endometrioid lab_ids are CG identifiers", {
  data(idat.endometrioid, package = "ovarian.subtypes")
  expect_true(all(grepl("^CG", idat.endometrioid$lab_id)))
})

test_that("idat.mucinous lab_ids are CG identifiers", {
  data(idat.mucinous, package = "ovarian.subtypes")
  expect_true(all(grepl("^CG", idat.mucinous$lab_id)))
})

test_that("idat.endometrioid has expected columns", {
  data(idat.endometrioid, package = "ovarian.subtypes")
  cols <- names(idat.endometrioid)
  expect_true(all(c("lab_id", "gene", "tumor_type", "alteration") %in% cols))
})

# ── hypermut ─────────────────────────────────────────────────────────────────

test_that("hypermut lab_ids are CG identifiers", {
  data(hypermut, package = "ovarian.subtypes")
  expect_true(all(grepl("^CG", hypermut$lab_id)))
})

test_that("hypermut samples are present in manifest or discordant", {
  data(manifest, package = "ovarian.subtypes")
  data(discordant, package = "ovarian.subtypes")
  data(hypermut, package = "ovarian.subtypes")
  known_ids <- c(manifest$lab_id, discordant$lab_id)
  expect_true(all(hypermut$lab_id %in% known_ids))
})

# ── clinical ─────────────────────────────────────────────────────────────────

test_that("clinical has expected columns", {
  data(clinical, package = "ovarian.subtypes")
  expect_true(all(c("lab_id", "subject_id", "tumor_type") %in% names(clinical)))
})

test_that("clinical lab_ids are CG identifiers", {
  data(clinical, package = "ovarian.subtypes")
  expect_true(all(grepl("^CG", clinical$lab_id)))
})

test_that("clinical lab_ids are a subset of manifest lab_ids", {
  data(manifest, package = "ovarian.subtypes")
  data(clinical, package = "ovarian.subtypes")
  expect_true(all(clinical$lab_id %in% manifest$lab_id))
})
