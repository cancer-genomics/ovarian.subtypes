## Unit tests for the mutational-signature port (card BASE-02).
##
## parse_mutation_coords()/signature_input_table() are pure R (no
## deconstructSigs/BSgenome dependency) and are tested directly here.
## fit_mutational_signatures()/endo_signature_matrix()/muc_signature_matrix()
## require `deconstructSigs` and `BSgenome.Hsapiens.UCSC.hg18` -- both are
## skipped when those packages are unavailable (see this card's Notes/log for
## why they could not be installed in the environment this port was written
## in) rather than faked.

library(ovarian.subtypes)
library(tibble)
library(dplyr)

# ── parse_mutation_coords ───────────────────────────────────────────────────

test_that("parse_mutation_coords parses the chrN_start-end_REF_ALT format", {
  out <- parse_mutation_coords("chr16_82782707-82782707_G_A")
  expect_equal(out$chr, "chr16")
  expect_equal(out$pos, 82782707L)
  expect_equal(out$ref, "G")
  expect_equal(out$alt, "A")
})

test_that("parse_mutation_coords parses the chrN:pos_REF/ALT format", {
  out <- parse_mutation_coords("chr2:231904809_G/A")
  expect_equal(out$chr, "chr2")
  expect_equal(out$pos, 231904809L)
  expect_equal(out$ref, "G")
  expect_equal(out$alt, "A")
})

test_that("parse_mutation_coords returns NA rows for unparseable strings", {
  out <- parse_mutation_coords(c("FAILED", "chr17_C_T"))
  expect_true(all(is.na(out$chr)))
  expect_true(all(is.na(out$pos)))
})

test_that("parse_mutation_coords handles a mix of formats and failures", {
  out <- parse_mutation_coords(c(
    "chr16_82782707-82782707_G_A",
    "chr2:231904809_G/A",
    "FAILED"
  ))
  expect_equal(out$chr, c("chr16", "chr2", NA_character_))
  expect_equal(out$pos, c(82782707L, 231904809L, NA_integer_))
})

# ── signature_input_table ───────────────────────────────────────────────────

make_mutations <- function() {
  tibble(
    lab_id = c("CGOV1T", "CGOV1T", "CGOV2T", "CGOV3T", "CGOV3T"),
    mutation = c(
      "chr16_82782707-82782707_G_A",
      "chr2:231904809_G/A",
      "chr1_100-105_AT_G",          # Deletion -- excluded by type filter
      "FAILED",                      # unparseable -- excluded by coord filter
      "chr7_42916463-42916463_C_T"
    ),
    type = c("Substitution", "Substitution", "Deletion", "Substitution", "Substitution")
  )
}

test_that("signature_input_table keeps only substitutions in lab_ids", {
  out <- signature_input_table(make_mutations(), c("CGOV1T", "CGOV3T"))
  expect_setequal(unique(out$lab_id), c("CGOV1T", "CGOV3T"))
  expect_equal(nrow(out), 3L)
})

test_that("signature_input_table excludes samples outside lab_ids", {
  out <- signature_input_table(make_mutations(), "CGOV1T")
  expect_false("CGOV2T" %in% out$lab_id)
  expect_equal(nrow(out), 2L)
})

test_that("signature_input_table drops unparseable mutation strings", {
  out <- signature_input_table(make_mutations(), "CGOV3T")
  ## CGOV3T has one FAILED (dropped) row and one parseable substitution row.
  expect_equal(nrow(out), 1L)
  expect_equal(out$chr, "chr7")
})

test_that("signature_input_table returns chr/pos/ref/alt columns", {
  out <- signature_input_table(make_mutations(), "CGOV1T")
  expect_true(all(c("lab_id", "chr", "pos", "ref", "alt") %in% names(out)))
})

# ── fit_mutational_signatures / endo_signature_matrix / muc_signature_matrix ─
# Require deconstructSigs + BSgenome.Hsapiens.UCSC.hg18. Skipped (not faked)
# when unavailable -- see cards/BASE-02-port-mutational-signatures.md.

test_that("fit_mutational_signatures errors clearly without deconstructSigs", {
  skip_if(
    requireNamespace("deconstructSigs", quietly = TRUE),
    "deconstructSigs is installed; the informative-error path is not exercised"
  )
  expect_error(
    fit_mutational_signatures(signature_input_table(make_mutations(), "CGOV1T")),
    "deconstructSigs"
  )
})

test_that("endo/muc signature matrices reproduce committed extdata files", {
  skip_if_not_installed("deconstructSigs")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg18")
  skip_if_not(
    file.exists(system.file("..", "..", "extdata", "mutations.tsv", package = "ovarian.subtypes")),
    "extdata/mutations.tsv not found relative to package root"
  )
  ## Intentionally left as a skip-only placeholder: a full run + numeric
  ## comparison against the committed .rds files is this card's Do step 4/6,
  ## and belongs in an interactive/manual comparison (see Notes/log), not a
  ## unit test that would need the committed files bundled into the package.
  skip("full-pipeline comparison is manual; see cards/BASE-02 Notes/log")
})
