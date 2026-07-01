## Unit tests for pure utility functions (no external data required).

library(ovarian.subtypes)
library(tibble)
library(dplyr)

# ── format_number (type, x) ────────────────────────────────────────────────
# format_number(type, x) — looks up x$n for rows matching tumor_type %in% type,
# sums them, and returns a named "(n=N)" string.

test_that("format_number returns (n=N) label for a single tumor type", {
  df <- tibble(tumor_type = c("ovarian endometrioid", "colorectal"), n = c(40L, 20L))
  result <- format_number("ovarian endometrioid", df)
  expect_equal(unname(result), "(n=40)")
  expect_equal(names(result), "ovarian endometrioid")
})

test_that("format_number sums across multiple matching tumor types", {
  df <- tibble(tumor_type = c("colorectal", "pancreas", "stomach"), n = c(20L, 10L, 15L))
  result <- format_number(c("colorectal", "pancreas", "stomach"), df)
  expect_equal(unname(result), "(n=45)")
})

test_that("format_number returns (n=0) when type is absent", {
  df <- tibble(tumor_type = "ovarian mucinous", n = 50L)
  result <- format_number("colorectal", df)
  expect_equal(unname(result), "(n=0)")
})

# ── tumor_colors ─────────────────────────────────────────────────────────────

test_that("tumor_colors returns a named character vector", {
  tc <- tumor_colors()
  expect_type(tc, "character")
  expect_true(!is.null(names(tc)))
})

test_that("tumor_colors covers the six primary tumor types", {
  tc <- tumor_colors()
  expected <- c(
    "Ovarian endometrioid", "Uterine endometrioid",
    "Ovarian mucinous", "Colorectal mucinous",
    "Pancreatic mucinous", "Stomach mucinous"
  )
  expect_true(all(expected %in% names(tc)))
})

# ── purity_filter ────────────────────────────────────────────────────────────
# purity_filter keeps rows where is.na(is_na_purity) OR purity > threshold.
# Rows with is_na_purity = TRUE (confirmed NA purity) are EXCLUDED.

make_manifest <- function() {
  tibble(
    lab_id       = c("CGOV1T", "CGOV2T", "CGOV3T", "CGOV4T"),
    purity       = c(0.05, 0.21, NA, 0.80),
    is_na_purity = c(FALSE, FALSE, TRUE, FALSE)
  )
}

test_that("purity_filter excludes low-purity samples", {
  out <- purity_filter(make_manifest())
  expect_false("CGOV1T" %in% out$lab_id)
})

test_that("purity_filter keeps samples above threshold", {
  out <- purity_filter(make_manifest())
  expect_true("CGOV2T" %in% out$lab_id)
  expect_true("CGOV4T" %in% out$lab_id)
})

test_that("purity_filter excludes samples flagged is_na_purity = TRUE", {
  out <- purity_filter(make_manifest())
  expect_false("CGOV3T" %in% out$lab_id)
})

# ── purity_exclusion ─────────────────────────────────────────────────────────
# purity_exclusion returns the complement of purity_filter.

test_that("purity_exclusion returns samples with purity <= threshold", {
  out <- purity_exclusion(make_manifest())
  expect_true("CGOV1T" %in% out$lab_id)
  expect_false("CGOV2T" %in% out$lab_id)
  expect_false("CGOV4T" %in% out$lab_id)
})

test_that("purity_exclusion returns samples flagged is_na_purity = TRUE", {
  out <- purity_exclusion(make_manifest())
  expect_true("CGOV3T" %in% out$lab_id)
})

test_that("purity_filter and purity_exclusion partition the manifest", {
  m <- make_manifest()
  kept <- purity_filter(m)
  excluded <- purity_exclusion(m)
  expect_setequal(m$lab_id, c(kept$lab_id, excluded$lab_id))
})

# ── collapse_gi (data frame version, from mutations.R) ───────────────────────

test_that("collapse_gi relabels GI mucinous subtypes to 'GI mucinous'", {
  dat <- tibble(tumor_type = c(
    "Pancreas mucinous",
    "Stomach mucinous",
    "Colorectal mucinous",
    "Ovarian mucinous"
  ))
  out <- collapse_gi(dat)
  expect_true(all(out$tumor_type[1:3] == "GI mucinous"))
  expect_equal(out$tumor_type[4], "Ovarian mucinous")
})

test_that("collapse_gi preserves non-GI rows unchanged", {
  dat <- tibble(tumor_type = c(
    "Ovarian endometrioid",
    "Uterine endometrioid",
    "Ovarian mucinous"
  ))
  out <- collapse_gi(dat)
  expect_equal(out$tumor_type, dat$tumor_type)
})

# ── collapse_gi_counts (named-vector version, from manuscript.R) ──────────────

test_that("collapse_gi_counts sums GI cancers into a single 'GI' entry", {
  types2 <- c(
    "ovarian endometrioid" = 40L,
    "ovarian mucinous" = 50L,
    "uterine endometrioid" = 30L,
    "colorectal" = 20L,
    "pancreas" = 10L,
    "stomach" = 15L
  )
  out <- collapse_gi_counts(types2)
  expect_equal(unname(out["GI"]), 45L)
  expect_equal(unname(out["OE"]), 40L)
  expect_equal(unname(out["OM"]), 50L)
  expect_equal(unname(out["UE"]), 30L)
  expect_equal(length(out), 4L)
})

test_that("collapse_gi_counts returns vector named OE/OM/UE/GI", {
  types2 <- c(
    "ovarian endometrioid" = 1L, "ovarian mucinous" = 2L,
    "uterine endometrioid" = 3L, "colorectal" = 4L,
    "pancreas" = 5L, "stomach" = 6L
  )
  out <- collapse_gi_counts(types2)
  expect_named(out, c("OE", "OM", "UE", "GI"))
})

# ── cancer_names ─────────────────────────────────────────────────────────────

test_that("cancer_names adds a 'tumor' column with display names", {
  dat <- tibble(tumor_type = c(
    "ovarian endometrioid",
    "ovarian mucinous",
    "uterine endometrioid",
    "colorectal",
    "pancreas",
    "stomach"
  ))
  out <- cancer_names(dat)
  expect_true("tumor" %in% names(out))
})

test_that("cancer_names output has one row per input row", {
  dat <- tibble(tumor_type = c("ovarian endometrioid", "colorectal"))
  out <- cancer_names(dat)
  expect_equal(nrow(out), nrow(dat))
})
