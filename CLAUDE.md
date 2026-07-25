# ovarian.subtypes — the shared base package

**This is the one shared-context module.** Both companion data packages depend on it; it
imports neither. Anything you change here can invalidate the sequencing arm, the
methylation arm, and the manuscript pipeline at once — which is why `BASE` cards are
marked `exclusive: true` and run alone.

It has its own git repository (remote `cancer-genomics/ovarian.subtypes`) and is gitignored
by the outer repo. Commits land in both.

## The package pipeline

`ovarian.subtypes/_targets.R` (58 targets) builds the analytic datasets from raw inputs in
`inst/extdata/` (`manifest.rds`, `sdat.rds`, `center_info.csv`,
`diagnosis_surgery_dates.csv`) plus `../output/facets*` text files, producing `manifest`,
`clinical`, `discordant`, `idat.*`, `methylation_se`. Run it **from inside this
directory**. It writes `.rda` via `save_object()` with `PHI_COLS` stripped.

This is distinct from the *manuscript* pipeline (`_targets.R` at the repo root), which
consumes this package's data plus one `OsSeqExpData` accessor.

## `data/` vs `inst/extdata/` — the contract

| | Contents | Published? |
|---|---|---|
| `inst/extdata/` | private raw inputs — BAM paths, DOB, PGDX IDs, surgery dates | **never** |
| `data/` | the `.rda` objects consumers load, PHI columns already stripped | yes |

**Critical note:** the upstream generating scripts in `code/` and this `_targets.R` will
reproduce the **unstripped** originals if re-run. The `.rda` files currently committed are
the authoritative public versions. If a rebuild reintroduces `pgdx_id`, `bamfile`,
`bam_local`, `genotype_id`, `facet_id`, `size`, or a `basename`/`Basename` column, that is
the bug — not the committed file. The stripped-column table is in the root `CLAUDE.md`.

## `R/` modules

| Module | Lines | Covers |
|---|---:|---|
| `manuscript.R` | 1059 | `\Sexpr{}` value helpers, absorbed from `manuscript/functions.R` |
| `visualization.R` | 867 | plotting; includes `gg_blank`, `plot_circos_grob` |
| `mutations.R` | 650 | mutation calls, spectra, drivers |
| `manifest.R` | 599 | manifest assembly and PHI stripping |
| `methylation.R` | 568 | methylation helpers |
| `clinical.R` | 266 | clinical/survival |
| `assay_functions.R` | 253 | SE-building; `segs_to_granges`, `load_tx` |
| `utils.R` | 210 | shared helpers |
| `coverage.R` | 96 | coverage stats |
| `colors.R` | 41 | palette |
| `functions.R` | 8 | vestigial stub from the Phase-2 split |

`ampliconGraph`, `gg_blank`, `load_tx`, `plot_circos_grob`, `segs_to_granges` are the five
functions `OsSeqExpData/data-raw` reaches into. Changing any of their signatures breaks
the sequencing arm silently — it calls them with `ovarian.subtypes::`, so there is no
import-time check.

## Known defect

`_targets.R:36` sets `facets.dir` from `output/facets-trellis/jhpce_directories`, which
**does not exist**. A cold rebuild of `manifest` therefore fails, and `WGS_SAMPLES`
derives from `manifest`. Card `BASE-01` owns this; it is the highest-severity item in the
base package.

## Gate

`Rscript tests/verify_snapshot.R` from the repo root, as always. Because everything
depends on this package, treat a green gate here as necessary but not sufficient — also
run `devtools::test()` in this directory.
