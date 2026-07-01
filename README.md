# ovarian.subtypes

An R package providing functions and data for reproducing the numeric results and
figures in:

> Hallberg D, Eastman AC, Koul S, Bruhm DC, Papp E, et al.
> "Genomic Landscapes of Endometrioid and Mucinous Ovarian Cancers and
> Morphologically Similar Tumor Types."
> *Cancer Research Communications* 5(11):1952–1966 (2025).

The companion analysis repository — including the `targets` pipeline that
reproduces all manuscript statistics, workflowr analysis pages, and the
LaTeX manuscript source — is at
[cancer-genomics/Hallberg2025](https://github.com/cancer-genomics/Hallberg2025).

## Installation

```r
remotes::install_github("cancer-genomics/ovarian.subtypes")
```

## Data

The package ships PHI-free sample-level datasets:

| Object | Description |
|--------|-------------|
| `manifest` | 220 tumor/normal samples; CG lab IDs only |
| `clinical` | Overall survival and clinical covariates |
| `idat.endometrioid` | Integrated somatic alterations — endometrioid subtypes |
| `idat.mucinous` | Integrated somatic alterations — mucinous subtypes |
| `idat.gi` | Integrated somatic alterations — GI mucinous subtypes |
| `methylation` | CpG methylation fractions (JHU cohort) |
| `methylation_se` | SummarizedExperiment: JHU + TCGA methylation β-values |
| `hypermut` | Hypermutated samples (excluded from mutation analyses) |
| `gene.pathway` | Gene–pathway annotations |

All identifiers are CG lab IDs (e.g., `CGOV167T`, `CGCRC245T`).

## License

GPL-3
