clean_names <- function(x) {
  x <- str_replace_all(x, " ", "_") %>%
    tolower()
  x
}

#' @export
clean_colnames3 <- function(x) {
  nms <- colnames(x)
  nms2 <- clean_names(nms)
  x2 <- set_colnames(x, nms2)
  x2
}

clean_colnames2 <- function(x) {
  cnms <- colnames(x) %>%
    str_replace_all(" ", "_") %>%
    tolower()
  colnames(x) <- cnms
  x
}

format_number <- function(type, x) {
  n <- x$n[x$tumor_type %in% type] %>%
    sum()
  Ns <- paste0("(n=", n, ")")
  labels <- paste(type, collapse = ",")
  names(Ns) <- labels
  Ns
}

remove_author <- function(x) {
  i <- grep("^\\\\author\\{\\}$", x)
  if (length(i) != 0 && grepl("^\\\\date\\{", x[i + 1])) x <- x[-i]
  i <- grep("^\\\\title", x)
  line <- paste0(
    "\\title{Genomic landscapes of endometrioid and mucinous ovarian cancers",
    " and morphologically similar tumor types}"
  )
  if (length(i) != 0) x[i] <- line
  i <- grep("Velculescu", x)
  line <- x[i]
  thanks <- paste0(
    "Velculescu \\\\thanks{To whom correspondence should be addressed: ",
    "velculescu@jhmi.edu (V.E.V.) and rscharpf@jhu.edu (R.B.S.)}"
  )
  line <- stringr::str_replace(line, "Velculescu", thanks)
  x[i] <- line
  x
}

saveit <- function(..., string, file) {
  x <- list(...)
  names(x) <- string
  save(list = names(x), file = file, envir = list2env(x))
}

save_object <- function(object, nm) {
  filepath <- file.path("data", paste0(nm, ".rda"))
  saveit(object, string = nm, file = filepath)
  filepath
}

load2 <- function(nm) {
  load(file.path("data", nm), envir = parent.frame(2))
}

key <- function(manifest.list) {
  manifest.list[["key"]]
}

man <- function(manifest.list) {
  manifest.list[["manifest"]]
}

#' Accessor for beta values of methylation SummarizedExperiment
#' @export
beta <- function(se) assays(se)[["beta"]]

pick_notna <- function(x) {
  x2 <- x[!is.na(x)]
  x2[1]
}

pathway_levels <- function() {
  pathway.levels <- c(
    "PI3K", "Ras and TK receptors",
    "Chromatin Regulating",
    "Cell cycle",
    "Notch", "DNA repair",
    "Mismatch repair",
    "WNT", "TGFBR", "JAK/STAT",
    "Other", "Large gene"
  )
  pathway.levels
}

mucinous_pathways <- function(mucinous.levels.file) {
  muc.pathways <- c(
    "Ras and TK receptors",
    "PI3K",
    "Chromatin Regulating",
    "Cell cycle",
    "Notch",
    "DNA repair",
    "Mismatch repair",
    "WNT",
    "TGFBR", "JAK/STAT",
    "Other", "Large gene"
  )
  levels <- readRDS(mucinous.levels.file)
  levels$pathway <- muc.pathways
  levels
}

compare <- function(obj1, obj2) {
  all(obj1$lab_id %in% obj2$lab_id)
  all(obj1$gene %in% obj2$gene)
  identical(
    levels(obj1$lab_id),
    levels(obj2$lab_id)
  )
  identical(
    levels(obj1$gene),
    levels(obj2$gene)
  )
  identical(
    levels(obj1$pathway),
    levels(obj2$pathway)
  )
}

suppl5_varnames <- function(s5) {
  orignames <- colnames(s5)
  varnames <- clean_names(orignames)
  varnames[c(
    8, 9, 13, 14, 15,
    16:19
  )] <- c(
    "olap_snps",
    "olap_het_snps",
    "cnv_type",
    "loh",
    "focal",
    "clin_gene",
    "biol_gene",
    "olap_gene",
    "olap_tx"
  )
  varnames
}

#' Convert a list of FACETS segment tibbles to a list of GRanges
#'
#' Converts each segment data frame (with columns seqnames, start, end, strand,

#' Convert a list of FACETS segment GRanges or tibbles to a list of GRanges
#'
#' Accepts each segment element as either a \code{GRanges} (as returned
#' directly by \code{segmentBins} / \code{read_one} in the trellis pipeline)
#' or a tibble/data.frame with columns \code{seqnames}, \code{start},
#' \code{end}, \code{strand}, \code{seg.mean}, \code{acn}.
#'
#' Seqinfo is extracted from the first element of \code{deletions}.
#' If the deletion is a \code{StructuralVariant} S4 object (as stored by
#' the trellis pipeline), seqinfo is taken from \code{trellis::variant()};
#' otherwise \code{GenomeInfoDb::seqinfo()} is called directly.
#'
#' @param segments  Named list of segment GRanges or tibbles.
#' @param deletions Named list of deletion objects (GRanges or
#'   StructuralVariant); provides seqinfo.
#' @return Named list of GRanges, same length and names as \code{segments}.
#' @export
segs_to_granges <- function(segments, deletions) {
    ## Extract seqinfo from the first deletion — handle both StructuralVariant
    ## S4 objects (access via trellis::variant()) and plain GRanges.
    del1 <- deletions[[1]]
    si <- if (is(del1, "StructuralVariant")) {
        GenomeInfoDb::seqinfo(trellis::variant(del1))
    } else {
        GenomeInfoDb::seqinfo(del1)
    }
    purrr::map(segments, function(x) {
        ## If already a GRanges, just ensure seqinfo is attached.
        if (is(x, "GRanges")) {
            GenomeInfoDb::seqinfo(x) <- si
            return(x)
        }
        ## Otherwise construct from tibble/data.frame columns.
        g <- GenomicRanges::GRanges(x$seqnames,
                                    IRanges::IRanges(x$start, x$end),
                                    strand = x$strand)
        g$seg.mean <- x$seg.mean
        g$acn      <- x$acn
        GenomeInfoDb::seqinfo(g) <- si
        g
    })
}

load_tx <- function(build = "hg18") {
    tx <- trellis:::loadTx(build)
    GenomeInfoDb::keepSeqlevels(tx, paste0("chr", c(1:22, "X")),
                                pruning.mode = "coarse")
}
