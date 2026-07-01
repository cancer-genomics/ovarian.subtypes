#' Consolidate per-sample coverage statistics from JHPCE array job output
#'
#' Reads per-sample reportRAW and reportMAPQ30 files from the coverage stats
#' array job and per-sample read-length files from the read-length array job.
#' Joins through bam_lookup to replace BAM basenames (which may contain PGDX
#' IDs) with CG lab IDs, then joins with the manifest for sample metadata and
#' computes average coverage.
#'
#' @param mapq30_dir Path to directory containing reportRAW.* and reportMAPQ30.*
#'   files produced by the coverage stats array job.
#' @param read_length_dir Path to directory containing per-sample read-length
#'   files (one integer per file, filename = bam_basename + ".txt").
#' @param bam_lookup Data frame with columns \code{bam_basename} and
#'   \code{lab_id}.
#' @param manifest Data frame containing at minimum \code{lab_id},
#'   \code{subject_id}, \code{platform}, \code{tumor.normal}, \code{purity},
#'   and \code{is_na_purity}.
#' @return A data frame with columns \code{lab_id}, \code{subject_id},
#'   \code{tumor_normal}, \code{platform}, \code{reference_genome},
#'   \code{mapq30_reads}, \code{read_length}, \code{average_coverage},
#'   \code{purity}, \code{purity_status}.
#' @export
consolidate_coverage_stats <- function(mapq30_dir, read_length_dir,
                                       bam_lookup, manifest) {
    # Read reportMAPQ30 files (tab-separated, no header):
    #   col1: bam_basename  col2: mapq30_reads  col3: mapq30_bases
    #   col4: mapq30_expected_coverage
    mapq30_files <- list.files(mapq30_dir, pattern = "^reportMAPQ30\\.",
                               full.names = TRUE)
    m30 <- do.call(rbind, lapply(mapq30_files, function(f) {
        read.delim(f, header = FALSE, col.names = c(
            "bamfile", "mapq30_reads", "mapq30_bases", "mapq30_expected_coverage"
        ))
    }))

    # Read reportRAW files (same format)
    raw_files <- list.files(mapq30_dir, pattern = "^reportRAW\\.",
                            full.names = TRUE)
    raw <- do.call(rbind, lapply(raw_files, function(f) {
        read.delim(f, header = FALSE, col.names = c(
            "bamfile", "reads", "bases", "expected_coverage"
        ))
    }))
    raw <- raw[!is.na(raw$bamfile), ]

    # Read read-length files (single integer per file;
    #   filename pattern: <bam_basename>.txt)
    rl_files <- list.files(read_length_dir, full.names = TRUE)
    rl <- data.frame(
        bamfile     = sub("\\.txt$", "", basename(rl_files)),
        read_length = vapply(rl_files, function(f) {
            as.integer(readLines(f, n = 1L))
        }, integer(1L)),
        stringsAsFactors = FALSE
    )

    # Manifest columns needed
    mfest <- manifest[, c("lab_id", "subject_id", "platform",
                          "tumor.normal", "purity", "is_na_purity")]

    # Join raw -> bam_lookup to get lab_id; drop rows with no match
    raw_labelled <- merge(raw, bam_lookup,
                          by.x = "bamfile", by.y = "bam_basename",
                          all.x = FALSE)
    raw_labelled <- raw_labelled[!is.na(raw_labelled$lab_id), ]

    # Join m30 and read lengths by bamfile, then manifest by lab_id
    combined <- merge(raw_labelled, m30, by = "bamfile", all.x = TRUE)
    combined  <- merge(combined,    rl,  by = "bamfile", all.x = TRUE)
    combined  <- combined[combined$lab_id %in% mfest$lab_id, ]
    combined  <- merge(combined, mfest, by = "lab_id", all.x = TRUE)

    # Compute average coverage
    # WGS denominator: ~3.1 Gb (hg18 mappable); WES: ~30 Mb
    combined$bases           <- as.numeric(combined$mapq30_reads) * combined$read_length
    combined$total_bases     <- ifelse(combined$platform == "WGS",
                                       3095693983, 30e6)
    combined$average_coverage <- combined$bases / combined$total_bases

    # Derive purity status from is_na_purity flag
    combined$purity_status <- dplyr::case_when(
        is.na(combined$is_na_purity) ~ "Not applicable",
        combined$is_na_purity        ~ "Not estimable",
        !combined$is_na_purity       ~ "Processed successfully"
    )

    # Return tidy output with R-friendly column names
    out <- combined[, c("lab_id", "subject_id", "tumor.normal", "platform",
                        "mapq30_reads", "read_length", "average_coverage",
                        "purity", "purity_status")]
    names(out)[names(out) == "tumor.normal"] <- "tumor_normal"
    out$reference_genome <- "hg18"
    out[, c("lab_id", "subject_id", "tumor_normal", "platform",
            "reference_genome", "mapq30_reads", "read_length",
            "average_coverage", "purity", "purity_status")]
}
