#' TCGA LDA scatter plot
#'
#' Builds a ggplot of TCGA samples projected onto LD1/LD2 with coloured
#' ellipses for the three training groups and decision-boundary line.
#'
#' @param obs Tibble of TCGA LD scores from \code{project_tcga_to_lda}.
#' @param ell Tibble of ellipse polygons from \code{project_tcga_to_lda}.
#' @param alpha_el Transparency for ellipse fill (default 0.4).
#' @return Named list: \code{grob} (ggplotGrob without shape legend) and
#'   \code{legend} (extracted shape legend).
#' @export
plot_tcga_lda <- function(obs, ell, alpha_el = 0.4) {
  dx.colors <- tumor_colors()
  groups <- c("Uterine endometrial", "Stomach mucinous", "Colorectal mucinous")
  fig <- obs %>%
    ggplot(aes(LD1, LD2)) +
    xlab("LD1") + ylab("LD2") +
    theme_bw(base_size = 20) +
    geom_point(data = obs,
               aes(color = Groups, shape = tumor.normal),
               size = 2, alpha = 1) +
    geom_polygon(data = ell, aes(group = Groups, fill = Groups),
                 alpha = alpha_el) +
    scale_fill_manual(values = dx.colors[groups]) +
    scale_color_manual(values = dx.colors[groups]) +
    scale_shape_manual(values = c("Tumor" = 19, "Normal" = 17)) +
    geom_vline(xintercept = -0.77, linetype = "dashed", color = "gray") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 22),
          panel.grid = element_blank()) +
    guides(color = "none", fill = "none", shape = guide_legend(title = "")) +
    annotate("text", 3, 5, label = "Stomach\nmucinous", size = 6) +
    annotate("text", 3, -3, label = "Colorectal\nmucinous", size = 6) +
    annotate("text", -3.6, 1.5, label = "Uterine\nendometrial", size = 6) +
    ggtitle("TCGA")
  ldafig.pre <- ggplotGrob(fig)
  legend <- ldafig.pre$grobs[[which(
    sapply(ldafig.pre$grobs, function(x) x$name) == "guide-box"
  )]]
  list(grob = ggplotGrob(fig + guides(shape = "none")), legend = legend)
}

#' Faceted JHU LDA scatter plot with TCGA ellipse background
#'
#' Filters \code{all2} to the five JHU cancer types, combines Pancreatic and
#' Stomach mucinous into a single facet panel, overlays TCGA ellipses as
#' background, and places each JHU cancer type in its own facet panel.
#'
#' @param all2 Tibble from \code{project_jhu_to_lda}.
#' @param ell Tibble of TCGA ellipse polygons from \code{project_tcga_to_lda}.
#' @param manifest Manifest data frame (for subject_id lookup).
#' @param alpha_el Transparency for ellipse fill (default 0.2).
#' @return Named list: \code{grob} (ggplotGrob) and \code{legend}
#'   (extracted shape legend).
#' @export
plot_jhu_lda <- function(all2, ell, manifest, alpha_el = 0.2) {
  dx.colors <- tumor_colors()
  border <- -0.77
  facet_levels <- c("Colorectal mucinous",
                    "Pancreatic / Stomach\nmucinous",
                    "Ovarian mucinous",
                    "Ovarian endometrioid")

  manifest2 <- select(manifest, subject_id, lab_id) %>%
    filter(lab_id %in% all2$lab_id) %>%
    distinct()

  all3 <- filter(all2, !Groups %in% c("TCGA", "Uterine endometrial")) %>%
    mutate(Groups = factor(Groups, c("Colorectal mucinous",
                                     "Pancreatic mucinous",
                                     "Stomach mucinous",
                                     "Ovarian endometrioid",
                                     "Ovarian mucinous")),
           group2 = Groups) %>%
    left_join(manifest2, by = "lab_id") %>%
    mutate(
      subject_id = ifelse(is.na(subject_id), lab_id, subject_id),
      group3 = ifelse(group2 %in% c("Pancreatic mucinous", "Stomach mucinous"),
                      "Pancreatic / Stomach\nmucinous", as.character(group2)),
      group3 = factor(group3, levels = facet_levels)
    )

  pan_stc_colors <- filter(all3, lab == "JHU",
                            Groups %in% c("Pancreatic mucinous", "Stomach mucinous")) %>%
    mutate(colors = sapply(as.character(Groups), function(x) dx.colors[x])) %>%
    pull(colors)

  fig <- all3 %>%
    ggplot(aes(LD1, LD2)) +
    xlab("LD1") + ylab("LD2") +
    theme_bw(base_size = 20) +
    geom_polygon(data = ell, aes(group = Groups, fill = Groups), alpha = alpha_el) +
    scale_fill_manual(
      values = dx.colors[c("Uterine endometrial", "Stomach mucinous", "Colorectal mucinous")]
    ) +
    facet_wrap(~group3, ncol = 2) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 22),
          axis.title.y = element_text(hjust = 1),
          panel.grid = element_blank()) +
    guides(color = "none", fill = "none") +
    geom_point(
      data = filter(all3, lab == "JHU", Groups == "Colorectal mucinous"),
      aes(shape = tumor.normal),
      color = dx.colors["Colorectal mucinous"], size = 3, alpha = 1
    ) +
    geom_point(
      data = filter(all3, lab == "JHU",
                    Groups %in% c("Pancreatic mucinous", "Stomach mucinous")),
      aes(shape = tumor.normal),
      color = pan_stc_colors, size = 3, alpha = 1
    ) +
    geom_point(
      data = filter(all3, lab == "JHU", Groups == "Ovarian mucinous"),
      aes(shape = tumor.normal),
      color = dx.colors["Ovarian mucinous"], size = 3, alpha = 1
    ) +
    geom_jitter(
      data = filter(all3, lab == "JHU", Groups == "Ovarian endometrioid"),
      aes(shape = tumor.normal),
      width = 0.8, height = 0.8,
      color = dx.colors["Ovarian endometrioid"], size = 3, alpha = 1
    ) +
    geom_vline(xintercept = border, linetype = "dashed", color = "gray") +
    guides(shape = guide_legend(title = ""))

  list(
    legend = cowplot::get_legend(fig),
    grob   = ggplotGrob(fig + guides(shape = "none"))
  )
}

#' @export update.reargraph
update.reargraph <- function(object, show_legend = TRUE, size = 15) {
  grobs <- object$grobs
  a <- grobs[["a"]] +
    theme(
      plot.title = element_text(size = size),
      panel.background = element_rect(fill = "gray95")
    )
  b <- grobs[["b"]] +
    theme(
      plot.title = element_text(size = size),
      panel.background = element_rect(fill = "gray95")
    )
  ag <- ggplotGrob(a)
  bg <- ggplotGrob(b)
  bg$widths <- ag$widths
  widths <- c(0.5, 0.5) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  heights <- c(0.95, 0.05) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  mat <- matrix(c(
    1, 2,
    3, 3
  ), byrow = TRUE, ncol = 2, nrow = 2)
  if (show_legend) {
    legend.grob <- grobs[["legend"]]
  } else {
    legend.grob <- grid::nullGrob()
    heights[2] <- unit(0, "npc")
  }
  gobj <- gridExtra::arrangeGrob(ag, bg,
    legend.grob,
    layout_matrix = mat,
    widths = widths,
    heights = heights
  )
  return(gobj)
}

#' @export
axis.labels <- function(ord_in, signif.digits = 1) {
  exp_var <- 100 * ord_in$svd^2 / sum(ord_in$svd^2)
  axes <- paste0("LD", 1:2)
  axes <- paste0(axes, " (", round(exp_var, signif.digits), "%)")
  axes
}

#' @export
my.ggord.lda <- function(ord_in, grp_in = NULL,
                         axes = c("1", "2"), ...) {
  obs <- data.frame(predict(ord_in)$x[, c("LD1", "LD2")]) %>%
    as_tibble() %>%
    mutate(lab = "TCGA")
  obs$Groups <- as.character(grp_in)
  obs
}

#' @export
my.ellipse <- function(obs) {
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ellipse_pro <- 0.95
  ell <- plyr::ddply(obs, "Groups", function(x) {
    if (nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$LD1, x$LD2))
    mu <- c(mean(x$LD1), mean(x$LD2))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = "+")) %>%
      as_tibble()
  })
  names(ell)[2:3] <- names(obs)[1:2]
  . <- plyr::.
  ell <- plyr::ddply(ell, .(Groups), function(x) x[chull(x$LD1, x$LD2), ])
  ell <- as_tibble(ell)
}

#' @export
fig1_ggplatform <- function(i, x, base_size = 15, colors, cancer) {
  dat <- x %>%
    mutate(
      platform = case_when(
        platform == "Methylation" ~ "Me",
        TRUE ~ platform
      ),
      platform = factor(platform,
        levels = c(
          "WGS", "WES",
          "Me"
        )
      )
    )
  orderby <- dat %>%
    group_by(subject_id) %>%
    summarize(
      nplatform = length(unique(platform)),
      tumor.normal = sum(tumor.normal == "normal,tumor"),
      nwgs = sum(platform == "WGS")
    ) %>%
    arrange(nplatform, tumor.normal, nwgs)
  dat2 <- dat %>%
    mutate(subject_id = factor(subject_id, orderby$subject_id))
  fig <- dat2 %>%
    ggplot(aes(platform, subject_id)) +
    geom_point(aes(fill = matched),
      pch = 21, size = 4
    ) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_x_discrete(drop = FALSE) +
    xlab("") +
    ylab("") +
    guides(fill = guide_legend(title = "")) +
    scale_fill_manual(
      values = colors,
      drop = FALSE
    ) +
    ggtitle(cancer[i])
  leg_raw <- cowplot::get_plot_component(fig, "guide-box",
    return_all = TRUE
  )[[1]]
  leg <- if (grid::is.grob(leg_raw)) leg_raw else leg_raw[[1]]
  fig2 <- fig + guides(fill = "none")
  result <- list(figure = fig2, legend = leg)
  result
}

find_group_samples <- function(sample_id, connection_matrix) {
  group_samples <- c()
  check_samples <- sample_id

  while (length(check_samples) != 0) {
    for (check_sample_id in check_samples) {
      if (!(check_sample_id %in% group_samples)) {
        group_samples <- c(group_samples, check_sample_id)
        local_samples <- names(which(connection_matrix[check_sample_id, ] == 1))
        check_samples <- setdiff(c(check_samples, local_samples), group_samples)
      }
    }
  }
  return(group_samples)
}

.ggRearrange2 <- function(df, ylabel = "Read pair index",
                          basepairs = 400, num.ticks = 5) {
  colors <- trellis:::readColors()[unique(df$read_type)]
  colors["splitread"] <- "black"
  nms <- names(trellis:::readColors())
  df$read_type <- factor(df$read_type, levels = nms)
  region <- read_type <- tagid <- NULL
  df1 <- filter(df, region == levels(region)[1])
  df2 <- filter(df, region == levels(region)[2])
  limits <- axis_limits(df, basepairs)
  gene1 <- levels(df$region)[1]
  gene2 <- levels(df$region)[2]
  xlim1 <- limits[[gene1]]
  xlim2 <- limits[[gene2]]
  labs1 <- trellis:::axis_labels5p(df1, xlim1, num.ticks)
  labs2 <- trellis:::axis_labels3p(df2, xlim2, num.ticks)
  a <- ggplot(df1, aes(
    ymin = tagid - 0.2,
    ymax = tagid + 0.2,
    xmin = start,
    xmax = end,
    color = read_type,
    fill = read_type, group = tagid
  )) +
    geom_rect() +
    ylab(ylabel) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_x_continuous(
      breaks = labs1[["breaks"]],
      labels = labs1[["labels"]]
    ) +
    coord_cartesian(xlim = xlim1) +
    xlab("") +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 5)
    ) +
    guides(color = "none", fill = "none") +
    geom_vline(xintercept = df$junction_5p[1], linetype = "dashed") +
    ggtitle(paste0(df1$region[1], " (", df1$seqnames[1], ")"))
  if (df1$reverse[1]) {
    a <- a + scale_x_reverse(
      breaks = labs1[["breaks"]],
      labels = labs1[["labels"]]
    )
  }
  b <- ggplot(df2, aes(
    ymin = tagid - 0.2,
    ymax = tagid + 0.2,
    xmin = start,
    xmax = end,
    color = read_type,
    fill = read_type, group = tagid
  )) +
    geom_rect() +
    ylab("read pair index") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_x_continuous(
      breaks = labs2[["breaks"]],
      labels = labs2[["labels"]]
    ) +
    coord_cartesian(xlim = xlim2) +
    xlab("") +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(size = 5)
    ) +
    guides(color = "none", fill = "none") +
    geom_vline(xintercept = df$junction_3p[1], linetype = "dashed") +
    ylab("") +
    ggtitle(paste0(df2$region[1], " (", df2$seqnames[1], ")"))
  if (df2$reverse[1]) {
    b <- b + scale_x_reverse(
      breaks = labs2[["breaks"]],
      labels = labs2[["labels"]]
    )
  }
  d <- ggplot(df, aes(
    ymin = tagid - 0.2,
    ymax = tagid + 0.2,
    xmin = start,
    xmax = end,
    color = read_type,
    fill = read_type, group = tagid
  )) +
    geom_rect() +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(color = guide_legend(title = ""), fill = guide_legend(title = ""))
  legend.grob <- trellis:::peelLegend(d)[[2]]
  agrob <- ggplotGrob(a)
  bgrob <- ggplotGrob(b)
  bgrob$widths <- agrob$widths
  list(
    a = a,
    b = b,
    `5p` = agrob,
    `3p` = bgrob,
    legend = legend.grob
  )
}

#' Used in 06.1-figure3.rmd to create amplicon graphs
#'
#' This is a modification to functions with same name (but without '2' postfix)
#' in svplots package. Instead of returning just a list of grobs, the
#' modification returns both the grobs and the ggplot objects, allowing further
#' modification of the ggplot objects prior to creating grobs for the figure.
#' @export
ggRearrange2 <- function(df, ylab = "Read pair index",
                         basepairs = 400, num.ticks = 5) {
  . <- NULL
  grobs <- .ggRearrange2(df, ylabel = ylab, basepairs, num.ticks)
  widths <- c(0.5, 0.5) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  heights <- c(0.95, 0.05) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  mat <- matrix(c(
    1, 2,
    3, 3
  ), byrow = TRUE, ncol = 2, nrow = 2)
  agrob <- grobs[["5p"]]
  bgrob <- grobs[["3p"]]
  legend.grob <- grobs[["legend"]]
  gobj <- gridExtra::arrangeGrob(agrob, bgrob,
    legend.grob,
    layout_matrix = mat,
    widths = widths,
    heights = heights
  )
  list(
    arranged.grobs = gobj,
    grobs = grobs
  )
}

#' @export
draw_heatmap <- function(heatmap.components) {
  cluster_data <- heatmap.components$cluster_data
  ha_rows <- heatmap.components$ha_rows
  hegt <- Heatmap(cluster_data,
    col = colorRamp2(
      c(0, 0.25, 0.5),
      c("#00FFCC", "#FFFFFF", "#0099FF")
    ),
    show_column_names = FALSE,
    left_annotation = ha_rows,
    show_column_dend = FALSE,
    show_heatmap_legend = FALSE,
    row_title_rot = 0,
    clustering_distance_rows = "euclidean"
  )
  return(hegt)
}

#' @export
heatmap_setup <- function(methylation_se) {
  metadata <- as_tibble(colData(methylation_se))
  df <- metadata %>%
    mutate(t.n = ifelse(t.n == "T", "Tumor", "Normal")) %>%
    rename(`Tissue type` = t.n) %>%
    mutate(endo.muc = ifelse(grepl("end", diagnosis),
      "Endometrial",
      "Mucinous"
    ))
  colnames(df) <- Hmisc::capitalize(colnames(df))
  bvals <- beta(methylation_se)
  df$Diagnosis <- factor(df$Diagnosis,
    levels = c(
      "Uterine endometrial",
      "Ovarian endometrioid",
      "Ovarian mucinous",
      "Colorectal mucinous",
      "Pancreatic mucinous",
      "Stomach mucinous"
    )
  )
  df$Tissue <- str_extract(df$Diagnosis, ".*(?= )")
  dx.colors <- tumor_colors()
  dx.colors <- dx.colors[-match(c("Uterine endometrioid", "Ovarian endometrioid"), names(dx.colors))]
  names(dx.colors) <- str_extract(names(dx.colors), ".*(?= )")
  study.colors <- c(
    "JHU" = "#002d72",
    "TCGA" = "gray90"
  )
  tn.colors <- c(
    "Normal" = "gray",
    "Tumor" = "black"
  )
  tn.lgd <- Legend(
    labels = names(tn.colors),
    legend_gp = gpar(fill = tn.colors),
    labels_gp = gpar(fontsize = 18),
    title_gp = gpar(fontsize = 22),
    title = " Tumor/Normal"
  )
  dx.lgd <- Legend(
    labels = names(dx.colors),
    legend_gp = gpar(fill = dx.colors),
    labels_gp = gpar(fontsize = 18),
    title_gp = gpar(fontsize = 22),
    title = " Tissue type"
  )
  histology.colors <- c("orange3", "steel blue")
  names(histology.colors) <- c("Endometrial", "Mucinous")
  histology.lgd <- Legend(
    labels = names(histology.colors),
    legend_gp = gpar(fill = histology.colors),
    labels_gp = gpar(fontsize = 18),
    title_gp = gpar(fontsize = 22),
    title = " Histology"
  )
  col_fun <- colorRamp2(
    c(0, 0.25, 0.5),
    c("#00FFCC", "#FFFFFF", "#0099FF")
  )
  beta.lgd <- Legend(
    col_fun = col_fun, title = expression(beta),
    labels_gp = gpar(fontsize = 18),
    title_gp = gpar(fontsize = 22)
  )
  horiz.legends <- packLegend(dx.lgd, histology.lgd,
    direction = "horizontal"
  )
  vert.legends <- packLegend(tn.lgd, beta.lgd, direction = "vertical")
  df <- df[, c("Tissue", "Tissue type", "Study", "Endo.muc")]
  colnames(df) <- c(" Tissue type", " Tumor/Normal", " Study", " Histology")
  ha_rows <- rowAnnotation(
    df = df[, c(" Tissue type", " Tumor/Normal", " Histology")],
    col = list(
      ` Tissue type` = dx.colors,
      ` Tumor/Normal` = tn.colors,
      ` Histology` = histology.colors
    ),
    show_legend = FALSE,
    annotation_name_rot = 45,
    annotation_name_side = "top"
  )
  ix <- seq_len(nrow(methylation_se))
  cluster_data <- t(bvals[ix, ])
  rownames(cluster_data) <- NULL
  result <- list(
    cluster_data = cluster_data,
    ha_rows = ha_rows,
    vert.legends = vert.legends,
    horiz.legends = horiz.legends
  )
  result
}

#' Build an amplicon network graph
#'
#' Plots a chromosomal network graph from an amplicon graph object, using a
#' seeded random colour palette per chromosome.  Requires `trellis`, `svplots`,
#' `GGally`, `network`, and `sna` to be installed.
#'
#' @param ag An amplicon graph object (from the trellis package).
#' @param tx A GRanges of transcript annotations (e.g. from `load_tx()`).
#' @param max_size Maximum node size passed to `GGally::ggnet2`.
#' @return A ggplot, or `NULL` if the graph has too few nodes.
#' @export
ampliconGraph <- function(ag, tx, max_size = 5) {
    if (length(trellis::graph(ag)@nodes) <= 1) return(NULL)
    B <- svplots:::plot_amplicons(ag)
    if (is.null(B)) return(NULL)
    B1 <- methods::as(B, "graphAM")
    am <- B1@adjMat
    net <- network::network(am, directed = FALSE)
    chroms <- sapply(strsplit(colnames(am), ":"), "[", 1)
    ar <- trellis::ampliconRanges(ag)
    hits <- GenomicRanges::findOverlaps(ar, tx, maxgap = 5000)
    cancer.con <- split(tx$cancer_connection[S4Vectors::subjectHits(hits)],
                        S4Vectors::queryHits(hits))
    is.driver <- sapply(cancer.con, any)
    is.driver2 <- rep(FALSE, ncol(am))
    is.driver2[as.integer(names(is.driver))] <- is.driver
    set.seed(123)
    ix <- grep("gr(a|e)y", grDevices::colors(), invert = TRUE)
    palette <- setNames(sample(grDevices::colors()[ix], 23),
                        paste0("chr", c(1:22, "X")))
    net <- network::set.vertex.attribute(net, "chrom", chroms)
    net <- network::set.vertex.attribute(net, "driver", is.driver2)
    GGally::ggnet2(net,
                   color = "chrom",
                   palette = palette,
                   shape = "driver",
                   size = "degree",
                   max_size = max_size) +
        ggplot2::guides(size = "none", shape = "none",
                        color = ggplot2::guide_legend(title = "")) +
        ggplot2::ggtitle("")
}

#' Build a ggbio circos plot as a grob
#'
#' Constructs a circular genome plot from a named list of structural variant
#' tracks and converts it to a grid grob.  Requires `svplots`, `ggbio`, and
#' `GenomeInfoDb` to be installed.
#'
#' @param all.sv Named list with elements `segments`, `deletions`, `amplicons`,
#'   `rears` (as returned by subsetting the `circos_data.rds` lists).
#' @param size Text size for chromosome labels.
#' @param radius Radius for the chromosome label track.
#' @return A grob (gtable), or `NULL` if the ggbio plot cannot be built.
#' @export
#\'  Build a circos plot grob for one sample
#\'
#\' Constructs a ggbio circos plot from a named list of structural variant
#\' tracks and returns a \code{ggplotGrob}.  The rearrangement link circle is
#\' only added when the sample has at least one rearrangement, matching the
#\' original \code{ext-figure4-7.rmd} behaviour.
#\'
#\' @param all.sv  Named list with elements \code{segments} (GRanges),
#\'   \code{deletions}, \code{amplicons}, \code{rears}; passed to
#\'   \code{svplots::circosTracks()}.
#\' @param size    Font size for chromosome labels (default 5).
#\' @param radius  Outer radius for chromosome labels (default 39).
#\' @param legend  If \code{TRUE}, return a list with elements \code{grob}
#\'   (the circos grob) and \code{legend} (the chromosome colour legend).
#\'   Default \code{FALSE}: return only the grob.
#\' @return A \code{gtable} grob, or \code{NULL} if rendering fails.
#\'   If \code{legend = TRUE}, a list with elements \code{grob} and
#\'   \code{legend} (the legend may be \code{NULL} for samples without
#\'   rearrangements).
#\' @export
plot_circos_grob <- function(all.sv, size = 5, radius = 39, legend = FALSE) {
    set.seed(123)
    ix <- grep("gr(a|e)y", grDevices::colors(), invert = TRUE)
    ## seqnames from circosTracks() are bare integers ("1".."22", "X"), not
    ## "chr1".."chr22".  Name palette2 to match what ggbio actually receives.
    palette2 <- setNames(sample(grDevices::colors()[ix], 23),
                         c(as.character(1:22), "X"))
    tracks <- svplots::circosTracks(all.sv)
    r      <- tracks[["r"]]
    GenomeInfoDb::seqinfo(r) <- GenomeInfoDb::seqinfo(tracks[["hg"]])
    GenomeInfoDb::seqinfo(r$linked.to) <- GenomeInfoDb::seqinfo(r)
    cnvs   <- tracks[["cnvs"]]
    has_rearrangements <- length(r) > 0
    ## Base plot: chromosome labels, CN segments, ideogram
    A <- ggbio::ggbio(buffer = 0) +
        ggbio::circle(tracks[["hg"]], geom = "text",
                      ggplot2::aes(label = seqnames),
                      vjust = 0, size = size, radius = radius)
    ## Rearrangement links — only when present
    if (has_rearrangements) {
        A <- A +
            ggbio::circle(r, geom = "link", linked.to = "linked.to",
                          radius = 25, color = "steelblue")
        if (length(cnvs) > 0)
            A <- A +
                ggbio::circle(cnvs, geom = "rect",
                              ggplot2::aes(color = type),
                              fill = "transparent", radius = 28) +
                ggplot2::scale_color_manual(values = c("red", "steelblue"))
    }
    A <- A +
        ggbio::circle(tracks[["gr.cn"]], geom = "segment",
                      color = "black",
                      ggplot2::aes(y = cn),
                      grid = FALSE, size = 0.5, radius = 31) +
        ggbio::circle(tracks[["hg"]], geom = "ideo",
                      ggplot2::aes(fill = seqnames),
                      color = "gray",
                      radius = 35) +
        ggplot2::scale_fill_manual(values = palette2) +
        ggplot2::guides(fill = "none", color = "none")
    grob <- tryCatch(ggplot2::ggplotGrob(A@ggplot), error = function(e) NULL)
    if (!legend) return(grob)
    ## Extract chromosome legend from samples with rearrangements
    leg <- if (has_rearrangements && !is.null(grob)) {
        fig <- A@ggplot +
            ggplot2::guides(fill = ggplot2::guide_legend(title = "Chromosome"))
        tryCatch(cowplot::get_legend(fig), error = function(e) NULL)
    } else NULL
    list(grob = grob, legend = leg)
}

#' Blank white ggplot panel
#'
#' Returns a ggplot with an all-white panel — used as a placeholder when a
#' sample has no amplicon graph to display.
#'
#' @return A ggplot object.
#' @export
gg_blank <- function() {
    ggplot2::ggplot() +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white",
                                                     color = "white")
        )
}

#' Build figure 3 plot data (circos and amplicon graph grobs)
#'
#' Runs the per-sample loop that produces circos grobs (via
#' \code{plot_circos_grob}) and amplicon network graphs (via
#' \code{ampliconGraph}) for all samples in \code{circos_data}.
#' Requires \pkg{ggbio}, \pkg{trellis}, and \pkg{svplots} to be installed.
#'
#' @param circos_data Named list from \code{circos_data.rds}; must contain
#'   elements \code{ids}, \code{segments}, \code{deletions},
#'   \code{amplicon_graphs}, \code{rlist}.
#' @param segments2 Named list of GRanges from \code{segs_to_granges}.
#' @param tx Transcript GRanges from \code{load_tx}.
#' @return Named list with elements \code{circos_list} and \code{ag_figs},
#'   each a named list of grobs (or \code{NULL} entries for samples that failed).
#' @export
build_figure3_plot_data <- function(circos_data, segments2, tx) {
    ids <- circos_data[["ids"]]
    track.names <- c("segments", "deletions", "amplicons", "rears")
    circos_list <- ag_figs <- vector("list", length(ids))
    names(circos_list) <- names(ag_figs) <- ids
    for (i in seq_along(ids)) {
        all.sv <- setNames(
            list(segments2[[i]],
                 circos_data[["deletions"]][[i]],
                 circos_data[["amplicon_graphs"]][[i]],
                 circos_data[["rlist"]][[i]]),
            track.names
        )
        circos_list[[i]] <- plot_circos_grob(all.sv)
        if (is.null(circos_list[[i]])) next
        B <- ampliconGraph(circos_data[["amplicon_graphs"]][[i]], tx)
        ## Convert to grob here, stripping the chromosome colour legend,
        ## matching the original 06.1-figure3.rmd pipeline which called
        ## ggplotGrob(g + guides(color = "none")) before saving ag.figs.rds.
        ag_figs[[i]] <- if (is.null(B)) gg_blank() else {
            ggplot2::ggplotGrob(B + ggplot2::guides(color = "none"))
        }
    }
    list(circos_list = circos_list, ag_figs = ag_figs)
}

#' Build the driver-gene legend for figure 3
#'
#' Creates a named list of gtable legends (one per sample) colouring driver
#' genes by chromosome.  Requires \pkg{trellis} and \pkg{ggpubr}.
#'
#' @param aglist Named list of amplicon graph objects for the figure 3 samples.
#' @param lab_id Character vector of sample IDs (same names as \code{aglist}).
#' @return Named list of gtable legend objects.
#' @export
build_amplicon_legend <- function(aglist, lab_id) {
    set.seed(123)
    ix <- grep("gr(a|e)y", grDevices::colors(), invert = TRUE)
    palette <- setNames(sample(grDevices::colors()[ix], 23),
                        paste0("chr", c(1:22, "X")))
    amps <- purrr::map(aglist, trellis::ampliconRanges) %>%
        purrr::map(function(x) x[!is.na(x$driver)])
    amps2 <- BiocGenerics::unlist(GenomicRanges::GRangesList(amps)) %>%
        tibble::as_tibble()
    ids <- rep(names(amps), purrr::map_int(amps, length))
    pal <- palette[c("chr12", "chr8")]
    dat <- amps2 %>%
        dplyr::mutate(lab_id = ids) %>%
        dplyr::select(lab_id, driver) %>%
        dplyr::distinct() %>%
        dplyr::mutate(driver = c("KRAS", "MYC"), chrom = c("chr12", "chr8"))
    names(pal) <- as.character(setNames(dat$driver, dat$chrom)[names(pal)])
    split(dat, dat$lab_id) %>%
        purrr::map(function(g, pal) {
            pal2 <- pal[dat$driver]
            fig <- ggplot2::ggplot(g, ggplot2::aes(lab_id, driver)) +
                ggplot2::geom_point(shape = 17, size = 5,
                                    ggplot2::aes(color = driver)) +
                ggplot2::scale_color_manual(values = pal2) +
                ggplot2::facet_wrap(~lab_id) +
                ggplot2::guides(color = ggplot2::guide_legend(title = "")) +
                ggplot2::theme(legend.key = ggplot2::element_rect(fill = "white"))
            ggpubr::get_legend(fig)
        }, pal = pal)
}

#' Build the rearrangement grobs for figure 3
#'
#' Extracts the specific chr19-chr22 rearrangement from CGOV161T and the
#' chr10-chr10 rearrangement from CGOV172T and returns an assembled
#' \code{arrangeGrob}.  Requires \pkg{trellis} and \pkg{gridExtra}.
#'
#' @param rlist Named list of rearrangement objects (from \code{circos_data}).
#' @return A gtable produced by \code{gridExtra::arrangeGrob}.
#' @export
build_figure3_rearrangements <- function(rlist) {
    update_galignment <- function(r) {
        for (i in seq_along(r)) {
            rr <- r[[i]]
            rr@improper <- BiocGenerics::updateObject(trellis::improper(rr))
            r[[i]] <- rr
        }
        r
    }
    rlist2 <- purrr::map(rlist[c("CGOV161T", "CGOV172T")], update_galignment)

    r.161t <- rlist2[["CGOV161T"]]
    lb <- trellis::linkedBins(r.161t)
    is_19.22 <- as.character(GenomicRanges::seqnames(lb)) == "chr19" &
        as.character(GenomicRanges::seqnames(trellis::linkedTo(lb))) == "chr22"
    r161.grob <- r.161t[is_19.22][1] %>%
        trellis::fiveTo3List() %>% `[[`(1) %>%
        trellis::rearDataFrame(build = "hg18") %>%
        ggRearrange2() %>% update.reargraph()

    r.172t <- rlist2[["CGOV172T"]]
    lb <- trellis::linkedBins(r.172t)
    is_10.10 <- as.character(GenomicRanges::seqnames(lb)) == "chr10" &
        as.character(GenomicRanges::seqnames(trellis::linkedTo(lb))) == "chr10"
    r172.grob <- r.172t[is_10.10][3] %>%
        trellis::fiveTo3List() %>% `[[`(1) %>%
        trellis::rearDataFrame(build = "hg18") %>%
        ggRearrange2() %>% update.reargraph()

    gridExtra::arrangeGrob(
        grobs   = list(grid::nullGrob(), r161.grob, grid::nullGrob(), r172.grob),
        heights = c(0.1, 1, 0.2, 1)
    )
}

#' Filter a RearrangementList to a specific breakpoint pair
#'
#' Searches both 5'→3' orientations so the caller does not need to know
#' which end is \code{linkedBins} and which is \code{linkedTo}.  Exactly
#' one matching rearrangement is expected; the function stops if zero or
#' more than one is found.
#'
#' After filtering, \code{fiveTo3List()} is called on the single matching
#' rearrangement and element \code{element} is passed to
#' \code{rearDataFrame(build = "hg18")} to produce the data frame consumed
#' by \code{\link{ggRearrange2}}.
#'
#' @param rlist   A \code{RearrangementList} (one sample).
#' @param chrom   Chromosome name as it appears in \code{seqlevels(rlist)},
#'   e.g. \code{"chr5"}.
#' @param pos1    Integer genomic position of the first breakpoint.
#' @param pos2    Integer genomic position of the second breakpoint.
#' @param maxgap  Maximum distance (bp) between a query position and a
#'   breakpoint to count as a match.  Default \code{5000}.
#' @param element Which element of \code{fiveTo3List()} to pass to
#'   \code{rearDataFrame()}.  Default \code{1L}; set to \code{2L} when
#'   manual inspection of the raw rearrangement indicates the second
#'   haplotype is the relevant one (as for CGOV161T chr5 breakpoint).
#'
#' @return A data frame returned by \code{trellis::rearDataFrame()}, ready
#'   for \code{\link{ggRearrange2}}.
#' @export
filter_rearrangement <- function(rlist, chrom, pos1, pos2,
                                 maxgap = 5000L, element = 1L) {
    si <- GenomicRanges::seqinfo(trellis::linkedBins(rlist))
    g1 <- GenomicRanges::GRanges(chrom,
                                  IRanges::IRanges(pos1, width = 1L),
                                  seqinfo = si)
    g2 <- GenomicRanges::GRanges(chrom,
                                  IRanges::IRanges(pos2, width = 1L),
                                  seqinfo = si)
    ## Search both orientations: linkedBins~pos1/linkedTo~pos2 and vice versa
    idxA <- which(
        IRanges::overlapsAny(trellis::linkedBins(rlist), g1, maxgap = maxgap) &
        IRanges::overlapsAny(trellis::linkedTo(rlist),   g2, maxgap = maxgap))
    idxB <- which(
        IRanges::overlapsAny(trellis::linkedTo(rlist),   g1, maxgap = maxgap) &
        IRanges::overlapsAny(trellis::linkedBins(rlist), g2, maxgap = maxgap))
    idx <- unique(c(idxA, idxB))
    if (length(idx) == 0L)
        stop("filter_rearrangement: no rearrangement found near ",
             chrom, ":", pos1, "-", pos2)
    if (length(idx) > 1L)
        stop("filter_rearrangement: ", length(idx),
             " rearrangements found near ",
             chrom, ":", pos1, "-", pos2,
             "; expected exactly 1. Consider reducing maxgap.")
    trellis::rearDataFrame(
        trellis::fiveTo3List(rlist[idx])[[element]],
        build = "hg18"
    )
}
