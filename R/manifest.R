subject_id <- function(x) {
  x %>%
    str_replace_all("_Ex", "") %>%
    str_replace_all("_WGS", "") %>%
    str_replace_all("T$", "") %>%
    str_replace_all("T_hg18_A", "") %>%
    str_replace_all("_WGS", "") %>%
    str_replace_all("([0-9])(T_[12])", "\\1") %>%
    str_replace_all("([0-9])(T[12])", "\\1") %>%
    str_replace_all(".mkdup", "") %>%
    str_replace_all("^[tn]_", "") %>%
    str_replace_all(".bam$", "") %>%
    str_replace_all("_eland", "") %>%
    str_replace_all(".final$", "") %>%
    str_replace_all("[TN]$", "") %>%
    str_replace_all("_Rep$", "") %>%
    str_replace_all("_Rpt$", "") %>%
    str_replace_all("T[be]$", "")
}

find_bams <- function(pgdx_id, bamfiles) {
  full.id <- pgdx_id
  abbrv.id <- full.id %>%
    strsplit("_") %>%
    "[["(1) %>%
    "["(1)
  index <- grep(abbrv.id, bamfiles)
  if (length(index) == 1) {
    bamfile <- bamfiles[index]
    return(bamfile)
  }
  lastchar <- substr(
    abbrv.id, nchar(abbrv.id),
    nchar(abbrv.id)
  ) %>%
    tolower()
  nm.split <- strsplit(abbrv.id, "[NT]")[[1]]
  abbrv.id2 <- nm.split[1]
  index <- grep(abbrv.id2, bamfiles)
  if (length(index) == 1) {
    bamfile <- bamfiles[index]
    return(bamfile)
  }
  if (length(index) == 0) {
    return(NA)
  }
  is_normal <- grepl("N", abbrv.id)
  tmp <- bamfiles[index]
  if (is_normal) {
    x <- paste0("n_", abbrv.id2)
    bamfile <- tmp[grep(x, tmp)]
  } else {
    bamfile <- tmp[grep("t_", tmp)]
  }
  if (length(bamfile) == 1) {
    return(bamfile)
  }
  if (length(bamfile) > 1) {
    bamfile <- bamfile[grep(abbrv.id, bamfile)]
  }
  if (length(bamfile) == 1) {
    return(bamfile)
  }
  return(NA)
}

find_bam2 <- function(pgdx_id, bamfile) {
  ix <- grep(pgdx_id, bamfile)
  if (length(ix) == 1) {
    return(bamfile[ix])
  }
  abbrv <- strsplit(pgdx_id, "_")[[1]][1]
  ix <- grep(abbrv, bamfile)
  if (length(ix) == 1) {
    return(bamfile[ix])
  }
  abbrv2 <- str_replace_all(abbrv, "[TN]$", "")
  ix <- grep(abbrv2, bamfile)
  if (length(ix) == 1) {
    return(bamfile[ix])
  }
  if (length(ix) == 2) {
    bams <- bamfile[ix]
    is_normal <- substr(abbrv, nchar(abbrv), nchar(abbrv)) == "N"
    if (is_normal) {
      ix <- grep("n_", basename(bams))
    } else {
      ix <- grep("t_", basename(bams))
    }
    bamfile <- bams[ix]
    return(bamfile)
  }
  NA
}

grep_bamfile <- function(i, manifest, tofix) {
  id <- tofix$prev_id[i]
  index <- grep(id, manifest$stripped_name)
  if (length(index) == 0) {
    return(manifest[index, ])
  }
  if (length(index) > 1) stop("Multiple BAM matches for id: ", id)
  manifest <- manifest[index, ] %>%
    mutate(prev_id = id) %>%
    ungroup() %>%
    select(lab_id, prev_id)
  return(manifest)
}

#' Make repairs to lab ids
#'
#' @param tofix:  a single column labeled lab_id
#' @param manifest: patient manifest
#' @return a two column tibble with lab_id and prev_id
repair_lab_id <- function(tofix, manifest, strip_Ex = FALSE) {
  alt0 <- filter(tofix, lab_id %in% manifest$lab_id)
  alt1 <- filter(tofix, !lab_id %in% manifest$lab_id)
  if (!strip_Ex) {
    tofix <- tibble(prev_id = unique(alt1$lab_id))
  } else {
    tofix <- tibble(prev_id = unique(alt1$lab_id)) %>%
      mutate(prev_id = str_replace_all(prev_id, "_Ex$", ""))
  }
  stripped.manifest <- manifest %>%
    mutate(
      x = str_replace_all(
        basename(bam_local),
        ".bam", ""
      ),
      x = str_replace_all(x, ".clean", ""),
      x = str_replace_all(x, ".mkdup", ""),
      x = str_replace_all(x, ".fxmt", "")
    ) %>%
    rename(stripped_name = x) %>%
    mutate(bam_local = basename(bam_local)) %>%
    select(
      subject_id, lab_id, stripped_name,
      bam_local, tumor.normal
    ) %>%
    filter(tumor.normal == "tumor")
  possible_matches <- seq_len(nrow(tofix)) %>%
    map_dfr(grep_bamfile,
      manifest = stripped.manifest,
      tofix = tofix
    )
  corrected <- tofix %>%
    left_join(possible_matches, by = "prev_id") %>%
    select(lab_id, prev_id)
  alt1.updated <- tofix %>%
    left_join(corrected, by = "prev_id") %>%
    select(lab_id, prev_id)
  if (strip_Ex) {
    alt1.updated$prev_id <- paste0(alt1.updated$prev_id, "_Ex")
  }
  alt0$prev_id <- alt0$lab_id
  alt3 <- bind_rows(alt0, alt1.updated)
  alt3
}

manifest_tumors <- function(manifest, tumor.levels) {
  manifest <- manifest %>%
    ungroup() %>%
    mutate(
      tumor_type = Hmisc::capitalize(tumor_type),
      tumor_type = case_when(
        tumor_type == "Colorectal" ~ "Colorectal mucinous",
        tumor_type == "Pancreas" ~ "Pancreas mucinous",
        tumor_type == "Stomach" ~ "Stomach mucinous",
        TRUE ~ tumor_type
      )
    ) %>%
    filter(
      tumor_type %in% tumor.levels,
      tumor.normal == "tumor"
    )
  manifest
}

#' Standardize pgdx identifiers
#'
#' @export
#' @param dat:  a tibble with pgdx_id and subject_id field
standardize_pgdx <- function(dat) {
  two.ids <- dat %>%
    mutate(two.ids = grepl("_Ex_PGDX", pgdx_id)) %>%
    filter(two.ids) %>%
    mutate(orig_id = pgdx_id) %>%
    mutate(
      pgdx_id2 = str_replace_all(pgdx_id, "^[tn]_", ""),
      pgdx_id2 = str_replace_all(pgdx_id2, "_Ex.mkdup.bam", ""),
      pgdx_id2 = str_replace_all(pgdx_id2, "_Ex", "")
    ) %>%
    mutate(
      id1 = sapply(strsplit(pgdx_id2, "_"), "[", 1),
      id2 = sapply(strsplit(pgdx_id2, "_"), "[", 2)
    ) %>%
    mutate(pgdx_id = ifelse(tumor.normal == "tumor", id2, id1)) %>%
    select(-c(orig_id, id1, id2, pgdx_id2))
  dat2 <- filter(dat, !grepl("_Ex_PGDX", pgdx_id)) %>%
    bind_rows(two.ids) %>%
    mutate(
      pgdx_id = str_replace_all(pgdx_id, "_WGS_Ex", ""),
      pgdx_id = str_replace_all(pgdx_id, "_$", ""),
      pgdx_id = str_replace_all(pgdx_id, "_Ex$", ""),
      pgdx_id = str_replace_all(pgdx_id, "_Ex_hg19", "")
    )
  temp <- filter(dat2, !is.na(pgdx_id)) %>%
    mutate(
      pgdx_id = str_replace_all(pgdx_id, "^[tn]_", ""),
      pgdx_id = str_replace_all(pgdx_id, ".mkdup.bam", ""),
      pgdx_id = str_replace_all(pgdx_id, "_Ex", ""),
      pgdx_id = paste0(pgdx_id, "_", platform),
      pgdx_id = str_replace_all(pgdx_id, "_WES", "_Ex")
    )
  dat3 <- filter(dat2, is.na(pgdx_id)) %>%
    bind_rows(temp)
  return(dat3)
}

#' Standardize pgdx identifiers
#'
#' @param dat:  a tibble with pgdx_id and subject_id field
standardize_pgdx2 <- function(dat) {
  two.ids <- dat %>%
    mutate(two.ids = grepl("_Ex_PGDX", pgdx_id)) %>%
    filter(two.ids) %>%
    mutate(orig_id = pgdx_id) %>%
    mutate(
      pgdx_id2 = str_replace_all(pgdx_id, "^[tn]_", ""),
      pgdx_id2 = str_replace_all(pgdx_id2, "_Ex.mkdup.bam", ""),
      pgdx_id2 = str_replace_all(pgdx_id2, "_Ex", "")
    ) %>%
    mutate(
      id1 = sapply(strsplit(pgdx_id2, "_"), "[", 1),
      id2 = sapply(strsplit(pgdx_id2, "_"), "[", 2)
    ) %>%
    mutate(pgdx_id = ifelse(tumor.normal == "tumor", id2, id1)) %>%
    select(-c(orig_id, id1, id2, pgdx_id2))
  dat2 <- filter(dat, !grepl("_Ex_PGDX", pgdx_id)) %>%
    bind_rows(two.ids) %>%
    mutate(
      pgdx_id = str_replace_all(pgdx_id, "_WGS_Ex", ""),
      pgdx_id = str_replace_all(pgdx_id, "_$", ""),
      pgdx_id = str_replace_all(pgdx_id, "_Ex$", ""),
      pgdx_id = str_replace_all(pgdx_id, "_Ex_hg19", ""),
      pgdx_id = str_replace_all(pgdx_id, " $", "")
    ) %>%
    select(-two.ids) %>%
    unite(pgdx_id2, c(pgdx_id, platform), sep = "_") %>%
    mutate(pgdx_id = str_replace_all(pgdx_id2, "WES$", "Ex")) %>%
    select(-pgdx_id2)
  return(dat2)
}

subject_id2 <- function(manifest) {
  one.to.many <- select(manifest, subject_id, genotype_id, platform) %>%
    group_by(genotype_id, platform) %>%
    summarize(
      one_to_many = length(unique(subject_id)) > 1,
      subject_ids = paste(unique(subject_id), collapse = ","),
      .groups = "drop"
    ) %>%
    filter(one_to_many)
  discordant <- unique(unlist(strsplit(one.to.many$subject_ids, ",")))
  manifest$subject_id2 <- ifelse(manifest$subject_id %in% discordant,
    manifest$genotype_id,
    manifest$subject_id
  )
  manifest
}

select_manifest_columns <- function(manifest) {
  manifest3 <- select(
    manifest, pgdx_id, lab_id,
    subject_id, bamfile, tumor.normal,
    bam_local,
    size, platform, tumor_type,
    genotype_id, subject_id, subject_id2
  )
  manifest3
}

update_manifest_ids <- function(manifest) {
  manifest2 <- manifest %>%
    mutate(
      temp = subject_id,
      subject_id = subject_id2,
      subject_id2 = temp
    ) %>%
    select(-temp)
}

discordant_tumor_type <- function(manifest4) {
  normals <- filter(manifest4, tumor.normal == "normal") %>%
    arrange(subject_id)
  tumors <- filter(manifest4, tumor.normal == "tumor") %>%
    arrange(subject_id)
  discord_label <- tumors$discordant_tumor_type
  names(discord_label) <- tumors$subject_id
  normals$discordant_tumor_type <- discord_label[normals$subject_id]
  manifest4 <- bind_rows(normals, tumors) %>% arrange(lab_id)
  manifest4
}

remove_any_duplicates <- function(manifest4) {
  dups <- manifest4$lab_id[duplicated(manifest4$lab_id)]
  manifest5 <- manifest4 %>%
    unite(uid, c(lab_id, platform),
      sep = "_",
      remove = FALSE
    ) %>%
    mutate(lab_id = ifelse(lab_id %in% dups,
      uid,
      lab_id
    )) %>%
    select(-uid) %>%
    filter(!duplicated(lab_id)) %>%
    distinct()
  manifest5
}

#' @export
read_facets <- function(file) {
  facets <- read_tsv(file, show_col_types = FALSE) %>%
    select(Sample) %>%
    distinct()
  facets
}

join_facets_to_manifest1 <- function(facets, manifest.list) {
  rename <- dplyr::rename
  manifest6 <- man(manifest.list)
  facets.matched <- inner_join(facets, key(manifest.list),
    by = c("Sample" = "pgdx_id")
  ) %>%
    rename(facet_id = Sample)
  facets.nomatch <- filter(facets, !Sample %in% manifest6$pgdx_id) %>%
    mutate(
      tmpid = str_replace_all(Sample, "t_", ""),
      tmpid = str_replace_all(tmpid, "_eland", ""),
      tmpid = str_replace_all(tmpid, "_hg18_", ""),
      tmpid = str_replace_all(tmpid, "_ExA", "_Ex")
    )
  facets.matched2 <- inner_join(facets.nomatch, key(manifest.list),
    by = c("tmpid" = "pgdx_id")
  ) %>%
    rename(facet_id = Sample)
  facets.nomatch <- filter(facets.nomatch, !tmpid %in% manifest6$pgdx_id)
  uid <- tibble(
    uid = unique(facets.nomatch$tmpid),
    uid2 = NA
  )
  for (i in seq_len(nrow(uid))) {
    id <- uid$uid[i]
    stripid <- strsplit(id, "_Ex_")
    if (length(stripid[[1]]) == 2) {
      uid$uid2[i] <- stripid[[1]][2]
      if (!uid$uid2[i] %in% manifest6$pgdx_id) stop()
      next()
    }
    stripid <- strsplit(id, "_")[[1]][1]
    ix <- grep(stripid, manifest6$pgdx_id)
    stripid <- manifest6$pgdx_id[ix]
    if (!stripid %in% manifest6$pgdx_id) stop()
    uid$uid2[i] <- stripid
  }
  facets.nomatch2 <- left_join(facets.nomatch, uid,
    by = c("tmpid" = "uid")
  )
  facets.matched3 <- inner_join(facets.nomatch2, key(manifest.list),
    by = c("uid2" = "pgdx_id")
  ) %>%
    rename(facet_id = Sample)
  facets2 <- bind_rows(
    facets.matched,
    facets.matched2,
    facets.matched3
  ) %>%
    select(facet_id, subject_id, lab_id)
  manifest7 <- left_join(manifest6, facets2, by = c("subject_id", "lab_id"))
  ix <- grep("177T", manifest7$lab_id)
  A <- manifest7[ix, ]
  manifest7 <- manifest7[-ix[2], ]
}

join_facets_to_manifest2 <- function(manifest7, directory.listing) {
  rename <- dplyr::rename
  facets <- tibble(Sample = directory.listing) %>%
    mutate(Sample = basename(Sample))
  man <- select(manifest7, subject_id, lab_id, pgdx_id)
  matched0 <- inner_join(facets, man,
    by = c("Sample" = "lab_id")
  ) %>%
    mutate(facet_id = Sample) %>%
    rename(lab_id = Sample)
  tmp <- filter(facets, !Sample %in% man$lab_id)
  matched1 <- inner_join(tmp, man,
    by = c("Sample" = "pgdx_id")
  ) %>%
    mutate(
      pgdx_id = Sample,
      facet_id = pgdx_id
    ) %>%
    select(-Sample)
  notmatched <- filter(tmp, !Sample %in% matched1$pgdx_id)
  tmp <- select(notmatched, Sample) %>%
    distinct() %>%
    mutate(
      uid = NA,
      lab_id = NA
    )
  for (i in seq_len(nrow(tmp))) {
    id <- tmp$Sample[i]
    if (id == "CGOV359T") {
      tmp$lab_id[i] <- id
      next()
    }
    if (id == "CGOV482") {
      id <- "CGOV482T"
      tmp$lab_id[i] <- tmp$uid[i] <- id
      next()
    }
    if (grepl("^LP", id)) {
      id2 <- paste0(id, "_WGS") %>%
        str_replace("LP6", "LP")
      ix <- match(id2, man$pgdx_id)
      tmp$uid[i] <- man$pgdx_id[ix]
      tmp$lab_id[i] <- man$lab_id[ix]
      next()
    }
    id2 <- str_replace(id, "_Ex", "")
    ix <- match(id2, man$pgdx_id)
    tmp$uid[i] <- man$pgdx_id[ix]
    tmp$lab_id[i] <- man$lab_id[ix]
    next()
  }
  tmp3 <- filter(tmp, !is.na(uid))
  notmatched2 <- filter(notmatched, Sample %in% tmp3$Sample)
  matched2 <- left_join(notmatched2, tmp3, by = "Sample") %>%
    rename(facet_id = Sample) %>%
    left_join(man, by = "lab_id") %>%
    select(-uid)
  tmp2 <- filter(tmp, is.na(uid)) %>%
    select(-uid) %>%
    rename(old_directory = Sample) %>%
    mutate(
      facet_id = NA,
      lab_id = NA
    )
  man <- select(manifest7, subject_id, lab_id, pgdx_id, subject_id2, tumor.normal) %>%
    filter(tumor.normal == "tumor")
  for (i in seq_len(nrow(tmp2))) {
    id <- tmp2$old_directory[i]
    ix <- grep(id, man$subject_id2)
    tmp2$facet_id[i] <- man$subject_id[ix]
    tmp2$lab_id[i] <- man$lab_id[ix]
  }
  matched3 <- tmp2 %>%
    mutate(subject_id = facet_id, pgdx_id = NA) %>%
    select(-old_directory)
  m <- bind_rows(
    matched0,
    matched1,
    matched2,
    matched3
  )
  manifest.notwgs <- filter(manifest7, is.na(facet_id))
  manifest.wes <- filter(manifest7, !is.na(facet_id))
  manifest.notwgs2 <- left_join(
    select(
      manifest.notwgs,
      -facet_id
    ),
    select(m, lab_id, facet_id),
    by = "lab_id"
  )
  manifest8 <- bind_rows(manifest.wes, manifest.notwgs2) %>%
    arrange(subject_id)
  manifest8
}

check_ids <- function(clinical, manifest) {
  stopifnot(all(clinical$lab_id %in% manifest$lab_id))
  stopifnot(all(clinical$subject_id %in% manifest$subject_id))
  stopifnot(all(manifest$subject_id %in% clinical$subject_id))
  tmp <- filter(manifest, lab_id == "CGOV177T_2")
  should.be.true <- nrow(tmp) == 1
  should.be.true
}

clean_manifest <- function(manifest2, cdat2) {
  manifest3 <- select_manifest_columns(manifest2)
  manifest4 <- update_manifest_ids(manifest3) %>%
    select(-tumor_type) %>%
    left_join(
      select(
        cdat2, lab_id, tumor_type,
        discordant_tumor_type
      ),
      by = "lab_id"
    )
  manifest4 <- discordant_tumor_type(manifest4)
  manifest5 <- remove_any_duplicates(manifest4)
  stopifnot(!any(duplicated(manifest5$lab_id)))
  manifest6 <- filter(manifest5, tumor.normal == "tumor") %>%
    select(subject_id, lab_id, pgdx_id)
  list(key = manifest6, manifest = manifest5)
}

join_with_facets <- function(manifest6, facets, directory.listing) {
  manifest7 <- join_facets_to_manifest1(facets, manifest6)
  manifest8 <- join_facets_to_manifest2(manifest7, directory.listing)
  manifest8
}

filter_discordant_tumors <- function(manifest8) {
  manifest <- filter(manifest8, !discordant_tumor_type) %>%
    filter(subject_id != "CGOV359") %>%
    mutate(lab_id = ifelse(lab_id == "CGOV151Tb_WES", "CGOV151Tb", lab_id))
  manifest
}

subset_by_tumors <- function(manifest, tumor.levels) {
  manifest <- manifest %>%
    cancer_names() %>%
    filter(tumor %in% tumor.levels) %>%
    filter(tumor.normal == "tumor")
  manifest
}

subset_endo <- function(manifest) {
  tumor.levels <- c("Ovarian endometrioid", "Uterine endometrioid")
  manifest2 <- subset_by_tumors(manifest, tumor.levels)
  manifest2
}

subset_gi <- function(manifest) {
  tumor.levels <- c("Colorectal mucinous", "Stomach mucinous", "Pancreas mucinous")
  manifest <- manifest %>%
    cancer_names() %>%
    filter(tumor %in% tumor.levels) %>%
    filter(tumor.normal == "tumor")
  manifest
}

subset_mucinous <- function(manifest) {
  tumor.levels <- c("Ovarian mucinous", "Colorectal mucinous")
  manifest <- manifest %>%
    cancer_names() %>%
    filter(tumor %in% tumor.levels) %>%
    filter(tumor.normal == "tumor")
  manifest
}

read_facets2 <- function(...) {
  dots <- list(...)
  facets_purity <- lapply(dots, read.delim) %>%
    bind_rows() %>%
    select(facet_id = Sample, purity = Purity, genotype_id = Genotype.ID) %>%
    mutate(is_na_purity = is.na(purity))
  return(facets_purity)
}

add_facets_purity <- function(manifest9, facets_purity) {
  manifest10 <- manifest9[!is.na(manifest9$facet_id), ]
  manifest11 <- left_join(manifest10, select(facets_purity, -genotype_id), by = "facet_id") %>%
    left_join(select(facets_purity, -facet_id), by = join_by("facet_id" == "genotype_id")) %>%
    mutate(
      purity = case_when(is.na(purity.x) & is.na(purity.y) ~ NA_real_,
        is.na(purity.x) & !is.na(purity.y) ~ purity.y,
        !is.na(purity.x) & is.na(purity.y) ~ purity.x,
        .default = NA_real_
      ),
      is_na_purity = case_when(is.na(is_na_purity.x) & is.na(is_na_purity.y) ~ NA,
        is.na(is_na_purity.x) & !is.na(is_na_purity.y) ~ is_na_purity.y,
        !is.na(is_na_purity.x) & is.na(is_na_purity.y) ~ is_na_purity.x,
        .default = NA
      )
    ) %>%
    replace_na(list(is_na_purity = FALSE)) %>%
    select(-c(purity.x, purity.y, is_na_purity.x, is_na_purity.y))
  manifest <- left_join(manifest9, select(manifest11, lab_id, purity, is_na_purity))
  return(manifest)
}

#' Drop all the samples with purity below threshold
#'
#' @export
purity_filter <- function(manifest, threshold = 0.2) {
  m <- manifest %>%
    filter(is.na(is_na_purity) | purity > threshold)
  m
}

#' @export
get_manifest <- function() {
  data(manifest, package = "ovarian.subtypes", envir = environment())
  manifest
}
