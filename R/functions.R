update.reargraph <- function(object, show_legend=TRUE, size=15){
    grobs <- object$grobs
    a <- grobs[["a"]] +
        theme(plot.title=element_text(size=size),
              panel.background=element_rect(fill="gray95"))
    b <- grobs[["b"]] +
        theme(plot.title=element_text(size=size),
              panel.background=element_rect(fill="gray95"))
    ag <- ggplotGrob(a)
    bg <- ggplotGrob(b)
    bg$widths <- ag$widths
    widths <- c(0.5, 0.5) %>%
        "/"(sum(.)) %>%
        unit(., "npc")
    heights <- c(0.95, 0.05) %>%
        "/"(sum(.)) %>%
        unit(., "npc")
    mat <- matrix(c(1, 2,
                    3, 3), byrow=TRUE, ncol=2, nrow=2)
    ##agrob <- grobs[["5p"]]
    ##bgrob <- grobs[["3p"]]
    if(show_legend){
        legend.grob <- grobs[["legend"]]
    } else{
        legend.grob <- nullGrob()
        heights[2] <- unit(0, "npc")
    }
    gobj <- arrangeGrob(ag, bg,
                         legend.grob,
                         layout_matrix=mat,
                         widths=widths,
                         heights=heights)
    return(gobj)
}

jags_data <- function(gdat){
    gdat <- gdat %>%
        readRDS() %>%
        as_tibble %>%
        set_colnames(gsub(" ", "_", colnames(.))) %>%
        set_colnames(gsub("\\.", "_", colnames(.))) %>%
        set_colnames(tolower(colnames(.)))
    ngenes <- length(levels(gdat$gene_symbol))
    nid <- length(levels(gdat$internal_id))
    ntumors <- gdat %>%
        group_by(tumor_type) %>%
        summarize(n=length(unique(internal_id)))
    gdat2 <- gdat %>%
        mutate(mutation=ifelse(!is.na(mutation), 1L, 0L),
               methylation=ifelse(!is.na(methylation), 1L, 0L),
               fusion=ifelse(!is.na(fusion), 1L, 0L),
               copynumber=ifelse(!is.na(copynumber), 1L, 0L)) %>%
        mutate(total=mutation+methylation+fusion+copynumber,
               is_altered=ifelse(total > 0, 1L, 0L))%>%
        select(gene_symbol, internal_id, tumor_type, is_altered) %>%
        mutate(tumor_type=ifelse(tumor_type=="Ovarian endometrioid",
                                 "ovarian",
                                 "uterine"),
               internal_id=as.character(internal_id),
               gene_symbol=as.character(gene_symbol)) %>%
        group_by(internal_id, gene_symbol) %>%
        summarize(tumor_type=unique(tumor_type),
                  number_altered=sum(is_altered),
                  is_altered=as.integer(number_altered > 0),
                  .groups="drop")
    tumortype <- group_by(gdat2, internal_id) %>%
        summarize(tumor_type=unique(tumor_type),
                  .groups="drop")
    gdat3 <- gdat2 %>%
        select(gene_symbol, internal_id, is_altered) %>%
        spread(gene_symbol, is_altered) %>%
        left_join(tumortype, by="internal_id")
    X <- select(ungroup(gdat3), -internal_id) %>%
        select(-tumor_type) %>%
        as.matrix()
    Y <- ifelse(gdat3$tumor_type=="ovarian", 1, 0)
    ##
    ## exclude hypermutated samples
    ##
    is_hypermutator <- X[, "hypermutator"] == 1
    X <- X[, -match("hypermutator", colnames(X))]
    Y <- Y[ !is_hypermutator ]
    X <- X[ !is_hypermutator, ]
    ## require at least 3 variants
    X <- X [, colSums(X) >= 3]
    X <- cbind(1, X)
    colnames(X)[1] <- "intercept"
    list(X=X, Y=Y)
}


contingency_table <- function(x, Ns){
    if(sum(x$number_altered) == 0) return(NULL)
    m <- x %>%
        mutate(tumor_type=as.character(tumor_type)) %>%
        group_by(tumor_type, is_altered) %>%
        summarize(n=n(),
                  .groups="drop") %>%
        pivot_wider(names_from="is_altered", values_from="n") %>%
        set_colnames(c("tumor_type", "wt", "mt")) %>%
        mutate(n=mt + wt)
    m2 <- left_join(Ns, m, by=c("tumor_type", "n")) %>%
        select(tumor_type, mt, n)
    m2[is.na(m2)] <- 0
    m <- select(m2, -tumor_type) %>%
        as.matrix()
    rownames(m) <- m2$tumor_type
    m
}

stan_inputs <- function(ct){
    ##ct = contingency table
    y <- ct[, 1]
    n <- ct[, 2]
    J <- length(y)
    ##
    ## With 4 cancers, we need 3 dummy variables
    ##  - first cancer will be the reference
    X <- diag(4)
    X <- X[, -1]
    list(y=y, n=n, J=J, K=6, X=X) ##, X=0:3)
}

#' Create list of y, n, J, and x for stan
#'
#' @export
inputs_endo <- function(ct){
    ##ct = contingency table
    y <- ct[1:2, 1]
    n <- ct[1:2, 2]
    J <- length(y)
    ##
    ## With 2 cancers, need 1 dummy variable
    x <- c(0, 1)
    list(y=y, n=n, J=J, x=x) ##, X=0:3)
}

#' Create list of y, n, J, and x for stan
#'
#' @export
inputs_mucinous <- function(ct){
    ##ct = contingency table
    y <- ct[3:4, 1]
    n <- ct[3:4, 2]
    J <- length(y)
    ##
    ## With 2 cancers, need 1 dummy variable
    x <- c(0, 1)
    list(y=y, n=n, J=J, x=x)
}

sampling2 <- function(data, model, params, ...){
    sampling(model,
             data=data,
             iter=params$iter,
             thin=params$thin,
             chains=params$chains,
             warmup=params$warmup,
             control=params$control, ...)
}

slice_params <- function(x){
    nms <- rownames(x)
    x %>%
        as_tibble() %>%
        mutate(parameter=nms) %>%
        filter(grepl("^beta", parameter) |
               grepl("^theta", parameter)) %>%
        select(parameter, mean, se_mean, sd, `2.5%`,
               `5%`, `50%`, `95%`, `97.5%`, `n_eff`, Rhat)
}

clean_names <- function(x){
    x <- str_replace_all(x, " ", "_") %>%
        tolower()
    x
}

clean_colnames3 <- function(x){
    nms <- colnames(x)
    nms2 <- clean_names(nms)
    x2 <- set_colnames(x, nms2)
    x2
}

subject_id <- function(x){
    x %>%
        str_replace_all("_Ex", "") %>%
        str_replace_all("_WGS", "") %>%
        str_replace_all("T$", "") %>%
        str_replace_all("T_hg18_A", "") %>%
        str_replace_all("_WGS", "") %>%
        ## CGST1 should stay CGST1
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

subject_id2 <- function(x){
    x %>%
        str_replace_all("N$", "") %>%
        str_replace_all("T$", "") %>%
        ## ST1 should stay ST1
        str_replace_all("([0-9])(T[12])", "\\1")
}

find_bams <- function(pgdx_id, bamfiles){
    ##for(i in seq_len(nrow(missing.bamfile))){
    ##full.id <- missing.bamfile$pgdx_id[i]
    full.id <- pgdx_id
    abbrv.id <- full.id %>%
        strsplit("_") %>%
        "[["(1) %>%
        "["(1)
    index <- grep(abbrv.id, bamfiles)
    if(length(index) == 1){
        bamfile <- bamfiles[index]
        return(bamfile)
    }
    ##missing.bamfile$bamfile[i] <- target[index]
    ##      next()
    ##}
    ## remove suffix and search for a n_ or t_
    lastchar <- substr(abbrv.id, nchar(abbrv.id),
                       nchar(abbrv.id)) %>%
        tolower()
    nm.split <- strsplit(abbrv.id, "[NT]")[[1]]
    abbrv.id2 <- nm.split[1]
    ##numbers_only <- str_replace_all(abbrv.id, "[NT]^", "")
    ##abbrv.id2 <- paste0(lastchar, "_", numbers_only)
    index <- grep(abbrv.id2, bamfiles)
    if(length(index) == 1){
        bamfile <- bamfiles[index]
        return(bamfile)
    }
    if(length(index)==0) return(NA)
    ## more than one hit
    is_normal <- grepl("N", abbrv.id)
    tmp <- bamfiles[index]
    if(is_normal){
        x <- paste0("n_", abbrv.id2)
        bamfile <- tmp[grep(x, tmp)]
    } else {
        ##x <- paste0("t_", abbrv.id2)
        bamfile <- tmp[grep("t_", tmp)]
    }
    if(length(bamfile) == 1) return(bamfile)
    if(length(bamfile) > 1){
        bamfile <- bamfile[grep(abbrv.id, bamfile)]
    }
    if(length(bamfile)==1)  return(bamfile)
    return(NA)
}

find_bam2 <- function(pgdx_id, bamfile){
    ix <- grep(pgdx_id, bamfile)
    if(length(ix) == 1){
        return(bamfile[ix])
    }
    abbrv <- strsplit(pgdx_id, "_")[[1]][1]
    ix <- grep(abbrv, bamfile)
    if(length(ix) == 1){
        return(bamfile[ix])
    }
    abbrv2 <- str_replace_all(abbrv, "[TN]$", "")
    ix <- grep(abbrv2, bamfile)
    if(length(ix) == 1){
        return(bamfile[ix])
    }
    if(length(ix)==2){
        bams <- bamfile[ix]
        is_normal <- substr(abbrv, nchar(abbrv), nchar(abbrv)) == "N"
        if(is_normal){
            ix <- grep("n_", basename(bams))
        } else {
            ix <- grep("t_", basename(bams))
        }
        bamfile <- bams[ix]
        return(bamfile)
    }
    NA
}

clean_colnames <- function(x){
    rename <- dplyr::rename
    nms <- colnames(x)
    nms2 <- nms %>%
        tolower(.) %>%
        str_replace_all(., "(\\.)\\1+", "_") %>%
        str_replace_all(., "\\.", "_") %>%
        str_replace_all(., "_$", "")
    colnames(x) <- nms2
    x  <- x %>%
        rename(years=years_from_diagnosis,
               overall_survival=overall_survival_status_0_alive_1_dead,
               histology=histological_tumor_type) %>%
        mutate(histology=str_replace_all(histology, "ovarian", "Ovarian"))
    x
}

clean_colnames2 <- function(x){
    cnms <- colnames(x) %>%
        str_replace_all(" ", "_") %>%
        tolower()
    colnames(x) <- cnms
    x
}

format_number <- function(type, x){
    n <- x$n[x$tumor_type %in% type] %>%
        sum()
    Ns <- paste0("(n=", n, ")")
    labels <- paste(type, collapse=",")
    names(Ns) <- labels
    Ns
}

remove_author <- function(x) {
    ## identify empty author line
    i <- grep("^\\\\author\\{\\}$", x)
    ## be sure it is the one pandoc inserts
    if(length(i) != 0 && grepl('^\\\\date\\{', x[i+1])) x <- x[-i]

    ## default puts thanks on the title
    i <- grep("^\\\\title", x)
    line <- "\\title{Genomic landscapes of endometrioid and mucinous ovarian cancers and morphologically similar tumor types}"
    if(length(i) != 0 ) x[i] <- line
    ## put thanks on Velculescu and Scharpf
    i <- grep("Velculescu", x)
    line <- x[i]
    line <- stringr::str_replace(line, "Velculescu", "Velculescu \\\\thanks{To whom correspondence should be addressed: velculescu@jhmi.edu (V.E.V.) and rscharpf@jhu.edu (R.B.S.)}")
    x[i] <- line
    ##i <- grep("\dagger", x)
    ##line <- x[i]
    ##line <- str_replace(line, "\\\\dagger", "*")
    ##x[i] <- line
    x
}

grep_bamfile <- function(i, manifest, tofix){
    ##id <- tofix$stripped_name[i]
    id <- tofix$prev_id[i]
    index <- grep(id, manifest$stripped_name)
    if(length(index)==0) return(manifest[index, ])
    if(length(index) > 1) browser()
    manifest <- manifest[index, ] %>%
        mutate(prev_id=id) %>%
        ungroup() %>%
        select(lab_id, prev_id)
    return(manifest)
}

#' Make repairs to lab ids
#'
#' @param tofix:  a single column labeled lab_id
#' @param manifest: patient manifest
#' @return a two column tibble with lab_id and prev_id
repair_lab_id <- function(tofix, manifest, strip_Ex=FALSE){
    alt0 <- filter(tofix, lab_id %in% manifest$lab_id)
    alt1 <- filter(tofix, !lab_id %in% manifest$lab_id)
    if(!strip_Ex){
        tofix <- tibble(prev_id=unique(alt1$lab_id))
    } else{
        tofix <- tibble(prev_id=unique(alt1$lab_id)) %>%
            mutate(prev_id=str_replace_all(prev_id, "_Ex$", ""))
    }
    ##mutate(stripped_name=str_replace_all(id, "_WGS_Ex", ""),
    ##stripped_name=str_replace_all(stripped_name, "_WGS", ""))
    stripped.manifest <- manifest %>%
        mutate(x=str_replace_all(basename(bam_local),
                                 ".bam", ""),
               x=str_replace_all(x, ".clean", ""),
               x=str_replace_all(x, ".mkdup", ""),
               x=str_replace_all(x, ".fxmt", "")) %>%
        rename(stripped_name=x) %>%
        mutate(bam_local=basename(bam_local)) %>%
        select(subject_id, lab_id, stripped_name,
               bam_local, tumor.normal) %>%
        filter(tumor.normal=="tumor")
    possible_matches <- seq_len(nrow(tofix)) %>%
        map_dfr(grep_bamfile, manifest=stripped.manifest,
                tofix=tofix)
    corrected <- tofix %>%
        ##rename(prev_id=id) %>%
        left_join(possible_matches, by="prev_id") %>%
        ##rename(prev_id=id) %>%
        select(lab_id, prev_id)
    ##alt1.updated <- alt1 %>%
    alt1.updated <- tofix %>%
        ##rename(prev_id=lab_id) %>%
        left_join(corrected, by="prev_id") %>%
        select(lab_id, prev_id)
    if(strip_Ex){
        alt1.updated$prev_id <- paste0(alt1.updated$prev_id, "_Ex")
    }
    alt0$prev_id <- alt0$lab_id
    alt3 <- bind_rows(alt0, alt1.updated)
    alt3
}

tile_theme <- function(){
    theme(axis.title=element_blank(),
          strip.placement="outside",
          panel.grid=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.background=element_rect(fill="white"),
          legend.background=element_rect(fill="white"),
          axis.line = element_blank(),
          strip.text=element_text(size=9,
                                  hjust=0.5,
                                  vjust=0.5),
          strip.text.y.left=element_text(angle=0,
                                         size=11),
          strip.background=element_rect(fill="grey88",
                                        color = "grey88"),
          panel.spacing.x = unit(1,"lines"))
}


suppl5_varnames <- function(s5){
    orignames <- colnames(s5)
    varnames <- clean_names(orignames)
    varnames[c(8, 9, 13, 14, 15,
               16:19)] <- c("olap_snps",
                            "olap_het_snps",
                            "cnv_type",
                            "loh",
                            "focal",
                            "clin_gene",
                            "biol_gene",
                            "olap_gene",
                            "olap_tx")
    varnames
}

order_samples <- function(x, gene.levels){
    x2 <- select(x, lab_id, gene_symbol) %>%
        mutate(hypermutator=gene_symbol=="hypermutator",
               gene1.alt=gene_symbol==gene.levels[1],
               gene2.alt=gene_symbol==gene.levels[2],
               gene3.alt=gene_symbol==gene.levels[3],
               gene4.alt=gene_symbol==gene.levels[4],
               gene5.alt=gene_symbol==gene.levels[5]) %>%
        group_by(lab_id) %>%
        summarize(hypermut=any(hypermutator),
                  gene1=any(gene1.alt),
                  gene2=any(gene2.alt),
                  gene3=any(gene3.alt),
                  gene4=any(gene4.alt),
                  gene5=any(gene5.alt))
    x3 <- x2 %>%
        arrange(hypermut,
                -gene1,
                -gene2,
                -gene3,
                -gene4,
                -gene5)
    x3
}

manifest_tumors <- function(manifest, tumor.levels){
    ##manifest <- here("output", "03-manifest.rmd",
    ##                 "manifest.rds") %>%
    ##    readRDS() %>%
    manifest <- manifest %>%
        ungroup() %>%
        mutate(tumor_type=Hmisc::capitalize(tumor_type),
               tumor_type=case_when(tumor_type=="Colorectal"~"Colorectal mucinous",
                                    tumor_type=="Pancreas"~"Pancreas mucinous",
                                    tumor_type=="Stomach"~"Stomach mucinous",
                                    TRUE~tumor_type)) %>%
        filter(tumor_type %in% tumor.levels,
               tumor.normal=="tumor")
    manifest
}

read_pathways <- function(pathway.file, tumor.type){
    rename <- dplyr::rename
    pathways <- read_csv(pathway.file, show_col_types=FALSE) %>%
        rename(gene_symbol=gene.symbol) %>%
        filter(tumor_type==tumor.type) %>%
        select(-tumor_type) %>%
        mutate(pathway=str_replace_all(pathway,
                                       "TGFBR pathway",
                                       "TGFBR"),
               pathway=str_replace_all(pathway, "BRCA",
                                       "DNA repair"))
    if(tumor.type=="mucinous"){
        ix <- which(pathways$gene_symbol == "JAK1" & pathways$pathway == "Cell cycle")
        ix2 <- which(pathways$gene_symbol == "MED1-STAT5B" & pathways$pathway == "Other")
        pathways2 <- pathways[-c(ix, ix2), ]
    } else pathways2 <- pathways
    pathways2
}

read_integrated_data <- function(manifest, pathways){
    ##manifest <- read_manifest(tumor.levels)
    tumortypes <- select(manifest, lab_id, tumor_type) %>%
        ungroup() %>%
        distinct()
    idat <- here("output", "01-data_integration.rmd",
                 "integrated_data.rds") %>%
        readRDS() %>%
        mutate(gene_symbol=gene) %>%
        filter(lab_id %in% manifest$lab_id) %>%
        left_join(pathways, by="gene_symbol") %>%
        left_join(tumortypes, by="lab_id") %>%
        distinct()
    idat
}

axis.labels <- function(ord_in, signif.digits=1){
    exp_var <- 100 * ord_in$svd^2 / sum(ord_in$svd^2)
    axes <- paste0("LD", 1:2)
    axes <- paste0(axes, ' (', round(exp_var, signif.digits), '%)')
    axes
}

my.ggord.lda <- function(ord_in, grp_in = NULL,
                         axes = c('1', '2'), ...){
    obs <- data.frame(predict(ord_in)$x[, c("LD1", "LD2")]) %>%
        as_tibble() %>%
        mutate(lab="TCGA")
    obs$Groups <- as.character(grp_in)
    obs
}

my.ellipse <- function(obs){
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ellipse_pro <- 0.95
    ell <- plyr::ddply(obs, 'Groups', function(x) {
        if(nrow(x) <= 2) {
            return(NULL)
        }
        sigma <- var(cbind(x$LD1, x$LD2))
        mu <- c(mean(x$LD1), mean(x$LD2))
        ed <- sqrt(qchisq(ellipse_pro, df = 2))
        data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+')) %>%
            as_tibble()
    })
    names(ell)[2:3] <- names(obs)[1:2]
    ## get convex hull for ell object, this is a hack to make it work with geom_polygon
    . <- plyr::.
    ell <- plyr::ddply(ell, .(Groups), function(x) x[chull(x$LD1, x$LD2), ])
    ell <- as_tibble(ell)
}

project.cancer <- function(se.jhu, pc.tcga, ld.tcga, cancertype, use_pcs=1:5){
    se.jhu2 <- se.jhu[, se.jhu$diagnosis==cancertype]
    jhu.meth <- assays(se.jhu2)[[1]] %>%
        t()
    jhu.pcs <- predict(pc.tcga, newdata=jhu.meth) %>%
        as_tibble() %>%
        select(paste0("PC", use_pcs)) %>%
        mutate(dx=se.jhu2$diagnosis)
    ## 2. Using the LDA classifier, make predictions on the JHU study
    jhu.class.predictions <- predict(ld.tcga, newdata=jhu.pcs)
    jhu.x <- jhu.class.predictions$x[, c("LD1", "LD2")] %>%
        as_tibble() %>%
        mutate(dx=as.character(se.jhu2$diagnosis),
               tumor=factor(se.jhu2$tumor, c("Normal", "Tumor"))) %>%
        rename(Groups=dx,
               tumor.normal=tumor) %>%
        mutate(lab="JHU")
    jhu.x
}

#' Wrapper for principal component analysis of SummarizedExperiment object
#'
#' @export
mypca <- function(se, scale=FALSE, center=TRUE, rk){
    x <- t(assays(se)[[1]])
    prcomp(x, scale=scale, center=center, rank.=rk)
}

## temporary -- to delete
collect_sampleinfo <- function(manifest, path){
    ##samps <- samps[c(1:8,11:15,18:45)]
    ##tmp <- list.files(dbruhm)[c(1:8,11:15,18:45)]
    ids <- manifest %>%
        mutate(path=path) %>%
        mutate(path=file.path(path, ".temp", vendor_id),
               rearfile=file.path(path,
                                  paste0(vendor_id, ".rearlist-final.rds")),
               ampfile=file.path(path,
                                 paste0(".unfiltered_amplicons.rds")),
               delfile=file.path(path,
                                 paste0(".unfiltered_deletions.rds")),
               segfile=file.path(path,
                                 paste0(".segments.rds")))
    ids
}


get_manifest <- function(samps){
    ids <- tibble(vendor_id=basename(samps),
                  lab_id=c("CGOV353T","CGOV354T","CGOV358T","CGOV359T",
                           "CGOV362T","CGOV365T","CGOV369T",
                           "CGOV375T","CGOV291T","CGOV292T",
                           "CGOV293T","CGOV295T","CGOV296T","CGOV127T",
                           "CGOV131T","CGOV136T","CGOV138T",
                           "CGOV139T","CGOV140T","CGOV141T","CGOV142T",
                           "CGOV170T","CGOV172T","CGOV173T","CGOV174T",
                           "CGOV176T","CGOV155T","CGOV159T",
                           "CGOV159T_3","CGOV144T","CGOV145T","CGOV145T_1",
                           "CGOV147T","CGOV148T","CGOV154T",
                           "CGOV157T","CGOV160T","CGOV160T_1","CGOV160T_2",
                           "CGOV161T","CGOV162T"
                           ),
                  subtypes=c("ovarian endometrioid",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian mucinous",
                             "ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian mucinous","ovarian endometrioid",
                             "ovarian mucinous",
                             "ovarian mucinous","ovarian endometrioid",
                             "colorectal","colorectal",
                             "colorectal","ovarian mucinous",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian mucinous","ovarian mucinous",
                             "ovarian endometrioid",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian endometrioid",
                             "ovarian endometrioid","ovarian endometrioid",
                             "ovarian endometrioid"))
    ids
}

#' Standardize pgdx identifiers
#'
#' @export
#' @param dat:  a tibble with pgdx_id and subject_id field
standardize_pgdx <- function(dat){
    two.ids <- dat %>%
        mutate(two.ids=grepl("_Ex_PGDX", pgdx_id)) %>%
        filter(two.ids) %>%
        mutate(orig_id=pgdx_id) %>%
        mutate(pgdx_id2=str_replace_all(pgdx_id, "^[tn]_", ""),
               pgdx_id2=str_replace_all(pgdx_id2, "_Ex.mkdup.bam", ""),
               pgdx_id2=str_replace_all(pgdx_id2, "_Ex", "")) %>%
        mutate(id1=sapply(strsplit(pgdx_id2, "_"), "[", 1),
               id2=sapply(strsplit(pgdx_id2, "_"), "[", 2)) %>%
        mutate(pgdx_id=ifelse(tumor.normal=="tumor", id2, id1)) %>%
        select(-c(orig_id, id1, id2, pgdx_id2))
    dat2 <- filter(dat, !grepl("_Ex_PGDX", pgdx_id)) %>%
        bind_rows(two.ids) %>%
        mutate(pgdx_id=str_replace_all(pgdx_id, "_WGS_Ex", ""),
               pgdx_id=str_replace_all(pgdx_id, "_$", ""),
               pgdx_id=str_replace_all(pgdx_id, "_Ex$", ""),
               pgdx_id=str_replace_all(pgdx_id, "_Ex_hg19", ""))
    ## Remove n_ t_
    temp <- filter(dat2, !is.na(pgdx_id)) %>%
        mutate(pgdx_id=str_replace_all(pgdx_id, "^[tn]_", ""),
               pgdx_id=str_replace_all(pgdx_id, ".mkdup.bam", ""),
               pgdx_id=str_replace_all(pgdx_id, "_Ex", ""),
               pgdx_id=paste0(pgdx_id, "_", platform),
               pgdx_id=str_replace_all(pgdx_id, "_WES", "_Ex"))
    dat3 <- filter(dat2, is.na(pgdx_id)) %>%
        bind_rows(temp)
    return(dat3)
}

#' Standardize pgdx identifiers
#'
#' @param dat:  a tibble with pgdx_id and subject_id field
standardize_pgdx2 <- function(dat){
    two.ids <- dat %>%
        mutate(two.ids=grepl("_Ex_PGDX", pgdx_id)) %>%
        filter(two.ids) %>%
        mutate(orig_id=pgdx_id) %>%
        mutate(pgdx_id2=str_replace_all(pgdx_id, "^[tn]_", ""),
               pgdx_id2=str_replace_all(pgdx_id2, "_Ex.mkdup.bam", ""),
               pgdx_id2=str_replace_all(pgdx_id2, "_Ex", "")) %>%
        mutate(id1=sapply(strsplit(pgdx_id2, "_"), "[", 1),
               id2=sapply(strsplit(pgdx_id2, "_"), "[", 2)) %>%
        mutate(pgdx_id=ifelse(tumor.normal=="tumor", id2, id1)) %>%
        select(-c(orig_id, id1, id2, pgdx_id2))
    dat2 <- filter(dat, !grepl("_Ex_PGDX", pgdx_id)) %>%
        bind_rows(two.ids) %>%
        mutate(pgdx_id=str_replace_all(pgdx_id, "_WGS_Ex", ""),
               pgdx_id=str_replace_all(pgdx_id, "_$", ""),
               pgdx_id=str_replace_all(pgdx_id, "_Ex$", ""),
               pgdx_id=str_replace_all(pgdx_id, "_Ex_hg19", ""),
               pgdx_id=str_replace_all(pgdx_id, " $", "")) %>%
        select(-two.ids) %>%
        unite(pgdx_id2, c(pgdx_id, platform), sep="_") %>%
        mutate(pgdx_id=str_replace_all(pgdx_id2, "WES$", "Ex")) %>%
        select(-pgdx_id2)
    return(dat2)
}

#' @export
fig1_ggplatform <- function(i, x, base_size=15, colors, cancer){
    dat <- x %>%
        mutate(platform=case_when(platform=="Methylation"~"Me",
                                  ##platform=="Survival"~"Surv",
                                  TRUE~platform),
               platform=factor(platform,
                               levels=c("WGS", "WES",
                                        "Me"))) ##, "Surv")))
    orderby <- dat %>%
        group_by(subject_id) %>%
        ##group_by(genotype_id) %>%
        summarize(nplatform=length(unique(platform)),
                  tumor.normal=sum(tumor.normal=="normal,tumor"),
                  nwgs=sum(platform=="WGS")) %>%
        arrange(nplatform, tumor.normal, nwgs)
    dat2  <- dat %>%
        mutate(subject_id=factor(subject_id, orderby$subject_id))
    fig <- dat2 %>%
        ggplot(aes(platform, subject_id)) +
        geom_point(aes(fill=matched),
                   pch=21, size=4) +
        theme_bw(base_size=base_size) +
        theme(panel.grid=element_blank(),
              axis.text.x=element_text(angle=45, hjust=1)) +
        scale_x_discrete(drop=FALSE) +
        xlab("") +
        ylab("") +
        guides(fill=guide_legend(title="")) +
        scale_fill_manual(values=colors,
                          drop=FALSE) +
        ggtitle(cancer[i])
    leg <- cowplot::get_legend(fig)
    fig2 <- fig + guides(fill="none")
    result <- list(figure=fig2, legend=leg)
    result
}

find_group_samples <- function(sample_id, connection_matrix) {

    group_samples <- c()
    check_samples <- sample_id

    while(length(check_samples) != 0) {
        for (check_sample_id in check_samples) {
            if(!(check_sample_id %in% group_samples)) {
                group_samples <- c(group_samples, check_sample_id)
                local_samples <- names(which(connection_matrix[check_sample_id,] == 1))
                check_samples <- setdiff(c(check_samples, local_samples), group_samples)
            }
        }
    }
    return(group_samples)
}


.ggRearrange2 <- function(df, ylabel="Read pair index",
                         basepairs=400, num.ticks=5){
  colors <- trellis:::readColors()[unique(df$read_type)]
  colors["splitread"] <- "black"
  nms <- names(trellis:::readColors())
  df$read_type <- factor(df$read_type, levels=nms)
  region <- read_type <- tagid <- NULL
  df1 <- filter(df, region==levels(region)[1])
  df2 <- filter(df, region==levels(region)[2])
  limits <- axis_limits(df, basepairs)
  gene1 <- levels(df$region)[1]
  gene2 <- levels(df$region)[2]
  xlim1 <- limits[[gene1]]
  xlim2 <- limits[[gene2]]
  labs1 <- trellis:::axis_labels5p(df1, xlim1, num.ticks)
  labs2 <- trellis:::axis_labels3p(df2, xlim2, num.ticks)
  a <- ggplot(df1, aes(ymin=tagid-0.2,
                       ymax=tagid+0.2,
                       xmin=start,
                       xmax=end,
                       color=read_type,
                       fill=read_type, group=tagid)) +
    geom_rect() +
    ylab(ylabel) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks=labs1[["breaks"]],
                       labels=labs1[["labels"]]) +
    coord_cartesian(xlim=xlim1) +
    xlab("") +
    theme(axis.text.x=element_text(size=7, angle=45, hjust=1),
          axis.text.y=element_blank(),
          plot.title=element_text(size=5)) +
    guides(color=FALSE, fill=FALSE) +
    geom_vline(xintercept=df$junction_5p[1], linetype="dashed") +
    ggtitle(paste0(df1$region[1], " (", df1$seqnames[1], ")"))
  if(df1$reverse[1]){
    a <- a + scale_x_reverse(breaks=labs1[["breaks"]],
                             labels=labs1[["labels"]])
  }
  b <- ggplot(df2, aes(ymin=tagid-0.2,
                       ymax=tagid+0.2,
                       xmin=start,
                       xmax=end,
                       color=read_type,
                       fill=read_type, group=tagid)) +
    geom_rect() +
    ylab("read pair index") +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks=labs2[["breaks"]],
                       labels=labs2[["labels"]]) +
    coord_cartesian(xlim=xlim2) +
    xlab("") +
    theme(axis.text.x=element_text(size=7, angle=45, hjust=1),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.direction="horizontal",
          plot.title=element_text(size=5)) +
    guides(color=FALSE, fill=FALSE) +
    geom_vline(xintercept=df$junction_3p[1], linetype="dashed") +
    ylab("") +
    ggtitle(paste0(df2$region[1], " (", df2$seqnames[1], ")"))
  if(df2$reverse[1]){
    b <- b + scale_x_reverse(breaks=labs2[["breaks"]],
                             labels=labs2[["labels"]])
  }
  ##
  ## plot both panels just to get the legend
  d <- ggplot(df, aes(ymin=tagid-0.2,
                      ymax=tagid+0.2,
                      xmin=start,
                      xmax=end,
                      color=read_type,
                      fill=read_type, group=tagid)) +
    geom_rect() +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    theme(legend.position="bottom", legend.direction="horizontal") +
    guides(color=guide_legend(title=""), fill=guide_legend(title=""))
  legend.grob <- trellis:::peelLegend(d)[[2]]
  agrob <- ggplotGrob(a)
  bgrob <- ggplotGrob(b)
  ##legend.grob <- gg.objs[[2]]
  bgrob$widths <- agrob$widths
  list(a=a,
       b=b,
       `5p`=agrob,
       `3p`=bgrob,
       legend=legend.grob)
}

#' Used in 06.1-figure3.rmd to create amplicon graphs
#'
#' This is a modification to functions with same name (but without '2' postfix) in svplots package.  Instead of returning just a list of grobs, the modification returns both the grobs and the ggplot objects.  This allows for further modification of the ggplot objects prior to creating grobs for the manuscript figure.
#' @export
ggRearrange2 <- function(df, ylab="Read pair index",
                        basepairs=400, num.ticks=5){
    . <- NULL
    grobs <- .ggRearrange2(df, ylabel=ylab,  basepairs, num.ticks)
    widths <- c(0.5, 0.5) %>%
        "/"(sum(.)) %>%
        unit(., "npc")
    heights <- c(0.95, 0.05) %>%
        "/"(sum(.)) %>%
        unit(., "npc")
    mat <- matrix(c(1, 2,
                    3, 3), byrow=TRUE, ncol=2, nrow=2)
    agrob <- grobs[["5p"]]
    bgrob <- grobs[["3p"]]
    legend.grob <- grobs[["legend"]]
    gobj <- arrangeGrob(agrob, bgrob,
                        legend.grob,
                        layout_matrix=mat,
                        widths=widths,
                        heights=heights)
    list(arranged.grobs=gobj,
         grobs=grobs)
}

#' Provides full names for cancer subtypes
#'
#' @export
cancer_names <- function(x){
    x2 <- x %>%
        mutate(tumor=Hmisc::capitalize(tumor_type))
    x3 <- x2 %>%
        mutate(tumor=case_when(tumor=="Colorectal"~"Colorectal mucinous",
                               tumor=="Pancreas"~"Pancreas mucinous",
                               tumor=="Stomach"~"Stomach mucinous",
                               TRUE~tumor))
    x3
}

#' Provide mucinous cancers
#'
#' @export
muc <- function() c("Colorectal mucinous", "Ovarian mucinous", "Pancreaas mucinous", "Stomach mucinous")

#' Provide endometrioid/endometrial cancers
#'
#' @export
endo <- function() c("Ovarian endometrioid", "Uterine endometrioid")

tumor_normal_matrix <- function(x){
    x.nested <- x %>%
        group_by(subject_id) %>%
        nest()
    nr <- map_int(x.nested$data, nrow)
    if(length(nr) < 4) return(NULL)
    x.nested2 <- x.nested[nr==2, ]
    x.nested2$data %>%
        map_dfr(function(x) x) %>%
        pull(propmeth) %>%
        matrix(nc=2, byrow=TRUE)
}

#' @export
pairedMeth <- function(methprop, manifest){
    ## comparison to matched normal for each tumor type
    manifest2 <- manifest %>%
        select(subject_id, lab_id, tumor_type,
               tumor.normal) %>%
        distinct()
    tumors <- filter(manifest, tumor.normal=="tumor")
    tumortypes <- tumors %>%
        select(subject_id, lab_id, tumor_type) %>%
        ungroup() %>%
        distinct()
    tt <- select(tumortypes, -lab_id) %>%
        distinct()
    meth2 <- methprop %>%
        select(Sample_Name, propmeth) %>%
        rename(lab_id=Sample_Name) %>%
        left_join(manifest2,
                  join_by(lab_id)) %>%
        select(-tumor_type) %>%
        left_join(tt, by="subject_id") %>%
        group_by(tumor_type) %>%
        nest()
    meth.matrix.list <- meth2$data %>%
        map(tumor_normal_matrix)
    nr <- sapply(meth.matrix.list, length)
    meth.matrix.list2 <- meth.matrix.list[nr > 0]
    meth.matrix.list2
}

#' Accessor for beta values of methylation SummarizedExperiment
#' @export
beta <- function(se) assays(se)[["beta"]]


#' Use consistent color scheme for cancer subtypes
#'
#' @export
tumor_colors <- function(){
    dx.colors <- c("Uterine endometrioid" = "#DDCC7F",
                   "Uterine endometrial" = "#DDCC7F",
                   "Ovarian endometrioid" = "#0F7554",
                   "Ovarian mucinous" = "#44AA99",
                   "Colorectal mucinous" = "#882255",
                   "Pancreatic mucinous" = "#AA4499",
                   "Stomach mucinous" = "#D695D0")
}

#' Create contingency table for comparing differences in mutation rates
#'
#' @export
complete_table <- function(x, Ns, tumor_order){
    x2 <- full_join(Ns, x, by=c("tumor_type", "n"))
    x2[is.na(x2)] <- 0L
    x3 <- select(x2, tumor_type, mt, n)
    x4 <- left_join(tumor_order, x3, by="tumor_type")
    x5 <- as.matrix(x4[, 2:3])
    rownames(x5) <- x4$tumor_type
    x5
}

#' Stan output
#'
#' @export
stan_output <- function(data.list, model,
                        params,
                        summaryfun,
                        probs=c(0.025, 0.05, 0.1,
                                0.5, 0.9, 0.95,
                                0.975)){
    tmp <- data.list %>%
        map(sampling2, model, params) %>%
        map(summaryfun, probs=probs) %>%
        map(1) %>%
        map(slice_params)
    ##%>%
    ##map(function(x) x[1, ])
    return(tmp)
}

#' Label Pancreas, Stomach, and Colorectal mucinous cancers as GI mucinous
#'
#' @export
collapse_gi <- function(dat){
    dat2  <- dat %>%
        mutate(tumor_type=as.character(tumor_type)) %>%
        mutate(tumor_type=case_when(tumor_type=="Pancreas mucinous"~"GI mucinous",
                                    tumor_type=="Stomach mucinous"~"GI mucinous",
                                    tumor_type=="Colorectal mucinous"~"GI mucinous",
                                    TRUE~tumor_type))
    dat2
}

#' x contains lab_id
#' manifest contains lab_id and lab_id2
#' replace lab_id with lab_id2 when not equal
#' @export
swap_lab_id <- function(x, y){
    y <- select(y, lab_id, lab_id2)
    isf <- is.factor(x$lab_id)
    if(isf){
        levs <- tibble(lab_id=levels(x$lab_id)) %>%
            inner_join(y, by="lab_id")
        levs2 <- levs$lab_id2
    }
    x.y <- inner_join(x, y, by="lab_id") %>%
        mutate(lab_id=ifelse(lab_id==lab_id2, lab_id, lab_id2)) %>%
        select(-lab_id2)
    if(isf){
        x.y$lab_id <- factor(x.y$lab_id, levs2)
    }
    return(x.y)
}


subject_id2 <- function(manifest){
    one.to.many <- select(manifest, subject_id, genotype_id, platform) %>%
        group_by(genotype_id, platform) %>%
        summarize(one_to_many=length(unique(subject_id)) > 1,
                  subject_ids=paste(unique(subject_id), collapse=","),
                  .groups="drop") %>%
        filter(one_to_many)
    discordant <- unique(unlist(strsplit(one.to.many$subject_ids, ",")))
    manifest$subject_id2 <- ifelse(manifest$subject_id %in% discordant,
                                   manifest$genotype_id,
                                   manifest$subject_id)
    manifest
}

.tumor_clinical_data <- function(obj){
    if(nrow(obj) == 1) return(obj)
    ntypes <- length(unique(obj$tumor_type))
    if(ntypes > 1){
        obj$discordant_tumor_type <- TRUE
    }
    obj2 <- filter(obj, tumor.normal=="tumor") %>%
        distinct()
    if(nrow(obj2)==1) return(obj2)
    obj2
}

clinical_data <- function(x){
    x2 <- select(x,
                 lab_id,
                 subject_id2,
                 tumor.normal,
                 tumor_type,
                 sex,
                 age_at_diagnosis_or_surgery,
                 stage_at_first_diagnosis,
                 smoker,
                 discordant_tumor_type)
    x3 <- ungroup(x2) %>%
        group_by(subject_id2) %>%
        nest()
    x3$data <- x3$data %>%
        map(.tumor_clinical_data)
    x4 <- unnest(x3, "data")
    x4
}

clean_sdata <- function(sdat){
    ##varnames <- tolower(colnames(sdat)) %>%
    ##    str_replace_all("([\\.])\\1+", ".") %>%
    ##    str_replace_all("\\.$", "") %>%
    ##    str_replace_all("\\.[01]\\.", ".") %>%
    ##    str_replace("internal.id", "lab_id") %>%
    ##    str_replace("source", "source_contact") %>%
    ##    str_replace("gender", "sex") %>%    
    ##    str_replace_all("\\.", "_")
    molecular <- sdat[, 17:26]
    sdat2 <- sdat[, 1:16]
    ## phi
    sdat3 <- sdat2[, -c(3, 4, 6, 11, 16)]
    varnames <- c("lab_id", "contact", 
                  "age_dx",
                  "age_surgery",
                  "sex",
                  "tumor_type",
                  "stage",
                  "pfs",
                  "is_alive",
                  "os",
                  "days_dx")
    colnames(sdat3) <- varnames
    ## clean up
    age.sx <- sdat3$age_surgery
    age.sx[age.sx==""] <- NA
    age.sx <- sapply(strsplit(sdat3$age_surgery, "/"), "[", 1)
    age.sx <- as.integer(age.sx)
    sdat3$age_surgery <- age.sx
    sdat3
}

clindata_description <- function(varnames){
    clindata.descr <- tibble(varname=varnames,
                             description=c("Unique lab identifier for the sample",
                                           "Contact PI providing sample",
                                           "Age at diagnosis (years)",
                                           "Age at surgery (years)",
                                           "Sex",
                                           "Histological tumor type",
                                           "FIGO stage (1988)",
                                           "Progression-free survival from diagnosis (days)",
                                           "Is alive (TRUE, FALSE)",
                                           "Overall survival from diagnosis (days)",
                                           "Days from diagnosis"))
    clindata.descr
}

join_clinical_data <- function(clinical.data, sdat3){
    cdat <- clinical.data %>%
        left_join(sdat3, join_by(lab_id, sex, tumor_type)) %>%
        mutate(age=ifelse(age_at_diagnosis_or_surgery=="NA", NA,
                          as.integer(age_at_diagnosis_or_surgery))) %>%
        select(-age_at_diagnosis_or_surgery) %>%
        mutate(stage=ifelse(is.na(stage), stage_at_first_diagnosis, stage)) %>%
        select(-stage_at_first_diagnosis) %>%
        mutate(age_surgery=ifelse(is.na(age_surgery), age, age_surgery)) %>%
        select(-age) %>%
        mutate(is_alive=ifelse(is_alive==0, TRUE, FALSE),
               pfs=ifelse(pfs=="-", NA, as.numeric(pfs)))
    cdat
}

overall_survival <- function(cdat){
    os <- cdat$os
    months <- os[grepl("m", os)] %>%
        strsplit("m") %>%
        sapply("[", 1) %>%
        as.numeric()
    days <- months * 30
    os[grepl("m", os)] <- days
    date.range <- os[grepl("/", os)] %>%
        str_split("-") %>%
        unlist %>%
        matrix(3, 2, byrow=TRUE) %>%
        set_colnames(c("start", "end")) %>%
        as_tibble() %>%
        mutate(start=ymd(start),
               end=ymd(end)) %>%
        mutate(days=end-start) %>%
        mutate(days=as.numeric(days))
    os[grepl("/", os)] <- date.range$days
    os <- as.numeric(os)
    return(os)
}

format_cancer_stage <- function(cdat){
    stage <- cdat$stage
    stage[grepl("n/a", stage)] <- NA
    stage <- ifelse(stage=="NA", NA, stage) %>%
        ifelse(.=="", NA, .) %>%
        str_replace_all("Ⅱ", "II") %>%
        str_replace_all("Ⅰ", "I") %>%
        str_replace_all("1", "I") %>%
        str_replace_all("2", "II") %>%
        str_replace_all("3", "III") %>%
        str_replace_all("4", "IV") %>%
        toupper()
    stage
}

check_lab_ids <- function(sdat3, manifest2){
    stopifnot(all(sdat3$lab_id %in% manifest2$lab_id))
    test <- select(manifest2, lab_id, stage_at_first_diagnosis,
                   age_at_diagnosis_or_surgery, sex, tumor_type) %>%
        mutate(age.y=age_at_diagnosis_or_surgery,
               age.y=ifelse(age.y=="NA", NA, age.y),
               age.y=as.integer(age.y)) %>%
        select(-age_at_diagnosis_or_surgery)
    check <- left_join(sdat3, test, by="lab_id") %>%
        select(lab_id, age_dx, age.y, age_surgery,
               sex.x, sex.y,
               tumor_type.x, tumor_type.y, stage, stage_at_first_diagnosis)
    ##    select(check, lab_id, age_dx, age_surgery, age.y) %>%
    ##        as.data.frame()
    ##    table(check$tumor_type.x, check$tumor_type.y)
    ##    select(check, stage, stage_at_first_diagnosis) %>%
    ##        as.data.frame()
    TRUE
}

select_clinical_columns <- function(cdat){
    cdat2 <- cdat %>%
        select(subject_id, lab_id, sex,
               tumor.normal, tumor_type,
               stage, smoker, age_dx,
               age_surgery, days_dx,
               is_alive, pfs, os, contact,
               discordant_tumor_type) %>%
        mutate(sex=ifelse(is.na(sex) & grepl("ovarian", tumor_type),
                          "Female", sex)) %>%
        ungroup()
    cdat2
}

select_manifest_columns <- function(manifest){
    manifest3 <- select(manifest, pgdx_id, lab_id,
                        subject_id, bamfile, tumor.normal,
                        bam_local,
                        size, platform, tumor_type, 
                        genotype_id, subject_id, subject_id2)
    manifest3
}

update_manifest_ids <- function(manifest){
    manifest2 <- manifest %>%
        mutate(temp=subject_id,
               subject_id=subject_id2,
               subject_id2=temp) %>%
        select(-temp)
}

discordant_tumor_type <- function(manifest4){
    ##
    ## update discordant_tumor_type for the matched normal sample according to
    ## the status of this variable for the cancer
    normals <- filter(manifest4, tumor.normal=="normal") %>%
        arrange(subject_id)
    tumors <- filter(manifest4, tumor.normal=="tumor") %>%
        arrange(subject_id)
    discord_label <- tumors$discordant_tumor_type
    names(discord_label) <- tumors$subject_id
    normals$discordant_tumor_type <- discord_label[ normals$subject_id ]
    manifest4 <- bind_rows(normals, tumors) %>% arrange(lab_id)
    manifest4
}
remove_any_duplicates <- function(manifest4){
    ##
    ## Make sure the lab ids are unique
    ##
    dups <- manifest4$lab_id[duplicated(manifest4$lab_id)]
    manifest5 <- manifest4 %>%
        unite(uid, c(lab_id, platform), sep="_",
              remove=FALSE) %>%
        mutate(lab_id=ifelse(lab_id %in% dups,
                             uid,
                             lab_id)) %>%
        select(-uid) %>%
        filter(!duplicated(lab_id)) %>%
        distinct()
    manifest5
}

read_facets <- function(file){
    facets <- read_tsv(file, show_col_types=FALSE) %>%
        select(Sample) %>%
        distinct()
    facets
}

join_facets_to_manifest1 <- function(facets, manifest.list){
    rename <- dplyr::rename
    manifest6 <- man(manifest.list)
    facets.matched <- inner_join(facets, key(manifest.list),
                                 by=c("Sample"="pgdx_id")) %>%
        rename(facet_id=Sample)
    facets.nomatch <- filter(facets, !Sample %in% manifest6$pgdx_id) %>%
        mutate(tmpid=str_replace_all(Sample, "t_", ""),
               tmpid=str_replace_all(tmpid, "_eland", ""),
               tmpid=str_replace_all(tmpid, "_hg18_", ""),
               tmpid=str_replace_all(tmpid, "_ExA", "_Ex"))
    facets.matched2 <- inner_join(facets.nomatch, key(manifest.list),
                                  by=c("tmpid"="pgdx_id")) %>%
        rename(facet_id=Sample)
    facets.nomatch <- filter(facets.nomatch, !tmpid %in% manifest6$pgdx_id) 
    uid <- tibble(uid=unique(facets.nomatch$tmpid),
                  uid2=NA)
    for(i in seq_len(nrow(uid))){
        id <- uid$uid[i]
        stripid <- strsplit(id, "_Ex_")
        if(length(stripid[[1]])==2){
            uid$uid2[i] <- stripid[[1]][2]
            if(!uid$uid2[i] %in% manifest6$pgdx_id) stop()
            next()
        }
        stripid <- strsplit(id, "_")[[1]][1]
        ix <- grep(stripid, manifest6$pgdx_id)
        stripid <- manifest6$pgdx_id[ix]
        if(!stripid %in% manifest6$pgdx_id) stop()
        uid$uid2[i] <- stripid
    }
    facets.nomatch2 <- left_join(facets.nomatch, uid,
                                 by=c("tmpid"="uid"))
    facets.matched3 <- inner_join(facets.nomatch2, key(manifest.list),
                                  by=c("uid2"="pgdx_id")) %>%
        rename(facet_id=Sample)
    facets2 <- bind_rows(facets.matched,
                         facets.matched2,
                         facets.matched3) %>%
        select(facet_id, subject_id, lab_id)
    ##manifest7 <- left_join(manifest5, facets2, by=c("subject_id", "lab_id"))
    manifest7 <- left_join(manifest6, facets2, by=c("subject_id", "lab_id"))
    ## The CGOV177 samples were mapped back to the same pgdx id
    ix <- grep("177T", manifest7$lab_id)
    A <- manifest7[ix, ]
    ## drop the second instance
    manifest7 <- manifest7[-ix[2], ]
}

join_facets_to_manifest2 <- function(manifest7, directory.listing){
    rename <- dplyr::rename
    ##
    ## Add FACETS id for WGS samples to facilitate merging FACET processed data
    ##
    facets <- tibble(Sample=directory.listing) %>% 
        mutate(Sample=basename(Sample))
    man <- select(manifest7, subject_id, lab_id, pgdx_id)
    matched0 <- inner_join(facets, man,
                           by=c("Sample"="lab_id")) %>%
        mutate(facet_id=Sample) %>%
        rename(lab_id=Sample)
    tmp <- filter(facets, !Sample %in% man$lab_id)
    matched1 <- inner_join(tmp, man,
                           by=c("Sample"="pgdx_id")) %>%
        mutate(pgdx_id=Sample,
               facet_id=pgdx_id) %>%
        select(-Sample)
    notmatched <- filter(tmp, !Sample %in% matched1$pgdx_id)
    tmp <- select(notmatched, Sample) %>%
        distinct() %>%
        mutate(uid=NA,
               lab_id=NA)
    for(i in seq_len(nrow(tmp))){
        id <- tmp$Sample[i]
        if(id == "CGOV359T"){
            tmp$lab_id[i] <- id
            next()
        }
        if(id == "CGOV482") {
            id <- "CGOV482T"
            tmp$lab_id[i] <- tmp$uid[i] <- id
            next()
        }
        if(grepl("^LP", id)){
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
    matched2 <- left_join(notmatched2, tmp3, by="Sample") %>%
        rename(facet_id=Sample) %>%
        left_join(man, by="lab_id") %>%
        select(-uid)
    ##
    ## These samples were rerun and have directories with labels GT_
    ##
    tmp2 <- filter(tmp, is.na(uid)) %>%
        select(-uid) %>%
        rename(old_directory=Sample) %>%
        mutate(facet_id=NA,
               lab_id=NA)
    man <- select(manifest7, subject_id, lab_id, pgdx_id, subject_id2, tumor.normal) %>%
        filter(tumor.normal=="tumor")
    for(i in seq_len(nrow(tmp2))){
        id <- tmp2$old_directory[i]
        ix <- grep(id, man$subject_id2)
        tmp2$facet_id[i] <- man$subject_id[ix]
        tmp2$lab_id[i] <- man$lab_id[ix]
    }
    matched3 <- tmp2 %>%
        mutate(subject_id=facet_id, pgdx_id=NA) %>%
        select(-old_directory)
    m <- bind_rows(matched0,
                   matched1,
                   matched2,
                   matched3)
    manifest.notwgs <- filter(manifest7, is.na(facet_id))
    manifest.wes <- filter(manifest7, !is.na(facet_id))
    manifest.notwgs2 <- left_join(select(manifest.notwgs,
                                         -facet_id),
                                  select(m, lab_id, facet_id),
                                  by="lab_id")
    manifest8 <- bind_rows(manifest.wes, manifest.notwgs2) %>%
        arrange(subject_id)
    manifest8
}

check_ids <- function(clinical, manifest){
    stopifnot(all(clinical$lab_id %in% manifest$lab_id))
    stopifnot(all(clinical$subject_id %in% manifest$subject_id))
    stopifnot(all(manifest$subject_id %in% clinical$subject_id))
    tmp <- filter(manifest, lab_id=="CGOV177T_2")
    should.be.true <- nrow(tmp)==1
    should.be.true
}

clean_clinical_data <- function(sdat, manifest2){
    rename <- dplyr::rename
    clinical.data <- manifest2 %>%
        mutate(discordant_tumor_type=FALSE) %>%
        clinical_data() %>%
        rename(subject_id=subject_id2)
    
    sdat3 <- clean_sdata(sdat)
    clindata.descr <- clindata_description(colnames(sdat3))
    stopifnot(check_lab_ids(sdat3, manifest2))
    ##
    ## 'stage' is more complete, but additional entries in manifest
    ##
    ## Merge with clinical data and keep sample identifier as some covariates are time-dependent
    ##
    clinical.data$sex[clinical.data$lab_id=="CGOV104T_Rep"] <- "Female"
    stopifnot(all(sdat3$lab_id %in% clinical.data$lab_id))
    cdat <- join_clinical_data(clinical.data, sdat3)
    cdat$os <- overall_survival(cdat)
    cdat$stage <- format_cancer_stage(cdat)
    cdat[cdat == "NA"] <- NA
    cdat2 <- select_clinical_columns(cdat)
    cdat2
}

key <- function(manifest.list){
    manifest.list[["key"]]
}

man <- function(manifest.list){
    manifest.list[["manifest"]]
}

clean_manifest <- function(manifest2, cdat2){
    manifest3 <- select_manifest_columns(manifest2)
    ##
    ## Use the tumor_type in the clinical.data
    ##
    manifest4 <- update_manifest_ids(manifest3) %>%
        select(-tumor_type) %>%
        left_join(select(cdat2, lab_id, tumor_type,
                         discordant_tumor_type),
                  by="lab_id")
    manifest4 <- discordant_tumor_type(manifest4)
    manifest5 <- remove_any_duplicates(manifest4)
    stopifnot(!any(duplicated(manifest5$lab_id)))
    ##
    ## Attach ids from FACETS copy number analysis
    ##
    manifest6 <- filter(manifest5, tumor.normal=="tumor") %>%
        select(subject_id, lab_id, pgdx_id)
    list(key=manifest6, manifest=manifest5)
}

join_with_facets <- function(manifest6, facets, directory.listing){
    manifest7 <- join_facets_to_manifest1(facets, manifest6)
    manifest8 <- join_facets_to_manifest2(manifest7, directory.listing)
    manifest8
}

filter_discordant_tumors <- function(manifest8){
    manifest <- filter(manifest8, !discordant_tumor_type) %>%
        filter(subject_id != "CGOV359") %>%
        mutate(lab_id=ifelse(lab_id == "CGOV151Tb_WES", "CGOV151Tb", lab_id))
    manifest
    
}

read_sdata <- function(sfile){
    sdat <- readRDS(sfile) %>%
        as_tibble()
    sdat
}

save_object <- function(object, nm){
    filepath <- file.path("data", paste0(nm, ".rda"))
    assign(nm, value=object)
    save(object, file=filepath, list=nm)
    filepath
}

load2 <- function(nm){
    load(file.path("data", nm), envir=parent.frame(2))
}

subset_by_tumors <- function(manifest, tumor.levels){
    manifest <- manifest %>%
        cancer_names() %>%
        filter(tumor %in% tumor.levels) %>%
        filter(tumor.normal=="tumor")
    manifest
}

pathway_levels <- function(){
    pathway.levels <- c("PI3K", "Ras and TK receptors",
                        "Chromatin Regulating",
                        "Cell cycle",
                        "Notch", "DNA repair",
                        "Mismatch repair",
                        "WNT", "TGFBR", "JAK/STAT",
                        "Other", "Large gene")
    pathway.levels
}

mucinous_pathways <- function(mucinous.levels.file){
    muc.pathways <- c("Ras and TK receptors",
                    "PI3K",
                    "Chromatin Regulating",
                    "Cell cycle",
                    "Notch",
                    "DNA repair",
                    "Mismatch repair",
                    "WNT",
                    "TGFBR", "JAK/STAT",
      "Other", "Large gene")
    levels <- readRDS(mucinous.levels.file)
    levels$pathway <- muc.pathways
    levels
}

read_idat <- function(idat.file, manifest, pathways){
    rename <- dplyr::rename
    tumortypes <- select(manifest, lab_id, tumor_type) %>%
        ungroup() %>%
        distinct()    
    idat <- idat.file %>%
        readRDS() %>%
        filter(lab_id %in% manifest$lab_id) %>%
        left_join(pathways, by=c("gene"="gene_symbol")) %>%
        left_join(tumortypes, by="lab_id") %>%
        distinct() %>%
        filter(!is.na(pathway)) %>%
        mutate(alteration=ifelse(type=="mutation",
                                 "mutation", alteration)) %>%
        cancer_names() %>%
        select(-tumor_type) %>%
        rename(tumor_type=tumor)

    idat
}

gene_list <- function(idat, pathway.levels2, tumor.levels){
    gene.list <- idat %>%
        ##filter(pathway != "Hypermutator") %>%
        mutate(pathway=factor(pathway, pathway.levels2),
               tumor_type=factor(tumor_type, tumor.levels)) %>%
        group_by(gene, pathway, tumor_type) %>%
        summarize(n=length(unique(lab_id)),
                  .groups="drop") %>%
        arrange(pathway, tumor_type, n) %>%
        group_by(pathway) %>%
        nest()
    gene.list
}

remove_duplicated_genes <- function(gene.list){
    gl <- filter(gene.list, pathway != "Hypermutator") %>%
        pull(data) %>%
        map(function(x){
            filter(x, !duplicated(gene)) %>%
                arrange(n)
        })
    gl
}

gene_levels <- function(gene.list){
    gene.levels <- unnest(gene.list, "data") %>%
        pull(gene) %>%
        unique()
}

endo_order <- function(idat){
    genes.for.sample.order <- idat %>%
        filter(pathway=="PI3K",
               tumor_type=="Ovarian endometrioid") %>%
        group_by(gene) %>%
        summarize(n=length(unique(lab_id)),
                  .groups="drop") %>%
        arrange(-n)
    genes.for.sample.order
}

muc_order <- function(idat){
    genes.for.sample.order <- idat %>%
        filter(pathway=="Ras and TK receptors",
               tumor_type=="Ovarian mucinous") %>%
        group_by(gene) %>%
        summarize(n=length(unique(lab_id)),
                  .groups="drop") %>%
        arrange(-n)
}

endo_id_levels <- function(idat, genes.for.sample.order){
    ovarian.order <- filter(idat, tumor_type=="Ovarian endometrioid") %>%
        mutate(gene_symbol=gene) %>%
        order_samples(gene.levels=genes.for.sample.order$gene)
    uterine.order <- filter(idat, tumor_type=="Uterine endometrial") %>%
        mutate(gene_symbol=gene) %>%
        order_samples(gene.levels=genes.for.sample.order$gene)
    id.levels <- c(ovarian.order$lab_id, uterine.order$lab_id)
}


muc_id_levels <- function(idat, genes.for.sample.order){
    ovarian.order <- filter(idat, tumor_type=="Ovarian mucinous") %>%
        mutate(gene_symbol=gene) %>%
        order_samples(gene.levels=genes.for.sample.order$gene)
    crc.order <- filter(idat, tumor_type=="Colorectal mucinous") %>%
        mutate(gene_symbol=gene) %>%
        order_samples(gene.levels=genes.for.sample.order$gene)
    id.levels <- c(ovarian.order$lab_id, crc.order$lab_id)
}

order_idat_endo <- function(idat, tumor.levels, pathway.levels, sample.order){
    ## **Genes:**
    ## - order by frequency within each pathway

    ## **Samples:**
    ##Order samples by
    ##  *  (1) hypermutator status (no hypermutator first)
    ##  *  (2) mutation status of most mutated gene in first pathway
    ##  *  (3) mutation status of second most comply mutated gene in first pathway

    ##  We can not determine the gene ordering solely by ovarian endometrioid samples as some genes are only altered in uterine endometrial.
    ## For each pathway, sort by tumor type and then frequency
    ##     - drop genes in uterine endometrial that are already in ovarian endometrioid
    ## decide gene order from ovarian endometrioid
    gene.list <- gene_list(idat, pathway.levels, tumor.levels)
    ## For each pathway, drop uterine endometrial genes that are already included in ovarian endometrioid
    gl <- remove_duplicated_genes(gene.list)
    gene.list$data[gene.list$pathway != "Hypermutator"] <- gl
    gene.levels <- gene_levels(gene.list)
    ##order_samples <- ovarian.subtypes:::order_samples
    ##genes.for.sample.order <- sample_order(idat)
    id.levels <- endo_id_levels(idat, sample.order)
    plevels <- pathway.levels
    gene.levels2 <- gene.levels[gene.levels != 'hypermutator']
    idat2 <- idat %>%
        filter(gene != "hypermutator") %>%
        mutate(lab_id=factor(lab_id, id.levels),
               gene=factor(gene, gene.levels2),
               pathway=factor(pathway, plevels),
               tumor_type=factor(tumor_type, tumor.levels))
    idat2
}

order_idat_mucinous <- function(idat, tumor.levels, pathway.levels, sample.order){
    ## **Genes:**
    ## - order by frequency within each pathway

    ## **Samples:**
    ##Order samples by
    ##  *  (1) hypermutator status (no hypermutator first)
    ##  *  (2) mutation status of most mutated gene in first pathway
    ##  *  (3) mutation status of second most comply mutated gene in first pathway

    ##  We can not determine the gene ordering solely by ovarian endometrioid samples as some genes are only altered in uterine endometrial.
    ## For each pathway, sort by tumor type and then frequency
    ##     - drop genes in uterine endometrial that are already in ovarian endometrioid
    ## decide gene order from ovarian endometrioid
    gene.list <- gene_list(idat, pathway.levels, tumor.levels)
    ## For each pathway, drop uterine endometrial genes that are already included in ovarian endometrioid
    gl <- remove_duplicated_genes(gene.list)
    gene.list$data[gene.list$pathway != "Hypermutator"] <- gl
    gene.levels <- gene_levels(gene.list)
    ##order_samples <- ovarian.subtypes:::order_samples
    ##genes.for.sample.order <- sample_order(idat)
    id.levels <- muc_id_levels(idat, sample.order)
    plevels <- pathway.levels
    gene.levels2 <- gene.levels[gene.levels != 'hypermutator']
    idat2 <- idat %>%
        filter(gene != "hypermutator") %>%
        mutate(lab_id=factor(lab_id, id.levels),
               gene=factor(gene, gene.levels2),
               pathway=factor(pathway, plevels),
               tumor_type=factor(tumor_type, tumor.levels))
    idat2
}

order_endo <- function(idat, manifest){
    tumor.levels <- c("Ovarian endometrioid", "Uterine endometrial")
    pathway.levels <- pathway_levels()
    pathway.levels2 <- c("Hypermutator", pathway.levels)
    sample.order <- endo_order(idat)
    idat2 <- order_idat_endo(idat, tumor.levels, pathway.levels2, sample.order)
    idat3 <- remove_duplicate_samples(idat2, manifest)
    idat3$pathway <- droplevels(idat3$pathway)
    idat3
}

exclude_genes_only_mutated_in_cgcrc254 <- function(idat){
    genes.to.drop <- idat %>% group_by(gene) %>%
        summarize(ids=paste(unique(lab_id), collapse=",")) %>%
        filter(ids=="CGCRC254T")
    idat2 <- filter(idat, !gene %in% genes.to.drop$gene)
    idat2
}

order_mucinous <- function(idat, manifest, pathway.levels){
    tumor.levels <- c("Ovarian mucinous", "Colorectal mucinous")
    sample.order <- muc_order(idat)
    ##trace(order_idat_mucinous, browser)
    pathway.levels2 <- c("Hypermutator", pathway.levels)
    idat2 <- order_idat_mucinous(idat, tumor.levels, pathway.levels2, sample.order)
    idat3 <- remove_duplicate_samples(idat2, manifest)
    ## Exclude rows that only appear in extreme hypermutator CGCRC254T
    idat4 <- exclude_genes_only_mutated_in_cgcrc254(idat3)
    idat4$pathway <- droplevels(idat4$pathway)
    idat4$gene <- droplevels(idat4$gene)
    idat4
}

order_gi <- function(idat, manifest){
    tumor.levels <- c("Colorectal mucinous", "Stomach mucinous", "Pancreas mucinous")
    idat2 <- order_idat(idat, tumor.levels)
    ##idat3 <- remove_duplicate_samples(idat2, manifest)
    idat2
}


remove_duplicate_samples <- function(idat2, manifest){
    ##
    ## For patients with multiple samples, only include a single sample
    ##
    dup.samples <- idat2 %>%
        select(lab_id) %>%
        distinct() %>%
        left_join(select(manifest, subject_id, lab_id), by="lab_id") %>%
        group_by(subject_id) %>%
        nest()
    nr <- map_dbl(dup.samples$data, nrow)
    dup.samples2  <- dup.samples[nr > 1, ]
    drop.samples <- dup.samples$data %>% map_dfr(function(x) x[-1, ])
    if(any("CGOV141T_1" %in% drop.samples$lab_id))
        drop.samples$lab_id[match("CGOV141T_1", drop.samples$lab_id)] <- "CGOV141T"
    idat2 <- filter(idat2, !lab_id %in% drop.samples$lab_id)
    id.levels <- levels(idat2$lab_id)
    id.levels <- id.levels[!id.levels %in% drop.samples$lab_id]
    idat2$lab_id <- factor(idat2$lab_id, id.levels)
    idat2
}

subset_endo <- function(manifest){
    tumor.levels <- c("Ovarian endometrioid", "Uterine endometrial")
    manifest2 <- subset_by_tumors(manifest, tumor.levels)
    manifest2
}


subset_gi <- function(manifest){
    tumor.levels <- c("Colorectal mucinous", "Stomach mucinous", "Pancreas mucinous")
    manifest <- manifest %>%
        cancer_names() %>%
        filter(tumor %in% tumor.levels) %>%
        filter(tumor.normal=="tumor")
    manifest
}

subset_mucinous <- function(manifest){
    tumor.levels <- c("Ovarian mucinous", "Colorectal mucinous")
    manifest <- manifest %>%
        cancer_names() %>%
        filter(tumor %in% tumor.levels) %>%
        filter(tumor.normal=="tumor")
    manifest
}

compare <- function(obj1, obj2){
    all(obj1$lab_id %in% obj2$lab_id)
    all(obj1$gene %in% obj2$gene)
    identical(levels(obj1$lab_id),
              levels(obj2$lab_id))
    ## levels are identical after dropping unused hypermutator level
    identical(levels(obj1$gene),
              levels(obj2$gene))
    ## 70 levels in idat.muc.ordered
    ## only 52 levels in idat.mucinous
    identical(levels(obj1$pathway),
              levels(obj2$pathway))
}

read_methylation_se <- function(file, manifest, discordant){
    rename <- dplyr::rename
    se <- readRDS(file)
    meth <- tibble(lab_id=colnames(se)) %>%
        mutate(platform="methylation") %>%
        filter(!lab_id %in% discordant$lab_id)
    meth2 <- select(manifest, -platform) %>%
        left_join(meth, by="lab_id") %>%
        filter(!is.na(platform))
    any(!meth$lab_id %in% meth2$lab_id)
    meth.fuzzymatch <- filter(meth, !lab_id %in% meth2$lab_id) %>%
        rename(alt_id=lab_id) %>%
        mutate(lab_id=c("CGCRC330N_1",
                        "CGCRC330T_1",
                        "CGCRC330T1_1",
                        "CGOV177T_2",
                        "CGOV179T_Rpt",
                        "CGOV186T_2",
                        "CGOV188N_2",
                        "CGOV188T_2",
                        "CGST1N_1",
                        "CGST1T_2",
                        "CGST2T_2")) %>%
        select(-alt_id)
    meth3 <- select(manifest, -platform) %>%
        left_join(meth.fuzzymatch, by="lab_id") %>%
        filter(!is.na(platform))
    meth4 <- bind_rows(meth2, meth3)
    methylation <- meth4
    methylation
}

read_methylation_data <- function(file, tcga.file){
    rename <- dplyr::rename
    se <- readRDS(tcga.file)
    metadata <- readRDS(file) %>%
        as_tibble() %>%
        filter(grepl("^C", Sample_Name)) %>%
        rename(lab_id=Sample_Name)
    metadata2 <- colData(se) %>%
        as_tibble()
    md <- left_join(metadata2, metadata,
                    join_by(lab_id)) %>%
        select(-c(Diagnosis, sampletype, Sample)) %>%
        set_colnames(tolower(colnames(.))) %>%
        mutate(t.n=substr(tumor, 1, 1))
    md2 <- as(md, "DataFrame")
    colData(se) <- md2
    colnames(se) <- se$lab_id
    se
}

check_against_manifest <- function(se, manifest, discordant){
    ## check that all jhu samples are in the manifest
    is_jhu <- se$study == "JHU"
    jhu <- se[, is_jhu]
    in_manifest <- colnames(jhu) %in% manifest$lab_id
    notin_manifest <- colnames(jhu)[!in_manifest]
    ## We only care about resolving matches for the samples that were not discordant
    notin_manifest2 <- notin_manifest[!notin_manifest %in% discordant$lab_id]
    dat <- tibble(to_map=notin_manifest2,
                  lab_id=NA)
    for(i in 1:nrow(dat)){
        id <- dat$to_map[i]
        if(id=="CGCRC330T"){
            dat$lab_id[i] <- "CGCRC330T_1"
            next()
        }
        ix <- grep(id, manifest$lab_id)
        if(length(ix)==1){
            dat$lab_id[i] <- manifest$lab_id[ix]
            next()
        }
        stop()
    }
    ix <- match(dat$to_map, colnames(jhu))
    colnames(jhu)[ix] <- dat$lab_id
    jhu$lab_id <- colnames(jhu)
    jhu2 <- jhu[, !colnames(jhu) %in% discordant$lab_id]
    stopifnot(all(colnames(jhu2) %in% manifest$lab_id))

    tcga <- se[, se$study=="TCGA"]
    methylation_se <- cbind(jhu2, tcga)
    methylation_se
}

get_manifest <- function(){
    data(manifest,package="ovarian.subtypes", envir=environment())
    manifest
}

#' @export
draw_heatmap <- function(heatmap.components){
    cluster_data <- heatmap.components$cluster_data
    ha_rows <- heatmap.components$ha_rows
    hegt <- Heatmap(cluster_data,
                    col = colorRamp2(c(0,0.25,0.5),
                                     c("#00FFCC","#FFFFFF","#0099FF")),
                    show_column_names = FALSE,
                    left_annotation = ha_rows,
                    show_column_dend = FALSE,
                    show_heatmap_legend=FALSE,
                    row_title_rot=0,
                    clustering_distance_rows = "euclidean")
    return(hegt)
}

#' @export
heatmap_setup <- function(methylation_se){
    metadata <- as_tibble(colData(methylation_se))
    df <- metadata %>%
        mutate(t.n=ifelse(t.n=="T", "Tumor", "Normal")) %>%
        rename(`Tissue type`=t.n) %>%
        mutate(endo.muc=ifelse(grepl("end", diagnosis),
                               "Endometrial",
                               "Mucinous"))
    colnames(df) <- Hmisc::capitalize(colnames(df))
    bvals <- beta(methylation_se)
    df$Diagnosis <- factor(df$Diagnosis,
                           levels = c("Uterine endometrial",
                                      "Ovarian endometrioid",
                                      "Ovarian mucinous",
                                      "Colorectal mucinous",
                                      "Pancreatic mucinous",
                                      "Stomach mucinous"))
    df$Tissue <- str_extract(df$Diagnosis, ".*(?= )")
    dx.colors <- tumor_colors()
    # use Ovarian mucinous colors as Ovarian tissue for both endo and muc
    # histology will differentiate between endo and muc 
    dx.colors <- dx.colors[-match(c("Uterine endometrioid", "Ovarian endometrioid"), names(dx.colors))]
    names(dx.colors) <- str_extract(names(dx.colors), ".*(?= )")
    study.colors <- c("JHU"="#002d72",
                      "TCGA"="gray90")
    tn.colors <- c("Normal" = "gray",
                   "Tumor" = "black")
    tn.lgd <- Legend(labels=names(tn.colors),
                      legend_gp=gpar(fill=tn.colors),
                      labels_gp=gpar(fontsize=18),
                      title_gp=gpar(fontsize=22),
                     title=" Tumor/Normal")
    dx.lgd <- Legend(labels=names(dx.colors),
                     legend_gp=gpar(fill=dx.colors),
                     labels_gp=gpar(fontsize=18),
                     title_gp=gpar(fontsize=22),
                     title=" Tissue type")
    histology.colors <- c("orange3", "steel blue")
    names(histology.colors) <- c("Endometrial", "Mucinous")
    histology.lgd <- Legend(labels=names(histology.colors),
                            legend_gp=gpar(fill=histology.colors),
                            labels_gp=gpar(fontsize=18),
                            title_gp=gpar(fontsize=22),
                            title=" Histology")
    col_fun <- colorRamp2(c(0,0.25,0.5),
                          c("#00FFCC","#FFFFFF","#0099FF"))
    beta.lgd <- Legend(col_fun = col_fun, title = expression(beta),
                       labels_gp=gpar(fontsize=18),
                       title_gp=gpar(fontsize=22))
    horiz.legends <- packLegend(dx.lgd, histology.lgd,
                                direction="horizontal")
    vert.legends <- packLegend(tn.lgd, beta.lgd, direction="vertical")
    df <- df[, c("Tissue", "Tissue type", "Study", "Endo.muc")]
    colnames(df) <- c(" Tissue type", " Tumor/Normal", " Study", " Histology")
    ha_rows  <-  rowAnnotation(df = df[, c(" Tissue type", " Tumor/Normal", " Histology")],
                               col = list(` Tissue type`=dx.colors,
                                          ` Tumor/Normal`=tn.colors,
                                          ` Histology`=histology.colors),
                               show_legend=FALSE,
                               annotation_name_rot=45,
                               annotation_name_side="top")
    ix <- seq_len(nrow(methylation_se))
    cluster_data <- t(bvals[ix, ])
    rownames(cluster_data) <- NULL
    result <- list(cluster_data=cluster_data,
                   ha_rows=ha_rows,
                   vert.legends=vert.legends,
                   horiz.legends=horiz.legends)
    result
}

read_facets2 <- function(...) {
    dots <- list(...)
    facets_purity <- lapply(dots, read.delim) %>%
        bind_rows() %>%
        select(facet_id = Sample, purity = Purity, genotype_id = Genotype.ID)
    return(facets_purity)
}

add_facets_purity <- function(manifest, facets_purity) {
    facets_purity <- mutate(facets_purity,
                            is_na_purity = is.na(purity))
    manifest2 <- manifest[!is.na(manifest$facet_id), ]
    manifest3 <- left_join(manifest2, select(facets_purity, facet_id, purity), by = "facet_id") %>%
        left_join(select(facets_purity, genotype_id, purity), by = join_by("facet_id" == "genotype_id")) %>%
        mutate(purity = case_when(is.na(purity.x) & is.na(purity.y) ~ NA,
                                  is.na(purity.x) & !is.na(purity.y) ~ purity.y,
                                  !is.na(purity.x) & is.na(purity.y) ~ purity.x,
                                  .default = NA)) %>%
        replace_na(list(is_na_purity = FALSE)) %>%
        select(-c(purity.x, purity.y))
    manifest4 <- left_join(manifest, select(manifest3, lab_id, purity))
    return(manifest4)
}

purity_filter <- function(manifest, threshold = 0.2) {
    manifest %>%
        mutate(purity = as.numeric(purity)) %>%
        filter(!(is.na(purity) | purity <= 0.2))
}

update_tcga_barcodes <- function(methylation_se, match_table) {
    # match_table was downloaded from the cluster
    barcodes <- match_table$Barcode.ID
    sum(is.na(barcodes)) # Only one is missing - check!
    # [1] 1
    # TCGA failed to update this sample
    # See https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/ for
    # more information
    barcodes[is.na(barcodes)] <- "TCGA-D5-6930-01A-11D-1926-05"
    names(barcodes) <- c(1:164)
    tissue_source <- match_table$Project.ID
    for (i in seq_along(tissue_source)) {
        if (is.na(tissue_source[i])) {
            # This is the NA observation
            tissue_source[i] <- "Colorectal mucinous"
        } else if (tissue_source[i] == "TCGA-COAD") {
            tissue_source[i] <- "Colorectal mucinous"
        } else if (tissue_source[i] == "TCGA-PAAD") {
            tissue_source[i] <- "Pancreatic mucinous"
        } else if (tissue_source[i] == "TCGA-STAD") {
            tissue_source[i] <- "Stomach mucinous"
        } else if (tissue_source[i] == "TCGA-UCEC") {
            tissue_source[i] <- "Uterine endometrial"
        } else {
            stop("Error! Wrong tissue source.")
            print(tissue_source[i])
        }
    }
    tissue_type <- match_table$Sample.Type
    tissue_type <- gsub("Primary ", "", tissue_type)
    tissue_type <- gsub("Solid Tissue ", "", tissue_type)
    tissue_type[is.na(tissue_type)] <- "Tumor"
    tissue_type_short <- sapply(tissue_type, function(x) substr(x, 1, 1))

    # Update the colData and colnames for the TCGA samples
    # Columns to update
    # lab_id, diagnosis, tumor, t.n
    # COAD = Colorectal mucinous
    # PAAD = Pancreas mucinous
    # STAD = Stomach mucinous
    # UCEC = Uterine endometrial
    allcols <- colnames(methylation_se)
    dx_levels <- levels(colData(methylation_se)$diagnosis)
    new_dx_levels <- dx_levels
    new_dx_levels[match("Pancreas mucinous", new_dx_levels)] <- "Pancreatic mucinous"
    colData(methylation_se)$diagnosis <- as.character(colData(methylation_se)$diagnosis)
    is_panc_muc <- colData(methylation_se)$diagnosis == "Pancreas mucinous"
    colData(methylation_se)$diagnosis[is_panc_muc] <- "Pancreatic mucinous"
    colData(methylation_se)$diagnosis[match(c(1:164), allcols)] <- tissue_source
    colData(methylation_se)$diagnosis <- factor(colData(methylation_se)$diagnosis,
                                                levels = new_dx_levels)
    colData(methylation_se)$tumor[match(c(1:164), allcols)] <- tissue_type
    colData(methylation_se)$t.n[match(c(1:164), allcols)] <- tissue_type_short
    colData(methylation_se)$lab_id[match(c(1:164), allcols)] <- barcodes
    colnames(methylation_se)[match(c(1:164), allcols)] <- barcodes
    return(methylation_se)
}
