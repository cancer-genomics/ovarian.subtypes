#knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
library(here)
library(magrittr)
library(ovarian.subtypes)
rename <- dplyr::rename
## Whether to include Hypermutator status as row in tile plot
hypermutator_as_gene <- FALSE
data(manifest, package="ovarian.subtypes")
tumor.levels <- c("Ovarian endometrioid", "Uterine endometrial")
manifest <- manifest %>%
    cancer_names() %>%
    filter(tumor %in% tumor.levels) %>%
    filter(tumor.normal=="tumor")
gene.pathway <- here("output", "gene_pathway.rmd",
               "gene_pathway.csv") %>%
    read_csv(show_col_types=FALSE)
pathway.levels <- c("PI3K", "Ras and TK receptors",
                    "Chromatin Regulating",
                    "Cell cycle",
                    "Notch", "DNA repair",
                    "Mismatch repair",
                    "WNT", "TGFBR", "JAK/STAT",
                    "Other", "Large gene")
pathways <- ovarian.subtypes:::read_pathways("endometrioid")
tumortypes <- select(manifest, lab_id, tumor_type) %>%
    ungroup() %>%
    distinct()
idat <- here("output", "05-data_integration.rmd",
             "integrated_data.rds") %>%
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
pathway.levels2 <- c("Hypermutator", pathway.levels)
gene.list <- idat %>%
    ##filter(pathway != "Hypermutator") %>%
    mutate(pathway=factor(pathway, pathway.levels2),
           tumor_type=factor(tumor_type, c("Ovarian endometrioid",
                                           "Uterine endometrial"))) %>%
    group_by(gene, pathway, tumor_type) %>%
    summarize(n=length(unique(lab_id)),
              .groups="drop") %>%
    arrange(pathway, tumor_type, n) %>%
    group_by(pathway) %>%
    nest()
## For each pathway, drop uterine endometrial genes that are already included in ovarian endometrioid
gl <- filter(gene.list, pathway != "Hypermutator") %>%
    pull(data) %>%
    map(function(x){
        filter(x, !duplicated(gene)) %>%
            arrange(n)
    })
gene.list$data[gene.list$pathway != "Hypermutator"] <- gl
gene.levels <- unnest(gene.list, "data") %>%
    pull(gene) %>%
    unique()
genes.for.sample.order <- idat %>%
    filter(pathway=="PI3K",
           tumor_type=="Ovarian endometrioid") %>%
    group_by(gene) %>%
    summarize(n=length(unique(lab_id)),
              .groups="drop") %>%
    arrange(-n)
order_samples <- ovarian.subtypes:::order_samples
ovarian.order <- filter(idat, tumor_type=="Ovarian endometrioid") %>%
    mutate(gene_symbol=gene) %>%
    order_samples(gene.levels=genes.for.sample.order$gene)
uterine.order <- filter(idat, tumor_type=="Uterine endometrial") %>%
    mutate(gene_symbol=gene) %>%
    order_samples(gene.levels=genes.for.sample.order$gene)
id.levels <- c(ovarian.order$lab_id, uterine.order$lab_id)
plevels <- pathway.levels
gene.levels2 <- gene.levels[gene.levels != 'hypermutator']
idat2 <- idat %>%
    filter(gene != "hypermutator") %>%
    mutate(lab_id=factor(lab_id, id.levels),
           gene=factor(gene, gene.levels2),
           pathway=factor(pathway, plevels),
           tumor_type=factor(tumor_type, tumor.levels))

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
##dup.samples3 <- unnest(dup.samples2, "data")
##stop('dups')
## keep the first
drop.samples <- dup.samples$data %>% map_dfr(function(x) x[-1, ])
drop.samples$lab_id[match("CGOV141T_1", drop.samples$lab_id)] <- "CGOV141T"

idat2 <- filter(idat2, !lab_id %in% drop.samples$lab_id)
id.levels <- levels(idat2$lab_id)
id.levels <- id.levels[!id.levels %in% drop.samples$lab_id]
idat2$lab_id <- factor(idat2$lab_id, id.levels)

idat.endometrioid <- idat2
save(idat.endometrioid,
     file=here("Rpackage",
               "ovarian.subtypes",
               "data",
               "idat.endometrioid.rda"))
