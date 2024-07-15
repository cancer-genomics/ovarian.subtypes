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
tumor.levels <- c("Ovarian mucinous", "Colorectal mucinous")
manifest <- manifest %>%
    cancer_names() %>%
    filter(tumor %in% tumor.levels) %>%
    filter(tumor.normal=="tumor")
gene.pathway <- here("output", "gene_pathway.rmd",
               "gene_pathway.csv") %>%
    read_csv(show_col_types=FALSE)
levels <- here("output", "05-data_integration.rmd",
               "mucinous_factor_levels.rds") %>%
    readRDS()
pathway.levels <- c("Ras and TK receptors",
                    "PI3K",
                    "Chromatin Regulating",
                    "Cell cycle",
                    "Notch",
                    "DNA repair",
                    "Mismatch repair",
                    "WNT",
                    "TGFBR", "JAK/STAT",
                    "Other", "Large gene")
levels$pathway <- pathway.levels
##pathway.levels <- levels$pathway
##pathway.levels <- c(pathway.levels[1:5], "Notch",
##                    pathway.levels[6:11])
pathways <- ovarian.subtypes:::read_pathways("mucinous") %>%
    filter(!(gene_symbol=="JAK1" & pathway=="Cell cycle")) %>%
    filter(!(gene_symbol=="MED1-STAT5B" & pathway=="Other"))
tmp <- pathways %>%
    group_by(gene_symbol) %>%
    summarize(n=length(unique(pathway)),
              pway=paste(pathway, collapse=",")) %>%
    filter(n > 1)
stopifnot(nrow(tmp) == 0)
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
##pathway.levels2 <- c("Hypermutator", pathway.levels)
gene.list <- idat %>%
    filter(gene!="hypermutator") %>%
    mutate(pathway=factor(pathway, pathway.levels),
           tumor_type=factor(tumor_type, tumor.levels)) %>%
    group_by(gene, pathway, tumor_type) %>%
    summarize(n=length(unique(lab_id)),
              .groups="drop") %>%
    arrange(pathway, tumor_type, n) %>%
    group_by(pathway) %>%
    nest()
fx <- function(x){
    filter(x, !duplicated(gene)) %>%
        arrange(n)
}
gl <- gene.list %>%
    filter(pathway != "Hypermutator") %>%
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
    filter(pathway=="Ras and TK receptors",
           tumor_type=="Ovarian mucinous") %>%
    group_by(gene) %>%
    summarize(n=length(unique(lab_id)),
              .groups="drop") %>%
    arrange(-n)
## Where is CGOV163T and CGOV473T have no events in genes and are
## not included
order_samples <- ovarian.subtypes:::order_samples
ovarian.order <- filter(idat, tumor_type=="Ovarian mucinous") %>%
    mutate(gene_symbol=gene) %>%
    order_samples(gene.levels=genes.for.sample.order$gene)
crc.order <- filter(idat, tumor_type=="Colorectal mucinous") %>%
    mutate(gene_symbol=gene) %>%
    order_samples(gene.levels=genes.for.sample.order$gene)
id.levels <- c(ovarian.order$lab_id, crc.order$lab_id)
gene.levels2 <- gene.levels[gene.levels != 'hypermutator']
plevels <- pathway.levels
idat2 <- idat %>%
    filter(gene != "hypermutator") %>%
    mutate(lab_id=factor(lab_id, id.levels),
           gene=factor(gene, gene.levels2),
           pathway=factor(pathway, plevels),
           tumor_type=factor(tumor_type, tumor.levels))

##
## For patients with multiple samples, include only a single sample
##
id.levels <- levels(idat2$lab_id)
dup.samples <- idat2 %>%
    select(lab_id) %>%
    distinct() %>%
    left_join(select(manifest, subject_id, lab_id), by="lab_id") %>%
    group_by(subject_id) %>%
    nest()
nr <- map_dbl(dup.samples$data, nrow)
dup.samples2  <- dup.samples[nr > 1, ]
## keep the first
drop.samples <- dup.samples$data %>% map_dfr(function(x) x[-1, ])
idat2 <- filter(idat2, !lab_id %in% drop.samples$lab_id)
id.levels <- id.levels[!id.levels %in% drop.samples$lab_id]
idat2$lab_id <- factor(idat2$lab_id, id.levels)

##
## Exclude genes that are only altered in CGCRC254T
##
exclude.genes <- idat2 %>%
    group_by(gene) %>%
    summarize(only.254=all(grepl("CGCRC254T", lab_id))) %>%
    filter(only.254)
idat2 <- filter(idat2, !gene %in% exclude.genes$gene)
gene.levels <- levels(idat2$gene)
gene.levels <- gene.levels[!gene.levels %in% exclude.genes$gene &
                           gene.levels %in% as.character(idat2$gene)]
idat2$gene <- factor(idat2$gene, gene.levels)

idat.mucinous <- idat2
save(idat.mucinous,
     file=here("../ovarian.subtypes/data/idat.mucinous.rda"))
