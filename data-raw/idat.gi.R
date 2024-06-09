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
tumor.levels <- c("Colorectal mucinous", "Stomach mucinous", "Pancreas mucinous")
manifest <- manifest %>%
    cancer_names() %>%
    filter(tumor %in% tumor.levels) %>%
    filter(tumor.normal=="tumor")
gene.pathway <- here("output", "gene_pathway.rmd",
               "gene_pathway.csv") %>%
    read_csv(show_col_types=FALSE)
##levels <- here("output", "05-data_integration.rmd",
##               "mucinous_factor_levels.rds") %>%
##    readRDS()
##pathway.levels <- levels$pathway
##pathway.levels <- c(pathway.levels[1:5], "Notch",
##                    pathway.levels[6:11])
pathways <- ovarian.subtypes:::read_pathways("mucinous")
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
idat.gi <- idat
save(idat.gi,
     file=here("Rpackage/ovarian.subtypes/data/idat.gi.rda"))
