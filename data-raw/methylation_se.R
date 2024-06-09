library(SummarizedExperiment)
library(here)
library(tidyverse)
library(magrittr)
load(here("Rpackage", "ovarian.subtypes",
          "data", "manifest.rda"))
load(here("Rpackage", "ovarian.subtypes",
          "data", "discordant.rda"))
rename <- dplyr::rename
select <- dplyr::select
se.lab.tcga <- here("extdata",
                    "se_lab_tcga.rds") %>%
    readRDS()
methylation_se <- se.lab.tcga

## merge with other metadata on the samples
metadata2 <- colData(methylation_se) %>%
    as_tibble()
metadata <- readRDS(here("extdata", "combmetadata.rds")) %>%
    as_tibble() %>%
    filter(grepl("^C", Sample_Name)) %>%
    rename(lab_id=Sample_Name)
md <- left_join(metadata2, metadata,
                join_by(lab_id)) %>%
    select(-c(Diagnosis, sampletype, Sample)) %>%
    set_colnames(tolower(colnames(.))) %>%
    mutate(t.n=substr(tumor, 1, 1))
md2 <- as(md, "DataFrame")
colData(methylation_se) <- md2
colnames(methylation_se) <- methylation_se$lab_id

## check that all jhu samples are in the manifest
is_jhu <- methylation_se$study == "JHU"
jhu <- methylation_se[, is_jhu]
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
tcga <- methylation_se[, methylation_se$study=="TCGA"]
methylation_se <- cbind(jhu2, tcga)
fname <- here("Rpackage", "ovarian.subtypes",
              "data",
              "methylation_se.rda")
save(methylation_se, file=fname)
