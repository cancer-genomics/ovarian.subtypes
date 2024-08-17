library(here)
library(tidyverse)
library(lubridate)
library(magrittr)
manifest <- here("inst", "extdata", 
                 "manifest.rds") %>%
    readRDS()
##
## subject_id2:  Define a group_id such that:
##    subject_id2 = subject_id when grouping by subject id and genotype id are the same
##                = genotype_id when grouping is discordant
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
##
## ASSUME THE TUMOR SAMPLE IS ANNOTATED CORRECTLY
##
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
manifest2 <- subject_id2(manifest)
clinical.data <- manifest2 %>%
    mutate(discordant_tumor_type=FALSE) %>%
    clinical_data() %>%
    rename(subject_id=subject_id2)
sdat <- readRDS("inst/extdata/sdat.rds") %>%
    as_tibble()
varnames <- tolower(colnames(sdat)) %>%
    str_replace_all("([\\.])\\1+", ".") %>%
    str_replace_all("\\.$", "") %>%
    str_replace_all("\\.[01]\\.", ".") %>%
    str_replace("internal.id", "lab_id") %>%
    str_replace("source", "source_contact") %>%
    str_replace("gender", "sex") %>%    
    str_replace_all("\\.", "_")
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
select(check, lab_id, age_dx, age_surgery, age.y) %>%
    as.data.frame()
table(check$tumor_type.x, check$tumor_type.y)
select(check, stage, stage_at_first_diagnosis) %>%
    as.data.frame()
##
## 'stage' is more complete, but additional entries in manifest
##
## Merge with clinical data and keep sample identifier as some covariates are time-dependent
##
clinical.data$sex[clinical.data$lab_id=="CGOV104T_Rep"] <- "Female"
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
cdat$os <- as.numeric(os)
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
cdat$stage <- stage
cdat[cdat == "NA"] <- NA
all(sdat3$lab_id %in% clinical.data$lab_id)
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
manifest3 <- select(manifest2, pgdx_id, lab_id,
                    subject_id, bamfile, tumor.normal,
                    bam_local,
                    size, platform, tumor_type, 
                    genotype_id, subject_id, subject_id2)
update_ids <- function(manifest){
    manifest2 <- manifest %>%
        mutate(temp=subject_id,
               subject_id=subject_id2,
               subject_id2=temp) %>%
        select(-temp)
}
##
## Use the tumor_type in the clinical.data
##
manifest4 <- update_ids(manifest3) %>%
    select(-tumor_type) %>%
    left_join(select(cdat2, lab_id, tumor_type,
                     discordant_tumor_type),
              by="lab_id")
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
stopifnot(!any(duplicated(manifest5$lab_id)))
##
## Attach ids from FACETS copy number analysis
##
manifest6 <- filter(manifest5, tumor.normal=="tumor") %>%
    select(subject_id, lab_id, pgdx_id)
facets <- read_tsv(here("output", "facets",
                        "merge-facets-tables.R",
                        "all-segments.txt"),
                   show_col_types=FALSE) %>%
    select(Sample) %>%
    distinct()
facets.matched <- inner_join(facets, manifest6,
                             by=c("Sample"="pgdx_id")) %>%
    rename(facet_id=Sample)
facets.nomatch <- filter(facets, !Sample %in% manifest6$pgdx_id) %>%
    mutate(tmpid=str_replace_all(Sample, "t_", ""),
           tmpid=str_replace_all(tmpid, "_eland", ""),
           tmpid=str_replace_all(tmpid, "_hg18_", ""),
           tmpid=str_replace_all(tmpid, "_ExA", "_Ex"))
facets.matched2 <- inner_join(facets.nomatch, manifest6, by=c("tmpid"="pgdx_id")) %>%
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
facets.matched3 <- inner_join(facets.nomatch2, manifest6, by=c("uid2"="pgdx_id")) %>%
    rename(facet_id=Sample)
facets2 <- bind_rows(facets.matched,
                     facets.matched2,
                     facets.matched3) %>%
    select(facet_id, subject_id, lab_id)
manifest7 <- left_join(manifest5, facets2, by=c("subject_id", "lab_id"))
## The CGOV177 samples were mapped back to the same pgdx id
ix <- grep("177T", manifest7$lab_id)
A <- manifest7[ix, ]
## drop the second instance
manifest7 <- manifest7[-ix[2], ]
##
## Add FACETS id for WGS samples to facilitate merging FACET processed data
##
facets <- tibble(Sample=fs::dir_ls(here("output", "facets-trellis/jhpce_directories"),
                                   type="directory")) %>%
    mutate(Sample=basename(Sample))
##facets <- read_tsv(here("output", "facets-trellis",
##                        "combined-deletion-table.txt"),
##                   show_col_types=FALSE) %>%
##    select(Sample) %>%
##    distinct()
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
discordant <- filter(manifest8, discordant_tumor_type)
save(discordant, file=here("Rpackage",
                           "ovarian.subtypes",
                           "data", "discordant.rda"))
manifest <- filter(manifest8, !discordant_tumor_type) %>%
    filter(subject_id != "CGOV359") %>%
    mutate(lab_id=ifelse(lab_id == "CGOV151Tb_WES", "CGOV151Tb", lab_id))
save(manifest, file=here("Rpackage",
                         "ovarian.subtypes",
                         "data", "manifest.rda"))
clinical <- filter(cdat2, lab_id %in% manifest$lab_id)
save(clinical, file=here("Rpackage",
                         "ovarian.subtypes",
                         "data", "clinical.rda"))
stopifnot(all(clinical$lab_id %in% manifest$lab_id))
stopifnot(all(clinical$subject_id %in% manifest$subject_id))
stopifnot(all(manifest$subject_id %in% clinical$subject_id))
tmp <- filter(manifest, lab_id=="CGOV177T_2")
stopifnot(nrow(tmp)==1)
