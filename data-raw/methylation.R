library(here)
library(tidyverse)
load(here("Rpackage", "ovarian.subtypes", "data", "manifest.rda"))
load(here("Rpackage", "ovarian.subtypes", "data", "discordant.rda"))
mdir <- here("code", "methylation.Rmd", "se.rds") %>%
    readRDS()
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
save(methylation, file=here("Rpackage",
                       "ovarian.subtypes",
                       "data",
                       "methylation.rda"))
