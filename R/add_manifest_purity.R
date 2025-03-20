load("data/manifest.rda")
source("R/functions.R")
library(tidyverse)
facets.file1 <- file.path("..", "output", "facets",
                          "merge-facets-tables.R",
                          "summary-stats.txt")
facets.file2 <- file.path("..", "output", "facets-trellis", "summary-stats.txt")
facets_purity <- read_facets2(facets.file1, facets.file2)
manifest9 <- add_facets_purity(manifest, facets_purity)
