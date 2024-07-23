library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(tidyverse)
library(dplyr)

source("cyemapplot.R")


# OPTION 1 - no gene expression data
edo.sim <- pairwise_termsim(edo)
cy.emapplot(edo.sim, analysis.name = "emapplot-nodata")

# OPTION 2 - with gene expression data (differentially expressed genes)
cy.emapplot(edo.sim, analysis.name = "emapplot-withdata", degs_data = degs)

