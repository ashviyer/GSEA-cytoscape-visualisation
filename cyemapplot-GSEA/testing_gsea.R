library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(tidyverse)
library(dplyr)
library(DOSE)

#setting the working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#call the gsea vizualisation function
source("cyemapplot_gsea.R")

#get deg data 
data(geneList, package="DOSE")

##deg data
degs_data<-as.data.frame(geneList)
degs_data$ID <- rownames(degs_data)
colnames(degs_data)<- c("log2FC", "ID")

#gsea 
ego_go<-gseGO(geneList,
              OrgDb = "org.Hs.eg.db",
              ont = "BP",
              pvalueCutoff = 0.01,
              maxGSSize = 100,
              minGSSize = 50)

#pairwise 
edo.sim <- pairwise_termsim(ego_go)

# OPTION 2 - with gene expression data (differentially expressed genes)
cy.emapplot(edo.sim, analysis.name = "emapplot-gsea-withdata",show_category=89, degs_data = degs_data)

