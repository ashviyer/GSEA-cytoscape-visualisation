library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(tidyverse)
library(dplyr)
library(DOSE)

#setting the working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#call the ora vizualisation function
source("cyemapplot_ora.R")

#get deg data
data(geneList, package="DOSE")

##deg data
degs_data<-as.data.frame(geneList)
degs_data$ID <- rownames(degs_data)
colnames(degs_data)<- c("log2FC", "ID")

##genes
sig.deg<- subset(degs_data, abs(log2FC) > 3 )


#enrichment analysis
edo_go <- enrichGO(sig.deg$ID, 
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   ont = "BP",
                   pvalueCutoff = 0.01,
                   minGSSize = 50,
                   universe = names(geneList))


# OPTION 1 - no gene expression data
edo.sim <- pairwise_termsim(edo_go)
cy.emapplot(edo.sim, analysis.name = "emapplot-ora-nodata", min_edge = 0.5)

# OPTION 2 - with gene expression data (differentially expressed genes)
cy.emapplot(edo.sim, analysis.name = "emapplot-ora-withdata", degs_data = degs_data)



