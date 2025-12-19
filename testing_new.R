library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(tidyverse)
library(dplyr)
library(DOSE)

#setting the working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#call the vizualisation function
source("cyemapplot.R")

#get deg data 
data(geneList, package="DOSE")

##deg data
degs_data<-as.data.frame(geneList)
degs_data$ID <- rownames(degs_data)
colnames(degs_data)<- c("log2FC", "ID")

##genes
sig.deg<- subset(degs_data, abs(log2FC) > 1 )

#enrichment analysis
edo_go <- enrichGO(sig.deg$ID, 
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   ont = "BP",
                   pvalueCutoff = 0.01,
                   minGSSize = 50,
                   universe = names(geneList))

edo.ora.sim <- pairwise_termsim(edo_go)


#gsea 
ego_go<-gseGO(geneList,
              OrgDb = "org.Hs.eg.db",
              ont = "BP",
              pvalueCutoff = 0.01,
              minGSSize = 50)

#pairwise 
edo.gsea.sim <- pairwise_termsim(ego_go)
y <- as.data.frame(ego_go)
source("cyemapplot.R")
cyemapplot(edo.ora.sim, show_category = 15, min_edge = 0.8, visualization="deg", degs_data = sig.deg)
cyemapplot(edo.gsea.sim, show_category = 15, min_edge = 0.8, visualization="deg", degs_data = sig.deg)

