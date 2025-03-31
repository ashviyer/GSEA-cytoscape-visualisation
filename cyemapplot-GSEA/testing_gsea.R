library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(tidyverse)
library(dplyr)

setwd("C:/Users/aishw/Documents/git-repositories/GSEA-cytoscape-visualisation/cyemapplot-GSEA")

source("cyemapplot_gsea.R")

library(DOSE)
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

# OPTION 1 - no gene expression data
edo.sim <- pairwise_termsim(ego_go)

cy.emapplot(edo.sim, analysis.name = "emapplot-nodata", min_edge = 0.5)

# OPTION 2 - with gene expression data (differentially expressed genes)
cy.emapplot(edo.sim, analysis.name = "emapplot-withdata",show_category=93, degs_data = degs_data)

##Adding the pie chart to pathway nodes to visualise the down-regulated, not -changed and the up-regulated genes.
setNodeCustomPieChart(c("down","up"), colors = c("#67A9CF","#EF8A62"), style.name = "cyemapplot.gsea.data" )
setEdgeColorDefault("#D7DBDD", style.name = "cyemapplot.gsea.data")
#to do 
#1.documentation
#2. test on gsea output
#3. r.project for deg and enrichment output or use th eones from clusterprofiler
