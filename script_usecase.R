library(clusterProfiler)
library(enrichplot)
library(RCy3)
library(tidyverse)
library(dplyr)
library(DOSE)
library(igraph)
library(renv)

##create output directory
outdir <- "Output-plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


#setting the working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#call the vizualisation function
source("cyemapplot.R")

#import usecase deg data 
ibm.data<- read.delim("res_sIBM_AMP_main_updated.txt", sep = " ")
ibm.data$ID<- rownames(ibm.data)
ibm.data <- ibm.data %>%
  separate(ID, into = c("ensembl_id", "symbol"), sep = ";", remove = FALSE)
ibm.data$ensembl_id <- sub("\\..*$", "", ibm.data$ensembl_id)
ibm.data<- ibm.data[complete.cases(ibm.data),]

ibm.data<- subset(ibm.data, select = c("ensembl_id", "log2FoldChange", "pvalue", "padj"))
rownames(ibm.data)<- ibm.data$ensembl_id
colnames(ibm.data)<- c("ID", "log2FC", "pvalue", "padj")

##Over -representation analysis
##select significant genes
sig.deg<- subset(ibm.data, pvalue < 0.05 & abs(log2FC) > 1 )

##basic visualisation using the Disease ontology database
#convert ids from ENSEMBL to ENTREZID
input_genes = bitr(sig.deg$ID, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
background_genes = bitr(rownames(ibm.data), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ora_do <- enrichDO(input_genes$ENTREZID, ont = "HDO",
                   organism = "hsa",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   universe = background_genes$ENTREZID,
                   maxGSSize = 500,
                   qvalueCutoff = 0.05,
                   readable = FALSE
)

#pairwise term calculculation
do.ora.sim <- pairwise_termsim(ora_do, showCategory = 328)

#Advanced cytoscape visualisation
cyemapplot(do.ora.sim, 
           show_category = 328, 
           min_edge = 0.4, 
           visualization="basic", 
           ig_layout = igraph::layout_with_kk, 
           layout_scale = 800, 
           filt_components = 8, 
           plot_components = TRUE, 
           top_components = 5, 
           analysis_name = "Disease-enrichment")


##pie chart visualisation for pathway database
ora_wp<- enrichWP(gene = input_genes$ENTREZID,
                  organism = "Homo sapiens",
                  maxGSSize = 500,
                  pAdjustMethod = "BH",
                  universe = background_genes$ENTREZID)

#pairwise term calculculation
wp.ora.sim <- pairwise_termsim(ora_wp, showCategory = 76)

#Advanced cytoscape visualisation
cyemapplot(wp.ora.sim, 
           show_category = 76, 
           min_edge = 0.2, 
           visualization="pie", 
           ig_layout = igraph::layout_with_kk, 
           layout_scale = 800, 
           filt_components = 1, 
           plot_components = TRUE, 
           top_components = 5, 
           analysis_name = "WikiPathway-Enrichment")

##deg visualisation using gsea using GO:Bp database
#ranking on product of logFC*-log10(adjusted p-value)
ibm.data$ranking<-ibm.data$log2FC*-log10(ibm.data$pvalue)

#select the ranking: logFC*-log10(adjusted p-value)
geneList= ibm.data$ranking
#add entrex gene id as gene names
names(geneList) = as.character(ibm.data$ID)
#sort on decreasing order
geneList = sort(geneList,decreasing = TRUE)
#gsea 
gsea_go_bp<-gseGO(geneList,
              OrgDb = "org.Hs.eg.db",
              keyType = "ENSEMBL",
              ont = "BP",
              pvalueCutoff = 0.01,
              pAdjustMethod = "fdr"
)

#pairwise 
gsea.go.bp.sim <- pairwise_termsim(gsea_go_bp, showCategory = 464)

##gsea visualisation
cyemapplot(gsea.go.bp.sim, 
           show_category = 464, 
           min_edge = 0.4, 
           visualization="deg", 
           degs_data = ibm.data, 
           ig_layout = igraph::layout_with_kk, 
           layout_scale = 800, 
           filt_components = 3, 
           plot_components = TRUE, 
           top_components = 10, 
           analysis_name = "GeneOntology-Enrichment")
