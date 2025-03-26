##' emapplot visualization in Cytoscape
##'
##' This function creates an emapplot in Cytoscape
##' @title cy.emapplot
##' @param gse output of enrichment (pairwise_termsim result)
##' @param show_category number of enriched terms to display
##' @param min_edge minimum percentage of overlap genes to display the edge, should between 0 and 1, default value is 0.5
##' @param degs_data differentially expressed genes with ID, log2FC and p-value columns
##' @param ... additional parameters, see also the parameters supported by the enricher() function
##' @return Location of Cytoscape session file
##' @export
##' @author ashviyer, mkutmon
##' 
cy.emapplot <- function(gse, analysis.name = "cy_emapplot", show_category=30, min_edge=0.2, degs_data = NULL,...) {
  res.gs.df <- as.data.frame(gse) 
  gs.info <- gs.info.basic(gse, res.gs.df, show_category)
  if(!is.null(degs_data)) { 
    gs.info <- gs.add.degs(gse, gs.info, degs_data) 
  }
  
  gse.sim.filt <- gse@termsim
  colnames(gse.sim.filt) <- res.gs.df$ID
  rownames(gse.sim.filt) <- res.gs.df$ID
  gse.sim.filt[gse.sim.filt < min_edge] <- 0
  
  graph <- igraph::graph_from_adjacency_matrix(gse.sim.filt, weighted=TRUE, mode="undirected", diag = FALSE)
  RCy3::createNetworkFromIgraph(graph, title = paste0(analysis.name, " - emapplot"), collection = analysis.name)
  RCy3::loadTableData(gs.info, data.key.column = "ID", table.key.column = "id")
  style.network(gs.info, !is.null(degs_data))
  
  output.file <- paste0(analysis.name,".cys")
  RCy3::saveSession(output.file)
  return(output.file)
}

gs.add.degs <- function(gse, gs.info, degs_data) {
  degs_data$condition <- 0
  degs_data$condition[degs_data$log2FC > 0] <- 1
  degs_data$condition[degs_data$log2FC < 0] <- 2
  degs_data$condition[is.na(degs_data$condition)]<- 0
  
  gs.info<-merge(pathway_data, psvdvshnl, by.x = "genes", by.y = "names", all.x = TRUE)
  return (gs.info)
}

gs.info.basic <- function(gse, res.gs.df, show_category) {
  if(nrow(res.gs.df) > show_category) {
    enriched.pathways <- res.gs.df[c(1:show_category),]
  } else {
    enriched.pathways <- res.gs.df
  }
  complete=list()
  for (pathways in res.gs.df$ID){
    complete[[pathways]]<-data.frame(genes = edo.sim@geneSets[[pathways]],
                                     pathway_id = rep(c(pathways), each = nrow(complete[[pathways]])))
    
  }
  ##Create a dataframe with pathways and genesets
  pathway_data <- bind_rows(complete)
  enriched.pathways <- enriched.pathways[,c(1,2,3,4,5,6)]
  enriched.pathways <- enriched.pathways %>%
    separate(GeneRatio, into = c("genesPathway", "genesAllPathways"), sep = "/") %>%
    mutate(genesPathway = as.numeric(genesPathway), genesAllPathways = as.numeric(genesAllPathways))
  enriched.pathways <- enriched.pathways %>%
    separate(BgRatio, into = c("bgPathway", "bgAllPathways"), sep = "/") %>%
    mutate(bgPathway = as.numeric(bgPathway), bgAllPathways = as.numeric(bgAllPathways))
  
  enriched.pathways$otherGenes <- enriched.pathways$bgPathway - enriched.pathways$genesPathway
  #enriched.pathways <- enriched.pathways[,c(1,4,2,8)]
  return(enriched.pathways)
}

style.network <- function(gs.info, withdata = FALSE) {
  RCy3::importVisualStyles("styles.xml")
  if(withdata) { style.name <- "cy.emapplot.data" } else { style.name <- "cy.emapplot.nodata" }
  RCy3::setVisualStyle(style.name)
}