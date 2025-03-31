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
  up <- degs_data[degs_data$log2FC > 0,]
  down <- degs_data[degs_data$log2FC < 0,]
  
  # Extract all gene sets from the enrichment result object
  gene_sets <- gse@geneSets
  # Initialize a list to store the counts
  gene_counts_list_up <- list()
  # Initialize a list to store the counts
  gene_counts_list_down <- list()
  
  # Iterate over each gene set and count how many genes from my_genes are in the set
  for (set_name in names(gene_sets)) {
    genes_in_set <- gene_sets[[set_name]]
    gene_counts_list_up[[set_name]] <- sum(genes_in_set %in% up$ID)
    gene_counts_list_down[[set_name]] <- sum(genes_in_set %in% down$ID)
  }
  
  # Convert the list to a data frame for easier viewing
  gene_counts_df_up <- data.frame(
    ID = names(gene_counts_list_up),
    up = unlist(gene_counts_list_up)
  )
  
  # Convert the list to a data frame for easier viewing
  gene_counts_df_down <- data.frame(
    ID = names(gene_counts_list_down),
    down = unlist(gene_counts_list_down)
  )
  
  gs.info <- merge(gs.info, gene_counts_df_up, by = "ID", keep.x = TRUE)
  gs.info <- merge(gs.info, gene_counts_df_down, by = "ID", keep.x = TRUE)
  return (gs.info)
}

gs.info.basic <- function(gse, res.gs.df, show_category) {
  if(nrow(res.gs.df) > show_category) {
    enriched.pathways <- res.gs.df[c(1:show_category),]
  } else {
    enriched.pathways <- res.gs.df
  }
  enriched.pathways <- enriched.pathways[,c(1,2,3,4,5,6,7)]
  return(enriched.pathways)
}

style.network <- function(gs.info, withdata = FALSE) {
  RCy3::importVisualStyles("styles.xml")
  if(withdata) { style.name <- "cyemapplot.gsea.data" } else { style.name <- "cyemapplot.gsea.data" }
  RCy3::setVisualStyle(style.name)
}