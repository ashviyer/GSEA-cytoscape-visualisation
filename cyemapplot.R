cyemapplot <- function(ea_sim, analysis_name = "Enrichment", show_category=30, min_edge=0.2, degs_data = NULL, visualization = "basic", ...) {
  ea.df <- as.data.frame(ea_sim) 
  if(nrow(ea.df) > show_category) {
    ea.df.filt <- ea.df[c(1:show_category),]
  } else {
    ea.df.filt <- ea.df
  }

  if(inherits(ea_sim, "gseaResult")) {
    gs.info <- gsea.info.basic(ea.df.filt)
    method <- "GSEA"
    if(analysis_name == "Enrichment") analysis_name <- "Enrichment_GSEA"
  } else if(inherits(ea_sim, "enrichResult")) {
    gs.info <- ora.info.basic(ea.df.filt)
    method <- "ORA"
    if(analysis_name == "Enrichment") analysis_name <- "Enrichment_ORA"
  }
  
  ea_sim_matrix <- filter.sim.matrix(ea_sim, ea.df.filt, min_edge)

  graph <- igraph::graph_from_adjacency_matrix(ea_sim_matrix, weighted=TRUE, mode="undirected", diag = FALSE)
  RCy3::createNetworkFromIgraph(graph, title = paste0(analysis_name, " - cyemapplot"), collection = analysis_name)
  RCy3::loadTableData(gs.info, data.key.column = "Description", table.key.column = "id")
  if(visualization == "basic") { basic_viz() }
  else if (visualization == "pie") { pie_chart_viz() }
  else if (visualization == "deg") { deg_viz(degs_data, gs.info) } # TODO: check if deg_data is added for deg visualization
}

basic_viz <- function() {
  if(!("cyemapplot_basic" %in% RCy3::getVisualStyleNames())) {
    RCy3::copyVisualStyle(from.style = "default", to.style = "cyemapplot_basic")
    RCy3::lockNodeDimensions(TRUE, "cyemapplot_basic")
    RCy3::setNodeSizeMapping(table.column = "setSize", table.column.values = c(0, 300), sizes = c(10,60), mapping.type = "c", style.name = "cyemapplot_basic")
    RCy3::setNodeShapeDefault("ELLIPSE", "cyemapplot_basic")
    RCy3::setNodeBorderWidthDefault(3, "cyemapplot_basic")
    RCy3::setNodeBorderColorMapping("NES_cat",table.column.values = c("up", "down"), colors = c("#D6604D","#4393C3"), mapping.type = "d", default.color = "#CCCCCC", style.name = "cyemapplot_basic")
    RCy3::setNodeColorDefault("white", "cyemapplot_basic")
    RCy3::setNodeLabelPositionDefault(new.nodeAnchor = "S",new.graphicAnchor = "N", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = "cyemapplot_basic")
    RCy3::setEdgeColorDefault("#CCCCCC", style.name = "cyemapplot_basic")
  }
  RCy3::setVisualStyle("cyemapplot_basic")
}

pie_chart_viz <- function() {
  if(!("cyemapplot_basic" %in% RCy3::getVisualStyleNames())) {
    basic_viz()
  }
  if(!("cyemapplot_piechart" %in% RCy3::getVisualStyleNames())) {
    RCy3::copyVisualStyle(from.style = "cyemapplot_basic", to.style = "cyemapplot_piechart")
    RCy3::setNodeCustomPieChart(columns = c("rest", "input"), colors = c("#FFFFFF","#998EC3"), slot = 1, startAngle = 90, style.name = "cyemapplot_piechart")
  }
  RCy3::setVisualStyle("cyemapplot_piechart")
}

deg_viz <- function(degs_data, gs.info) {
  up <- degs_data[degs_data$log2FC > 0,]
  down <- degs_data[degs_data$log2FC < 0,]
  
  gs.info$n_up <- sapply(strsplit(gs.info$geneID, "/"), function(genes) {
    sum(genes %in% up$ID)
  })
  
  # Add column for number of downregulated genes
  gs.info$n_down <- sapply(strsplit(gs.info$geneID, "/"), function(genes) {
    sum(genes %in% down$ID)
  })
  
  RCy3::loadTableData(gs.info, data.key.column = "ID", table.key.column = "id")
  
  
  if(!("cyemapplot_basic" %in% RCy3::getVisualStyleNames())) {
    basic_viz()
  }
  if(!("cyemapplot_degs" %in% RCy3::getVisualStyleNames())) {
    #TODO: fix for GSEA - should not use core_enrichment but all genes in the pathway
    RCy3::copyVisualStyle(from.style = "cyemapplot_basic", to.style = "cyemapplot_degs")
    RCy3::setNodeCustomPieChart(columns = c("rest", "n_up","n_down"), colors = c("#FFFFFF","#D6604D","#4393C3"), slot = 2, startAngle = 90, style.name = "cyemapplot_degs")
  }
  RCy3::setVisualStyle("cyemapplot_degs")
}

filter.sim.matrix <- function(ea_sim, ea.df.filt, min_edge) {
  ea_sim_matrix <- ea_sim@termsim
  ea_sim_matrix.filt <- ea_sim_matrix[rownames(ea_sim_matrix) %in% ea.df.filt$Description, 
                                      colnames(ea_sim_matrix) %in% ea.df.filt$Description]
  
  ea_sim_matrix.filt[ea_sim_matrix.filt < min_edge] <- 0
  return(ea_sim_matrix.filt)
}

ora.info.basic <- function(ea.df.filt) {
  ea.df.info <- ea.df.filt[,c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","geneID")]

  ea.df.info <- ea.df.info %>%
    separate(GeneRatio, into = c("degs.pathway", "degs"), sep = "/") %>%
    mutate(degs.pathway = as.integer(degs.pathway), degs = as.integer(degs))
  ea.df.info <- ea.df.info %>%
    separate(BgRatio, into = c("setSize", "genes"), sep = "/") %>%
    mutate(setSize = as.integer(setSize), genes = as.integer(genes))
  
  ea.df.info$notdegs.pathway <- ea.df.info$setSize - ea.df.info$degs.pathway
  ea.df.info$input <- ea.df.info$degs.pathway
  ea.df.info$rest <- ea.df.info$notdegs.pathway
  
  ea.df.info <- ea.df.info %>% select(-degs, -genes)
  return(ea.df.info)
}

gsea.info.basic <- function(ea.df.filt, show_category) {
  ea.df.info <- ea.df.filt[,c("ID","Description","setSize","NES","pvalue","p.adjust","core_enrichment")]
  ea.df.info$NES_cat <- ifelse(ea.df.info$NES > 0, "up", "down")
  ea.df.info$input <- str_count(ea.df.info$core_enrichment, "/") + 1 
  ea.df.info$geneID <- ea.df.info$core_enrichment
  ea.df.info$rest <- ea.df.info$setSize - ea.df.info$input
  return(ea.df.info)
}