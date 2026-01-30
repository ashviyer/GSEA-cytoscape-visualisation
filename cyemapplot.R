cyemapplot <- function(ea_sim, analysis_name = "Enrichment", show_category=30, min_edge=0.2,
                       degs_data = NULL, visualization = "basic",
                       ig_layout = igraph::layout_with_kk, layout_scale = 500,
                       filt_components = 1,
                       plot_components = TRUE,          # <--- new
                       top_components = 5,   # <--- add this (min component size)
                       ...) {
  
  ea.df <- as.data.frame(ea_sim)
  ea.df.filt <- if (nrow(ea.df) > show_category) ea.df[1:show_category, ] else ea.df
  
  is_gsea <- inherits(ea_sim, "gseaResult")
  if (is_gsea) {
    gs.info <- gsea.info.basic(ea_sim, ea.df.filt)
    if (analysis_name == "Enrichment") analysis_name <- "Enrichment_GSEA"
  } else if (inherits(ea_sim, "enrichResult")) {
    gs.info <- ora.info.basic(ea.df.filt)
    if (analysis_name == "Enrichment") analysis_name <- "Enrichment_ORA"
  }
  
  ea_sim_matrix <- filter.sim.matrix(ea_sim, ea.df.filt, min_edge)
  
  graph <- igraph::graph_from_adjacency_matrix(ea_sim_matrix, weighted=TRUE, mode="undirected", diag=FALSE)
  
  # ---- filter components by size ----
  comp <- igraph::components(graph)
  keep_comps <- which(comp$csize >= filt_components)
  
  if (length(keep_comps) == 0) {
    stop("No components with size >= filt_components = ", filt_components)
  }
  
  nodes_to_keep <- which(comp$membership %in% keep_comps)
  g_filtered <- igraph::induced_subgraph(graph, nodes_to_keep)
  # -----------------------------------
  
  viz_label <- switch(visualization, basic="basic", pie="pie", deg="deg", visualization)
  net_title <- paste0("cyemapplot-", viz_label, "-main")
  RCy3::createNetworkFromIgraph(g_filtered, title = net_title, collection = analysis_name)
  
  # layout main
  xy <- ig_layout(g_filtered)
  rescale_to <- function(v, to = c(-layout_scale, layout_scale)) {
    r <- range(v, finite = TRUE)
    if (diff(r) == 0) return(rep(mean(to), length(v)))
    (v - r[1]) / diff(r) * diff(to) + to[1]
  }
  nodes <- igraph::V(g_filtered)$name
  RCy3::setNodePositionBypass(node.names = nodes,
                              new.x.locations = rescale_to(xy[,1]),
                              new.y.locations = rescale_to(xy[,2]))
  RCy3::fitContent()
  
  # load attributes + visualize main
  node.cols <- RCy3::getTableColumnNames("node")
  key <- if ("name" %in% node.cols) "name" else if ("shared name" %in% node.cols) "shared name" else stop("No node key column found")
  RCy3::loadTableData(gs.info, data.key.column = "Description", table.key.column = key)
  
  if (visualization == "basic") basic_viz(is_gsea = is_gsea)
  else if (visualization == "pie") pie_chart_viz()
  else if (visualization == "deg") deg_viz(degs_data, gs.info)
  
  # ---- COMPONENT SUBNETWORKS ----
  if (plot_components) {
    comp2 <- igraph::components(g_filtered)
    
    # order components by size (largest first)
    ord <- order(comp2$csize, decreasing = TRUE)
    ord <- ord[seq_len(min(top_components, length(ord)))]
    
    for (i in seq_along(ord)) {
      comp_id <- ord[i]
      v_idx <- which(comp2$membership == comp_id)
      g_sub <- igraph::induced_subgraph(g_filtered, v_idx)
      
      sub_title <- paste0("cyemapplot-", viz_label, "-top", i, "cluster")
      RCy3::createNetworkFromIgraph(g_sub, title = sub_title, collection = analysis_name)
      
      # layout subnetwork
      xy2 <- ig_layout(g_sub)
      nodes2 <- igraph::V(g_sub)$name
      RCy3::setNodePositionBypass(node.names = nodes2,
                                  new.x.locations = rescale_to(xy2[,1]),
                                  new.y.locations = rescale_to(xy2[,2]))
      RCy3::fitContent()
      
      # load attributes + apply same visualization to subnetwork
      node.cols2 <- RCy3::getTableColumnNames("node")
      key2 <- if ("name" %in% node.cols2) "name" else if ("shared name" %in% node.cols2) "shared name" else key
      RCy3::loadTableData(gs.info, data.key.column = "Description", table.key.column = key2)
      
      if (visualization == "basic") basic_viz(is_gsea = is_gsea)
      else if (visualization == "pie") pie_chart_viz()
      else if (visualization == "deg") deg_viz(degs_data, gs.info)
    }
  }
}

basic_viz <- function(is_gsea = FALSE) {
  if(!("cyemapplot_basic" %in% RCy3::getVisualStyleNames())) {
    RCy3::copyVisualStyle(from.style = "default", to.style = "cyemapplot_basic")
    RCy3::lockNodeDimensions(TRUE, "cyemapplot_basic")
    RCy3::setNodeSizeMapping(table.column = "setSize", table.column.values = c(0, 300), sizes = c(10,60), mapping.type = "c", style.name = "cyemapplot_basic")
    RCy3::setNodeShapeDefault("ELLIPSE", "cyemapplot_basic")
    #RCy3::setNodeBorderWidthDefault(3, "cyemapplot_basic")
    RCy3::setNodeColorDefault("#F0F0F0", "cyemapplot_basic")
    RCy3::setNodeLabelPositionDefault(new.nodeAnchor = "S",new.graphicAnchor = "N", new.justification = "c", new.xOffset = 0, new.yOffset = 0, style.name = "cyemapplot_basic")
    RCy3::setEdgeColorDefault("#969696", style.name = "cyemapplot_basic")
    
    node.cols <- RCy3::getTableColumnNames("node")
    # Only for GSEA (i.e., only when NES_cat exists)
    if (is_gsea && ("NES_cat" %in% node.cols)) {
      RCy3::setNodeBorderColorMapping(
        "NES_cat",
        table.column.values = c("up", "down"),
        colors = c("#D6604D", "#4393C3"),
        mapping.type = "d",
        default.color = "#CCCCCC",
        style.name = "cyemapplot_basic"
      )
      RCy3::setNodeBorderWidthDefault(9, "cyemapplot_basic") 
    } else {
      RCy3::setNodeBorderColorDefault("#CCCCCC", style.name = "cyemapplot_basic")
      RCy3::setNodeBorderWidthDefault(3, "cyemapplot_basic")
    }
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
  
  node.cols <- RCy3::getTableColumnNames("node")
  key <- if ("name" %in% node.cols) "name" else if ("shared name" %in% node.cols) "shared name" else stop("No node key column found")
  
  RCy3::loadTableData(gs.info, data.key.column = "Description", table.key.column = key)
  
  
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
    separate(GeneRatio, into = c("query.genes.in.term", "query.genes"), sep = "/") %>%
    mutate(query.genes.in.term = as.integer(query.genes.in.term), query.genes = as.integer(query.genes))
  ea.df.info <- ea.df.info %>%
    separate(BgRatio, into = c("setSize", "background.genes"), sep = "/") %>%
    mutate(setSize = as.integer(setSize), background.genes = as.integer(background.genes))
  
  ea.df.info$rest <- ea.df.info$setSize - ea.df.info$query.genes.in.term
  ea.df.info$input <- ea.df.info$query.genes.in.term ## double check this
  
  
  ea.df.info <- ea.df.info %>% select(-query.genes, -background.genes)
  return(ea.df.info)
}

gsea.info.basic <- function(ea_sim, ea.df.filt) {
  ea.df.info <- ea.df.filt[, c("ID","Description","setSize","NES","pvalue","p.adjust","core_enrichment")]
  
  ea.df.info$NES_cat <- ifelse(ea.df.info$NES > 0, "up", "down")
  ea.df.info$input <- stringr::str_count(ea.df.info$core_enrichment, "/") + 1
  ea.df.info$geneID_leading_edge <- ea.df.info$core_enrichment
  
  ea.df.info$geneID <- vapply(ea.df.info$ID, function(id) {
    gs <- ea_sim@geneSets[[id]]
    if (is.null(gs)) NA_character_ else paste(gs, collapse = "/")
  }, character(1))
  
  ea.df.info$setSize <- ea.df.info$setSize
  ea.df.info$rest <- ea.df.info$setSize - ea.df.info$input
  ea.df.info
}