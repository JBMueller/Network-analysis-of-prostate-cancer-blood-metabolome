library(readxl) # version: 1.4.3
library(RCy3) # version: 2.24.0
library(corrr) # version: 0.4.4
library(ggplot2) # version: 3.5.1
library(igraph) # version: 2.0.3
library(rstudioapi) # version: 0.16.0
library(ComplexHeatmap) # version: 2.20.0
library(RColorBrewer) # version: 1.1.3

# setting working directory to current file location
wd <- getActiveDocumentContext()$path
setwd(dirname(wd))

# import metabolomics measurements and meta data
met_matrix <- read.csv("Data/MSB1014_metabolomics_preprocessed.csv", header = T, row.names = 1)
colnames(met_matrix) <- gsub("^X", "", colnames(met_matrix))
meta_data <- read_xlsx("Data/EDRN_ClinicalData.xlsx", sheet = 2)

# identifying samples without metabolite measurements
unmeasured_sample_id <- meta_data$SAMPLE_ID[!meta_data$SAMPLE_ID %in% rownames(met_matrix)]

# removing samples without metabolite measurements from meta_data
meta_data <- meta_data[!meta_data$SAMPLE_ID %in% unmeasured_sample_id,]

# subset metabolite measurements into patients and control
if (all(rownames(met_matrix) == meta_data$SAMPLE_ID)){
  met_cancer <- met_matrix[meta_data$Diagnosis == 1, ]
  met_control <- met_matrix[meta_data$Diagnosis == 0, ]
} else {
  print("metabolomics samples and meta_data sample IDs are not in the same order!")
}
# create a function to quickly infer networks from metabolite data
make_network <- function(met_matrix, group, corr_cutoff, component_size){
# calculating pearson correlation of metabolites
cancer.cor <- as.data.frame(correlate(met_matrix, diagonal = 0, method = "pearson"))
row.names(cancer.cor) <- cancer.cor$term
cancer.cor[1] <- NULL

is_symmetric <- all(cancer.cor == t(cancer.cor))

cancer.cor.vector <- as.matrix(cancer.cor)
cancer.cor.vector <- as.vector(cancer.cor.vector)
cancer.cor.filtered <- cancer.cor
cancer.cor.filtered[cancer.cor.filtered < corr_cutoff] <- 0
is_symmetric <- all(cancer.cor.filtered == t(cancer.cor.filtered))


graph <- igraph::graph_from_adjacency_matrix(as.matrix(cancer.cor.filtered), weighted=TRUE, mode = "undirected")
# delete unconnected nodes
isolated <- which(degree(graph)==0)
graph <- delete_vertices(graph, isolated)

# identify number of components and remove small components
num_components <- components(graph, mode = 'weak')
remove <- vector()
for (i in 1:length(num_components$membership)){
  curr_node <- num_components$membership[i]
  if (num_components$csize[curr_node] < component_size) {
    remove[length(remove) + 1] <- i
  }
}
graph <- delete_vertices(graph, remove)
network_name <- paste("correlation_network", group)
createNetworkFromIgraph(graph,title=network_name,collection=group)
return(graph)
}

# calculate the log2FC to visualize on the networks
table <- data.frame(name = character(ncol(met_matrix)))
table$name <- colnames(met_matrix)
exponent_met_cancer <- exp(met_cancer)
exponent_met_control <- exp(met_control)
table$mean_abundance_cancer <- apply(exponent_met_cancer, 2, FUN = mean)
table$mean_abundance_control <- apply(exponent_met_control, 2, FUN = mean)
table$log2FC <- log(table$mean_abundance_cancer/table$mean_abundance_control, base = 2)

# create the correlation networks for non-cancer samples
control_network <- make_network(met_control, "control", corr_cutoff = 0.6, component_size = 10)
loadTableData(table, data.key.column="name")

# create a style (after loading the first network because the table data is needed)
style.name = "log2FC mapping"
defaults <- list(NODE_SHAPE="ellipse",
                 NODE_SIZE=40,
                 EDGE_TRANSPARENCY=120,
                 NODE_BORDER_WIDTH = 2,
                 NODE_LABEL_POSITION="N,S,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','log2FC','c',c(-0.3,0.0,0.3), c("blue","white","red"))
createVisualStyle(style.name, defaults, list(nodeFills, nodeLabels))
setVisualStyle(style.name)

# create the correlation networks for cancer samples
cancer_network <- make_network(met_cancer, "cancer", corr_cutoff = 0.6, component_size = 10)
loadTableData(table, data.key.column="name")
setVisualStyle(style.name)

# find the overlap of nodes between both networks
overlap <- intersect(V(control_network)$name, V(cancer_network)$name)

# highlight the nodes in the cytoscape view
networkSuid = getNetworkSuid()
clearSelection()
selectNodes(network=networkSuid, overlap, by.col="id")

# create a function to extract and visualize all components
visualize_components <- function(graph, group, visualize){
  subgraphs_list <- list()
  
  components <- components(graph, mode = 'weak')
  for (component in 1:components$no) {
    subgraph_name <- paste0("subgraph_", component)
    subgraphs_list[[subgraph_name]] <- subgraph(graph, which(components$membership == component))
    curr.subgraph <- subgraphs_list[[subgraph_name]]
    if (visualize == TRUE) {
      createNetworkFromIgraph(curr.subgraph ,title=paste(subgraph_name, group),collection=group)
      loadTableData(table, data.key.column="name")
      setVisualStyle(style.name)
    }
  }
  return(subgraphs_list)
}

# extract the individual components from the non-cancer network
control_components <- visualize_components(control_network, "control", visualize = F) # set visualize to TRUE if components should be visualized in cytoscape

# extract the individual components from the cancer network
cancer_components <- visualize_components(cancer_network, "cancer", visualize = F) # set visualize to TRUE if components should be visualized in cytoscape

# performing GN community detection on large components

GN_clustering <- function (g, cutoff) {
  GN_clusters <- cluster_edge_betweenness(g)
  
  community_subgraphs <- list()
  for (i in unique(membership(GN_clusters))) {
    community_nodes <- which(membership(GN_clusters) == i)
    if (length(community_nodes) >= cutoff) {
      community_subgraphs[[i]] <- induced_subgraph(g, community_nodes)
    }
  }
  community_subgraphs <- Filter(Negate(is.null), community_subgraphs)
  return(community_subgraphs)
}

# select the large components to cluster
GN_control_1 <- GN_clustering(control_components[["subgraph_1"]], cutoff = 6)
GN_cancer_1 <- GN_clustering(cancer_components[["subgraph_1"]], cutoff = 6)
GN_cancer_2 <- GN_clustering(cancer_components[["subgraph_2"]], cutoff = 6)

# remove the larger components that were clustered
control_components$subgraph_1 <- NULL
cancer_components$subgraph_1 <- NULL
cancer_components$subgraph_2 <- NULL

# combine the smaller components with the GN clusters of the large components 
all_control_subgraphs <- c(control_components, GN_control_1)
all_cancer_subgraphs <- c(cancer_components, GN_cancer_1, GN_cancer_2)

# Initialize a matrix to store the overlap scores
overlap_matrix <- matrix(0, nrow = length(all_control_subgraphs), ncol = length(all_cancer_subgraphs))

# Iterate over both lists of subgraphs
for (i in seq_along(all_control_subgraphs)) {
  for (j in seq_along(all_cancer_subgraphs)) {
    # Calculate overlap score
    overlap_matrix[i, j] <- length(intersect(names(V(all_control_subgraphs[[i]])), names(V(all_cancer_subgraphs[[j]]))))/min(length(all_control_subgraphs[[i]]), length(all_cancer_subgraphs[[j]]))
  }
}

# rename all subgraphs to assign a number to them (helps with knowing which to visualize later)
row_names <- paste("control", 1:length(all_control_subgraphs), sep = "_")
col_names <- paste("cancer", 1:length(all_cancer_subgraphs), sep = "_")


# define a color palette, make 1 stand out from the rest
color_palette <- c(colorRampPalette(c("white", "red"))(49), "#4f0202")

# Create custom breaks ensuring '1' is highlighted in dark red
breaks <- c(seq(0, 0.999, length.out = 50), 1)

# Create the heatmap with labels
overlap_heatmap <- pheatmap(overlap_matrix,
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            main = "Overlap Heatmap",
                            labels_row = row_names,  
                            labels_col = col_names,
                            legend = TRUE,
                            color = color_palette,  
                            breaks = breaks,
                            fontsize_row = 10,      
                            fontsize_col = 10        
)


names(all_cancer_subgraphs) <- col_names
names(all_control_subgraphs) <- row_names



# visualize the subgraphs that are unique to cancer or control
createNetworkFromIgraph(all_cancer_subgraphs[["cancer_22"]] ,title="cancer_22",collection="cancer")
loadTableData(table, data.key.column="name")
setVisualStyle(style.name)

createNetworkFromIgraph(all_cancer_subgraphs[["cancer_16"]] ,title="cancer_16",collection="cancer")
loadTableData(table, data.key.column="name")
setVisualStyle(style.name)

createNetworkFromIgraph(all_control_subgraphs[["control_3"]] ,title="control_3",collection="control")
loadTableData(table, data.key.column="name")
setVisualStyle(style.name)

createNetworkFromIgraph(all_control_subgraphs[["control_8"]] ,title="control_8",collection="control")
loadTableData(table, data.key.column="name")
setVisualStyle(style.name)

createNetworkFromIgraph(all_control_subgraphs[["control_25"]] ,title="control_25",collection="control")
loadTableData(table, data.key.column="name")
setVisualStyle(style.name)

createNetworkFromIgraph(all_control_subgraphs[["control_27"]] ,title="control_27",collection="control")
loadTableData(table, data.key.column="name")
setVisualStyle(style.name)

# sort the groups to divide them in the heatmap
meta_data_sorted <-  meta_data[order(meta_data$Diagnosis), c("SAMPLE_ID", "Diagnosis")]
met_matrix_sorted <- met_matrix[meta_data_sorted$SAMPLE_ID, ]

meta_data_sorted$Diagnosis <- ifelse(meta_data_sorted$Diagnosis == 0, "non-cancer", "cancer")
#meta_data_sorted$Diagnosis <- as.factor(meta_data_sorted$Diagnosis)

group <- meta_data_sorted$Diagnosis
group_colors <- c("non-cancer" = "orange", "cancer" = "purple")
col_annotation <- HeatmapAnnotation(Group = group, col = list(Group = group_colors))

# creating metabolite subsets of nodes present in the unique modules
control_3_subset <- met_matrix_sorted[, names(V(all_control_subgraphs[["control_3"]]))]
control_8_subset <- met_matrix_sorted[, names(V(all_control_subgraphs[["control_8"]]))]
control_25_subset <- met_matrix_sorted[, names(V(all_control_subgraphs[["control_25"]]))]
cancer_13_subset <- met_matrix_sorted[, names(V(all_cancer_subgraphs[["cancer_13"]]))]
cancer_22_subset <- met_matrix_sorted[, names(V(all_cancer_subgraphs[["cancer_22"]]))]

# pareto-scaling the data for better visualization on heatmap
scaled_data <- apply(cancer_22_subset, 2, function(col){
  mean_col <- mean(col, na.rm = T)
  sd_col <- sd(col, na.rm = T)
  (col-mean_col)/sqrt(sd_col)
})

# Create heatmap object to find clusters between cancer an non-cancer
# data not shown in report due to unclear results
Heatmap(t(scaled_data), 
        name = "metabolite abundance",
        top_annotation = col_annotation,
        column_split = group,
        column_gap = unit(2, "mm"),  # Add space between groups
        show_column_names = FALSE,
        show_row_names = TRUE,
        cluster_rows = FALSE,
        heatmap_legend_param = list(title = "Metabolite abundance"),
        border = TRUE,
        col = colorRampPalette(c("navy" ,"white", "firebrick"))(100))  # Add border around the heatmap
