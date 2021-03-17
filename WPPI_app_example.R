### PPI creation

source("WPPI_app.R")

# example gene set
genes.interest <- c("ERCC8","AKT3","NOL3","TTK","GFI1B","CDC25A","TPX2","SHE")

# example HPO annotations set 
HPO.interest <- filter(HPO.data,
                       grepl("Diabetes",HPO_Name)) %>%
  dplyr::select(HPO_Name) %>% unique() %>% pull(HPO_Name)

# Application of the main function
genes_ranked <- score_candidate_genes_from_PPI(genes_interest = genes.interest,
                                               HPO_interest = HPO.interest,
                                               percentage_output_genes = 10,
                                               graph_order = 1)

############################################################
# All available tools and steps needed for the methodology:

# igraph PPI object
graph_op <- graph_from_op(Omnipath.human.data)
edges_op <- E(graph_op)
vertices_op <- V(graph_op)

# genes mapped and not mapped in Omnipath 
genes_mapped <- isgene_omnipath(graph_op,genes.interest,1)
genes_notmapped <- isgene_omnipath(graph_op,genes.interest,0)

# Subgraphs
graph_op_0 <- subgraph_op(graph_op,genes.interest,0)
graph_op_1 <- subgraph_op(graph_op,genes.interest,1)

# Subnetworks 
net_op_1 <- as_data_frame(graph_op_1)

# Adjacency Matrix
adj_op_1 <- graph_to_adjacency(graph_op_1)
saveRDS(adj_op_1, file = "WPPI_Data/Adjacency_Matrix.RDS")

# Common neighbours
shared_neighbors <- common_neighbors(graph_op_1)

# Filter annotation datasets with graph object
GO_data_op <- filter_annot_with_network(GO.data,graph_op_1)
HPO.data.interest <- HPO.data %>% filter(HPO_Name %in% HPO.interest)
HPO_data_op <- filter_annot_with_network(HPO.data.interest,graph_op_1)

# Adjacency Matrix to Weighted Adjacency Matrix
weighted_adj_op_1 <- weighted_adj(graph_op_1,
                                  shared_neighbors,
                                  GO_data_op,
                                  HPO_data_op,
                                  nr_genes(GO_data_op),
                                  nr_genes(HPO_data_op))
saveRDS(weighted_adj_op_1, file = "WPPI_Data/Weighted_Matrix.RDS")

# Random Walk on weighted matrix
random_walk_op_1 <- random_walk(weighted_adj_op_1)

# Prioritization of genes based on functional scores 
genes_ranked <- prioritization_genes(graph_op_1,random_walk_op_1,genes.interest)