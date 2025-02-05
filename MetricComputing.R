# Metric calculation, file name input (e.g. 30_SSN.csv)

library(dplyr)
library(igraph)

##Single Sample Network Estimates 
SSN_files <- "/datos/ot/anakamura/Nuri/TALL_Ranked_NoNormal"
SSN_list <- list.files(SSN_files)

## New data frame 
metrics <- data.frame()

for (x in SSN_list) {
  Graph <- read.csv(file.path(SSN_files, x)) 
  Graph <- graph_from_data_frame(Graph[, c("Source", "Target")], directed = FALSE)
  
  # Example of file name input: "30_SSN.csv"
  
  patient <- sub("_.*", "", x)
  
  
  # Metrics not related to communities
  diameter <- diameter(Graph)  
  #edge_density <- edge_density(Graph)  # Redundant in a network with a set number of edges
  num_nodes <- vcount(Graph) 
  #mean_degree <- mean(degree(Graph))  # Redundant in a network with a set number of edges
  global_clustering <- transitivity(Graph, type = "global") 
  #avg_path_length <- mean_distance(Graph)  # Redundant in a network with a set number of edges
  efficiency <- global_efficiency(Graph)  
  
  # Identification and filtering of communities
  set.seed(123)
  louvain <- cluster_louvain(Graph)
  community_sizes <- sizes(louvain)
  
  # Communities bigger than 10 nodes
  filtered_community_sizes <- community_sizes[community_sizes >= 10]
  
  
  # Metrics related to communities 
  num_communities <- length(filtered_community_sizes) 
  mean_CommunitySize <- mean(filtered_community_sizes)  
  max_community_size <- max(filtered_community_sizes)  
  
  # Metrics vector
  metricsV <- c(patient, num_communities, diameter, num_nodes, modularity(louvain),
                global_clustering, efficiency, mean_CommunitySize, max_community_size)
  
  # rowbind
  metrics <- rbind(metrics, metricsV)
}


# Colname vector
cols <- c("patient", "num_communities", "diameter", "num_nodes", "modularity_score", 
          "global_clustering", "efficiency", "mean_CommunitySize", "max_community_size")
colnames(metrics) <- cols

#Arranged patients 
metrics <- metrics %>% mutate(patient = as.numeric(patient)) %>% arrange(patient)



write.csv(metrics, file.path("~/TALL/metrics123.csv") )
