library(igraph)

SSN_files <- "/datos/ot/anakamura/Nuri/TALL_Ranked_NoNormal"
SSN_list <- list.files(SSN_files)


#SSN_list <- SSN_list[1:30]

metrics <- data.frame()

for (x in SSN_list) {
  Graph <- read.csv(file.path(SSN_files, x)) 
  Graph <- graph_from_data_frame(Graph[, c("Source", "Target")], directed = FALSE)
  
  patient <- sub("_.*", "", x)
  
  # Métricas de la red completa
  diameter <- diameter(Graph)  
  num_nodes <- vcount(Graph)  
  global_clustering <- transitivity(Graph, type = "global")  
  efficiency <- global_efficiency(Graph)  
  
  # Detección de comunidades con Louvain
  set.seed(123)
  louvain <- cluster_louvain(Graph)
  community_sizes <- sizes(louvain)
  
  # Identificar comunidades con al menos 10 nodos
  large_communities <- which(community_sizes >= 10)
  
  # Crear un subgrafo con solo los nodos en comunidades grandes
  nodes_in_large_communities <- which(membership(louvain) %in% large_communities)
  subgraph <- induced_subgraph(Graph, nodes_in_large_communities)
  
  # Recalcular modularidad en la subred filtrada
  if (vcount(subgraph) > 0) {
    filtered_louvain <- cluster_louvain(subgraph)
    modularity_value <- modularity(filtered_louvain)
  } else {
    modularity_value <- NA
  }
  
  # Calcular métricas de comunidades grandes
  num_communities <- length(large_communities)
  mean_CommunitySize <- ifelse(length(large_communities) > 0, mean(community_sizes[large_communities]), NA)
  max_community_size <- ifelse(length(large_communities) > 0, max(community_sizes[large_communities]), NA)
  
  # Guardar métricas
  metricas <- c(patient, num_communities, diameter, num_nodes, modularity_value,
                global_clustering, efficiency, mean_CommunitySize, max_community_size)
  
  metrics <- rbind(metrics, metricas)
}



# Asignar nombres de columnas
cols <- c("patient", "num_communities", "diameter", "num_nodes", "modularity_score", 
          "global_clustering", "efficiency", "mean_CommunitySize", "max_community_size")
colnames(metrics) <- cols
metrics <- metrics %>% mutate(patient = as.numeric(patient)) %>% arrange(patient)


metricsold <- read.csv("~/TALL/Final/8metrics242.csv")
metrics[] <- lapply(metrics, as.numeric)

write.csv(metrics, file.path("~/TALL/metricsfinal.csv"), row.names = FALSE  )
