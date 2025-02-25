
#metrics matrix as input

metrics


####PLOT WSS

wss <- 0

# Look over 1 to 15 possible clusters
for (i in 1:15) {
  # Fit the model: km.out
  km.out <- kmeans(metrics, centers = i, nstart = 20, iter.max = 50)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
  
}
wss_df <- tibble(clusters = 1:15, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14)) +
  xlab('Number of clusters') +
  geom_hline(
    yintercept = wss, 
    linetype = 'dashed')

print(scree_plot)


#centroids identification

set.seed(123)
kmeans_result_2 <- kmeans(metrics, centers = 2)  # 2 clusters
kmeans_result_3 <- kmeans(metrics, centers = 3)  # 3 clusters
kmeans_result_4 <- kmeans(metrics, centers = 4)  # 4 clusters



kclusters <- data.frame(kmeans2 = kmeans_result_2$cluster, kmeans3 = kmeans_result_3$cluster, kmeans4 = kmeans_result_4$cluster)

write.csv(kclusters, file = ("~/TALL/Final/KmeansClustering/kmeanclustersfinal.csv"))
