library(tidyverse)
library(ggforce)
library(factoextra)
library(FactoMineR)
data.master <- read.csv("3_Outputs/sheets/RiskMasterSheet.csv")
clust_data = data.master[c("Livelihood","Demography","Cultural","Health","Learning","Assets","Flexibility","Agency","Organisation","ic2020","imp.ssp585.2050","Villages","ISO3")]

clust_data = clust_data[sample(1:nrow(clust_data)), ]

# Decide how many clusters to look at
n_clusters <- 20
# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)
# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(scale(clust_data[1:11]), centers = i, nstart = 25)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot #From the graph, you can see the optimal k is seven, where the curve is starting to have a diminishing return.
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)
(scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
    geom_point(size = 4)+
    geom_line() +
    scale_x_continuous(breaks = seq(1, n_clusters, 1)) +
    xlab('Number of clusters'))
#ggsave("optimalN_cluster.png", dpi = 1200, width = 5, height = 5)

# Execution of k-means with k=8
set.seed(123)
km.res <- kmeans(scale(clust_data[1:11]), 4, nstart = 25)
clust_data <- cbind(clust_data, cluster = km.res$cluster)
head(clust_data)

KMNC <- as.factor(clust_data$cluster)
res.pca.cluster <- prcomp(clust_data[1:11], scale. = TRUE, center = TRUE)
summary(res.pca.cluster)

(plt.bip.cluster <- fviz_pca_biplot(res.pca.cluster, 
                                    axes = c(1, 2),
                                    repel = TRUE,
                                    col.var = "black", # Variables color
                                    col.ind = KMNC, 
                                    label = "var", 
                                    labelsize = 3, 
                                    arrowsize = .5)+
    theme_bw(base_size = 10)+
    geom_mark_hull(concavity = 5, expand = 0, radius = 0, aes(fill=KMNC), alpha = 0.2, linewidth=.2)+
    scale_color_manual(name = "Clusters", 
                       values = c("1" = "darkred", 
                                  "2" = "#749B58FF",
                                  "3" = "#466983FF",
                                  "4" = "pink4"))+
    scale_fill_manual(name = "Clusters", 
                      values = c("1" = "darkred", 
                                 "2" = "#749B58FF",
                                 "3" = "#466983FF",
                                 "4" = "pink4"))+
    guides(colour = "none", shape = "none")+
    labs(title = "", 
         x = paste("PC-1", "(",round(get_eig(res.pca.cluster)[,2][1],1),"%",")"), 
         y = paste("PC-2", "(",round(get_eig(res.pca.cluster)[,2][2],1),"%",")"))+
    theme(legend.position = "right",
          #panel.border = element_blank(), 
          panel.grid = element_blank(), 
          legend.title = element_blank()))

ggsave(plot = plt.bip.cluster, "Clusters.png", dpi = 1200, height = 4, width = 4)
