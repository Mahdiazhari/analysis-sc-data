# Analyze hierarchical clustering of cell lines with each measure of association.

#rm(list = ls())
setwd("~/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis")


###########
clean_metric = function(vec) {
  as.character(fct_recode(vec,
                          "Spearman correlation" = "spearman",
                          "Pearson correlation" = "pearson",
                          "phi_s" = "phi_s",
                          "rho_p" = "rho_p"
  ))
}

clean_theme = theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        strip.text = element_text(size = 8),
        strip.background = element_rect(color = "grey50", size = 0),
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"), 
        axis.ticks = element_line(colour = "grey50"),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(size = 12, hjust = 0.5))

is_bounded = function(metric) {
  metric %in% c("pearson", "spearman", "rho_p")
}
##############################

options(stringsAsFactors = F)

pacman::p_load(tidyverse, magrittr, ClusterR, mclust, factoextra, R.utils) 
library(tidyverse)
library(magrittr)
library(ClusterR)
library(mclust)
#install.packages("factoextra")
library(factoextra)
library(R.utils)


# create results container
#results = data.frame()

# Analyze each matrix in turn 

#files = list.files("/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/results/cor", pattern = "*.Rdata", 
#                   full.names = T)
#for (file in files) {

  file = "/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/results/cor/saver/GSE81861-spearman.Rdata"
  #file = "/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/results/cor/saver/GSE81861-pearson.Rdata"
  #file = "/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/results/cor/saver/GSE81861-phi_s.Rdata"
  #file = "/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/results/cor/saver/GSE81861-rho_p.Rdata"

  coef = gsub("^.*-|\\.Rdata", "", basename(file))
  load(file, verbose = T) ## coexpr
  
  coexpr[1:5, 1:5]
  sum(is.na(coexpr)==TRUE)
  sum(is.infinite(coexpr)==TRUE)
  
  # scale metrics to [-1, 1] range
  if (is_bounded(coef)) {
    coexpr = scales::rescale(coexpr, to = c(-1, 1))
  }
  
  # hierarchically cluster: Base R: stats
  clust = hclust(as.dist(-coexpr))
  clusters = cutree(clust, k = 7)
  types = gsub("_.*$|\\-.*$", "", names(clusters)) #mars:these are the true label/cell types
  
  #Mars added plotting
  table(clusters)
  fviz_cluster(list(data = coexpr, cluster = clusters), main = "(imputed) rho_p Cluster Plot") #this is okay
  #fviz_dend(clust, cex = 0.5, k = 7, show_labels = FALSE, rect = TRUE, main = "rho_p Dendogram")
  
  # ++++++++++++++++++++++++
  # Hierarchical clustering using factoExtra ()
  
  # Use hcut() which compute hclust and cut the tree
  hc.cut <- hcut(coexpr, k = 7, hc_func = "hclust", hc_method = "complete")
  # Visualize dendrogram
  fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE, main = "(imputed) Spearman  Dendogram")
  # Visualize cluster
  fviz_cluster(hc.cut, ellipse.type = "convex", main = "(imputed) Spearman Cluster Plot") #other ellipse=norm,t
  # Visualize the silhouette
  fviz_silhouette(hc.cut)
  # Visualize clusters as scatter plots
  fviz_cluster(hc.cut)
  # Cluster assignments of observations
  hc.cut$cluster
  # Size of clusters
  hc.cut$size
  
  # ++++++++++++++++++++++++

  # calculate ARI
  ARI = mclust::adjustedRandIndex(clusters, types)
  #ARI = mclust::adjustedRandIndex(hc.cut$cluster, types)
  ARI
  # calculate NMI
  NMI = ClusterR::external_validation(as.integer(factor(types)), 
                                      clusters, method = 'nmi')
  
 