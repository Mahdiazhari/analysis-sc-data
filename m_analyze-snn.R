# Analyze shared-nearest-neighbor (snn) clustering of cell lines with each measure of association.

setwd("~/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis")

options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(scran)
library(igraph)
library(mclust)
library(bluster)
library(RColorBrewer)
library(cowplot)

#########
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
#########

# create results container
results = data.frame()

# set up parameter values
k_vals = c(2, 5, 10, 20, 50)
#k_vals =c(2)

# define SNN+louvain clustering function
snn_louvain = function(coexpr, k) {
  # build SNN graph
  ## adapted from scran::buildSNNGraph
  
  #knn = BiocNeighbors::findKNN(-coexpr, k = k, get.distance = F)
  #g = .Call(scran:::cxx_build_snn_rank, knn$index)
  #edges = g[[1]]
  #weights = g[[2]]
  #g = make_graph(edges, directed = F)
  #E(g)$weight = weights
  #g = simplify(g, edge.attr.comb = "first")
  g <- buildSNNGraph(-coexpr,k)
  # cluster with louvain
  cl = igraph::cluster_louvain(g)
  
  # pull out cluster membership
  clusters = membership(cl)
  return(clusters)
}

# analyze each matrix in turn 
files = list.files("results/cor", pattern = "*.Rdata", 
                   full.names = T)
for (file in files) {
  coef = gsub("^.*-|\\.Rdata", "", basename(file))
  load(file, verbose = T) ## coexpr
  
  # replace missing values with median (ZI kendall)
  coexpr[is.na(coexpr)] = median(coexpr, na.rm = T)
  # replace infinite values with the minimum (binomial)
  coexpr[is.infinite(coexpr)] = min(coexpr, na.rm = T)
  
  # scale metrics to [-1, 1] range
  #if (!dismay::is_bounded(coef)) {
  #  coexpr = scales::rescale(coexpr, to = c(-1, 1))
  #}
  
  if (is_bounded(coef)) {
    coexpr = scales::rescale(coexpr, to = c(-1, 1))
  }
  
  # analyze a range of k
  for (k in k_vals) {
    print(paste("Value tested:",k))
    # get cluster membership
    clusters = snn_louvain(-coexpr, k)
    # also get ground truth
    types = gsub("_.*$|\\-.*$", "", colnames(coexpr))
    
    # calculate ARI
    ARI = mclust::adjustedRandIndex(clusters, types)
    # calculate NMI
    NMI = ClusterR::external_validation(as.integer(factor(types)), 
                                        as.integer(clusters), method = 'nmi')
    
    # also test within batches
    batches = c("H1", "GM12878")
    patt = paste0(batches, collapse = "|")
    subset = coexpr[grepl(patt, rownames(coexpr)), 
                    grepl(patt, rownames(coexpr))]
    clusters = snn_louvain(-coexpr, k)
    types = gsub("_.*$|\\-.*$", "", colnames(coexpr))
    ARI_batch = mclust::adjustedRandIndex(clusters, types)
    NMI_batch = ClusterR::external_validation(
      as.integer(factor(types)), as.integer(clusters), method = 'nmi')
    
    # add to results
    results %<>% rbind(list(coefficient = coef, k = k, 
                            ARI = ARI, ARI_batch = ARI_batch,
                            NMI = NMI, NMI_batch = NMI_batch))
    gc()
  }
}

# write
write.csv(results, "results/clustering/sNN-louvain2.csv", row.names = F)

# print results
results %>% group_by(coefficient) %>% 
  summarise(ARI = median(ARI)) %>% arrange(-ARI)
results %>% group_by(coefficient) %>% 
  summarise(NMI = median(NMI)) %>% arrange(-NMI)

# plot
results$coefficient %<>% clean_metric()
labels = levels(with(results, reorder(coefficient, -ARI, median)))
p1 = ggplot(results, aes(x = reorder(coefficient, -ARI, median), y = ARI, 
                         fill = coefficient)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_boxplot() +
  #scale_x_discrete(labels = c(
  #  expression(rho[p]), labels[2], expression(phi[s]), labels[-c(1:3)])) +
  scale_x_discrete(labels = c(expression(phi[s]),
                              expression(rho[p]), labels[-c(1:2)])) +
  scale_y_continuous("ARI", limits = c(0, 1), 
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = "none") +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
p1

# plot NMI
labels = levels(with(results, reorder(coefficient, -NMI, median)))
p2 = ggplot(results, aes(x = reorder(coefficient, -NMI, median), 
                         y = NMI, fill = coefficient)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_x_discrete(labels = c(
  #  labels[1], expression(rho[p]), expression(phi[s]), labels[-c(1:3)])) +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]), labels[-c(1:2)])) +
  scale_y_continuous("NMI", limits = c(0, 1), 
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = "none") +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
p2

# save
s10 = plot_grid(p1, p2, labels = letters, label_size = 10)
ggsave("results/fig/round2/Figure_s10_Snn.pdf", s10, width = 17, height = 8, units = 'cm', 
       useDingbats = F)

