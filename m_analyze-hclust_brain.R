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

pacman::p_load(tidyverse, magrittr, ClusterR, mclust) 
#library(tidyverse)
#library(magrittr)


# create results container
results = data.frame()

# Analyze each matrix in turn 
files = list.files("/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/results/cor", pattern = "*.Rdata", 
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
  
  # hierarchically cluster
  clust = hclust(as.dist(-coexpr))
  clusters = cutree(clust, k = 7)
  types = gsub("_.*$|\\-.*$", "", names(clusters))
  
  # calculate ARI
  ARI = mclust::adjustedRandIndex(clusters, types)
  # calculate NMI
  NMI = ClusterR::external_validation(as.integer(factor(types)), 
                                      clusters, method = 'nmi')
  
  # also test within batches
  batches = c("H1", "GM12878")
  patt = paste0(batches, collapse = "|")
  subset = coexpr[grepl(patt, rownames(coexpr)), grepl(patt, rownames(coexpr))]
  clust = hclust(as.dist(-subset))
  clusters = cutree(clust, k = 2)
  types = gsub("_.*$|\\-.*$", "", names(clusters)[grepl(patt, names(clusters))])
  # calculate ARI
  ARI_batch = mclust::adjustedRandIndex(clusters, types)
  # calculate NMI
  NMI_batch = ClusterR::external_validation(as.integer(factor(types)), 
                                            clusters, method = 'nmi')
  
  # add to results
  results %<>% rbind(list(coefficient = coef, ARI = ARI, ARI_batch = ARI_batch,
                          NMI = NMI, NMI_batch = NMI_batch))
}

# Write
write.csv(results, "results/clustering/hclust.csv2", row.names = F)


# print results
results %>% arrange(-ARI)
results %>% arrange(-ARI_batch)
results %>% arrange(-NMI)
results %>% arrange(-NMI_batch)

labels

# 1. plot main figure for ARI
results$coefficient %<>% clean_metric()
labels = levels(with(results, reorder(coefficient, -ARI)))

pARI = ggplot(results, aes(x = reorder(coefficient, -ARI), y = ARI, 
                        fill = coefficient)) +
  geom_col() +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                              labels[-c(1:2)])) +
  #scale_x_discrete(labels = labels) +
  scale_y_continuous("adjusted Rand index", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = "none") +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
pARI
ggsave("results/fig/round2/Figure_4_ARI.pdf", pARI, width = 8.9, height = 8, units = 'cm', 
       useDingbats = F)

# 2. plot within-batch analyses for ARI
labels = levels(with(results, reorder(coefficient, -ARI_batch)))
pARIbatch = ggplot(results, aes(x = reorder(coefficient, -ARI_batch), y = ARI_batch,
                        fill = coefficient)) +
  geom_col() +
  #scale_x_discrete(labels = c(labels[1], expression(rho[p]), expression(phi[s]),
  #                            labels[-c(1:3)])) +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                              labels[-c(1:2)])) +
  scale_y_continuous("\nARI, H1/GM12878 cells only", limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = "none") +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
pARIbatch
ggsave("results/fig/round2/Figure_S8_ARIbatch.pdf", pARIbatch, width = 8.9, height = 8, units = 'cm', 
       useDingbats = F)

# 3. plot comparison between ARI and ARI within-batch analyses
ARIcompare = plot_grid(pARI, pARIbatch, labels = letters, label_size = 10)
ggsave("results/fig/round2/Figure_ARIbatchcompare.pdf", ARIcompare, width = 17, height = 8, units = 'cm', 
       useDingbats = F)

# 4. plot NMI figure 
labels = levels(with(results, reorder(coefficient, -NMI)))
pNMI = ggplot(results, aes(x = reorder(coefficient, -NMI), y = NMI, 
                         fill = coefficient)) +
  geom_col() +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                              labels[-c(1:2)])) +
  scale_y_continuous("NMI", limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = "none") +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
pNMI
ggsave("results/fig/round2/Figure_NMI.pdf", pNMI, width = 8.9, height = 8, units = 'cm', 
       useDingbats = F)

# 5. plot within-batch under NMI
labels = levels(with(results, reorder(coefficient, -NMI_batch)))
pNMIbatch = ggplot(results, aes(x = reorder(coefficient, -NMI_batch), y = NMI_batch,
                         fill = coefficient)) +
  geom_col() +
  #scale_x_discrete(labels = c(labels[1], expression(rho[p]), expression(phi[s]),
  #                            labels[-c(1:3)])) +
  scale_x_discrete(labels = c(expression(rho[p]), expression(phi[s]),
                              labels[-c(1:2)])) +
  scale_y_continuous("\nNMI, H1/GM12878 cells only", limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = colorRampPalette(rev(
    brewer.pal(n = 12, name = "Paired")))(17), guide = "none") +
  clean_theme + 
  theme(axis.title.x = element_blank(), 
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())
pNMIbatch
ggsave("results/fig/round2/Figure_NMIbatch.pdf", pNMIbatch, width = 8.9, height = 8, units = 'cm', 
       useDingbats = F)

# 6. plot comparison between NMI and NMI within-batch analyses
s9 = plot_grid(pNMI, pNMIbatch, labels = letters, label_size = 10)
ggsave("results/fig/round2/Figure_S9_NMIbatchcompare.pdf", s9, width = 17, height = 8, units = 'cm', 
       useDingbats = F)

# 7. plot comparison between ARI and NMI
ARINMI = plot_grid(pARI, pNMI, labels = letters, label_size = 10)
ggsave("results/fig/round2/Figure_ARINMIcompare.pdf", ARINMI, width = 17, height = 8, units = 'cm', 
       useDingbats = F)




