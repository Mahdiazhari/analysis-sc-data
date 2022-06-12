# Construct cell correlation matrices for RCA cell lines dataset.

#setwd("~/git/SCT-MoA")
setwd("~/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis")

library(methods)
#library(dismay)
library(propr)
library(stats)
library(WGCNA)
library(pcaPP)
library(anndata)

## Measure of Association Function
assoc = function(mat, metric = c(
  'pearson', 'spearman'), ...) {
  
  # first, convert to numeric, if needed 
  if (typeof(mat) == "integer") {
    mat[] <- as.numeric(mat)
  }
  
  # switch over metrics
  cor = NULL
  if (metric == 'pearson') {
    cor = WGCNA::cor(mat, method = 'pearson', ...)
  } else if (metric == 'spearman') {
    cor = stats::cor(mat, method = 'spearman', ...)
  } else if (metric == 'phi_s') {
    cor = -1.0 * propr::phis(mat, select = colnames(mat))@matrix
  } else if (metric == 'rho_p') {
    cor = propr::perb(mat, select = colnames(mat))@matrix
  } else{
    stop("invalid distance/similarity metric: ", metric)
  }
  
  return(cor)
}

## Metric Function
metrics = function() {
  c('pearson', 'spearman', 'phi_s', 'rho_p')
}

is_bounded = function(metric) {
  metric %in% c("pearson", "spearman", "rho_p")
}

## Read file
#brain1= read_h5ad("/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/data/brain_qc_threecells_extremefiltered.h5ad")
#X=brain1$X
#file2 = "data/brain_qc_threecells2_extremefiltered.h5ad"
brain2 = read_h5ad("/Users/mbmedrano/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis/data/brain_qc_threecells2_extremefiltered.h5ad")
X2 = brain$X
names = X[[1]]
#names = dat2[1]
cell = t(as.matrix(dat[, -1]))
colnames(cell) = paste0(names, '-', seq_len(ncol(cell)))

## calculate each correlation
#coefs = dismay::metrics()
coefs = metrics()
for (coef in coefs) {
  message("  calculating ", coef, " matrix ...")
  # create correlation matrix 
  output = file.path("results/cor", paste0(
    "GSE81861-", coef, ".Rdata"))
  if (!file.exists(output)) {
    #coexpr = dismay::dismay(cell, metric = coef)
    coexpr = assoc(cell, metric = coef)
    save(coexpr, file = output)
  }
}


rho_p <- load("results/cor/GSE81861-rho_p.RData")
write.csv(coexpr, "GSE81861-rho_p.csv")

rho_p <- read.csv("results/cor/GSE81861-rho_p.csv")
dim(rho_p)

phi_s<- load("results/cor/GSE81861-phi_s.RData")
write.csv(coexpr, "results/cor/GSE81861-phi_s.csv")
phi_s<- read.csv("results/cor/GSE81861-phi_s.csv")
dim(phi_s)

