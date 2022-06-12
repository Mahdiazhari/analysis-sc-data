# Construct cell correlation matrices for RCA cell lines dataset.
# This rev version converts 0 values in the count matrix to 1.

#setwd("~/git/SCT-MoA")
setwd("~/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis")

library(methods)
#library(dismay)
library(propr)
library(stats)
library(WGCNA)
library(pcaPP)
library(Seurat)
library(Matrix)
library(dplyr)
#install.packages('gdata')
library(gdata)
library(data.table)
install.packages('R.utils')


## Measure of Association Function
assoc = function(mat, metric = c(
  'pearson', 'spearman'), ...) {
  
  # first, convert to numeric, if needed 
  if (typeof(mat) == "integer") {
    mat[] <- as.numeric(mat)
  }
  
  #replace 0s with 1s
  
  
  # switch over metrics
  cor = NULL
  if (metric == 'pearson') {
    cor = WGCNA::cor(log(mat), method = 'pearson', ...)
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

options(stringsAsFactors = F)

## Read file

## Read file
file = "data/GSE81861.txt.gz"
dat = read.delim(file, row.names = NULL)
names = dat[[1]]
cell = t(as.matrix(dat[, -1]))
colnames(cell) = paste0(names, '-', seq_len(ncol(cell)))
cell[cell==0] <- 1
sum((cell==0)==TRUE)

## calculate each correlation
#coefs = dismay::metrics()
coefs = metrics()
for (coef in coefs) {
  message("  calculating ", coef, " matrix ...")
  # create correlation matrix 
  output = file.path("results/cor/rev", paste0(
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


library(ggplot2)
ggplot(dat) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  #xlab("Raw expression counts") +
  ylab("Raw expression counts") +
  #ylab("Number of genes")+
  xlab("Number of genes")







