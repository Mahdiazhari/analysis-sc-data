
setwd("~/Documents/BSE DSM/4 Master Thesis/Thesis Code/thesis")
pacman::p_load(tidyverse, magrittr, ClusterR, mclust) 
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2")
library(easyGgplot2)

file = "data/GSE81861.txt.gz"
dat <- read.delim(file, row.names = NULL)
dat[1:5, 1:5]
dat <- dat[,-1]
dat[1:5, 1:5]
sum((dat==0)==TRUE)
dat[dat==0] <- 1
sum((dat==0)==TRUE)


#our dat has cells in rows and genes as columns (cells x genes)

counts_per_cell <- Matrix::rowSums(dat)
hist(log10(counts_per_cell+1),main='reads per cell',col='wheat',breaks=100)
hist(counts_per_cell,main='reads per cell',col='wheat',breaks=100)
length(counts_per_cell)

counts_per_gene <- Matrix::colSums(dat)
hist(log10(counts_per_gene+1),col='wheat',main='reads per gene',breaks=100)
hist(counts_per_gene,col='wheat',main='reads per gene',breaks=100)
length(counts_per_gene)

## plots to use: Genes and Cells with Non-zero reads
genes_per_cell <- Matrix::rowSums(dat>0) # count gene only if it has non-zero reads mapped.
#h1 <- hist(log10(genes_per_cell+1),col='wheat',main='genes per cell')
hist(genes_per_cell,col='wheat',ylab="N cells", xlab="N genes", main='genes with non-zero reads', breaks=100)

cells_per_gene <- Matrix::colSums(dat>0) # only count cells where the gene is expressed
#hist(log10(cells_per_gene+1),col='wheat',main='cells per gene')
hist((cells_per_gene),col='wheat', ylab="N genes", xlab="N cells", main='cells per gene', breaks=100)
h = grid.arrange(h1, h2, labels = letters, label_size = 10)


#plot(counts_per_cell, genes_per_cell,log='xy',col='blue')
#plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)', col='red')

#Proportion of Zeroes/Dropouts
sum((dat==0)==TRUE)
zeros <- colSums(dat==0)/nrow(dat)*100
hist(zeros, xlab="percentage", ylab="N genes", main='percentage of zeros per gene',breaks = 100, col='wheat') #use
#hist(log(zeros), main='percentage of zeros per gene',breaks = 100, col='wheat')

zeros_cell <- rowSums(dat==0)/ncol(dat)*100
hist(zeros_cell, xlab="percentage", ylab="N cells", main='percentage of zeros per cell',breaks = 100, col='wheat') #use
#hist(log(zeros), main='percentage of zeros per gene',breaks = 100, col='wheat')

#genes_per_cell_zero <- Matrix::rowSums(dat==0) # count gene only if it has zero reads mapped.
#cells_per_gene_zero <- Matrix::colSums(dat==0) # only count cells where the gene is not expressed
#hist(genes_per_cell_zero,col='wheat',main='zero genes per cell', breaks=100)
#hist(log(genes_per_cell_zero),col='wheat',main='zero genes per cell', breaks=100)
#h2 <- hist(log10(cells_per_gene+1),col='wheat',main='cells per gene')

#not use: mean-variance relationship
logcnt <- log2(dat+1)
meanlogcntpergene <- colMeans(logcnt)
variabilitypergene <- apply(logcnt, 2, function(x) (sd(x))^2)
plot(meanlogcntpergene, variabilitypergene)

#logcnt <- log2(dat+1)
meancntpergene <- colMeans(dat)
variabilitypergene <- apply(dat, 2, function(x) sqrt(sd(x))) #2 here refer to column
plot(meancntpergene, variabilitypergene)
