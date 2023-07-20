library(survival)
library(Formula)
library(ggplot2)
library(lattice)
library(Hmisc)


setwd("C:/Users/Administrator/Desktop")
genus <- read.csv('genusA7.csv', row.name = 1, check.names = FALSE)
genus <- genus[which(rowSums(genus) >= 0.005), ]
genus1 <- genus
genus1[genus1>0] <- 1
genus <- genus[which(rowSums(genus1) >= 1), ]
genus_corr <- rcorr(t(genus), type = 'spearman')
r <- genus_corr$r
r[abs(r) < 0.7] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')
p[p>=0.01] <- -1
p[p<0.01 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0  
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'SB genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)



