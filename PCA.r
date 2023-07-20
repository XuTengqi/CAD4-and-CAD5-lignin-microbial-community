library(vegan)

##################################
#读取环境数据
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

##PCA 排序
#环境变量需要标准化，详情 ?rda
pca_env <- rda(env, scale = TRUE)

#I 型标尺
summary(pca_env, scaling = 1)
#II 型标尺
#summary(pca_env, scaling = 2)

##重要内容提取
#各主成分（PCA轴）特征根，特征根是固定的，和标尺选择无关
pca_eig <- pca_env$CA$eig
#除以特征根总和，即可得各主成分（PCA轴）解释量
pca_exp <- pca_env$CA$eig / sum(pca_env$CA$eig)

#提取对象排序坐标，以 I 型标尺为例，以前两轴为例
#scores() 提取样方坐标
site.scaling1 <- scores(pca_env, choices = 1:2, scaling = 1, display = 'site')
#或者在 summary 中提取
site.scaling1 <- summary(pca_env, scaling = 1)$sites[ ,1:2]
site.scaling1
#write.table(site.scaling1, 'site.scaling1.txt', col.names = NA, sep = '\t', quote = FALSE)

#提取变量排序坐标，以 II 型标尺为例，以前两轴为例
#scores() 提取环境变量坐标
env.scaling2 <- scores(pca_env, choices = 1:2, scaling = 2, display = 'sp')
#或者 summary 中提取
env.scaling2 <- summary(pca_env, scaling = 2)$species[ ,1:2]
env.scaling2
#write.table(env.scaling2, 'env.scaling2.txt', col.names = NA, sep = '\t', quote = FALSE)

##PCA 作图
#简洁版双序图，详情 ?biplot
pc1 <- paste('PC1:', round(pca_exp[1]*100, 2), '%')
pc2 <- paste('PC2:',round(pca_exp[2]*100, 2), '%')

par(mfrow = c(1, 2))
biplot(pca_env, choices = c(1, 2), scaling = 1, main = 'I 型标尺', col = c('red', 'blue'), xlab = pc1, ylab = pc2)
biplot(pca_env, choices = c(1, 2), scaling = 2, main = 'II 型标尺', col = c('red', 'blue'), xlab = pc1, ylab = pc2)

#ordiplot() 作图，I 型标尺为例，详情 ?ordiplot
env.scaling1 <- summary(pca_env, scaling = 1)$species[ ,1:2]

png('pca_env.png', width = 600, height = 600, res = 100, units = 'px')
ordiplot(pca_env, type = 'n', choices = c(1, 2), scaling = 1, main = 'I 型标尺', xlab = pc1, ylab = pc2)
points(pca_env, dis = 'site', choices = c(1, 2), scaling = 1, pch = 21, bg = c(rep('red', 4), rep('orange', 4), rep('green3', 4)), col = NA, cex = 1.2)
arrows(0, 0, env.scaling1[ ,1], env.scaling1[ ,2], length = 0.1, lty = 1, col = 'blue')
text(pca_env, dis = 'sp', choices = c(1, 2), scaling = 1, col = 'blue', cex = 0.8)
dev.off()

#ggplot2 作图，I 型标尺为例
#提取样方和环境变量排序坐标，前两轴
site.scaling1 <- data.frame(summary(pca_env, scaling = 1)$sites[ ,1:2])
env.scaling1 <- data.frame(summary(pca_env, scaling = 1)$species[ ,1:2])

#添加分组信息
site.scaling1$sample <- rownames(site.scaling1)
site.scaling1$group <- c(rep('A', 4), rep('B', 4), rep('C', 4))
env.scaling1$group <- rownames(env.scaling1)

#ggplot2 作图
library(ggplot2)

p <- ggplot(site.scaling1, aes(PC1, PC2)) +
geom_point(aes(color = group)) +
scale_color_manual(values = c('red', 'orange', 'green3')) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
labs(x = pc1, y = pc2) +
geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
geom_segment(data = env.scaling1, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.1, 'cm')), size = 0.3, color = 'blue') +
geom_text(data = env.scaling1, aes(PC1 * 1.1, PC2 * 1.1, label = group), color = 'blue', size = 3)

#ggsave('pca_env.pdf', p, width = 4.5, height = 4)
ggsave('pca_env.png', p, width = 4.5, height = 4)

##################################
#读取物种数据
species <- read.delim('phylum_table.txt', sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

#物种数据 Hellinger 预转化，详情 ?decostand
species_hel <- decostand(species, method = 'hellinger')

#使用转化后的物种数据执行 PCA，无需标准化
pca_sp <- rda(species_hel, scale = FALSE)

#summary(pca_sp, scaling = 1)
#summary(pca_sp, scaling = 2)

#这里的排序图就只简单展示样方，不展示物种了
par(mfrow = c(1, 2))
ordiplot(pca_sp, scaling = 1, display = 'site', type = 'text')
ordiplot(pca_sp, scaling = 1, display = 'site', type = 'points')

#其它细节可参考上文

##################################
##帮助评估有解读价值的 PCA 轴
#Kaiser-Guttman 准则
#计算特征值的平均值，可据此保留高于平均值的主成分轴
pca_eig[pca_eig > mean(pca_eig)]

#断棍模型
#将单位长度的“棍子”随机分成与 PCA 轴数一样多的几段，由高往低排序
#可据此选取特征根大于所对应断棍长度的轴，或者总和大于所对应断棍长度总和的前几轴
n <- length(pca_eig)
bsm <- data.frame(j=seq(1:n), p = 0)
bsm$p[1] <- 1/n
for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
bsm$p <- 100*bsm$p/n
bsm

#绘制每轴的特征根和方差百分比 
par(mfrow = c(2, 1))
barplot(pca_eig, main = '特征根', col = 'bisque', las = 2)
abline(h = mean(pca_eig), col = 'red')
legend('topright', '平均特征根', lwd = 1, col = 2, bty = 'n')
barplot(t(cbind(100 * pca_eig/sum(pca_eig), bsm$p[n:1])), beside = TRUE, main = '% 变差', col = c('bisque', 2), las = 2)
legend('topright', c('% 特征根', '断棍模型'), pch = 15, col = c('bisque', 2), bty = 'n')

##某大佬的函数
source('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/cleanplot.pca.R')
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pca_env <- rda(env, scale = TRUE)

par(mfrow = c(1, 2))
cleanplot.pca(pca_env, scaling = 1)	#I 型标尺中加入平衡贡献圈，评估重要的变量
cleanplot.pca(pca_env, scaling = 2)

##对比环境变量标准化前后的差异
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

pca_env1 <- rda(env, scale = FALSE)	#不执行标准化
pca_env2 <- rda(env, scale = TRUE)	#执行标准化

#I 型标尺
par(mfrow = c(1, 2))
cleanplot.pca(pca_env1, scaling = 1)
cleanplot.pca(pca_env2, scaling = 1)

round(pca_env1$CA$eig/sum(pca_env1$CA$eig)*100, 2)
round(pca_env2$CA$eig/sum(pca_env2$CA$eig)*100, 2)

##对比物种数据 Hellinger 转化前后的差异
otu <- read.delim('otu_table.txt', sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
otu_hel <- decostand(otu, method = 'hellinger')

pca_sp1 <- rda(otu, scale = FALSE)	#使用 Hellinger 转化前的数据
pca_sp2 <- rda(otu_hel, scale = FALSE)	#使用 Hellinger 转化后的数据
pca_sp3 <- rda(otu, scale = TRUE)	#使用 Hellinger 转化前的数据，但执行标准化

#特征值提取
pca_exp1 <- pca_sp1$CA$eig / sum(pca_sp1$CA$eig)
pc1_sp1 <- paste('PC1:', round(pca_exp1[1]*100, 2), '%')
pc2_sp1 <- paste('PC2:',round(pca_exp1[2]*100, 2), '%')

pca_exp2 <- pca_sp2$CA$eig / sum(pca_sp2$CA$eig)
pc1_sp2 <- paste('PC1:', round(pca_exp2[1]*100, 2), '%')
pc2_sp2 <- paste('PC2:',round(pca_exp2[2]*100, 2), '%')

pca_exp3 <- pca_sp3$CA$eig / sum(pca_sp3$CA$eig)
pc1_sp3 <- paste('PC1:', round(pca_exp3[1]*100, 2), '%')
pc2_sp3 <- paste('PC2:',round(pca_exp3[2]*100, 2), '%')

#I 型标尺
par(mfrow = c(2, 2))
#ordiplot(pca_sp1, scaling = 1, display = 'site', type = 'text')
#ordiplot(pca_sp2, scaling = 1, display = 'site', type = 'text')
#ordiplot(pca_sp3, scaling = 1, display = 'site', type = 'text')
ordiplot(pca_sp1, dis = 'site', type = 'n', choices = c(1, 2), scaling = 1, main = 'Hellinger 前，不标准化', xlab = pc1_sp1, ylab = pc2_sp1)
points(pca_sp1, dis = 'site', choices = c(1, 2), scaling = 1, pch = 21, bg = c(rep('red', 12), rep('orange', 12), rep('green3', 12)), col = NA, cex = 1.2)
ordiplot(pca_sp2, dis = 'site', type = 'n', choices = c(1, 2), scaling = 1, main = 'Hellinger 后，不标准化', xlab = pc1_sp2, ylab = pc2_sp2)
points(pca_sp2, dis = 'site', choices = c(1, 2), scaling = 1, pch = 21, bg = c(rep('red', 12), rep('orange', 12), rep('green3', 12)), col = NA, cex = 1.2)
ordiplot(pca_sp3, dis = 'site', type = 'n', choices = c(1, 2), scaling = 1, main = 'Hellinger 前，标准化', xlab = pc1_sp3, ylab = pc2_sp3)
points(pca_sp3, dis = 'site', choices = c(1, 2), scaling = 1, pch = 21, bg = c(rep('red', 12), rep('orange', 12), rep('green3', 12)), col = NA, cex = 1.2)

##排序图中添加聚类树
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pca_env <- rda(env, scale = TRUE)

#I 型标尺
p <- ordiplot(pca_env, dis = 'site', type = 'n', choices = c(1, 2), scaling = 1)
points(pca_env, dis = 'site', choices = c(1, 2), scaling = 1, pch = 21, bg = c(rep('red', 4), rep('orange', 4), rep('green3', 4)), col = NA, cex = 1.2)

#例如 UPGMA 聚类
env_upgma <- hclust(dist(scale(env)), method = 'average')
ordicluster(p, env_upgma, col = 'gray')
