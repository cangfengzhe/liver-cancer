load('./data/DCGL0905.rda')
dc_gene <- DCsum.res$DCGs$DCG %>% as.character()
hcv_data <- hcv_data[, 1:63]
dc_heat_map <- hcv_data[dc_gene,]


#全部基因
exprs.1 <- hcv_data[, 1:16]
exprs.2 <- hcv_data[, 17: 63]

# 差异共表达基因
exprs.1 <- dc_heat_map[, 1:16]
exprs.2 <- dc_heat_map[, 17: 63]
r.method = 'pear'
genes <- rownames(exprs.1)
exprs.1 <- as.matrix(exprs.1)
exprs.2 <- as.matrix(exprs.2)
cor.1 <- cor(t(exprs.1), method = r.method, use = "pairwise.complete.obs")
cor.2 <- cor(t(exprs.2), method = r.method, use = "pairwise.complete.obs")
cor.1 <- cor.1[lower.tri(cor.1, diag = F)]
cor.2 <- cor.2[lower.tri(cor.2, diag = F)]

name.row <- matrix(rep(genes, length(genes)), length(genes), 
                   length(genes))
name.col <- matrix(rep(genes, length(genes)), length(genes), 
                   length(genes), byrow = T)
name.pairs <- matrix(paste(name.row, name.col, sep = ","), 
                     length(genes), length(genes))

name.pairs <- name.pairs[lower.tri(name.pairs, diag = F)]

levels(sample$type) <- c('nohcc', 'hcc')

dc_heatmap <- rbind(data.frame(cor = cor.1, type = 'hcc'),data.frame(cor = cor.2, type = 'nocc'))

cor_diff <- data.frame(cor=cor.1-cor.2, type='diff')

library(ggplot2)
library(reshape2)
len <- nrow(dc_heatmap)/2 %>% floor
sample <- dc_heatmap[(len-300): (len+300),]



save.image('./data/heatmap.rda')

t.test(cor~type, data = dc_heatmap)


dc_heatmap$cor[dc_heatmap$type == 'hcc'] <- dc_heatmap$cor[dc_heatmap$type == 'hcc'] -0.07
dc_heatmap$type %>% levels()
library(ggplot2)
library(ggthemes)
library(scales)


color1 = rgb(255, 193, 193, maxColorValue = 255)
color2 = rgb(192, 194, 255, maxColorValue = 255)
color1 = '#f48c37'; color2 = '#7030a0'
g <- ggplot(data = dc_heatmap, aes(x = cor, fill = type)) +
  geom_density(alpha=.5)+
  scale_fill_manual(values = c(color2, color1), name = 'Type', labels = c('Cirrhosis with HCC','Cirrhosis without HCC'))+
  theme_bw()+
  theme( axis.title = element_text(size= 14),  # 设置轴标题，可以细致到 axis.title.x, axis.title.y
         axis.text = element_text(size = 10),      # 设置轴刻度
         legend.title = element_text(size = 12),   # 设置图例的标题
         legend.text = element_text(size = 12),
         panel.grid.major = element_blank()
         ) +
  scale_y_continuous( expand = c(0, 0),
                     limits= c(0, 1.3))+
  scale_x_continuous(expand = c(0.01,0))
g

ggsave('distribute.png', dpi = 600, width = 16, height = 13, unit = 'cm')

