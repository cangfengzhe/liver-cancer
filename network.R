load("~/baiduyun/work/RFile/liver_cancer/data/dc.rdata")
gene_node <- read_csv('./data/hcv_node.csv')
gene_node_dc <- gene_node %>% left_join(dc, by = c('name' = 'gene')) %>% 
  as.data.frame()
View(gene_node_dc)
colnames(gene_node_dc)
cor_data <- gene_node_dc[, c('Degree', 'adj.P.Val')] %>% na.omit()
cor(cor_data)
unique(cor_data$Degree)

cor_data_2 <- ldply(unique(cor_data$Degree), function(x){
  data_mean <- cor_data[which(cor_data$Degree == x), 2] %>% mean()
  c(x, data_mean)
})
View(cor_data_2)
cor(cor_data_2)
plot(cor_data_2)
