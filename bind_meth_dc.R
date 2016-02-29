# 甲基化差异基因相关性
meth_expr <- left_join(data.frame(gene = meth_gene, stringsAsFactors = F),
                       hcv_data, by = 'gene') %>% 
  na.omit()
rownames(meth_expr) <- meth_expr$gene
meth_expr <- meth_expr[, 2:64]
exprs.1 <- meth_expr[, 1:16]
exprs.2 <- meth_expr[, 17: 63]
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
rm(list = c("name.row", "name.col"))
name.pairs <- name.pairs[lower.tri(name.pairs, diag = F)]
names(cor.1) <- names(cor.2) <- name.pairs
gene_name <- ldply(1:length(cor.1), function(ii){
  tmp <- strsplit(names(cor.1)[ii], ',')
  do.call('rbind', tmp)
})

meth_dc <- data.frame(gene1 = gene_name[,1], gene2 = gene_name[,2],hcc = cor.1, nohcc = cor.2)
View(meth_dc)
meth_dc$diff <- meth_dc$hcc - meth_dc$nohcc
meth_dc <- meth_dc[meth_dc$diff %>% abs() > 0.3, ]
View(meth_dc)

library(sqldf)
meth_gene_df <- data.frame(gene = meth_gene, stringsAsFactors = F)
meth_gene_df <- meth_gene_df %>% filter(!grepl('(MME)|(MYB)', gene))

dcls <- DRsort.res$DCLs
meth_dco <- sqldf('select * from meth_gene_df left join dcls  on meth_gene_df.gene = dcls.`Gene.1` or meth_gene_df.gene = dcls.`Gene.2`')


meth_dco1 <- sqldf('select * from meth_gene_df left join dcls on meth_gene_df.gene = dcls.`DCG`')
# 去掉共有 MME MYB
meth_dco <- meth_dco[setdiff(1:254, c(55, 67, 122, 155, 227)), ]

meth_dco <- meth_dco[, c('Gene.1', 'Gene.2', 'cor.1', 'cor.2')]
meth_dco$diff <- meth_dco$cor.1 - meth_dco$cor.2
colnames(meth_dco) <- colnames(meth_dc)
meth_dc_full <-  rbind(meth_dco, meth_dc) %>% unique
meth_dc_full %>% View

dcg_link1 <- dcg_link[, c('Gene.1', 'Gene.2', 'cor.1', 'cor.2')]
dcg_link1$diff <- dcg_link1$cor.1 - dcg_link1$cor.2
colnames(dcg_link1) <- colnames(meth_dc)
dcg_link_full <- rbind(dcg_link1, meth_dc_full) %>% na.omit()
View(dcg_link_full)
write.csv(dcg_link_full, './data/dcg_link_full.csv', quote = F)

dcg_link_full$gene1 <- as.character(dcg_link_full$gene1)
dcg_link_full$gene2 <- as.character(dcg_link_full$gene2)
dc_meth_gene <- unique(c(dcg_link_full$gene1, dcg_link_full$gene2))



