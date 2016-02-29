# 重新进行计算
load('./data/eset_data.rda')
library(hgu133plus2.db)
symbol <-  as.data.frame(hgu133plus2SYMBOL)

hcv_data <- left_join(eset_data, symbol, by = c('gene' = 'probe_id')) %>% na.omit() %>% dplyr::select(c(1:63,65)) %>% 
  group_by(symbol) %>% 
  summarise_each(funs(mean)) %>% 
  as.data.frame()

rownames(hcv_data) <- hcv_data$symbol
hcv_data <- hcv_data[, 2:64]

# 采用DCGL package 进行计算

exprs_hcc <- hcv_data[, 1:16]
exprs_nohcc <- hcv_data[, 17: 63]

library(DCGL)
DCp.res <- DCp(exprs_hcc, exprs_nohcc,
               r.method = c("pearson", "spearman")[1],
               link.method = c("qth", "rth", "percent")[1],
               cutoff = 0.01,
               N=100,
               N.type = c("pooled", "gene_by_gene")[2],
               q.method = c("BH", "holm", "hochberg", "hommel", "bonferroni","BY", "fdr")[1])

save(DCp.res, file = './data/DCp_res_hcv.rdata')
DCe.res <- DCe(exprs_hcc, exprs_nohcc,
               link.method = c("qth", "rth", "percent")[1],
               cutoff = 0.01,
               r.method = c("pearson", "spearman")[1],
               q.method = c("BH", "holm", "hochberg", "hommel","bonferroni", "BY", "fdr")[1],
               nbins = 20, p = 0.1)


# colnames(DCp.res)[3:4] <-  colnames(DCp.res)[4:3]
# 修改包数据， 无需pvalue校正

DCp.res <- DCp.res[, c(1, 2, 4, 3)]
colnames(DCp.res)[3:4] <-  colnames(DCp.res)[4:3]
DCsum.res <- DCsum(DCp.res, DCe.res,
                   DCpcutoff = 0.01,
                   DCecutoff = 0.01)


# 差异共表达基因与差异表达基因合并
dce <- DCsum.res$DCGs
dce_pvalue <- dce %>% left_join(dc_fc, by = c('DCG' = 'gene'))


DCsum.res$DCGs %>% View
DCsum.res$DCLs %>% nrow

# colnames(DCsum.res$DCGs)


dcg <- DCsum.res$DCGs
dc_link <- DCsum.res$DCLs
# 合并DCp与差异表达
DCp.res$DCG <- rownames(DCp.res)
dcp_de <- left_join(DCp.res, dc_fc, by = c('DCG' = 'gene'))

# 寻找 TF 调控关系
data(tf2target)
DRsort.res <- DRsort(DCsum.res$DCGs, DCsum.res$DCLs, tf2target, D2)

# 导出差异共表达网络
dcn <- DRsort.res$DCLs
dcn$diff <- dcn$cor.2 - dcn$cor.1  
write.csv(dcn, file = './data/diff_co_express_network.csv')  

# 对转录因子进行排序
exp_gene <- rownames(D2)
DRrank.res <- DRrank(DCsum.res$DCGs, DCsum.res$DCLs, tf2target, exp_gene, Nperm = 10)

View(DRrank.res)

# 抽取 差异共表达基因与转录因子的联系
View(DRsort.res$DCLs)
View(DRsort.res$DCGs)

# 在DRsort.res$DRL中选择 2个基因都是差异共表达基因
dcg_link <- filter(DRsort.res$DRLs, grepl(';', DCG) ) %>% as.data.frame()

unique(c(dcg_link$Gene.1, dcg_link$Gene.2)) %>% length
dcg_link_1 <- DRsort.res$DCLs
dcg_link_1$diff <- dcg_link_1$cor.1 - dcg_link_1$cor.2
dcg_link_2 <- filter(dcg_link_1, diff<=-1 & diff)
summary(dcg_link_1$diff)
save.image('./data/DCGL0828.rda')



