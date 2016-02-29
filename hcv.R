# 丙肝 肝硬化 肝癌----
# 处理原始microarray 数据
library(affy)
# dir.create('hcv_data')
setwd('./hcv_data/')
raw_data <- ReadAffy(compress = T)
# exprs(raw_data) %>% View
# rma process ----
# eset <- rma(raw_data)
eset <- expresso(raw_data, normalize = F,
                 bgcorrect.method="rma",pmcorrect.method="pmonly",
                 summary.method="liwong")
# get the data ----
eset_data <- exprs(eset)
# normalize by median
library(limma)
eset_median <- normalizeMedianValues(eset_data)
eset_log <- log2(eset_median)

# apply(eset_log, 2, summary)
# boxplot(eset_log, boxwex=0.6, notch=T, outline=T, las=2)
eset_data <- eset_log

gsm_type <- read_csv('./raw_data/gsm_type.csv')
gsm_type <- as.data.frame(gsm_type)
gsm_type$type <- NA
gsm_type$type[grep('without', gsm_type$desc)] <- 'nohcc'
gsm_type$type[grep('without', gsm_type$desc, invert = T)] <- 'hcc'
col_name <- sapply(colnames(eset_data), function(x){
  substring(x, 1, 9)
})
colnames(eset_data) <- col_name
# rownames(eset_data) <- eset_data[,1]
eset_data <- as.data.frame(eset_data)
# eset_data <- eset_data[, -1]
# the first 16 are the hcc, the last 47 are no hcc ----
col_index <- c(dplyr::filter(gsm_type, type == 'hcc') %>% .[,1] , dplyr::filter(gsm_type, type == 'nohcc') %>% .[,1])

eset_data <- eset_data[, col_index]


# save.image(file = '../data/hcv.rdata')


#  DCGL differential co-expression ----
load('./data/hcv.rdata')
library(DCGL)
library(hgu133plus2.db)
symbol <-  as.data.frame(hgu133plus2SYMBOL)

# 方差对探针进行过滤
var_gene <- varianceBasedfilter(hcv_data, 0.05)
var_gene_name <- rownames(var_gene)
# 添加差异表达探针
library(limma)
# browseVignettes('limma')
tinyCond <- c(rep(1,16),rep(2,47)) # 1是肝癌，2是肝硬化
design <- model.matrix(~ 0+factor(c(rep(1,16),rep(2,47))))
colnames(design) <- c("hcc", "nohcc") 
eset_data_2 <- eset_data[, c(1:63)]
fit <- lmFit(eset_data_2, design)

contrast.matrix <- makeContrasts(hcc-nohcc, levels=design) 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
dc_raw <- topTable(fit2, adjust="BH", number = 100000)
dc_raw$gene <- rownames(dc_raw)
dc_raw_filter <- filter(dc_fc, adj.P.Val<0.01)
dc_raw_gene <- dc_raw_filter$gene
# 两者合并
raw_gene <- c(var_gene_name, dc_raw_gene, meth_gene) %>% unique
raw_gene_df <- data.frame(gene = raw_gene, stringsAsFactors = F)
# 整合芯片数据
eset_data$gene <- rownames(eset_data)
raw_sample <- raw_gene_df %>% left_join(hcv_data, by = 'gene')
save(raw_sample, eset_data, file = './data/eset_data.rda')

D1 <- raw_sample

D2 <- left_join(D1, symbol, by = c('gene' = 'probe_id')) %>% na.omit() %>% dplyr::select(2:65) %>% 
  group_by(symbol) %>% 
  summarise_each(funs(mean)) %>% 
  as.data.frame()
rownames(D2) <- D2$symbol
D2 <- D2[, 2:64]
 
# 转到 diff_co_express.r 进行共表达分析


# 下面的方法采用了DCGL包进行差异共表达分析
exprs_hcc <- D2[, 1:16]
exprs_nohcc <- D2[, 17: 63]


  DCp.res <- DCp(exprs_hcc, exprs_nohcc,
                 r.method = c("pearson", "spearman")[1],
                 link.method = c("qth", "rth", "percent")[1],
                 cutoff = 0.01,
                 N=100,
                 N.type = c("pooled", "gene_by_gene")[2],
                 q.method = c("BH", "holm", "hochberg", "hommel", "bonferroni","BY", "fdr")[1])
  
View(DCp.res)
save(DCp.res, file = './data/DCp_res_hcv.rdata')
DCe.res <- DCe(exprs_hcc, exprs_nohcc,
               link.method = c("qth", "rth", "percent")[1],
               cutoff = 0.01,
               r.method = c("pearson", "spearman")[1],
               q.method = c("BH", "holm", "hochberg", "hommel","bonferroni", "BY", "fdr")[1],
               nbins = 20, p = 0.1)


# 修改包数据， 无需pvalue校正
# DCp.res <- DCp.res[, c(1, 2, 4, 3)]
# colnames(DCp.res)[3:4] <-  colnames(DCp.res)[4:3]

DCsum.res <- DCsum(DCp.res, DCe.res,
                   DCpcutoff = 0.01,
                   DCecutoff = 0.01)


# 差异共表达基因与差异表达基因合并
dce <- DCsum.res$DCGs
dce_pvalue <- dce %>% left_join(dc_fc, by = c('DCG' = 'gene'))
write.csv(dce_pvalue, file = './data/dce_pvalue_0906.csv')

DCsum.res$DCGs %>% View
DCsum.res$DCLs %>% nrow

# colnames(DCsum.res$DCGs)


dcg <- DCsum.res$DCGs
dc_link <- DCsum.res$DCLs
# 合并DCp与差异表达

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

# 在DRsort.res$DCL中选择 2个基因都是差异共表达基因
# 加入甲基化的差异基因
dcg_link <- filter(DRsort.res$DCLs, grepl(';', DCG) ) %>% as.data.frame()

unique(c(dcg_link$Gene.1, dcg_link$Gene.2)) %>% length

dcg_link_1 <- DRsort.res$DCLs
dcg_link_1$diff <- dcg_link_1$cor.1 - dcg_link_1$cor.2
dcg_link_2 <- filter(dcg_link_1, diff<=-1 & diff)
summary(dcg_link_1$diff)
save.image('./data/DCGL0828.rda')
save.image('./data/DCGL0905.rda')
