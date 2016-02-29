library(GEOquery)

## Array quality in mice studies was assessed with the arrayQualityMetrics package of Bioconductor (http://www.bioconductor.org). Data were normalized and logged-2 transformed using the RMA algorithm.

hbv_raw_data <- getGEO(filename = './raw_data/GSE25097_family.soft')
# dir.create('./data')
# save.image('./data/data.rdata')
gsmlist <- GSMList(hbv_raw_data)
# rm(hbv_raw_data)
gc()
probe_ids <- as.vector(Table(gsmlist[[1]])$ID_REF)
sample_names <- sapply(1:length(gsmlist), function(ii){
  names(gsmlist[ii])
})
matrix_data <- ldply(1: length(gsmlist), .progress = 'text', .fun = function(ii){
  tab <- Table(gsmlist[[ii]])
  id_match <- match(tab$ID_REF, probe_ids)
  tab$VALUE[id_match] %>% as.numeric()
})

matrix_data <- t(matrix_data)
colnames(matrix_data) <- sample_names
rownames(matrix_data) <- probe_ids

# expriment design ----
experTable <- ldply(1: length(gsmlist), function(ii){
  c(names(gsmlist[ii]), gsmlist[[ii]]@header$source_name_ch1)
})

phenoData()
# mapping gene ----
matrix_data_df <- matrix_data %>% as.data.frame()
# matrix_data_df$probe_id <- rownames(matrix_data)
# hcc & nohcc
hbv_data <- matrix_data_df[, c(290:557, 7:46)]
colnames(matrix_data_df)
# 差异共表达分析  DCGL
hbv_data <- log2(hbv_data)

library(DCGL)

# 方差对探针进行过滤
var_gene <- varianceBasedfilter(hbv_data, 0.01)
var_gene_name <- rownames(var_gene)

hbv_hcc <- hbv_data[, 1:268]
hbv_nohcc <- hbv_data[, 269: 308]

# 差异表达
# 添加差异表达探针
library(limma)
# browseVignettes('limma')
tinyCond <- c(rep(1,268),rep(2,40)) # 1是肝癌，2是肝硬化
design <- model.matrix(~ 0+factor(c(rep(1, 268),rep(2, 40))))
colnames(design) <- c("hcc", "nohcc") 

fit <- lmFit(hbv_data, design)

contrast.matrix <- makeContrasts(hcc-nohcc, levels=design) 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
dc_raw <- topTable(fit2, adjust="BH", number = 100000)
dc_raw$gene <- rownames(dc_raw)
dc_raw_filter <- filter(dc_raw, adj.P.Val<0.01)
dc_raw_gene <- dc_raw_filter$gene
# 两者合并
raw_gene <- c(var_gene_name, dc_raw_gene) %>% unique
raw_gene_df <- data.frame(gene = raw_gene, stringsAsFactors = F)
hbv_data_1 <- hbv_data[raw_gene, ]

hbv_hcc <- hbv_data_1[, 1:268]
hbv_nohcc <- hbv_data_1[, 269: 308]
save.image('./data/hbv.rda')
DCp.res <- DCp(hbv_hcc, hbv_nohcc,
               r.method = c("pearson", "spearman")[1],
               link.method = c("qth", "rth", "percent")[1],
               cutoff = 0.01,
               N=100,
               N.type = c("pooled", "gene_by_gene")[2],
               q.method = c("BH", "holm", "hochberg", "hommel", "bonferroni","BY", "fdr")[1])

View(DCp.res)
save(DCp.res, file = './data/DCp_res_hbv.rdata')
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
