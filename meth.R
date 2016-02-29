library(limma)
library(GEOquery)
hcv_meth_raw_data <- getGEO(filename = './raw_data/GSE18081_family.soft')
# dir.create('./data')
# save.image('./data/data.rdata')
gsmlist <- GSMList(hcv_meth_raw_data)
# rm(hbv_raw_data)
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
  c(names(gsmlist[ii]), gsmlist[[ii]]@header$title)
})

# mapping gene ----
matrix_data_df <- matrix_data %>% as.data.frame()
# matrix_data_df$probe_id <- rownames(matrix_data)
# 20 个 hcc &  16个 nohcc
hcv_meth_nohcc_header <- experTable$V1[ grep('Cirrhosis non', experTable$V2)]
length(hcv_meth_nohcc_header)

hcv_meth_hcc_header <- experTable$V1[ grep('Cirrhosis$|Cirrhosis-', experTable$V2)]
length(hcv_meth_hcc_header)
hcv_meth_data <- matrix_data_df[, c(hcv_meth_hcc_header, hcv_meth_nohcc_header)]
# 标准化
hcv_meth_data <- normalizeMedianValues(hcv_meth_data)
#
boxplot(hcv_meth_data)

# 差异分析 ----
library(limma)
data_type <- factor(c(rep('hcc', 26),rep('nohcc',16)))
design <- model.matrix(~ 0 + data_type)
colnames(design) <- c('hcc', 'nohcc')
fit <- lmFit(hcv_meth_data, design)
contrast.matrix <- makeContrasts(hcc-nohcc, levels=design) 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
dc <- topTable(fit2, adjust="BH", number = 1000000)
dc$gene <- rownames(dc)
View(dc)
boxplot(hcv_meth_data)
meth_dc <- dc
save(meth_dc, file = './data/meth_dc.rda')

View(meth_dc)
meth_0.01 <- meth_dc[meth_dc$P.Value<0.01,]
meth_gene <- meth_0.01$gene
meth_gene <- sapply(meth_gene, function(x){
  
  xx <- strsplit(x, '_')
  xx[[1]][1]
})
meth_gene -> meth_0.01$gene
# remove the repeat
meth_0.01 <- meth_0.01[c(1:7, 9:35),]
View(meth_0.01)
# meth_0.01 用于作图， ggplot2.r
# 进入 fisher.r 进行miRNA富集

View(DCsum.res$DCLs)

#  fold change

meth_fc <-  apply(hcv_meth_data, 1, function(x){
  mean(x[1:26]) - mean(x[27:42])
})
head(meth_fc)

meth_fc <- data.frame(gene = names(meth_fc),fc = meth_fc)
meth_dc_fc <- meth_dc %>% left_join(meth_fc, by = 'gene')
View(meth_dc_fc)
