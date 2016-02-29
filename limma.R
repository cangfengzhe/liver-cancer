



# dc = function(df, design){

# hcv.R 中得到了标准化后的 eset_data,现在对eset_data做差异表达分析

# 基因的确定
library(hgu133plus2.db)
symbol <-  as.data.frame(hgu133plus2SYMBOL)  # 探针名字-基因名字

# 探针归并到基因----
eset_data$gene <- rownames(eset_data)
eset_data_gene <- left_join(eset_data, symbol, by = c('gene' = 'probe_id')) %>% na.omit() %>% 
  dplyr::select(-c(64)) %>% 
  group_by(symbol) %>% 
  summarise_each(funs(mean)) %>% 
  as.data.frame()

rownames(eset_data_gene) <- eset_data_gene$symbol
eset_data_gene <- eset_data_gene[, c(2:64)]

# 差异分析 ----
data_type <- factor(c(rep('hcc',16),rep('nohcc',47)))
design <- model.matrix(~ 0 + data_type)
colnames(design) <- c('hcc', 'nohcc')
fit <- lmFit(eset_data_gene, design)
contrast.matrix <- makeContrasts(hcc-nohcc, levels=design) 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
dc <- topTable(fit2, adjust="BH", number = 100000)
dc$gene <- rownames(dc)
View(dc)

# 计算 fold change

fc <- ldply(1: nrow(eset_data_gene), .progress = 'text', function(ii){
  hcc_mean <- mean(eset_data_gene[ii, 1:16] %>% as.numeric())
  nohcc_mean <-  mean(eset_data_gene[ii, 17:63] %>% as.numeric())
  fc <- hcc_mean-nohcc_mean
 c(hcc_mean, nohcc_mean, fc)
})

fc$symbol <- rownames(eset_data_gene)

dc_fc <- dc %>% left_join(fc, by = c('gene' = 'symbol'))
dc_fc$square <- 2^(dc_fc$logFC)
save(dc_fc, file = './data/diff_express.rda')
dc_filter <- filter(dc_fc, adj.P.Val<0.01)
save(dc_filter, file = './data/de_0.01.rda')
write_csv(dc_filter, path = './data/dc_p0.01.csv')
View(dc_fc)
