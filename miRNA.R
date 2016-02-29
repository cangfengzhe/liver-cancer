# 读取多个数据库 miRNA-target数据

# microRNA ------------------
microRNA_raw <- read_delim('./raw_data/miRNA_target/microRNA', delim = '\t')
microRNA <-  microRNA_raw[, c(2:3)]
colnames(microRNA) <- c('mirna', 'gene_id')

rm(microRNA_raw)
gc()

# miRDB---------------
miRDB_raw <- read_delim('./raw_data/miRNA_target/miRDB_v5.0_prediction_result.txt', delim = '\t', col_names = c('mirRNA', 'target_acc', 'score'))
miRDB_hsa <- miRDB_raw[grep('hsa', miRDB_raw$mirRNA), ]

acc2geneid <- read_csv('./raw_data/miRNA_target/geneid2refseq.csv')
miRDB <- left_join(miRDB_hsa, acc2geneid, by = c('target_acc' = 'accession'))
miRDB <- miRDB[, c(1, 4)]
colnames(miRDB) <- c('mirna', 'gene_id')

rm(miRDB_raw, acc2geneid, miRDB_hsa)

# microT ----------------------------------------------
microT <- read.delim('./raw_data/miRNA_target/microT_v3.0.txt', sep = '|', row.names = NULL)

microT <- read_delim('./raw_data/miRNA_target/microT_v3.0.txt', delim = '|' ) %>% select(c(2, 3))

microT_proc <- microT[, c(2,3)]
rm(microT); gc()
colnames(microT_proc) <- c('miRNA', 'ensg')
microT <- microT_proc[grep('hsa', microT_proc$miRNA), ]
# convert ENSG to gene id
library(org.Hs.eg.db)
ensg2geneid <- as.data.frame(org.Hs.egENSEMBL)

microT_geneid <- microT %>% 
  left_join(ensg2geneid, by = c('ensg' = 'ensembl_id')) %>% 
  dplyr::select(-2)

colnames(microT_geneid) <- c('mirna', 'gene_id')
rm(microT, microT_proc, ensg2geneid); gc()
# target_scan------------
target_scan_raw <- read_delim('./raw_data/miRNA_target/target_scan.txt', delim = '\t') 
target_scan <- target_scan_raw[, c(1,2)]
# 处理ensg的版本号
colnames(target_scan) <- c('miRNA', 'ensg')
target_scan <- unique(target_scan)
ensg <- strsplit(target_scan$ensg, '\\.')
ensg_df <- do.call('rbind', ensg)
target_scan$ensg_1 <- ensg_df[,1]
target_scan <- target_scan[, c(1, 3)]

library(splitstackshape)


# 将基因家族分开-----
target_scan_1 <- cSplit(target_scan %>% as.data.frame(), splitCols = 1, direction = 'long', sep = '/')
# 
target_scan_1$miRNA <- as.character(target_scan_1$miRNA)
target_scan_1$miRNA[grep('(miR)|(let)', target_scan_1$miRNA, invert = T)] <- 
  paste('miR-', target_scan_1$miRNA[grep('(miR)|(let)', target_scan_1$miRNA, invert = T)], sep = '')
# ensg mapping
target_scan_geneid <- target_scan_1 %>% as.data.frame() %>% 
  left_join(ensg2geneid, by = c('ensg_1' = 'ensembl_id'))

target_scan_geneid <- target_scan_geneid[, c(1, 3)]
target_scan_geneid$miRNA <- paste('hsa-', target_scan_geneid$miRNA, sep = '')
colnames(target_scan_geneid) <- c('mirna', 'gene_id')

rm(target_scan_1, target_scan_raw); gc()






# mirTar expriment data----------
mirTar <- read_csv('./raw_data/miRNA_target/mirTar_experiment_valid.csv')
colnames(mirTar) <- c('mirna', 'gene_id')

# 进行计数-------
microRNA <- unique(microRNA)
microT_geneid <- unique(microT_geneid)
miRDB <- unique(miRDB)
mirTar <- unique(mirTar)
target_scan_geneid <- unique(target_scan_geneid)

microRNA_bind <- paste(microRNA$mirna, microRNA$gene_id, sep = '=')
microT_bind <- paste(microT_geneid$mirna, microT_geneid$gene_id, sep = '=')
miRDB_bind <- paste(miRDB$mirna, miRDB$gene_id, sep = '=')
# 实验验证数据暂不处理
# mirTar_bind <- paste(mirTar$mirna, mirTar$gene_id, sep = '=')

target_scan_bind <- paste(target_scan_geneid$mirna, target_scan_geneid$gene_id, sep = '=')

bind <- c(microRNA_bind %>% as.character() %>% unique, microT_bind %>% as.character() %>% unique, miRDB_bind %>% as.character(), target_scan_bind %>% as.character() %>% unique)

bind_t <- table(bind)
length(bind_t)
bind_gt2 <- bind_t[bind_t>1]
bind_miR_target <- strsplit(rownames(bind_gt2), '=') %>% 
  do.call('rbind', .) 

bind_miR_target <- bind_miR_target %>% as.data.frame(stringsAsFactors = F)
bind_miR_target$count <- bind_gt2

# mapping gene name---------
gene_symbol <- as.data.frame(org.Hs.egSYMBOL)
colnames(bind_miR_target) <- c('mirna', 'gene_id', 'count')
# 合并mirTar试验数据
mirTar$count <- 5
bind_miR_target <- rbind(bind_miR_target, mirTar)

miR_target <- bind_miR_target %>% 
  na.omit() %>% 
  left_join(gene_symbol, by = c('gene_id' = 'gene_id'))

save.image('./data/miRNA_target.rda')

# 合并差异共表达数据 ------
load('./data/DCGL0828.rda')
# 使用 dcg_link
library(sqldf)
# dcg_mir <- sqldf('select * from dcg_link left join miR_target on dcg_link.`Gene.1` = miR_target.symbol and dcg_link.`Gene.2` = miR_target.symbol')

# 查找共有的miR
common_mir <- sapply(1:nrow(dcg_link), function(ii){
   intersect(miR_target[which(miR_target$symbol == dcg_link$Gene.1[ii]), 1],
             miR_target[which(miR_target$symbol == dcg_link$Gene.2[ii]), 1])
})
# 合并 dcg_link与 common_mir
dcg_link$common_tf_len <- NA
dcg_link$common_mir <- NA
dcg_link$common_mir_len <- NA
dcg_link$common.TF <- as.character(dcg_link$common.TF)

sapply(1:nrow(dcg_link), function(ii){
 
  dcg_link$common_mir_len[ii] <<- length(common_mir[ii][[1]])
  if(dcg_link$common_mir_len[ii]>0)  dcg_link$common_mir[ii] <<- paste(common_mir[ii][[1]], collapse = ';')
  
  dcg_tf <- strsplit(dcg_link$common.TF[ii], ';')
  dcg_link$common_tf_len[ii] <<- length(dcg_tf[[1]])
  length(dcg_tf[[1]])
})
dcg_link$diff <- dcg_link$cor.1- dcg_link$cor.2
write_csv(dcg_link, path = './data/dcg_link.csv')

View(dcg_link)

#合并差异表达基因
dcg_link$Gene.2 %>% class
dcg_link <- left_join(dcg_link, dc_fc, by = c('Gene.1' = 'gene')) %>% 
  left_join(dc_fc, by = c('Gene.2' = 'gene'))
dcg_link <- dcg_link[, c(1:14, 19, 24, 29, 34)]

mir_link <- data.frame(gene1 = dcg_mir$mirna, gene2 = dcg_mir$symbol, count = dcg_mir$count) %>% unique()
nrow(mir_link)
write_csv(dcg_link, path = './data/dcg_link.csv')
write_csv(mir_link, path = './data/mir_link.csv')
library(DCGL)
data("tf2target")
nrow(tf2target)
dcg <- c(dcg_link$Gene.1 %>% as.character(), dcg_link$Gene.2 %>% as.character()) %>% unique %>% as.data.frame(stringsAsFactors = F)

colnames(dcg) <- 'gene'
tf2target$gene <- as.character(tf2target$gene)
dcg_tf <- dcg %>% 
  left_join(tf2target, by = 'gene')

nrow(dcg_tf)
write_csv(dcg_tf, path = './data/dcg_tf.csv')

# 抽取 gene-TF, gene - miRNA
gene_tf_link <- ldply(1:nrow(dcg_link), function(ii){
  dcg_tf <- strsplit(dcg_link$common.TF[ii], ';')
  rbind(data.frame(obj1 = dcg_link$Gene.1[ii], obj2 = dcg_tf[[1]], type='gene2tf', link = 5, stringsAsFactors = F), 
        data.frame(obj1 = dcg_link$Gene.1[ii], obj2 = dcg_tf[[1]], type='gene2tf', link = 5,  stringsAsFactors = F) )
        
})
gene_tf_link <- unique(gene_tf_link)

gene_mir_link <- ldply(1:nrow(dcg_link), function(ii){
  dcg_mir <- strsplit(dcg_link$common_mir[ii], ';')
  rbind(data.frame(obj1 = dcg_link$Gene.1[ii], obj2 = dcg_mir[[1]], type='gene2mir', link = 5, stringsAsFactors = F), 
        data.frame(obj1 = dcg_link$Gene.1[ii], obj2 = dcg_mir[[1]], type='gene2mir', link = 5,  stringsAsFactors = F) )
  
})
gene_mir_link <- unique(gene_mir_link) %>% na.omit()

gene_gene <-  data.frame(obj1 = dcg_link$Gene.1, obj2 = dcg_link$Gene.2, type='gene2gene', link = dcg_link$diff,  stringsAsFactors = F) 

gene_link <- rbind(gene_gene, gene_tf_link, gene_mir_link)
write_csv(gene_link, path = './data/gene_link_cys.csv')


# GEO 中导入 miRNA microarray 数据----
nrow(dcg_link)
nrow(mir_link)
nrow(dcg_tf)

# 采用全部miRNA-----------
# miR_target_full <- strsplit(rownames(bind_t), '=') %>% 
#   do.call('rbind', .)
# miR_target_full <- miR_target_full %>% as.data.frame(stringsAsFactors = F)
# miR_target_full$count <- bind_t
# colnames(miR_target_full) <- c('mirna', 'gene_id', 'count')
# miR_target_full <- rbind(miR_target_full, mirTar)
# 
# miR_target_full <- miR_target_full %>% 
#   na.omit() %>% 
#   left_join(gene_symbol, by = c('gene_id' = 'gene_id'))
# 
# library(sqldf)
# dcg_mir_full<- sqldf('select * from dcg_link left join miR_target_full on dcg_link.`Gene.1` = miR_target_full.symbol and dcg_link.`Gene.2` = miR_target_full.symbol')
# View(dcg_mir_full)
# dcg_mir3 <- dcg_mir[dcg_mir[,'count']>2,]





library(GEOquery)

geo_data_raw <- getGEO(filename = './raw_data/miRNA_target/GSE40744_family.soft')

# dir.create('./data')
# save.image('./data/data.rdata')
gsmlist <- GSMList(geo_data_raw)
rm(GSE63898)
gc()
# 得到探针id
probe_ids <- as.vector(Table(gsmlist[[1]])$ID_REF)
sample_names <- ldply(1:length(gsmlist), function(ii){
  c(names(gsmlist[ii]), gsmlist[[ii]]@header$title)
})

# 抽出试验数据----
matrix_data <- ldply(1: length(gsmlist), .progress = 'text', .fun = function(ii){
  tab <- Table(gsmlist[[ii]])
  id_match <- match(tab$ID_REF, probe_ids)
  tab$VALUE[id_match] %>% as.numeric()
})

matrix_data <- t(matrix_data)
colnames(matrix_data) <- sample_names$V1
rownames(matrix_data) <- probe_ids
matrix_data <- as.data.frame(matrix_data)



# 探针名转换为 miRNA 名字 ----
# 去除_之前和之后的字符
split_str <- strsplit(probe_ids, '_')
hsa_str <- sapply(split_str, function(x){
 out <- x[grep('hsa', x)]
 if(grepl('star', out)){
   out <- out[grep('star', out)] %>% substr(1, nchar(out) - 5)
 }
 
 out
})

matrix_data$name <- tolower(hsa_str)

# 相同miRNA 取平均
matrix_data_mean <- group_by(matrix_data, name) %>% 
  summarise_each(funs(mean)) %>% 
  as.data.frame()

rownames(matrix_data_mean) <- matrix_data_mean$name
matrix_data_mean <- matrix_data_mean[, -1]
# boxplot(data_log2, boxwex=0.6, notch=T, outline=T, las=2)
# summary(data_log2[,1])


# log2 运算 ------
data_log2 <- log2(matrix_data_mean)
# View(data_log2)


# expriment design ----
# sample_names %>% View

# 区分肝硬化与肝癌-----
sample_names$hcc[grep(pattern = 'HCV-associated cirrhosis surrounding HCC', sample_names$V2)] <- 'hcc'
sample_names$hcc[grep(pattern = 'HCV-associated cirrhosis without HCC', sample_names$V2)] <- 'nohcc'

sample_names <- na.omit(sample_names)
# 抽出肝硬化 肝癌相关的样本
data_log2_hcc <- data_log2[, sample_names$V1]
# data_hcc <- matrix_data[, sample_names$V1]
# limma 进行差异表达分析
design <- model.matrix(~ 0 + sample_names$hcc)
colnames(design) <- c('hcc', 'nohcc')
fit <- lmFit(data_log2_hcc, design)
contrast.matrix <- makeContrasts(hcc-nohcc, levels=design) 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
dc <- topTable(fit2, adjust="BH", number = 100000)
dc$gene <- rownames(dc)
View(dc)
save.image('./data/miRNA_chip.rda')
write_csv(dc, path = './data/miRNA_diff_express.csv')

