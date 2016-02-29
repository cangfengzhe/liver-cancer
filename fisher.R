miR_target
library(DCGL)
data("tf2target")
load('./data/DCGL0828.rda')
load('./data/miRNA_target.rda')
load('./data/miRNA_chip.rda')
dc_link <- filter(DRsort.res$DCLs, grepl(';', DCG) & cor.diff>=0.5 ) %>% as.data.frame(stringsAsFactor = F)
dc_link$Gene.1 <-  as.character(dc_link$Gene.1)
dc_link$Gene.2 <-  as.character(dc_link$Gene.2)

#加入甲基化的基因
dcg <- dc_meth_gene
all_gene_len <- length(unique(c(dcg, tf2target$gene, miR_target$symbol)))
all_gene_len - 1088
# 对 TF， miRNA进行富集分析
mir <- miR_target$mirna %>% tolower() %>% unique
miR_target$mirna <- tolower(miR_target$mirna)
mir_enrich <-   ldply(1: length(mir), .progress = 'text', function(ii){
    mir2gene <- miR_target$symbol[miR_target$mirna == mir[ii]]
  
    mir_len <- length(mir2gene)
    mir_gene_len <- length(intersect(mir2gene, dcg))
    
    fisher_data <- fisher.test(matrix(c(mir_gene_len, mir_len - mir_gene_len, length(dcg) -  mir_gene_len, all_gene_len - mir_len + mir_gene_len - length(dcg)), 2, 2))
    c(mir[ii], fisher_data$p.value)
  })

mir_enrich$p_adjust <- p.adjust(mir_enrich$V2, method = 'BH')
# View(mir_enrich)
mir_enrich$V1 <- tolower(mir_enrich$V1)
mir_enrich_dc <- left_join(mir_enrich, dc, by = c('V1' = 'gene'))
View(mir_enrich_dc)
write.csv(mir_enrich_dc, file = './data/mir_enrich_dc.csv', quote = F)

# 对 TF进行富集分析 
tf <- tf2target$TF %>% unique %>% as.character()
tf2target$TF <- as.character(tf2target$TF)
tf2target$gene <- as.character(tf2target$gene)

tf_enrich <-   ldply(1: length(tf),.progress = 'text', function(ii){
  tf2gene <- tf2target$gene[tf2target$TF == tf[ii]]
  tf_len <- length(tf2gene)
  tf_gene_len <- length(intersect(tf2gene, dcg))
  
  fisher_data <- fisher.test(matrix(c(tf_gene_len, tf_len-tf_gene_len, length(dcg) -  tf_gene_len, all_gene_len - tf_len - length(dcg) + tf_gene_len), 2, 2))
  c(tf[ii], fisher_data$p.value)
})

tf_enrich$p_adjust <- p.adjust(tf_enrich$V2, method = 'BH')

tf_enrich_dc <- tf_enrich %>% left_join(dc_fc, by = c('V1' = 'gene'))
View(tf_enrich_dc)
save(mir_enrich_dc, tf_enrich_dc, file = './data/enrich.rda')


# 提取富集到的miRNA TF ------
miR_target$mirna <- miR_target$mirna %>% tolower()
mir_de <- filter(mir_enrich, p_adjust<0.01)
mir_gene <- ldply(mir_de$V1, function(x){
  mir2gene <- intersect(miR_target$symbol[miR_target$mirna == x],dcg)
  print(length(mir2gene))
  if(length(mir2gene)>0){
    out <- data.frame(mir = x, gene = mir2gene, stringsAsFactors = F)
  } else{
    out <- data.frame(mir = NA, gene = NA)
  }
  out  
})
write.csv(mir_gene, file = './data/mir2gene.csv', quote = F)

tf_de <- filter(tf_enrich_dc, p_adjust<0.01) %>% dplyr::select(c(1))

tf_gene <- ldply(tf_de$V1, function(x){
  tf2gene <- intersect(tf2target$gene[tf2target$TF == x], dcg)
  print(length(tf2gene))
  if(length(tf2gene)>0){
    out <- data.frame(tf = x, gene = tf2gene, stringsAsFactors = F)
  } else{
    out <- data.frame(tf = NA, gene = NA)
  }
  out  
})
write.csv(tf_gene, './data/tf_enrich_gene.csv', quote = F)

# DCGL 获得的转录因子
dcgl_tf <- DRrank.res
dcgl_tf$tf <- rownames(dcgl_tf)
dcgl_tf <- filter(dcgl_tf, p.adj<0.01)
dcgl_tf_gene <- ldply(dcgl_tf$tf, function(x){
  tf2gene <- intersect(tf2target$gene[tf2target$TF == x], dcg)
  print(length(tf2gene))
  if(length(tf2gene)>0){
    out <- data.frame(tf = x, gene = tf2gene, stringsAsFactors = F)
  } else{
    out <- data.frame(tf = NA, gene = NA)
  }
  out  
})

write.csv(dcgl_tf_gene, './data/dcgl_tf_enrich_gene.csv', quote = F)
write.csv(meth_gene_df, './data/meth_gene_dc.csv')
