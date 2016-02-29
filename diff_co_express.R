# 差异共表达基因

# source('http://bioconductor.org/biocLite.R')
# biocLite('EBcoexpress')
# library(EBcoexpress)
# vignette('EBcoexpress')
#
browseVignettes('EBcoexpress')


# -----------------
### R code from vignette source 'EBcoexpressVignette.Rnw'

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("EBcoexpress")


###################################################
### code chunk number 2: packageLoad
###################################################
library(EBcoexpress)
tinyCond <- c(rep(1,16),rep(2,47))
tinyPat <- ebPatterns(c("1,1","1,2"))

###################################################
### code chunk number 3: Dcreation
###################################################
D <- makeMyD(D2, tinyCond, useBWMC=F)

set.seed(3)
initHP <- initializeHP(D, tinyCond)

print(initHP)
save(initHP, file = './data/initHP.rda')
oout <- ebCoexpressOneStep(D, tinyCond, tinyPat, initHP)

result1 <- oout$POSTPROBS
priorDiagnostic(D, tinyCond, oout, 1)
priorDiagnostic(D, tinyCond, oout, 2)

ppbDC1 <- 1-result1[,1] 
crit_s <- crit.fun(result1[,1], 0.01) # soft thread
kept_s <- ppbDC1[ppbDC1 >= crit_s]

kept_h <- ppbDC1[ppbDC1 >= 0.99]


klabs_h <- names(kept_h)

gene_name <- strsplit(klabs_h, '~')

gene_name <- do.call('rbind', gene_name)
gene_name <- as.data.frame(gene_name)
gene_name$full_name <- klabs_h

D_cor <- as.data.frame(D)
D_cor$full_name <- rownames(D_cor)
hcv_link <- gene_name %>% left_join(D_cor, by = 'full_name')

hcv_link$diff <- hcv_link$Condition1-hcv_link$Condition2

median(hcv_link$diff %>% abs)

write_csv(hcv_link, path = './data/hcv_link.csv')



save.image('./data/hcv.rdata')




# gene_node ------

gene_node <- read_csv('./data/gene_node.csv')

# Fold Change
D_power <- 2^(D2)
D_power[1,1]
gene_fc <- apply(D_power, 1, function(x){
   mean(x[1:16])/mean(x[17:63])
})
gene_fc <- data.frame(gene = names(gene_fc), fc = gene_fc)
gene_fc$gene <- as.character(gene_fc$gene)
gene_node <- left_join(gene_node, gene_fc, by=c('name' = 'gene'))
cor(gene_node$Degree, gene_node$fc)
hub_gene <- filter(gene_node, Degree>=50 & BetweennessCentrality >= 0.05)

# write_csv(gene_node, path = './data/gene_node.csv')
hcv_link$lt0 <- NA
hcv_link$gt0 <- NA
hcv_link$gt0[which(hcv_link$diff>=0)] <- 1
hcv_link$lt0[which(hcv_link$diff<=0)] <- 1
gt <- hcv_link[, c(1, 7, 8)]
colnames(gt)[1] <- 'gene_name'
lt <- hcv_link[, c(2,7,8)]
colnames(lt)[1] <- 'gene_name'
gene_cor_change <- rbind(gt, lt)
View(gene_cor_change)


