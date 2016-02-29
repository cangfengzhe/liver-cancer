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

D_filter <- as.data.frame(D)
D_filter$diff <- D[,1] - D[,2]
D_filter$names <- rownames(D_filter)
D_median <- median(abs(D_filter$diff))
D_diff <-  dplyr::filter(D_filter, abs(diff)>D_median)
rownames(D_diff) <- D_diff$names
D_diff <- as.matrix(D_diff[,c(1,2)])

# aa <- cor(express_data_order[, 1:228] %>% t())

# cor(express_data_order[1, 1:228] %>% as.numeric(), express_data_order[2, 1:228] %>% as.numeric())

#cor(fiftyGenes[49,1:100], fiftyGenes[50, 1:100])
#cor(fiftyGenes[49,101:125], fiftyGenes[50, 101:125])


###################################################
### code chunk number 4: initialization
###################################################
set.seed(3)
initHP1 <- initializeHP(D_diff, tinyCond)
print(initHP)

save.image('./data/hcv.rdata')
oout1 <- ebCoexpressOneStep(D, tinyCond, tinyPat, initHP1)
result11 <- oout1$POSTPROBS
priorDiagnostic(D_diff, tinyCond, oout1, 1)
priorDiagnostic(D_diff, tinyCond, oout1, 2)

ppbDC11 <- 1-result11[,1] 
crit_s1 <- crit.fun(result11[,1], 0.01) # soft thread
kept_s1 <- ppbDC11[ppbDC11 >= crit_s1]

kept_h <- ppbDC11[ppbDC11 >= 0.99]


klabs_s <- names(kept_s)
klabs_s <- ldply(klabs_s, function(x){
  strsplit('sfd~dsf', '~')
})
soft_name <- strsplit(klabs_s, '~')

soft_name <- do.call('rbind', soft_name)
soft_name <- as.data.frame(soft_name)
soft_name$full_name <- klabs_s
D_cor <- as.data.frame(D)
D_cor$full_name <- rownames(D_cor)
hcv_link <- soft_name %>% left_join(D_cor, by = 'full_name')

hcv_link$diff <- hcv_link$Condition1-hcv_link$Condition2
summary(hcv_link$diff %>% abs)
write_csv(hcv_link, path = './data/hcv_link.csv')



save.image('./data/hcv.rdata')
# DC pair names, under soft thresholding
klabs_h <- names(kept_h)
# DC pair names, under hard thresholding
klabs_s %>% length
klabs_h %>% length
###################################################
### code chunk number 15: TPgeneration
