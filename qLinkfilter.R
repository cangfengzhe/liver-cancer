function (exprs.1, exprs.2, cutoff = 0.25, r.method = c("pearson", 
                                                        "spearman")[1], q.method = c("BH", "holm", "hochberg", "hommel", 
                                                                                     "bonferroni", "BY", "fdr")[1]) 
{
  degree.1 <- ncol(exprs.1) - 2
  degree.2 <- ncol(exprs.2) - 2
  genes <- rownames(exprs.1)
  exprs.1 <- as.matrix(exprs.1)
  exprs.2 <- as.matrix(exprs.2)
  cor.1 <- cor(t(exprs.1), method = r.method, use = "pairwise.complete.obs")
  cor.2 <- cor(t(exprs.2), method = r.method, use = "pairwise.complete.obs")
  cor.1 <- cor.1[lower.tri(cor.1, diag = F)]
  cor.2 <- cor.2[lower.tri(cor.2, diag = F)]
  rm(exprs.1)
  rm(exprs.2)
  t.1 <- cor.1 * sqrt(degree.1)/sqrt(1 - cor.1 * cor.1)
  t.2 <- cor.2 * sqrt(degree.2)/sqrt(1 - cor.2 * cor.2)
  p0.1 <- 2 * pt(-abs(t.1), degree.1, lower.tail = TRUE, log.p = FALSE)
  p0.2 <- 2 * pt(-abs(t.2), degree.2, lower.tail = TRUE, log.p = FALSE)
  rm(t.1)
  rm(t.2)
  q.1 <- p.adjust(p0.1, method = q.method) # 校正pvalue
  q.2 <- p.adjust(p0.2, method = q.method)
  rm(p0.1)
  rm(p0.2)
  rth.1 <- abs(cor.1[which.min(abs(q.1 - cutoff))])
  rth.2 <- abs(cor.2[which.min(abs(q.2 - cutoff))])
  cor.1[q.1 >= cutoff & q.2 >= cutoff] <- cor.2[q.1 >= cutoff & 
                                                  q.2 >= cutoff] <- 0 # 大于阈值的归零
  name.row <- matrix(rep(genes, length(genes)), length(genes), 
                     length(genes))
  name.col <- matrix(rep(genes, length(genes)), length(genes), 
                     length(genes), byrow = T)
  name.pairs <- matrix(paste(name.row, name.col, sep = ","), 
                       length(genes), length(genes))
  rm(list = c("name.row", "name.col"))
  name.pairs <- name.pairs[lower.tri(name.pairs, diag = F)]
  names(cor.1) <- names(cor.2) <- name.pairs
  cor.filtered <- list(rth.1 = rth.1, rth.2 = rth.2, cor.filtered.1 = cor.1, 
                       cor.filtered.2 = cor.2, genes = genes)
  return(cor.filtered)
}
