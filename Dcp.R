function (exprs.1, exprs.2, r.method = c("pearson", "spearman")[1], 
          link.method = c("qth", "rth", "percent")[1], cutoff = 0.25, 
          N = 0, N.type = c("pooled", "gene_by_gene")[1], q.method = c("BH", 
                                                                       "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr")[1]) 
{
  cor.filtered.1 <- cor.filtered.2 <- NULL
  if (nrow(exprs.1) != nrow(exprs.2)) 
    stop("the two expression matrices must have the same number of rows (genes).")
  if (length(rownames(exprs.1)) == 0 | length(rownames(exprs.2)) == 
      0) 
    stop("the expression matrices must have row names specifying the gene names.")
  if (min(ncol(exprs.1), ncol(exprs.2)) < 3) {
    stop("each expression matrix must have at least three or more columns.")
  }
  else if (min(ncol(exprs.1), ncol(exprs.2)) < 5) {
    warning("the minimum number of columns is less than five and the result may not be reliable.")
  }
  m <- nrow(exprs.1)
  if (m > 5000) 
    warning("the number of genes exceeds 5000 and the program may takes long time to run.")
  genes = rownames(exprs.1)
  # cor.filtered  过滤方法筛选之后的相关性
  cor.filtered = switch(link.method, rth = rLinkfilter(exprs.1, 
                                                       exprs.2, r.method = r.method, cutoff = cutoff), qth = qLinkfilter(exprs.1, 
                                                                                                                         exprs.2, r.method = r.method, cutoff = cutoff), percent = percentLinkfilter(exprs.1, 
                                                                                                                                                                                                     exprs.2, r.method = r.method, cutoff = cutoff))
  rth.1 = cor.filtered$rth.1
  rth.2 = cor.filtered$rth.2
  cor.filtered.1 = cor.filtered$cor.filtered.1 # 相关性
  cor.filtered.2 = cor.filtered$cor.filtered.2
  calc.dC <- function(cor.filtered.1, cor.filtered.2, genes) {
    nzero.vec <- (cor.filtered.1 != 0) | (cor.filtered.2 != 
                                            0)
    nzero.sm <- diag(rep(0, length(genes))) # 对角为1
    nzero.sm[lower.tri(nzero.sm, diag = F)] <- nzero.vec # 填充matrix下半部分
    nzero.sm = nzero.sm + t(nzero.sm) # 补全matrix
    number_uniq <- apply(nzero.sm, 1, sum) # 横行的和，每一个基因对应相关性大于阈值的基因的个数
    squares = (cor.filtered.1 - cor.filtered.2)^2 # 相关性的平方差
    squares.sm <- diag(rep(0, length(genes)))
    squares.sm[lower.tri(squares.sm, diag = F)] <- squares # 平方差填充
    squares.sm = squares.sm + t(squares.sm)
    ss = apply(squares.sm, 1, sum)
    LNED.result = as.vector(matrix(NA, length(genes), 1))
    LNED.result[number_uniq != 0] = sqrt(ss[number_uniq != 
                                              0])/sqrt(number_uniq[number_uniq != 0]) # 非0的number_uniq对应的相关性的平方差 开方 除以 对应基因数目开方
    names(LNED.result) <- genes
    list(dC = LNED.result, length = number_uniq)
  }
  dC.length = calc.dC(cor.filtered.1, cor.filtered.2, genes)
  dC = dC.length$dC
  number_uniq = dC.length$length
  if (N > 0) {
    dC0 <- matrix(nrow = length(genes), ncol = N)
    rownames(dC0) <- genes
    exprs <- cbind(exprs.1, exprs.2)
    expSamples <- colnames(exprs)
    n.1 = ncol(exprs.1)
    n.2 = ncol(exprs.2)
    cat.j = 0
    for (j in 1:N) { # 扰动的次数
      if ((j * 100/N)%/%10 > cat.j) {
        cat.j = cat.j + 1
        cat(cat.j * 10, "%", "\n")
      }
      seq <- sample(n.1 + n.2)
      exprs.1 <- exprs[, seq[1:n.1]]
      exprs.2 <- exprs[, seq[(n.1 + 1):(n.1 + n.2)]]
      rownames(exprs.1) <- rownames(exprs.2) <- genes
      cor.filtered = switch(link.method, rth = rLinkfilter(exprs.1, 
                                                           exprs.2, r.method = r.method, cutoff = cutoff), 
                            qth = qLinkfilter(exprs.1, exprs.2, r.method = r.method, 
                                              cutoff = cutoff), percent = percentLinkfilter(exprs.1, 
                                                                                            exprs.2, r.method = r.method, cutoff = cutoff))
      dC0.j <- calc.dC(cor.filtered$cor.filtered.1, cor.filtered$cor.filtered.2, 
                       cor.filtered$genes)
      dC0[, j] = dC0.j$dC
    }
    p.value = switch(N.type, gene_by_gene = apply(cbind(dC0, 
                                                        dC), 1, function(x) sum(x[1:(length(x) - 1)] > x[length(x)], # 随机扰动的值与真实值比较
                                                                                na.rm = T)/length(!is.na(x[1:(length(x) - 1)]))), 
                     pooled = sapply(dC, function(x) sum(as.vector(dC0) > 
                                                           x, na.rm = T)/length(!is.na(as.vector(dC0)))))
    q.value <- p.adjust(p.value, method = q.method)
    Result <- data.frame(dC = dC, links = number_uniq, p.value = p.value, 
                         q.value = q.value)
    row.names(Result) <- genes
  }
  else {
    Result <- data.frame(dC = dC, links = number_uniq, p.value = rep(NA, 
                                                                     length(dC)), q.value = rep(NA, length(dC)))
    row.names(Result) <- genes
  }
  return(Result)
}
