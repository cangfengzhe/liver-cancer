function (exprs.1, exprs.2, link.method = c("qth", "rth", "percent")[1], 
          cutoff = 0.25, r.method = c("pearson", "spearman")[1], q.method = c("BH", 
                                                                              "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr")[1], 
          nbins = 20, p = 0.1, figname = c("LFC.s.jpeg", "LFC.d.jpeg")) 
{
  cor.filtered.1 <- cor.filtered.2 <- rth.1 <- rth.2 <- NULL
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
  genes = rownames(exprs.1)
  cor.filtered = switch(link.method, rth = rLinkfilter(exprs.1, 
                                                       exprs.2, r.method = r.method, cutoff = cutoff), qth = qLinkfilter(exprs.1, 
                                                                                                                         exprs.2, r.method = r.method, cutoff = cutoff), percent = percentLinkfilter(exprs.1, 
                                                                                                                                                                                                     exprs.2, r.method = r.method, cutoff = cutoff))
  rth.1 = cor.filtered$rth.1
  rth.2 = cor.filtered$rth.2
  cor.filtered.1 = cor.filtered$cor.filtered.1
  cor.filtered.1[is.na(cor.filtered.1)] = 0
  cor.filtered.2 = cor.filtered$cor.filtered.2
  cor.filtered.2[is.na(cor.filtered.2)] = 0
  idx.same = (cor.filtered.1 * cor.filtered.2) > 0
  idx.same[is.na(idx.same)] <- TRUE
  idx.diff = (cor.filtered.1 * cor.filtered.2) < 0
  idx.diff[is.na(idx.diff)] <- FALSE
  idx.switched = (cor.filtered.1 * cor.filtered.2 < 0) & (abs(cor.filtered.1) >= 
                                                            rth.1 & abs(cor.filtered.2) >= rth.2)
  idx.switched[is.na(idx.switched)] <- FALSE
  cor.same = cbind(cor.filtered.1[idx.same], cor.filtered.2[idx.same])
  rownames(cor.same) <- names(cor.filtered.1)[idx.same]
  cor.switched = cbind(cor.filtered.1[idx.switched], cor.filtered.2[idx.switched])
  rownames(cor.switched) <- names(cor.filtered.1)[idx.switched]
  cor.diff = cbind(cor.filtered.1[idx.diff & (!idx.switched)], 
                   cor.filtered.2[idx.diff & (!idx.switched)])
  rownames(cor.diff) <- names(cor.filtered.1)[idx.diff & (!idx.switched)]
  n.switchedDCL = nrow(cor.switched)
  if (is.null(rownames(cor.same))) {
    name.same = NULL
  }
  if (!is.null(rownames(cor.same))) {
    name.same = strsplit(rownames(cor.same), ",")
    name.same = matrix(unlist(name.same), length(name.same), 
                       2, byrow = T)
  }
  if (is.null(rownames(cor.switched))) {
    name.switched = NULL
  }
  if (!is.null(rownames(cor.switched))) {
    name.switched = strsplit(rownames(cor.switched), ",")
    name.switched = matrix(unlist(name.switched), length(name.switched), 
                           2, byrow = T)
  }
  if (is.null(rownames(cor.diff))) {
    name.diff = NULL
  }
  if (!is.null(rownames(cor.diff))) {
    name.diff = strsplit(rownames(cor.diff), ",")
    name.diff = matrix(unlist(name.diff), length(name.diff), 
                       2, byrow = T)
  }
  name.all = rbind(name.same, name.switched, name.diff)
  if (nrow(cor.same) > 1) {
    de.s = LFC(cor.same, nbins, p, sign = "same", figname = figname[1])
    DCL.same = cor.same[de.s, ]
    name.same = name.same[de.s, ]
    n.sameDCL = nrow(DCL.same)
    DCL.same <- data.frame(name.same, DCL.same)
    colnames(DCL.same) <- c("Gene.1", "Gene.2", "cor.1", 
                            "cor.2")
  }
  else stop("only one or no same-signed pair in all!")
  if (nrow(cor.diff) > 1) {
    de.d = LFC(cor.diff, nbins, p, sign = "diff", figname = figname[2])
    DCL.diff = cor.diff[de.d, ]
    name.diff = name.diff[de.d, ]
    n.diffDCL = nrow(DCL.diff)
    DCL.diff <- data.frame(name.diff, DCL.diff)
    colnames(DCL.diff) <- c("Gene.1", "Gene.2", "cor.1", 
                            "cor.2")
  }
  else stop("only one or no differently-signed pair in all!")
  pairs = rbind(name.same, name.diff, name.switched)
  if (n.switchedDCL > 0) {
    DCL.switched <- data.frame(name.switched, cor.switched)
    colnames(DCL.switched) <- c("Gene.1", "Gene.2", "cor.1", 
                                "cor.2")
    cor.max <- apply(abs(cor.switched), 1, max)
    middle <- sort(cor.max, method = "quick", index.return = TRUE, 
                   decreasing = TRUE)$ix
    DCL.switched <- DCL.switched[middle, ]
  }
  g.all <- graph.data.frame(name.all)
  gene.all <- as.matrix(V(g.all)$name)
  de.all <- degree(g.all)
  g <- graph.data.frame(pairs)
  gene.1 <- as.matrix(V(g)$name)
  de <- degree(g)
  g.same <- graph.data.frame(name.same)
  g.same.name <- as.matrix(V(g.same)$name)
  degree.same <- as.matrix(degree(g.same))
  g.diff <- graph.data.frame(name.diff)
  g.diff.name <- as.matrix(V(g.diff)$name)
  degree.diff <- as.matrix(degree(g.diff))
  if (n.switchedDCL > 0) {
    g.switch <- graph.data.frame(name.switched)
    g.switch.name <- as.matrix(V(g.switch)$name)
    degree.switch <- as.matrix(degree(g.switch))
  }
  else {
    degree.switch = matrix(0, 1, 1)
    DCL.switched = matrix("NULL", 1, 1)
  }
  degree.bind <- matrix(0, m, 5)
  row.names(degree.bind) <- genes
  colnames(degree.bind) <- c("All.links", "DC.links", "DCL.same", 
                             "DCL.diff", "DCL.switched")
  degree.bind[gene.all, 1] = de.all
  degree.bind[gene.1, 2] = de
  degree.bind[g.same.name, 3] = degree.same
  degree.bind[g.diff.name, 4] = degree.diff
  if (n.switchedDCL > 0) {
    degree.bind[g.switch.name, 5] = degree.switch
  }
  prob <- nrow(pairs)/nrow(name.all)
  p.value <- pbinom(degree.bind[, "DC.links"] - 1, degree.bind[, 
                                                               "All.links"], prob, lower.tail = F, log.p = FALSE)
  q.value <- p.adjust(p.value, method = q.method)
  degree.bind <- cbind(degree.bind, p.value, q.value)
  colnames(degree.bind) <- c("All.links", "DC.links", "DCL_same", 
                             "DCL_diff", "DCL_switch", "p", "q")
  middle <- sort(as.numeric(degree.bind[, "q"]), method = "quick", 
                 decreasing = FALSE, index.return = TRUE)$ix
  DCGs <- degree.bind[middle, ]
  DCLs <- rbind(data.frame(DCL.same, type = "same signed"), 
                data.frame(DCL.diff, type = "diff signed"))
  if (n.switchedDCL > 0) 
    DCLs <- rbind(DCLs, data.frame(DCL.switched, type = "switched opposites"))
  DCLs <- data.frame(DCLs, cor.diff = FALSE)
  DCLs[, "cor.diff"] <- abs(DCLs[, "cor.1"] - DCLs[, "cor.2"])
  Result <- list(DCGs = DCGs, DCLs = DCLs)
  return(Result)
}
