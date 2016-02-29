head(bb)
View(head(bb))
library(graph)
library(igraph)
gg <- graph.adjacency(bb, mode = 'undirected', weighted = T)
save(gg, bb, file = './data/matrix_data.rdata')

gg1 <- graph.adjacency(bb, mode = 'undirected', weighted = T)

eb <- edge.betweenness.community(gg1, directed = F)
wc <- walktrap.community(gg1)
membership(eb)

fc <- fastgreedy.community(gg1) 
wc <- walktrap.community(gg1)
 membership(wc) %>% table
sc <- spinglass.community(gg1,  weights = E(gg1)$weight)

imc <- infomap.community(g)
membership(imc)
communities(imc)
# K = igraph.to.graphNEL(gg)
# save(K, file = 'K.rdata')
# maximalCliques = maxClique(K)
# igraph.ad
adjec <- get.adjacency()
