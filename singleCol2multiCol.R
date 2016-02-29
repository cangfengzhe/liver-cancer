# change one column to multiple column
raw_data <- read_csv('./data/hcv_node.csv') %>% as.data.frame()
head(raw_data)
hcv_node_cluster <- singleCol2multiCol(raw_data, '__glayCluster')
singleCol2multiCol <- function(raw_data, 1){
  
cluster_size <- unique(raw_data[,1])
multi_column <- sapply(1: length(cluster_size), function(ii){
  raw_data$name[which(raw_data[,1]==ii)]
})
length(multi_column)
max_length <- sapply(multi_column, length) %>% max
cluster_lt20 <- NA
n <- 0
complecation <- ldply(multi_column, function(x){
    if(length(x)<20) {
      n <<- n+1
      cluster_lt20[n] <<- n
    }
    if(length(x) < max_length){
      x <-  c(x, rep(NA, max_length - length(x)))
    
    return(x)
    }
  
})
complecation
}
hcv_node_cluster <- t(hcv_node_cluster)
hcv_node_cluster_filter <- sapply(1:ncol(hcv_node_cluster), function(ii){
  
})
write_csv(hcv_node_cluster %>% as.data.frame(), path='./data/hcv_node_cluster.csv')
hcv_node_cluster %>% class
