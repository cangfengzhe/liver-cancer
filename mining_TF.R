library(RCurl)
library(rvest)
library(RSQLite)

tf_list <- read.table('./raw_data/tf_list.txt', stringsAsFactors = F)
sql_db <- dbConnect(SQLite(), './data/sqlite.db')

url_str_pre <- 'https://cb.utdallas.edu/cgi-bin/TRED/tred.cgi?process=searchTFGene&sel_type=factor_name&factor_organism=any&tx_search_terms='
url_str_post <- '&target_organism=human&prom_quality=0&bind_quality=0&submit=SEARCH'
n <- 0
error_out <- NA
tf_genes <- ldply(tf_list[,1],.progress = 'text', function(x){
  tryCatch({
    print(x)
    url_str <- paste(c(url_str_pre, x, url_str_post), collapse = '')
    table_data <- html(url_str) %>% html_nodes('table')
    if(table_data[[3]] %>% 
       html_table() %>% length() == 1){
      return(data.frame(name=x, X2=NA, X6 = NA, X7 = NA))
    }
    gene_num <- table_data[[3]] %>% 
      html_table() %>% .[1,2] %>% strsplit(' ') %>% .[[1]] %>% .[5] %>% as.numeric()
    if(gene_num>20){
      gene_df <- gene_page(tf = x, gene_num)
    }else {
      gene_data <-  table_data[[4]] %>%  html_table(fill = F)
#       gene_name <- gene_data[2: nrow(gene_data), 2]
#       gene_pq <- gene_data[2: nrow(gene_data), 6]
#       gene_bq <- gene_data[2: nrow(gene_data), 7]
      gene_full <- gene_data[2: nrow(gene_data), c(2, 6, 7)]
      # gene_df <- data.frame(name=x, gene = gene_name, pq = gene_pq, bq = gene_bq)
      gene_df <- data.frame(name=x, gene_full )
      
    }
    gene_df
    dbWriteTable(sql_db, 'tf_gene', gene_df, append = T)
  }, error = function(e){
    n <<- n+1;
    error_out[n] <<- x 
    print(e)
  })
  
})




page_url_1 <- 'https://cb.utdallas.edu/cgi-bin/TRED/tred.cgi?hitCount='

page_url_2 <- '&process=searchTFGene&factor_organism=any&target_organism=human&sel_type=factor_name&prom_quality=0&bind_quality=0&tx_search_terms='

page_url_3 <- '&start='

page_url_4 <- '&anotherPage=GO'
gene_page <- function(tf, gene_num){
  tf_gene <- ldply(1:ceiling(gene_num/20), function(ii){
    start_num <- (ii-1)*20
    page_url <- paste(c(page_url_1, gene_num, page_url_2, tf, page_url_3, start_num, page_url_4), collapse = '')
    table_data <- html(page_url) %>% html_nodes('table')
    gene_data <- table_data[[4]] %>%  html_table(fill = F)
    out <- gene_data[2: nrow(gene_data), c(2, 6, 7)] 
    out <- data.frame(name = tf,  out)
    out
  })
  
  tf_gene
}




# tf_genes <- unique(tf_genes)
# save(tf_genes, file = './data/tf_genes.rdata')

