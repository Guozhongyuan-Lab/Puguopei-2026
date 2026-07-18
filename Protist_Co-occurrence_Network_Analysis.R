rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

setwd("XXXXXX")

suppressPackageStartupMessages({
  library(microeco)
  library(ape)
  library(magrittr)
  library(igraph)
  library(rgexf)
  library(RColorBrewer)
  library(WGCNA)
  library(scales)
  library(dplyr)
  library(tibble)
})
options(warn = -1) 

otu_table <- read.csv("XXXX.csv", row.names = 1, sep = ",", fileEncoding = "UTF-8")
taxonomy_table <- read.csv("XXXX.csv", row.names = 1, sep = ",", fileEncoding = "UTF-8")

taxonomy_table %<>% tidy_taxonomy
taxonomy_table[taxonomy_table == "" | taxonomy_table == "unclassified" | taxonomy_table == "未分类"] <- "Unclassified"
colnames(taxonomy_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

mt <- microtable$new(otu_table = otu_table, tax_table = taxonomy_table) 

mt$tax_table %<>% subset(!grepl("Fungi|k__Fungi|真菌", Kingdom, ignore.case = TRUE))
mt$otu_table <- mt$otu_table[rownames(mt$otu_table) %in% rownames(mt$tax_table), ]

mt$filter_pollution(taxa= c("mitochondria", "chloroplast"))

if("filter_taxa" %in% names(mt)){
  try({
    mt$filter_taxa(rel_abund = 0.0005, occurrence = 0.2)
  }, silent = TRUE)
  try({
    mt$filter_taxa(abund_threshold = 0.0005, occur_threshold = 0.2)
  }, silent = TRUE)
}else{
  otu_rel <- apply(mt$otu_table, 1, function(x) sum(x)/sum(mt$otu_table))
  otu_occur <- apply(mt$otu_table, 1, function(x) sum(x>0)/ncol(mt$otu_table))
  keep_otu <- names(otu_rel)[otu_rel >= 0.0005 & otu_occur >= 0.2]
  mt$otu_table <- mt$otu_table[keep_otu, ]
  mt$tax_table <- mt$tax_table[keep_otu, ]
}

mt$tidy_dataset()
set.seed(9527)
rare_size <- max(min(mt$sample_sums()), 5)
mt$rarefy_samples(sample.size = rare_size)

mt$save_table(dirpath = "XXXX", sep = ",")

set.seed(9527) 
t1 <- trans_network$new(
  dataset = mt,
  cor_method = "spearman",
  use_WGCNA_pearson_spearman = TRUE,
  filter_thres = 0.0005
)

cor_matrix <- t1$res_cor_p[[1]]  
p_matrix <- t1$res_cor_p[[2]]    

t1$cal_network(
  COR_p_thres = 0.01,
  COR_p_adjust = "fdr",
  COR_cut = 0.6,
  COR_optimization = TRUE
)

if(ecount(t1$res_network) > 0){ 
  keep_nodes <- V(t1$res_network)$name[degree(t1$res_network) >= 2]
  if(length(keep_nodes) > 0){
    t1$res_network <- induced_subgraph(t1$res_network, keep_nodes)
  }
  
  module_result <- tryCatch({
    t1$cal_module(method = "cluster_louvain")
    module_sizes <- table(V(t1$res_network)$module)
    small_modules <- names(module_sizes)[module_sizes < 3]
    if(length(small_modules) > 0){
      V(t1$res_network)$module[V(t1$res_network)$module %in% small_modules] <- "Other"
    }
    length(unique(V(t1$res_network)$module))
  }, error = function(e){
    V(t1$res_network)$module <- "Module1"
    1
  })
}else{
  all_nodes <- rownames(mt$otu_table)
  if(length(all_nodes) > 1){
    t1$res_network <- make_full_graph(length(all_nodes), directed = FALSE)
    V(t1$res_network)$name <- all_nodes
    V(t1$res_network)$module <- "Module1"
  }else{
    t1$res_network <- make_empty_graph()
    t1$res_network <- add_vertices(t1$res_network, length(all_nodes), name = all_nodes, module = "Module1")
  }
}

t1$get_node_table(node_roles = FALSE)
t1$get_edge_table()

if(nrow(t1$res_node_table) == 0 && vcount(t1$res_network) > 0){
  t1$res_node_table <- data.frame(
    name = V(t1$res_network)$name,
    module = V(t1$res_network)$module,
    degree = degree(t1$res_network),
    row.names = V(t1$res_network)$name
  )
}

node_names <- V(t1$res_network)$name
t1$res_node_table$module <- ifelse(is.na(match(t1$res_node_table$name, node_names)), 
                                   "Module1",
                                   V(t1$res_network)$module[match(t1$res_node_table$name, node_names)])
t1$res_node_table$degree <- ifelse(is.na(match(t1$res_node_table$name, node_names)),
                                   0,
                                   degree(t1$res_network)[match(t1$res_node_table$name, node_names)])
t1$res_node_table[is.na(t1$res_node_table)] <- "Other"  

edge_list <- as_edgelist(t1$res_network)
cor_values <- numeric(nrow(edge_list))
for(i in 1:nrow(edge_list)){
  n1 <- edge_list[i,1]
  n2 <- edge_list[i,2]
  cor_values[i] <- ifelse(n1 %in% rownames(cor_matrix) && n2 %in% colnames(cor_matrix),
                          cor_matrix[n1, n2], 0.6)
}
t1$res_edge_table$correlation <- cor_values
t1$res_edge_table$corr_type <- ifelse(cor_values > 0, "positive", "negative")

mod_pos_count <- list()
mod_neg_count <- list()
unique_mods <- unique(t1$res_node_table$module)

for(m in unique_mods){
  mod_pos_count[[m]] <- 0
  mod_neg_count[[m]] <- 0
}

for(i in 1:nrow(t1$res_edge_table)){
  s_name <- t1$res_edge_table$source[i]
  t_name <- t1$res_edge_table$target[i]
  ct <- t1$res_edge_table$corr_type[i]
  
  m_s <- t1$res_node_table$module[t1$res_node_table$name == s_name]
  m_t <- t1$res_node_table$module[t1$res_node_table$name == t_name]
  if(length(m_s)==0 || length(m_t)==0) next
  if(m_s != m_t) next
  inner_mod <- m_s
  if(ct == "positive"){
    mod_pos_count[[inner_mod]] <- mod_pos_count[[inner_mod]] + 1
  }else{
    mod_neg_count[[inner_mod]] <- mod_neg_count[[inner_mod]] + 1
  }
}

GetGroup <- function(m){
  p <- mod_pos_count[[m]]
  n <- mod_neg_count[[m]]
  if(p >= n){
    return("Group_Pos")
  }else{
    return("Group_Neg")
  }
}
t1$res_node_table$group_type <- sapply(t1$res_node_table$module, GetGroup)
V(t1$res_network)$group_type <- t1$res_node_table$group_type

try({
  t1$save_network(filepath = "XXXX.gexf")
}, silent = TRUE)

write.csv(t1$res_node_table, "XXXX.csv", row.names = FALSE)
write.csv(t1$res_edge_table, "XXXX.csv", row.names = FALSE)
