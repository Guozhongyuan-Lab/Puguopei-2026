rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ===================== 1. 环境配置 =====================
# 请确保已根据实际情况修改此路径
setwd("D:/A研究生各种安排/AAAAA小论文/A四周极限写作/JHM审后大修/重跑图/NPSS微生物网络图")

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

# ===================== 2. 细菌数据预处理 =====================
otu_table <- read.csv("细菌NPSS.csv", row.names = 1, sep = ",", fileEncoding = "UTF-8")
taxonomy_table <- read.csv("细菌注释表.csv", row.names = 1, sep = ",", fileEncoding = "UTF-8")

taxonomy_table %<>% tidy_taxonomy
taxonomy_table[taxonomy_table == "" | taxonomy_table == "unclassified" | taxonomy_table == "未分类"] <- "Unclassified"
colnames(taxonomy_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

mt <- microtable$new(otu_table = otu_table, tax_table = taxonomy_table) 

# 筛选细菌界
mt$tax_table %<>% subset(grepl("Bacteria", Kingdom, ignore.case = TRUE))
mt$otu_table <- mt$otu_table[rownames(mt$otu_table) %in% rownames(mt$tax_table), ]

# 过滤污染
mt$filter_pollution(taxa= c("mitochondria", "chloroplast"))

# 丰度&出现率过滤
if("filter_taxa" %in% names(mt)){
  try({ mt$filter_taxa(rel_abund = 0.0005, occurrence = 0.2) }, silent = TRUE)
} else {
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

mt$save_table(dirpath = "非塑料际组细菌_Basic_files", sep = ",")

# ===================== 3. 共现网络计算 =====================
set.seed(9527) 
t1 <- trans_network$new(dataset = mt, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0005)
cor_matrix <- t1$res_cor_p[[1]]  

t1$cal_network(COR_p_thres = 0.01, COR_p_adjust = "fdr", COR_cut = 0.6, COR_optimization = TRUE)

# 网络拓扑计算与模块划分
if(ecount(t1$res_network) > 0){ 
  keep_nodes <- V(t1$res_network)$name[degree(t1$res_network) >= 2]
  t1$res_network <- induced_subgraph(t1$res_network, keep_nodes)
  
  tryCatch({
    t1$cal_module(method = "cluster_louvain")
    module_sizes <- table(V(t1$res_network)$module)
    small_modules <- names(module_sizes)[module_sizes < 3]
    V(t1$res_network)$module[V(t1$res_network)$module %in% small_modules] <- "Other"
  }, error = function(e){ V(t1$res_network)$module <- "Module1" })
}

t1$get_node_table(node_roles = FALSE)
t1$get_edge_table()

# 边属性填充
edge_list <- as_edgelist(t1$res_network)
cor_values <- sapply(1:nrow(edge_list), function(i) cor_matrix[edge_list[i,1], edge_list[i,2]])
t1$res_edge_table$correlation <- cor_values
t1$res_edge_table$corr_type <- ifelse(cor_values > 0, "positive", "negative")

# 模块正负连线统计与分组
mod_pos_count <- setNames(numeric(length(unique(t1$res_node_table$module))), unique(t1$res_node_table$module))
mod_neg_count <- mod_pos_count

for(i in 1:nrow(t1$res_edge_table)){
  m_s <- t1$res_node_table$module[t1$res_node_table$name == t1$res_edge_table$source[i]]
  m_t <- t1$res_node_table$module[t1$res_node_table$name == t1$res_edge_table$target[i]]
  if(m_s == m_t) {
    if(t1$res_edge_table$corr_type[i] == "positive") mod_pos_count[m_s] %<>% +1
    else mod_neg_count[m_s] %<>% +1
  }
}

t1$res_node_table$group_type <- sapply(t1$res_node_table$module, function(m) {
  if(mod_pos_count[m] >= mod_neg_count[m]) "Group_Pos" else "Group_Neg"
})

# ===================== 输出文件 =====================
t1$save_network(filepath = "NPSS细菌_network.gexf")
write.csv(t1$res_node_table, "Gephi_节点表_NPSS细菌.csv", row.names = FALSE)
write.csv(t1$res_edge_table, "Gephi_边表_NPSS细菌.csv", row.names = FALSE)
