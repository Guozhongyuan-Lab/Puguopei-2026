rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ===================== 1. 环境配置 =====================
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

# ===================== 2. 真菌数据预处理 =====================
otu_table <- read.csv("真菌NPSS.csv", row.names = 1, sep = ",", fileEncoding = "UTF-8")
taxonomy_table <- read.csv("真菌注释表.csv", row.names = 1, sep = ",", fileEncoding = "UTF-8")

taxonomy_table %<>% tidy_taxonomy
taxonomy_table[taxonomy_table == "" | taxonomy_table == "unclassified" | taxonomy_table == "未分类"] <- "Unclassified"
colnames(taxonomy_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

mt <- microtable$new(otu_table = otu_table, tax_table = taxonomy_table) 

# 筛选真菌界
mt$tax_table %<>% subset(grepl("Fungi|真菌", Kingdom, ignore.case = TRUE))
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

# ===================== 3. 共现网络计算 =====================
set.seed(9527) 
t1 <- trans_network$new(dataset = mt, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE, filter_thres = 0.0005)
cor_matrix <- t1$res_cor_p[[1]] 

t1$cal_network(COR_p_thres = 0.01, COR_p_adjust = "fdr", COR_cut = 0.6, COR_optimization = TRUE)

# 拓扑过滤与模块划分
if(ecount(t1$res_network) > 0){ 
  keep_nodes <- V(t1$res_network)$name[degree(t1$res_network) >= 2]
  t1$res_network <- induced_subgraph(t1$res_network, keep_nodes)
  
  tryCatch({
    t1$cal_module(method = "cluster_louvain")
    module_sizes <- table(V(t1$res_network)$module)
    V(t1$res_network)$module[V(t1$res_network)$module %in% names(module_sizes)[module_sizes < 3]] <- "Other"
  }, error = function(e){ V(t1$res_network)$module <- "Module1" })
}

# 导出并处理节点属性
t1$get_node_table(node_roles = FALSE)
t1$get_edge_table()

# 完善节点表
t1$res_node_table$degree <- degree(t1$res_network)[match(t1$res_node_table$name, V(t1$res_network)$name)]
t1$res_node_table$module[is.na(t1$res_node_table$module)] <- "Other"

# 计算边相关性类型
cor_values <- sapply(1:nrow(t1$res_edge_table), function(i) {
  cor_matrix[t1$res_edge_table$source[i], t1$res_edge_table$target[i]]
})
t1$res_edge_table$correlation <- cor_values
t1$res_edge_table$corr_type <- ifelse(cor_values > 0, "positive", "negative")

# 自动划分 Group_Pos / Group_Neg
mod_stats <- t1$res_edge_table %>% 
  filter(source %in% t1$res_node_table$name & target %in% t1$res_node_table$name) %>%
  mutate(mod_s = t1$res_node_table$module[match(source, t1$res_node_table$name)],
         mod_t = t1$res_node_table$module[match(target, t1$res_node_table$name)]) %>%
  filter(mod_s == mod_t) %>%
  group_by(mod_s, corr_type) %>%
  summarise(count = n(), .groups = 'drop')

t1$res_node_table$group_type <- sapply(t1$res_node_table$module, function(m) {
  pos <- sum(mod_stats$count[mod_stats$mod_s == m & mod_stats$corr_type == "positive"], na.rm = TRUE)
  neg <- sum(mod_stats$count[mod_stats$mod_s == m & mod_stats$corr_type == "negative"], na.rm = TRUE)
  ifelse(pos >= neg, "Group_Pos", "Group_Neg")
})

# ===================== 输出文件 =====================
t1$save_network(filepath = "NPSS真菌_network.gexf")
write.csv(t1$res_node_table, "Gephi_节点表_NPSS真菌.csv", row.names = FALSE)
write.csv(t1$res_edge_table, "Gephi_边表_NPSS真菌.csv", row.names = FALSE)

cat("\n所有网络分析已完成。下一步请在 Gephi 中使用 ForceAtlas 2 布局，并通过 group_type 进行节点颜色分区。\n")
