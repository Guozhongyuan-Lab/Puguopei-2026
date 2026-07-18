rm(list=ls())

library(ggplot2)
library(ggExtra)
library(vegan)
library(ggthemes)
library(readxl)

setwd("XXXXXX")

data_raw <- read_excel("ASV_abundance.xlsx", sheet = 1)  

rownames(data_raw) <- data_raw[[1]]  
data <- data_raw[, -1]  
data <- as.matrix(data)  

data <- data / apply(data, 2, sum)
data <- t(data)  

group_raw <- read_excel("group_info.xlsx", sheet = 1)  

colnames(group_raw)[colnames(group_raw) == "ID"] <- "Sample_ID"  

rownames(group_raw) <- group_raw$Sample_ID
group <- group_raw  

bray <- vegdist(data, method = 'bray')
bray_matrix <- as.matrix(bray)
write.table(bray_matrix, "bray_curtis_matrix.txt", sep = "\t", quote = FALSE)

pcoa <- cmdscale(bray, k = 3, eig = TRUE)  

pcoa_data <- data.frame(pcoa$point)  
pcoa_data$Sample_ID <- rownames(pcoa_data)  
names(pcoa_data)[1:3] <- c("PCoA1", "PCoA2", "PCoA3")  

eig_total <- sum(pcoa$eig)  
eig_percent <- round((pcoa$eig / eig_total) * 100, 1)  

pcoa_result <- merge(pcoa_data, group, by = "Sample_ID")
head(pcoa_result)  

dune.div <- adonis2(data ~ group, data = group, permutations = 999, method = "bray")
dune.div  
dune_adonis <- paste0("Adonis R²: ", round(dune.div$R2[1], 2), "; P-value: ", round(dune.div$`Pr(>F)`[1], 3))
dune_adonis  

p <- ggplot(pcoa_result, aes(x = PCoA1, y = PCoA2, color = group)) +
  geom_point(aes(shape = group), size = 5) +
  labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
       y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
       caption = dune_adonis) +
  scale_colour_manual(values = c("#2f7eb0", "#d00000")) +
  theme(legend.position = c(0.9, 0.19),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", fill = "transparent"),
        axis.text = element_text(color = "black", size = 10)) +
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")

ggsave("pcoa_plot.png", plot = p, width = 5, height = 5, dpi = 600)  
ggsave("pcoa_plot.tiff", plot = p, width = 5, height = 5, dpi = 600)  
ggsave("pcoa_plot.pdf", plot = p, width = 5, height = 5)  

p_with_ellipse <- p + stat_ellipse(data = pcoa_result,
                                   geom = "polygon",
                                   level = 0.9,
                                   linetype = 2,
                                   linewidth = 0.5,
                                   aes(fill = group),
                                   alpha = 0.3,
                                   show.legend = TRUE) +
  scale_fill_manual(values = c("#2f7eb0", "#d00000"))
print(p_with_ellipse)

ggsave("pcoa_ellipse.png", plot = p_with_ellipse, width = 5, height = 5, dpi = 600)  
ggsave("pcoa_ellipse.tiff", plot = p_with_ellipse, width = 5, height = 5, dpi = 600)  
ggsave("pcoa_ellipse.pdf", plot = p_with_ellipse, width = 5, height = 5)  

png(file = "pcoa_marginal.png", width = 5, height = 5, res = 600, units = "in")  
p_marginal <- ggMarginal(
  p_with_ellipse,
  type = "density",
  margins = "both",
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE
)
print(p_marginal)
dev.off()  

tiff(file = "pcoa_marginal.tiff", width = 5, height = 5, res = 600, units = "in")  
p_marginal <- ggMarginal(
  p_with_ellipse,
  type = "density",
  margins = "both",
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE
)
print(p_marginal)
dev.off()  

pdf(file = "pcoa_marginal.pdf", width = 5, height = 5)  
p_marginal <- ggMarginal(
  p_with_ellipse,
  type = "density",
  margins = "both",
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE
)
print(p_marginal)
dev.off()
