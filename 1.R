# 统计扩增子测序结果
library(ggplot2)
setwd("D:/大亚湾/18S.paired.report/report/OTU_Cluster_Taxonomy/1.OTU_Cluster_Taxonomy")
otu_data = read.csv("OTU_taxonomy.csv", header = T,
                    check.names = F)
colnames(otu_data)
otu_state = data.frame()
j = 1
for (i in unique(otu_data$Phylum)) {
  otu_state[j, 1] = i
  otu_state[j, 2] = length(subset(otu_data, Phylum == i)$Phylum)
  otu_state[j, 3] = sum(subset(otu_data, Phylum == i)$Abundance)/23978326
  j = j + 1
}
x11()

colnames(otu_state) = c("Phylum", "OTU count", "OTU sum")
p = ggplot(data = otu_state,
           aes(x = reorder(Phylum, -`OTU count`),
               y = `OTU count`)) +
  geom_col() +
  labs(x = "Phylum", y = "OTU count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10))
p

ggsave("phylum_otu.tiff", dpi = 300)

sum(otu_data$Abundance)

p = ggplot(data = otu_state,
           aes(x = reorder(Phylum, -`OTU sum`),
               y = `OTU sum`)) +
  geom_col() +
  labs(x = "Phylum", y = "reads number") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10))
p

ggsave("phylum_reads_number.tiff", dpi = 300)
