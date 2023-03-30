setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")

library(stringr)
library(vegan)
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(lemon)
library(igraph)
library(ggrepel)
library(ggprism)
library(ggsignif)
library(ggpubr)


# 导入大亚湾的数据
daya_env = read.csv("env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann = read.csv("OTU_taxonomy.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)
colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")

# 选出属于浮游植物的子集
phy_ann_data = subset(daya_ann, Phylum %in% c("Dinoflagellata",
                                              "Ochrophyta",
                                              "Chlorophyta",
                                              "Haptophyta",
                                              "Cryptophyta"))

phy_ann_data$OTU_id = rownames(phy_ann_data)

daya_abu = read.csv("OTU_shared_final.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
colnames(daya_abu)[1] = "OTU_id"

# 整合丰度数据与注释结果
phy_ann_abu = merge(phy_ann_data[c("OTU_id", "Class", "Genus")], 
                    daya_abu,
                    by = "OTU_id")
colSums(phy_ann_abu[3:294])
# 计算不同样品的多样性
count_data = as.data.frame(t(phy_ann_abu[3:294]))
observed_species = estimateR(count_data)[1,]
chao1 = estimateR(count_data)[2,]
shannon <- diversity(count_data, index = 'shannon', base = exp(1))
simpson <- diversity(count_data, index = 'simpson', base = exp(1))
divs = cbind.data.frame(observed_species, chao1, shannon, simpson)

divs["Stations"] = gsub("_", "-", rownames(divs))
divs["Stations"] = gsub("-3um-[1,2,3]", "", divs$Stations)
divs["Stations"] = gsub("-0.2um-[1,2,3]", "", divs$Stations)

divs$Sample_id = rownames(divs)
divs["Fraction"] = ""

for (i in c(1:length(divs$Fraction))) {
  if (str_count(rownames(divs)[i], "3um") > 0) {
    divs$Fraction[i] = "3-200 μm"
  }
  if (str_count(rownames(divs)[i], "0.2um") > 0) {
    divs$Fraction[i] = "0.2-3 μm"
  }
}

divs = merge(divs, daya_env[c("Stations", "Temperature")], by = "Stations")

divs_long = melt(divs, measure.vars = c("observed_species", "chao1", 
                                        "shannon", "simpson"))

x11()
windowsFonts(Arial = windowsFont("Arial"))
p = ggplot(divs_long, aes(x = Temperature,
                          y = value,
                          col = variable)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F, size = 2) +
  facet_wrap(Fraction~variable, scales = "free_y", nrow = 2) +
  theme(legend.position = "None",
        strip.text = element_text(size = 12, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
  

p
ggsave("daya_divs.tiff", dpi = 300, height = 6, width = 12)


rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL

phy_ann_genus_abu = aggregate(phy_ann_abu[3:294],
                              list(Class = phy_ann_abu$Class,
                                   Genus = phy_ann_abu$Genus),
                                sum)

phy_ann_genus_abu$Class[which(phy_ann_genus_abu$Class == "")] = "Unknown"
phy_ann_genus_abu$Genus[which(phy_ann_genus_abu$Genus == "")] = "Unknown"


rownames(phy_ann_genus_abu) = paste(phy_ann_genus_abu$Class,
                                    phy_ann_genus_abu$Genus,
                                    sep = " ")

phy_ann_genus_abu$Class = NULL
phy_ann_genus_abu$Genus = NULL

LF_phy_ann_genus_abu = phy_ann_genus_abu[grepl(pattern = "3um", 
                                               colnames(phy_ann_genus_abu))]

SF_phy_ann_genus_abu = phy_ann_genus_abu[grepl(pattern = "0.2um", 
                                               colnames(phy_ann_genus_abu))]


LF_phy_ann_genus_abu = as.data.frame(t(LF_phy_ann_genus_abu))
SF_phy_ann_genus_abu = as.data.frame(t(SF_phy_ann_genus_abu))


LF_phy_ann_genus_abu = LF_phy_ann_genus_abu/rowSums(LF_phy_ann_genus_abu)
SF_phy_ann_genus_abu = SF_phy_ann_genus_abu/rowSums(SF_phy_ann_genus_abu)

LF_phy_ann_genus_abu = as.data.frame(t(LF_phy_ann_genus_abu))
SF_phy_ann_genus_abu = as.data.frame(t(SF_phy_ann_genus_abu))

LF_phy_ann_genus_abu$Sum = rowSums(LF_phy_ann_genus_abu)
SF_phy_ann_genus_abu$Sum = rowSums(SF_phy_ann_genus_abu)

attach(LF_phy_ann_genus_abu)
attach(SF_phy_ann_genus_abu)
LF_phy_ann_genus_abu = LF_phy_ann_genus_abu[order(LF_phy_ann_genus_abu$Sum,
                                                  decreasing = T),]
SF_phy_ann_genus_abu = SF_phy_ann_genus_abu[order(SF_phy_ann_genus_abu$Sum,
                                                  decreasing = T),]
detach(LF_phy_ann_genus_abu)
detach(SF_phy_ann_genus_abu)

LF_phy_ann_genus_abu$Sum = NULL
SF_phy_ann_genus_abu$Sum = NULL



LF_phy_ann_genus_abu = as.data.frame(t(LF_phy_ann_genus_abu))

LF_phy_ann_genus_abu_top = LF_phy_ann_genus_abu[c("Dinophyceae Heterocapsa",
                                                  "Dinophyceae Dinophyceae_XXX",
                                                  "Dinophyceae Pyrophacus",
                                                  "Dinophyceae Tripos",
                                                  "Bacillariophyta Thalassiosira",
                                                  "Bacillariophyta Polar-centric-Mediophyceae_X",
                                                  "Dinophyceae Fragilidium",
                                                  "Dinophyceae Gonyaulax",
                                                  "Noctilucophyceae Noctiluca")]


SF_phy_ann_genus_abu = as.data.frame(t(SF_phy_ann_genus_abu))
SF_phy_ann_genus_abu_top = SF_phy_ann_genus_abu[c("Mamiellophyceae Ostreococcus",
                                                  "Mamiellophyceae Micromonas",
                                                  "Mamiellophyceae Bathycoccus",
                                                  "Syndiniales Dino-Group-I-Clade-1_X",
                                                  "Syndiniales Dino-Group-II-Clade-6_X",
                                                  "Syndiniales Dino-Group-III_XX",
                                                  "Syndiniales Dino-Group-I-Clade-5_X",
                                                  "Syndiniales Dino-Group-II-Clade-14_X",
                                                  "Trebouxiophyceae Picochlorum")]
LF_phy_ann_genus_abu_top$Stations = gsub("_", "-", 
                                         rownames(LF_phy_ann_genus_abu_top))
LF_phy_ann_genus_abu_top$Stations = gsub("-3um-[1,2,3]", "",
                                         LF_phy_ann_genus_abu_top$Stations)

SF_phy_ann_genus_abu_top$Stations = gsub("_", "-", 
                                         rownames(SF_phy_ann_genus_abu_top))
SF_phy_ann_genus_abu_top$Stations = gsub("-0.2um-[1,2,3]", "",
                                         SF_phy_ann_genus_abu_top$Stations)


#### 
LF_phy_ann_genus_abu_top = merge(daya_env[c("Stations", "Temperature")],
                                 LF_phy_ann_genus_abu_top,
                                 by = "Stations")

LF_phy_ann_genus_abu_top_long = melt(LF_phy_ann_genus_abu_top,
                                     id.vars = c("Stations", "Temperature"))

LF_phy_ann_genus_abu_top_long[c("Phylum", "Genus")] = str_split_fixed(LF_phy_ann_genus_abu_top_long$variable,
                                                                      pattern = " ",
                                                                      2)

SF_phy_ann_genus_abu_top = merge(daya_env[c("Stations", "Temperature")],
                                 SF_phy_ann_genus_abu_top,
                                 by = "Stations")

SF_phy_ann_genus_abu_top_long = melt(SF_phy_ann_genus_abu_top,
                                     id.vars = c("Stations", "Temperature"))

SF_phy_ann_genus_abu_top_long[c("Phylum", "Genus")] = str_split_fixed(SF_phy_ann_genus_abu_top_long$variable,
                                                                      pattern = " ",
                                                                      2)

x11()
p = ggplot(LF_phy_ann_genus_abu_top_long,
           aes(x = Temperature,
               y = value,
               col = Phylum)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F, size = 2) +
  facet_wrap(~Genus, scales = "free_y") +
  theme(legend.position = "top",
        strip.text = element_text(size = 8, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 8))
p
p1 = ggplot(SF_phy_ann_genus_abu_top_long,
           aes(x = Temperature,
               y = value,
               col = Phylum)) +
  scale_color_brewer(palette = "Set1") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F, size = 2) +
  facet_wrap(~Genus, scales = "free_y") +
  theme(legend.position = "top",
        strip.text = element_text(size = 8, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 8))
p1

p2 = grid_arrange_shared_legend(p, p1, ncol = 2, 
                                position = "top",
                                plot = TRUE)
p2 = p + p1

p2

ggsave("daya_genus.tiff", dpi = 300, width = 16, height = 8)


# 导入Tara的数据
tara_env = read.csv("tara_env.txt", sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
tara_phy = read.csv("phy_data.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
tara_phy_abu = tara_phy[10:343]

tara_phy_abu = as.data.frame(t(tara_phy_abu))
observed_species = estimateR(tara_phy_abu)[1,]
chao1 = estimateR(tara_phy_abu)[2,]

tara_phy_abu[1:5, 1:5]

shannon = vegan::diversity(tara_phy_abu, index = 'shannon', base = exp(1))
simpson = vegan::diversity(tara_phy_abu, index = 'simpson', base = exp(1))
divs = cbind.data.frame(observed_species, chao1, shannon, simpson)
divs$Barcode = rownames(divs)

tara_env_divs = merge(tara_env[c("Barcode", "Temperature",
                                 "Depth_Nominal", "Fraction_Lower",
                                 "Fraction_Upper")],
                      divs,
                      by = "Barcode")


tara_env_divs$Fraction = paste(tara_env_divs$Fraction_Lower,
                               tara_env_divs$Fraction_Upper,
                              sep = "-")
tara_env_divs$Fraction = paste(tara_env_divs$Fraction,
                               "μm")
unique(tara_env_divs$Fraction)
tara_env_divs_2 = subset(tara_env_divs, Fraction %in% c("0.8-5.00 μm",
                                                        "5-20.00 μm",
                                                        "20-180.00 μm"))

tara_env_divs_2 = subset(tara_env_divs_2, Depth_Nominal == "SRF")

tara_divs_long = melt(tara_env_divs_2, measure.vars = c("observed_species", "chao1", 
                                        "shannon", "simpson"))

tara_divs_long = subset(tara_divs_long, Temperature > 23&Temperature < 30)

# x11()
tara_divs_long$Fraction = factor(tara_divs_long$Fraction,
                                 levels = c("0.8-5.00 μm",
                                            "5-20.00 μm",
                                            "20-180.00 μm"))

p = ggplot(tara_divs_long, aes(x = Temperature,
                               y = value,
                               col = variable)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F) +
  facet_wrap(Fraction~variable, scales = "free_y", nrow = 3) +
  theme(legend.position = "None",
        strip.text = element_text(size = 10, family = "Arial", face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


p
ggsave("tara_divs.tiff", dpi = 300, height = 6, width = 8)

# Tara Ocean 3um key genus
# 3um: "Dinophyceae Heterocapsa", "Dinophyceae Dinophyceae_XXX", "Dinophyceae Pyrophacus",
# 0.2um: "Mamiellophyceae Ostreococcus", "Mamiellophyceae Micromonas", "Mamiellophyceae Bathycoccus"

tara_phy[1:5, 1:11]

target_sp = c("Dinophyceae Heterocapsa", "Dinophyceae Dinophyceae_XXX", 
              "Dinophyceae Pyrophacus", "Mamiellophyceae Ostreococcus", 
              "Mamiellophyceae Micromonas", "Mamiellophyceae Bathycoccus")

tara_phy = read.csv("phy_data.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

"Dinophyceae" %in% tara_phy$Class
"Tripos" %in% tara_phy$Genus


tara_phy$Class_Genus = paste(tara_phy$Class, tara_phy$Genus)
tara_phy$Class_Genus[1:5]

tara_phy = subset(tara_phy, Class_Genus %in% target_sp)

length(tara_phy$Domain)
length(tara_phy)

tara_phy_sum = aggregate(tara_phy[10:343],
                         by = list(Genus = tara_phy$Class_Genus),
                         sum)

rownames(tara_phy_sum) = tara_phy_sum$Genus

tara_phy_sum$Genus = NULL

tara_phy_sum = as.data.frame(t(tara_phy_sum))



# 比较Tara的几个物种随温度的变化
setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S/atlas")
p <- dir(full.names = T)
file_names = data.frame(str_split_fixed(p, "_", 6))
file_names$filename = p
file_names = file_names[5:7]
colnames(file_names) = c("Class", "Genus", "file_path")



data = data.frame(matrix(ncol = 7))
colnames(data) = c("Barcode", 
                   "Depth_Nominal",
                   "Fraction_Lower (μm)",
                   "Fraction_Upper (μm)",
                   "Temperature (°C)",
                   "sum",
                   "Genus")


for (i in c(1:length(file_names$Class))) {
  title = paste(file_names$Class[i], file_names$Genus[i])
  abu_file_path = paste(file_names$file_path[i],
                        "abundance_matrix.csv", sep = "/")
  env_file_path = paste(file_names$file_path[i],
                        "environmental_parameters.csv", sep = "/")
  abu_data = read.csv(abu_file_path, sep = "\t", header = T,
                      check.names = F, encoding = "UTF-8",
                      row.names = 1)
  abu_data$taxonomy = NULL
  abu_data = data.frame(t(abu_data))
  abu_data$sum = rowSums(abu_data)
  abu_data$Barcode = rownames(abu_data)
  env_data = read.csv(env_file_path, sep = "\t", header = T,
                      check.names = F, encoding = "UTF-8")
  temp_data = merge(env_data[c("Barcode", 
                               "Depth_Nominal",
                               "Fraction_Lower (μm)",
                               "Fraction_Upper (μm)",
                               "Temperature (°C)")],
                    abu_data[c("Barcode", "sum")],
                    by = "Barcode")
  temp_data$Genus = title
  data = rbind(data, temp_data)
}

data$Fraction = paste(data$`Fraction_Lower (μm)`,
                      data$`Fraction_Upper (μm)`,
                      sep = "-")
data$Fraction = paste(data$Fraction,
                      "μm")


srf_data = subset(data, Depth_Nominal == "SRF"&
                    `Temperature (°C)`>23&
                    Fraction %in% c("0.8-5 μm",
                                    "5-20 μm",
                                    "20-180 μm"))

# X11()
srf_data[c("Class", "Genus")] = str_split_fixed(srf_data$Genus,
                                                pattern = " ",
                                                2)
srf_data$Fraction = factor(srf_data$Fraction, 
                           levels = c("0.8-5 μm",
                                      "5-20 μm",
                                      "20-180 μm"))

p = ggplot(srf_data,
           aes(x = `Temperature (°C)`,
               y = sum,
               col = Class)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F, size = 2) +
  facet_wrap(Genus~Fraction, scales = "free_y") +
  theme(legend.position = "top",
        strip.text = element_text(size = 10, family = "Arial", face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10))
p

ggsave("tara_genus.tiff", dpi = 300, width = 16, height = 8)

# 滑动分组求网络参数
setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")

daya_env = read.csv("env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann = read.csv("OTU_taxonomy.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)


colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")
# 选出属于浮游植物的子集
phy_ann_data = subset(daya_ann, Phylum %in% c("Dinoflagellata",
                                              "Ochrophyta",
                                              "Chlorophyta",
                                              "Haptophyta",
                                              "Cryptophyta"))

phy_ann_data$OTU_id = rownames(phy_ann_data)

daya_abu = read.csv("OTU_shared_final.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

# 整合丰度数据与注释结果
colnames(daya_abu)[1] = "OTU_id"
phy_ann_abu = merge(phy_ann_data[c("OTU_id", "Class", "Genus")], 
                    daya_abu,
                    by = "OTU_id")


rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL

phy_ann_genus_abu = aggregate(phy_ann_abu[3:294],
                              list(Class = phy_ann_abu$Class,
                                   Genus = phy_ann_abu$Genus),
                              sum)

phy_ann_genus_abu$Class[which(phy_ann_genus_abu$Class == "")] = "Unknown"
phy_ann_genus_abu$Genus[which(phy_ann_genus_abu$Genus == "")] = "Unknown"


rownames(phy_ann_genus_abu) = paste(phy_ann_genus_abu$Class,
                                    phy_ann_genus_abu$Genus,
                                    sep = " ")

phy_ann_genus_abu$Class = NULL
phy_ann_genus_abu$Genus = NULL

phy_ann_genus_abu = as.data.frame(t(phy_ann_genus_abu))

phy_ann_genus_abu = phy_ann_genus_abu/rowSums(phy_ann_genus_abu)

phy_ann_genus_abu = phy_ann_genus_abu[!grepl(pattern = "Unknown", 
                                             colnames(phy_ann_genus_abu))]



phy_ann_genus_abu["Stations"] = gsub("_", "-", rownames(phy_ann_genus_abu))
phy_ann_genus_abu["Stations"] = gsub("-3um-[1,2,3]", "", phy_ann_genus_abu$Stations)
phy_ann_genus_abu["Stations"] = gsub("-0.2um-[1,2,3]", "", phy_ann_genus_abu$Stations)


phy_ann_genus_abu$Sample_id = rownames(phy_ann_genus_abu)
phy_ann_genus_abu["Fraction"] = ""

for (i in c(1:length(phy_ann_genus_abu$Fraction))) {
  if (str_count(rownames(phy_ann_genus_abu)[i], "3um") > 0) {
    phy_ann_genus_abu$Fraction[i] = "3-200 μm"
  }
  if (str_count(rownames(phy_ann_genus_abu)[i], "0.2um") > 0) {
    phy_ann_genus_abu$Fraction[i] = "0.2-3 μm"
  }
}

phy_ann_genus_env = merge(daya_env[c("Stations", "Temperature")], 
                          phy_ann_genus_abu,
                          by = "Stations")


# 按Station和Fraction求均值

phy_ann_genus_env = aggregate(phy_ann_genus_env[2:232],
                              by = list(Stations = phy_ann_genus_env$Stations,
                                        Fraction = phy_ann_genus_env$Fraction),
                              mean)

phy_ann_genus_env_sum = aggregate(phy_ann_genus_env[4:233],
                                  by = list(Fraction = phy_ann_genus_env$Fraction),
                                  sum)

rownames(phy_ann_genus_env_sum) = phy_ann_genus_env_sum$Fraction
phy_ann_genus_env_sum$Fraction = NULL
phy_ann_genus_env_sum = as.data.frame(t(phy_ann_genus_env_sum))
phy_ann_genus_env_sum$ratio = phy_ann_genus_env_sum$`3-200 μm`/phy_ann_genus_env_sum$`0.2-3 μm`

LF_genus = rownames(subset(phy_ann_genus_env_sum, ratio >= 1))
SF_genus = rownames(subset(phy_ann_genus_env_sum, ratio < 1))

LF_genus = append(x = "Stations", LF_genus)
SF_genus = append(x = "Stations", SF_genus)

LF_phy_ann_genus_env = subset(phy_ann_genus_env, Fraction == "3-200 μm")
SF_phy_ann_genus_env = subset(phy_ann_genus_env, Fraction == "0.2-3 μm")

LF_phy_ann_genus_env = LF_phy_ann_genus_env[LF_genus]
SF_phy_ann_genus_env = SF_phy_ann_genus_env[SF_genus]

phy_ann_genus_env = merge(LF_phy_ann_genus_env,
                          SF_phy_ann_genus_env,
                          by = "Stations")


phy_ann_genus_env = merge(daya_env[c("Stations", "Temperature")],
                          phy_ann_genus_env,
                          by = "Stations")



# 移位分组
setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")
library(igraph)
library(NetSwan)
library(psych)
library(dplyr)



genus_env = phy_ann_genus_env
rownames(genus_env) = genus_env$Stations
genus_env$Stations = NULL

attach(genus_env)
genus_env = genus_env[order(genus_env$Temperature),]
detach(genus_env)


# 移位分组
i = 0
len = length(genus_env$Temperature) - 3
df = as.data.frame(matrix(nrow = len-3, ncol = 4))
colnames(df) = c("Temperature", "Edge number")
df$Temperature = genus_env$Temperature[4:len]


for (i in c(4:len)) {
  j = i - 3
  genus_env_part = genus_env[c(i-3,i-2,i-1,i, i+1, i+2,i+3),]
  
  for (k in colnames(genus_env_part)) {
    n = nrow(genus_env_part[genus_env_part[k]>0,])
    if (n < 3) {
      genus_env_part[k] = NULL
    }
  }
  genus_env_part$Temperature = NULL
  genus_env_part = scale(genus_env_part)
  corr = corr.test(genus_env_part, method = "pearson", 
                   ci = FALSE, adjust = "fdr")
  relation = corr$r
  pvalue = corr$p
  relation[abs(relation)<0.6|pvalue>0.05] = 0
  dr = as.data.frame(relation)
  dr$species1 = rownames(dr)
  edge = melt(dr, id.var = "species1")
  colnames(edge)[2:3] = c("species2", "correction")
  edge$color = ""
  edge$color = ifelse(edge$correction>0, "positive", "negative")
  edge$correction = abs(edge$correction)
  edge = edge[edge$species1 != edge$species2,]
  edge = subset(edge, abs(correction) > 0)
  colnames(edge) = c("Source", "Target", "weight", "color")
  edge$unique = ""
  for (i in c(1:length(edge$Source))) {
    if (as.character(edge$Source[i]) >= as.character(edge$Target[i])) edge$unique[i] = paste(edge$Source[i], edge$Target[i])
    if (as.character(edge$Source[i]) < as.character(edge$Target[i])) edge$unique[i] = paste(edge$Target[i], edge$Source[i])
  }
  edge = edge %>% distinct(unique, .keep_all =TRUE)
  edge$unique = NULL
  file_name = paste("edge_inf_pearson/", df$Temperature[j], ".csv", sep = "")
  write.csv(edge, file_name, row.names=FALSE, quote=FALSE)
  edge_number = length(edge$Source)
  df$`Edge number`[j] = edge_number
}


setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S/edge_inf_pearson")
edge_files = list.files()

daya_ann = read.csv("../OTU_taxonomy.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)
colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")
# 选出属于浮游植物的子集
phy_ann_data = subset(daya_ann, Phylum %in% c("Dinoflagellata",
                                              "Ochrophyta",
                                              "Chlorophyta",
                                              "Haptophyta",
                                              "Cryptophyta"))


dino_genus = subset(phy_ann_data, Phylum == "Dinoflagellata")
dino_genus = unique(paste(dino_genus$Class, dino_genus$Genus))
diat_genus = subset(phy_ann_data, Class == "Bacillariophyta")
diat_genus = unique(paste(diat_genus$Class, diat_genus$Genus))
chlo_genus = subset(phy_ann_data, Phylum == "Chlorophyta")
chlo_genus = unique(paste(chlo_genus$Class, chlo_genus$Genus))

dino_genus = as.data.frame(dino_genus)
diat_genus = as.data.frame(diat_genus)
chlo_genus = as.data.frame(chlo_genus)


df = as.data.frame(matrix(nrow = len-2, ncol = 12))

colnames(df) = c("Temperature", "Edge number", "Positive edge number", 
                 "Negative edge number", "Max remove nodes", 
                 "Node number", "Weight mean",
                 "Core node number", "Mean of edge number of nodes",
                 "Core node number in Dinoflagellata",
                 "Core node number in Bacillariophyta",
                 "Core node number in Chlorophyta")


k = 0

edge_files = edge_files[grep(".csv", edge_files)]

for (i in edge_files) {
  k = k + 1
  Temperature = gsub(".csv", "", i)
  temp_file = read.csv(i)
  edge_num = length(temp_file$color)
  pos_num = length(subset(temp_file, color == "positive")$color)
  neg_num = length(subset(temp_file, color == "negative")$color)
  nodes = unique(temp_file$Source)
  nodes = append(nodes, unique(temp_file$Target))
  nodes = unique(nodes)
  edge_num_set = c()
  core_node_num = 0
  dino_cores = 0
  diat_cores = 0
  chlo_cores = 0
  gra = graph.edgelist(as.matrix(temp_file[,c("Source", "Target")]), 
                       directed = FALSE)
  f4<-swan_combinatory(gra,10)
  f4 = as.data.frame(f4)
  colnames(f4) = c("remove", "betweeness", "degree", "cascading", "random")
  f5 = subset(f4, random<0.3)
  max_remove = max(f5$remove)*length(nodes)
  for (j in nodes) {
    nodes_edge_num = length(subset(temp_file, Source == j)$color)
    nodes_edge_num = nodes_edge_num + length(subset(temp_file, Target == j)$color)
    edge_num_set = append(edge_num_set, nodes_edge_num)
    if (nodes_edge_num > 10) {
      core_node_num = core_node_num + 1
      genu = gsub("_LF", "", j)
      genu = gsub("_SF", "", genu)
      if (genu %in% dino_genus$dino_genus) dino_cores = dino_cores + 1
      if (genu %in% diat_genus$diat_genus) diat_cores = diat_cores + 1
      if (genu %in% chlo_genus$chlo_genus) chlo_cores = chlo_cores + 1
    }
  }
  nodes_edge_num_mean = mean(edge_num_set)
  df$Temperature[k] = Temperature
  df$`Edge number`[k] = edge_num
  df$`Positive edge number`[k] = pos_num
  df$`Negative edge number`[k] = neg_num
  df$`Max remove nodes`[k] = max_remove
  df$`Node number`[k] = length(nodes)
  df$`Weight mean`[k] = mean(temp_file$weight)
  df$`Core node number`[k] = core_node_num
  df$`Mean of edge number of nodes`[k] = nodes_edge_num_mean
  df$`Core node number in Dinoflagellata`[k] = dino_cores
  df$`Core node number in Bacillariophyta`[k] = diat_cores
  df$`Core node number in Chlorophyta`[k] = chlo_cores
}


library(ggplot2)
library(reshape2)

env_data = read.csv("../env_data_final.csv", header = T)

df = merge(df, env_data[, c("Temperature", "TChla")], by = "Temperature")

colnames(df)

df$`Positive edge number` = NULL
df$`Negative edge number` = NULL
df$`Core node number in Dinoflagellata` = NULL
df$`Core node number in Bacillariophyta` = NULL
df$`Core node number in Chlorophyta` = NULL



df_long = melt(df, id.vars = c("Temperature", "TChla"))

df_long$Temperature = as.numeric(df_long$Temperature)


df_long$Region = ''
for (i in c(1:length(df_long$Temperature))) {
  if (df_long$Temperature[i] <= 26.6) df_long$Region[i] = "Low Temperature"
  if (df_long$Temperature[i] > 26.6) df_long$Region[i] = "High Temperature"
}


p = ggplot(df_long, aes(x = Temperature, y = value)) +
  labs(x = "Temperature [°C]",
       y = "Value") +
  theme_bw() +
#  geom_vline(xintercept = 26.6, linewidth = 2) +
  geom_point(alpha = 0.5) +
  labs(y = "") +
  facet_wrap(~variable, scales = "free_y") +
  geom_smooth(method = "glm", linewidth = 2, se = F, formula = y ~ x+I(x^2)) +
  theme(legend.position = "none",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 15))
# x11()
p

ggsave("../network_temperatrue_pearson.tiff", dpi = 300)


# 磷酸盐与色素含量
library(ggplot2)
library(lattice)
library(reshape2)
setwd("D:/论文/大亚湾/18S扩增子/色素")
daya_env = read.csv("../Daya Bay 18S/env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
pig_env = read.csv("pigment.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
pig_env = merge(daya_env[c("Stations", "PO43-", "NO3-", "SiO32-")],
                pig_env, by = "Stations")
colnames(pig_env)
pig_env$Station = NULL
pig_env$Longitude = NULL
pig_env$Latitude = NULL
pig_env$Stations = NULL

pig_cor = data.frame(cor(pig_env, method = "pearson"))
colnames(pig_cor)
pig_cor = pig_cor[c("PO43.", "NO3.", "SiO32.")]

pig_cor$pigments = rownames(pig_cor)
pig_cor_long = melt(pig_cor, id.vars = "pigments")

pig_cor_long$pigments = factor(pig_cor_long$pigments,
                               levels = pig_cor$pigments)

# X11()
p = ggplot(pig_cor_long, aes(x = variable, y = pigments)) +
  xlab("") + 
  ylab("") +
  scale_fill_gradient2(high="red", low = "blue", mid = "white") +
  theme_bw() +
#  theme(panel.grid.major = element_blank()) +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10)) +
  geom_tile(aes(fill = value))
p
ggsave(p, filename = "heatmap.tiff", width = 6, 
       height = 5, dpi = 300)

## CCA分析

setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")

daya_env = read.csv("env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann = read.csv("OTU_taxonomy.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)
colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")
# 选出属于浮游植物的子集
phy_ann_data = subset(daya_ann, Phylum %in% c("Dinoflagellata",
                                              "Ochrophyta",
                                              "Chlorophyta",
                                              "Haptophyta",
                                              "Cryptophyta"))

phy_ann_data$OTU_id = rownames(phy_ann_data)

daya_abu = read.csv("OTU_shared_final.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
colnames(daya_abu)[1] = "OTU_id"

# 整合丰度数据与注释结果
phy_ann_abu = merge(phy_ann_data[c("OTU_id", "Class", "Genus")], 
                    daya_abu,
                    by = "OTU_id")
rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL

phy_ann_genus_abu = aggregate(phy_ann_abu[3:294],
                              list(Class = phy_ann_abu$Class,
                                   Genus = phy_ann_abu$Genus),
                              sum)

phy_ann_genus_abu$Class[which(phy_ann_genus_abu$Class == "")] = "Unknown"
phy_ann_genus_abu$Genus[which(phy_ann_genus_abu$Genus == "")] = "Unknown"


rownames(phy_ann_genus_abu) = paste(phy_ann_genus_abu$Class,
                                    phy_ann_genus_abu$Genus,
                                    sep = " ")

phy_ann_genus_abu$Class = NULL
phy_ann_genus_abu$Genus = NULL

phy_ann_genus_abu = as.data.frame(t(phy_ann_genus_abu))

phy_ann_genus_abu = phy_ann_genus_abu/rowSums(phy_ann_genus_abu)

phy_ann_genus_abu = phy_ann_genus_abu[!grepl(pattern = "Unknown", 
                                             colnames(phy_ann_genus_abu))]



phy_ann_genus_abu["Stations"] = gsub("_", "-", rownames(phy_ann_genus_abu))
phy_ann_genus_abu["Stations"] = gsub("-3um-[1,2,3]", "", phy_ann_genus_abu$Stations)
phy_ann_genus_abu["Stations"] = gsub("-0.2um-[1,2,3]", "", phy_ann_genus_abu$Stations)


phy_ann_genus_abu$Sample_id = rownames(phy_ann_genus_abu)
phy_ann_genus_abu["Fraction"] = ""

for (i in c(1:length(phy_ann_genus_abu$Fraction))) {
  if (str_count(rownames(phy_ann_genus_abu)[i], "3um") > 0) {
    phy_ann_genus_abu$Fraction[i] = "3-200 μm"
  }
  if (str_count(rownames(phy_ann_genus_abu)[i], "0.2um") > 0) {
    phy_ann_genus_abu$Fraction[i] = "0.2-3 μm"
  }
}

phy_ann_genus_env = merge(daya_env[c("Stations", "Temperature")], 
                          phy_ann_genus_abu,
                          by = "Stations")


# 按Station和Fraction求均值

phy_ann_genus_env = aggregate(phy_ann_genus_env[2:232],
                              by = list(Stations = phy_ann_genus_env$Stations,
                                        Fraction = phy_ann_genus_env$Fraction),
                              mean)

phy_ann_genus_env_sum = aggregate(phy_ann_genus_env[4:233],
                                  by = list(Fraction = phy_ann_genus_env$Fraction),
                                  sum)

rownames(phy_ann_genus_env_sum) = phy_ann_genus_env_sum$Fraction
phy_ann_genus_env_sum$Fraction = NULL
phy_ann_genus_env_sum = as.data.frame(t(phy_ann_genus_env_sum))
phy_ann_genus_env_sum$ratio = phy_ann_genus_env_sum$`3-200 μm`/phy_ann_genus_env_sum$`0.2-3 μm`

LF_genus = rownames(subset(phy_ann_genus_env_sum, ratio >= 1))
SF_genus = rownames(subset(phy_ann_genus_env_sum, ratio < 1))

LF_genus = append(x = "Stations", LF_genus)
SF_genus = append(x = "Stations", SF_genus)

LF_phy_ann_genus_env = subset(phy_ann_genus_env, Fraction == "3-200 μm")
SF_phy_ann_genus_env = subset(phy_ann_genus_env, Fraction == "0.2-3 μm")

LF_phy_ann_genus_env = LF_phy_ann_genus_env[LF_genus]
SF_phy_ann_genus_env = SF_phy_ann_genus_env[SF_genus]

phy_ann_genus = merge(LF_phy_ann_genus_env,
                      SF_phy_ann_genus_env,
                      by = "Stations")

env_data = daya_env[c("Temperature", "SiO32-", "PO43-", "NO2-",
                      "NH4+", "NO3-", "Stations")]

env_data = subset(env_data, env_data$Stations != "M27-1")

rownames(phy_ann_genus) = phy_ann_genus$Stations
phy_ann_genus$Stations = NULL

group_data = env_data[c("Stations","Temperature")]

group_data["group"] = group_data["Stations"]

for (i in c(1:length(group_data$Stations))){
  if (group_data$Temperature[i]<25) {group_data$region[i] = "(,25)"}
  if (group_data$Temperature[i]>=28) {group_data$region[i] = "[28,)"}
  if (group_data$Temperature[i]>=25 & group_data$Temperature[i]<28) {group_data$region[i] = "[25,28)"}
}
group_data = group_data["region"]

rownames(env_data) = env_data$Stations
env_data$Stations = NULL

rownames(phy_ann_genus)
rownames(env_data)

setdiff(rownames(env_data), rownames(phy_ann_genus))




phy_ann_genus_2 = phy_ann_genus
env_data_2 = env_data

phy_ann_genus_2 = decostand(phy_ann_genus_2, method = "hellinger")


for (col in names(env_data_2)) {
  env_data_2[[col]] = (env_data_2[[col]] - min(env_data_2[[col]])) / (max(env_data_2[[col]]) - min(env_data_2[[col]]))
  env_data_2[[col]] = env_data_2[[col]] + 1
}

  
rda_result = rda(phy_ann_genus_2~., env_data_2)

rda_data = summary(rda_result)

sp = as.data.frame(rda_data$species[,1:2])*1
st = as.data.frame(rda_data$sites[,1:2])*1
yz = as.data.frame(rda_data$biplot[,1:2])*1

color = c("blue", "green", "red")

p = ggplot() +
  geom_point(data = st, aes(RDA1, RDA2, 
                            col = group_data$region),
             size = 3) +
  geom_point(data = sp, aes(RDA1*3, RDA2*3)) +
  scale_color_manual(values = color) +
  geom_segment(data = yz, aes(x = 0, y = 0, 
                              xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.35, "cm"),
                             type = "closed"),
               linetype = 1,
               linewidth = 1,
               color = "black") +
  geom_text_repel(data = yz,
                  aes(RDA1, RDA2, label = row.names(yz))) +
  labs(x = paste("RDA 1(", format(100 *rda_data$cont[[1]][2,1], digits = 4), "%)", sep = ""),
       y = paste("RDA 2(", format(100 *rda_data$cont[[1]][2,2], digits = 4), "%)", sep = ""))+
  theme_bw()
p

ggsave("RDA.pdf")


# 排序

#描述统计
data<-summary(rda_result)
#检验环境因子相关显著性（Monte Carlo permutation test

df_permutest <- permutest(rda_result, permu=999) # permu=999是表示置换循环的次数
#每个环境因子显著性检验
df_envfit <- envfit(rda_result, env_data_2, permu=999)
#数据处理
cor_data<-data.frame(data$constr.chi/data$tot.chi, data$unconst.chi/data$tot.chi)

cor_com <- data.frame(tax=colnames(env_data_2),r=df_envfit$vectors$r,p=df_envfit$vectors$pvals)

cor_com[1:6,3] = cor_com[,3]>0.05 # 将p<0.05标记为FALSE，p>0.05标记为TRUE，使用此数据绘制柱形图。

p3 <- ggplot(cor_com,aes(x =tax, y = r),size=2) +
  geom_bar(aes(fill=tax),stat = 'identity', width = 0.8)+
  geom_text(aes(y = r+0.05, label = ifelse(p==T,"","*")),size = 5, fontface = "bold") +
  labs(x = '', y = '')+
  xlab("Environmental factor")+
  ylab(expression(r^"2"))+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_prism(palette = "candy_bright") +
  theme(legend.position = "none")

p3
ggsave("RDA_env_r2.pdf")

# LASSO模型分析

setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S")

library(stringr)
library(vegan)
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(lemon)
library(igraph)
library(ggrepel)
library(ggprism)

library(glmnet) # LASSO


# 导入大亚湾的数据
daya_env = read.csv("env_data_final.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
daya_pig = read.csv("pigment.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")


env_data = daya_env[c("Stations", "Temperature", "NO3-", "NO2-", 
                      "PO43-", "SiO32-", "NH4+")]

pig_data = daya_pig
pig_data$Longitude = NULL
pig_data$Latitude = NULL


rownames(env_data) = env_data$Stations
env_data$Stations = NULL

length(pig_data)
pig_data = aggregate(pig_data[2:20],
                     by = list(Stations = pig_data$Stations),
                     mean)

rownames(pig_data) = pig_data$Stations
pig_data$Stations = NULL

setdiff(rownames(env_data), rownames(pig_data))
setdiff(rownames(pig_data), rownames(env_data))

env_data = subset(env_data, rownames(env_data) != "M27-1")
env_data$Tabs = abs(env_data$Temperature - 24.8805)

# 先画19种色素相对温度的变化图
pig_env = merge(env_data,
                pig_data,
                by = "row.names")

pig_env$Row.names = NULL
pig_env_long = melt(pig_env, id.vars = "Temperature")
p = ggplot(data = pig_env_long,
           aes(x = Temperature,
               y = value)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~variable, scale = "free_y") +
  geom_smooth(method = "lm", formula = y~x+I(x^2))
  
p

colnames(pig_env)

data = pig_env[1:8]
model = lm(`Chlorophyll c3`~`Tabs`+ `NO3-` + `NO2-` + `PO43-`+ `SiO32-` + `NH4+`,
           data = data)
a = summary(model)

a$fstatistic[1]

model = glm(pigment ~`Temperature`+ `NO3-` + `NO2-` + `PO43-`+ `SiO32-` + `NH4+`,
           data = data)
a = summary(model)


model = lm(`Chlorophyll c3`~`Temperature`+ `NO3-` + `NO2-` + `SiO32-` + `NH4+`,
           data = data)
a = summary(model)
a$aic

model = glm(`Chlorophyll c3`~`Tabs`+ `NO3-` + `NO2-` + `SiO32-` + `NH4+`,
           data = data)
a = summary(model)
a$fstatistic[1]


df = as.data.frame(matrix(nrow = 19, ncol = 3))
row = 0
for (i in colnames(pig_data)) {
  row = row + 1
  data = merge(env_data, pig_data[c(i)], by = "row.names")
  data$Row.names = NULL
  colnames(data)[8] = "pigment"
  model = lm(pigment ~  Temperature + `NO3-` + `NO2-` + `PO43-`+ `SiO32-` + `NH4+`,
             data = data)
  a = summary(model)
  df[row, 1] = as.numeric(a$fstatistic[1])
  model = lm(pigment ~ `Temperature`,
             data = data)
  a = summary(model)
  df[row, 2] = as.numeric(a$fstatistic[1])
  model = lm(pigment ~ `Temperature` + I(Temperature^2),
             data = data)
  a = summary(model)
  df[row, 3] = as.numeric(a$fstatistic[1])
}

rownames(df) = colnames(pig_data)
colnames(df) = c("T + Inorganic nutrient",
                  "T",
                  "T + T^2")

write.csv(df, "lm_model.csv")


# df = as.data.frame(matrix(nrow = 19, ncol = 3))
# row = 0
# for (i in colnames(pig_data)) {
#   row = row + 1
#   data = merge(env_data, pig_data[c(i)], by = "row.names")
#   data$Row.names = NULL
#   colnames(data)[8] = "pigment"
#   model = glm(pigment ~  Temperature + `NO3-` + `NO2-` + `PO43-`+ `SiO32-` + `NH4+`,
#              data = data)
#   a = summary(model)
#   df[row, 1] = as.numeric(a$aic)
#   model = glm(pigment ~ `Temperature`,
#              data = data)
#   a = summary(model)
#   df[row, 2] = as.numeric(a$aic)
#   model = glm(pigment ~ `Temperature` + I(Temperature^2),
#              data = data)
#   a = summary(model)
#   df[row, 3] = as.numeric(a$aic)
# }
# 
# df


df$pigment = rownames(df)
df_long = melt(df, id.vars = "pigment")
p = ggplot(df_long, aes(x = variable, 
                        y = value,
                        fill = variable)) +
  scale_fill_brewer(palette = "Dark2") +
  geom_col() +
  theme_bw() +
  labs(x = "Combination of variables",
       y = "F-statistic") +
  facet_wrap(~pigment, scale = "free_y") +
  theme(legend.position = c(0.9, 0.1),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 15))
p

ggsave("lm_model.pdf")

# 功能类群
phy_ann_abu[1:10, 1:10]
length(phy_ann_abu)
fun_abu = aggregate(phy_ann_abu[4:295],
                    by = list(phy_ann_abu$Class),
                    sum)
rownames(fun_abu) = fun_abu$Group.1

fun_abu = as.data.frame(t(fun_abu))

fun_abu$Stations = gsub("_", "-", rownames(fun_abu))
fun_abu$Stations = gsub("-3um-[1,2,3]", "", fun_abu$Stations)
fun_abu$Stations = gsub("-0.2um-[1,2,3]", "",fun_abu$Stations)



fun_abu$Fraction = rownames(fun_abu)
fun_abu$Fraction = ifelse(grepl("0.2um", fun_abu$Fraction), 
                          "0.2-3 μm", "3-200 μm")

colnames(daya_env)
fun_abu_env = merge(daya_env[c("Stations", "Temperature")], 
                    fun_abu,
                    by = "Stations")


lf_fun_abu_env = subset(fun_abu_env, Fraction == "3-200 μm")
Sf_fun_abu_env = subset(fun_abu_env, Fraction == "0.2-3 μm")

lf_fun_abu_env[3:42] <- apply(lf_fun_abu_env[3:42], 2,
                              as.numeric)
lf_fun_abu_env[3:42] = lf_fun_abu_env[3:42]/rowSums(lf_fun_abu_env[3:42])

colSums(lf_fun_abu_env[3:42])

data = lf_fun_abu_env[3:42]

lf_fun_abu_env[3:42] = data[, order(-colSums(data))]

colnames(lf_fun_abu_env)[3:42] = colnames(data[, order(-colSums(data))])

colSums(lf_fun_abu_env[3:42])

lf_fun_abu_env = lf_fun_abu_env[1:8]
lf_fun_abu_env_long = melt(lf_fun_abu_env, id.vars = c("Stations","Temperature"))



Sf_fun_abu_env[3:42] <- apply(Sf_fun_abu_env[3:42], 2, as.numeric)
Sf_fun_abu_env[3:42] = Sf_fun_abu_env[3:42]/rowSums(Sf_fun_abu_env[3:42])


data = Sf_fun_abu_env[3:42]

Sf_fun_abu_env[3:42] = data[, order(-colSums(data))]

colnames(Sf_fun_abu_env)[3:42] = colnames(data[, order(-colSums(data))])

colSums(Sf_fun_abu_env[3:42])

Sf_fun_abu_env = Sf_fun_abu_env[1:8]
Sf_fun_abu_env_long = melt(Sf_fun_abu_env, id.vars = c("Stations","Temperature"))

colnames(lf_fun_abu_env_long)
lf_fun_abu_env_long = aggregate(lf_fun_abu_env_long$value,
                                by = list(Stations = lf_fun_abu_env_long$Stations,
                                          Class = lf_fun_abu_env_long$variable,
                                          Temperature = lf_fun_abu_env_long$Temperature),
                                mean)
colnames(lf_fun_abu_env_long)
Sf_fun_abu_env_long = aggregate(Sf_fun_abu_env_long$value,
                                by = list(Stations = Sf_fun_abu_env_long$Stations,
                                          Class = Sf_fun_abu_env_long$variable,
                                          Temperature = Sf_fun_abu_env_long$Temperature),
                                mean)

colnames(Sf_fun_abu_env_long)


lf_fun_abu_env_long$Region = ''
for (i in c(1:length(lf_fun_abu_env_long$Temperature))) {
  if (lf_fun_abu_env_long$Temperature[i] <= 25) lf_fun_abu_env_long$Region[i] = "Low Temperature"
  if (lf_fun_abu_env_long$Temperature[i] > 28) lf_fun_abu_env_long$Region[i] = "High Temperature"
  if (lf_fun_abu_env_long$Temperature[i] > 25&lf_fun_abu_env_long$Temperature[i] <= 28) {
    lf_fun_abu_env_long$Region[i] = "Middle Temperature"
  }
}

lf_fun_abu_env_long$Region = factor(lf_fun_abu_env_long$Region,
                                    levels = c("Low Temperature",
                                               "Middle Temperature",
                                               "High Temperature"))


p1 = ggplot(lf_fun_abu_env_long,
            aes(x = Temperature,
                y = x*100,
                col = Region)) +
  geom_point(alpha = 0.5) +
#  geom_boxplot(aes(group = Region, y = x*100, col = Region)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Temperature [°C]",
       y = "Relative abundance (%)",
       title = "3-200 μm Fraction") +
  facet_wrap(~Class, scale = "free_y") +
  geom_smooth(method = "glm", se = F, linewidth = 2, col = "red", formula = y ~ x+I(x^2)) +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

#p1

Sf_fun_abu_env_long$Region = ''
for (i in c(1:length(Sf_fun_abu_env_long$Temperature))) {
  if (Sf_fun_abu_env_long$Temperature[i] <= 25) Sf_fun_abu_env_long$Region[i] = "Low Temperature"
  if (Sf_fun_abu_env_long$Temperature[i] > 28) Sf_fun_abu_env_long$Region[i] = "High Temperature"
  if (Sf_fun_abu_env_long$Temperature[i] > 25&Sf_fun_abu_env_long$Temperature[i] <= 28) {
    Sf_fun_abu_env_long$Region[i] = "Middle Temperature"
  }
}

p2 = ggplot(Sf_fun_abu_env_long,
           aes(x = Temperature,
               y = x*100,
               col = Region)) +
  geom_point(alpha = 0.5) +
  #  geom_boxplot(aes(group = Region, y = x*100, col = Region)) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Temperature [°C]",
       y = "Relative abundance (%)",
       title = "0.2-3 μm Fraction") +
  facet_wrap(~Class, scale = "free_y") +
  geom_smooth(method = "glm", se = F, linewidth = 2, col = "Blue", formula = y ~ x+I(x^2)) +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p = p1 + p2 
p

ggsave("Daya_class.tiff", dpi = 300)


# Tara

tara_env = read.csv("tara_env.txt", sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

tara_phy = read.csv("phy_data.csv",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")

tara_class = aggregate(tara_phy[10:343],
                       by = list(tara_phy$Class),
                       sum)
row.names(tara_class) = tara_class$Group.1

tara_class$Group.1 = NULL
tara_class = as.data.frame(t(tara_class))
tara_class$Barcode = row.names(tara_class)

tara_class_env = merge(tara_env[c("Barcode", "Temperature",
                                  "Depth_Nominal", "Fraction_Lower",
                                  "Fraction_Upper")],
                       tara_class,
                       by = "Barcode")


tara_class_env$Fraction = paste(tara_class_env$Fraction_Lower,
                                tara_class_env$Fraction_Upper,
                                sep = "-")
tara_class_env$Fraction = paste(tara_class_env$Fraction, "μm")

tara_class_env_2 = subset(tara_class_env, Fraction %in% c("0.8-5.00 μm",
                                                        "5-20.00 μm",
                                                        "20-180.00 μm"))

tara_class_env_2 = subset(tara_class_env_2, Depth_Nominal == "SRF")
colnames(tara_class_env_2)
tara_class_env_2$Fraction_Lower = NULL
tara_class_env_2$Fraction_Upper = NULL
tara_class_env_2$Depth_Nominal = NULL


tara_class_env_2 = subset(tara_class_env_2, !Barcode %in% c("TARA_X000001286",
                                                           "TARA_X000001288",
                                                           "TARA_X000001338"))


len = length(tara_class_env_2) - 1
tara_class_env_2[3:len] = apply(tara_class_env_2[3:len], 
                                  2, 
                                  as.numeric)
tara_class_env_2[3:len] = tara_class_env_2[3:len]/rowSums(tara_class_env_2[3:len])


data = tara_class_env_2[3:len]

tara_class_env_2[3:len] = data[, order(-colSums(data))]

colnames(tara_class_env_2)[3:len] = colnames(data[, order(-colSums(data))])

colSums(tara_class_env_2[3:len])

tara_class_env_3 = tara_class_env_2[1:7]

tara_class_env_3 = cbind(tara_class_env_3, tara_class_env_2$Fraction)

colnames(tara_class_env_3)[length(tara_class_env_3)] = "Fraction"
tara_class_env_2$Fraction
tara_class_env_long = melt(tara_class_env_3, id.vars = c("Barcode",
                                                         "Temperature",
                                                         "Fraction"))

tara_class_env_long = subset(tara_class_env_long, Temperature > 15)

# x11()
tara_class_env_long$Fraction = factor(tara_class_env_long$Fraction,
                                      levels = c("0.8-5.00 μm","5-20.00 μm","20-180.00 μm"))


tara_class_env_long$value = as.numeric(tara_class_env_long$value)


p = ggplot(tara_class_env_long, aes(x = Temperature,
                                    y = value,
                                    col = variable)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "") +
  geom_smooth(se = F, size = 1, method = "lm", formula = y~x+I(x^2)) +
#  geom_smooth(se = F, size = 1, method = "loess") +
  facet_wrap(Fraction~variable, scales = "free_y", nrow = 3) +
  theme(legend.position = "None",
        strip.text = element_text(size = 10, family = "Arial", face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
p
ggsave("tara_class.tiff", dpi = 300)




# 绘制网络图
library(igraph)
setwd("D:/论文/大亚湾/18S扩增子/Daya Bay 18S/edge_inf")
edge_data = read.csv("29.csv", sep = ",",  header = T)
node_data = rbind(edge_data$Source, edge_data$Target)
node_data = str_split_fixed(node_data, " ", 2)
node_data = as.data.frame(node_data)
node_data = unique(node_data)

colnames(node_data) = c("class", "genus")
node_data = node_data[c("genus", "class")]
node_data$genus = paste(node_data$class, node_data$genus)

edge_data$color[edge_data$color == "positive"] = "red"
edge_data$color[edge_data$color == "negative"] = "green"

node_data$class[!(node_data$class %in% c("Dinophyceae", "Bacillariophyta", "Mamiellophyceae"))] = "Others"




network = graph_from_data_frame(d=edge_data, vertices=node_data, directed=F) 

library(RColorBrewer)

coul  <- brewer.pal(length(unique(node_data$class)), "Dark2") 




# Create a vector of color
my_color <- coul[as.numeric(as.factor(V(network)$class))]

# Make the plot
plot(network, vertex.color = my_color, layout = layout_with_fr,
     vertex.size = 7, vertex.label = "")
legend("bottomleft", legend=levels(as.factor(V(network)$class)), 
       col = coul , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, 
       text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))


# 高温区是否有低温区没有的浮游植物OTU？


daya_env = read.csv("env_data_final.csv", header = T, 
                    check.names = F, encoding = "UTF-8")

daya_ann = read.csv("OTU_taxonomy.xls", sep = "\t", 
                    header = T, check.names = F, encoding = "UTF-8")

daya_ann[3:9] = str_split_fixed(daya_ann$Taxonomy, ";", 7)
colnames(daya_ann)[3:9] = c("Kingdom", "Phylum",
                            "Class", "Order", "Family", 
                            "Genus", "Species")

phy_ann_data = subset(daya_ann, Phylum %in% c("Dinoflagellata",
                                              "Ochrophyta",
                                              "Chlorophyta",
                                              "Haptophyta",
                                              "Cryptophyta"))

phy_ann_data$OTU_id = rownames(phy_ann_data)

daya_abu = read.csv("OTU_shared_final.xls",
                    sep = "\t",
                    header = T,
                    check.names = F,
                    encoding = "UTF-8")
colnames(daya_abu)[1] = "OTU_id"

# 整合丰度数据与注释结果
phy_ann_abu = merge(phy_ann_data[c("OTU_id")], 
                    daya_abu,
                    by = "OTU_id")


rownames(phy_ann_abu) = phy_ann_abu$OTU_id
phy_ann_abu$OTU_id = NULL
phy_ann_abu = as.data.frame(t(phy_ann_abu))



phy_ann_abu["Stations"] = gsub("_", "-", rownames(phy_ann_abu))
phy_ann_abu["Stations"] = gsub("-3um-[1,2,3]", "", phy_ann_abu$Stations)
phy_ann_abu["Stations"] = gsub("-0.2um-[1,2,3]", "", phy_ann_abu$Stations)

length(phy_ann_abu)
colnames(phy_ann_abu)[1510]


phy_ann_abu_sum = aggregate(phy_ann_abu[1:1509],
                            by = list(Stations = phy_ann_abu$Stations),
                            sum)


rownames(phy_ann_abu_sum) = phy_ann_abu_sum$Stations
phy_ann_abu_sum$Stations = NULL

data = ifelse(phy_ann_abu_sum > 1, 1, 0)
data = as.data.frame(data)
data$Stations = rownames(data)
data = merge(data, daya_env[c("Stations", "Temperature")],
             by = "Stations")
data = data[order(data$Temperature),]



rownames(data) = data$Stations
data$Stations = NULL
data$Temperature = NULL
data = as.data.frame(t(data))

df = as.data.frame(matrix(ncol = 2, nrow = 49))
colnames(df) = c("Stations", "Number of new OTUs")

new_otu_df = subset(data, `M24-2` == 0)

new_otu_df$`M24-2` = NULL

for (i in c(1:length(colnames(new_otu_df)))) {
  Station = colnames(new_otu_df)[i]
  len_1 = length(new_otu_df[new_otu_df[Station] == 1,]$T24)
  
  new_otu_df = new_otu_df[new_otu_df[Station] == 0,]

  df$Stations[i] = Station
  df$`Number of new OTUs`[i] = len_1
}

df = merge(df, daya_env[c("Stations", 
                          "Temperature")],
           by = "Stations")


p1 = ggplot(df,
           aes(x = Temperature,
               y = `Number of new OTUs`)) +
  geom_point() +
  theme_bw() +
  labs(x = "Temperature [°C]")
p1
ggsave("new OTUs and Temperature.tiff", dpi = 300)










