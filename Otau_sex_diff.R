##### load packages######
library(tidyverse)
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library
library(edgeR)
library(stats)
library(topGO)
library(RColorBrewer)

# I. RNAseq data analyses #############################################################################
##### Step 1: import RNAseq read abundance from salmon #####
# read in study design
Ot_sample_table <- read.delim("/Users/ericanadolski/GitHub/Otaurus_sexual_dimorphism/Ot_sample_info_RNA.txt") 
#Ot F11 E removed from Ot_sample_info2.txt
Ot_sample_table <- Ot_sample_table %>% mutate(Sex_Trait = paste0(Sex, "_", Trait)) 

# create file paths to the abundance files using the 'file.path' function
path <- file.path("/Volumes/T7_Drive_2/RNAseq/salmon-dsx", paste0(Ot_sample_table$Sample, "_quant/quant.sf")) 
file.exists(path)

# import counts with tximport
Ot_salmon_tx <- tximport(path, 
                         type = "salmon",
                         txOut = TRUE,
                         countsFromAbundance = "lengthScaledTPM",
                         ignoreTxVersion = TRUE)
class(Ot_salmon_tx)
names(Ot_salmon_tx)

Ot_counts <- as.data.frame(Ot_salmon_tx$counts)
colnames(Ot_counts) <- Ot_sample_table$Sample

##### Step 2: data filtering and visualization #####
# generate DESeq dataset for exploratory analysis
Ot_dds_exp <- DESeqDataSetFromTximport(Ot_salmon_tx, Ot_sample_table, design = ~ 1)

# pre-filtering data
nrow(Ot_dds_exp) # 19500
# filter out genes with very low read count across most samples
Ot_dds_exp <- Ot_dds_exp[rowSums(counts(Ot_dds_exp)) > 5, ]
nrow(Ot_dds_exp) # 15979

Ot_dds_rld <- rlog(Ot_dds_exp)
Ot_dds_rld <- as.data.frame(assay(Ot_dds_rld))
Ot_dds_rld <- as.data.frame(t(Ot_dds_rld))

Ot_sampleDists <- dist(Ot_dds_rld)
Ot_sampleDistMatrix <- as.matrix(Ot_sampleDists)

### heatmap all traits ####
rownames(Ot_sampleDistMatrix) <- paste(Ot_sample_table$Sex,Ot_sample_table$Replicate, Ot_sample_table$Trait,sep="-")
colnames(Ot_sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
Ot_heatmap <- pheatmap(Ot_sampleDistMatrix,
                       clustering_distance_rows=Ot_sampleDists,
                       clustering_distance_cols=Ot_sampleDists,
                       col=colors)

ggsave(Ot_heatmap, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_heatmap_salmon.pdf",
       width = 8, height = 8)

### principal component analysis (PCA) all traits ####
PCA_Ot_RNA <- prcomp(Ot_dds_rld, center = TRUE)
PCs_Ot_RNA <- as.data.frame(PCA_Ot_RNA$x)
PCs_Ot_RNA$Sample <- row.names(PCs_Ot_RNA)
PCs_Ot_RNA$Sample <- Ot_sample_table$Sample
PCs_Ot_RNA <- inner_join(PCs_Ot_RNA, Ot_sample_table, by='Sample') 
summary(PCA_Ot_RNA)

Ot_PCA <- ggplot(data = PCs_Ot_RNA, aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 3) +
  scale_fill_manual(values=c("#F8766D", "#4393C3")) + 
  scale_shape_manual(values = c(17, 3, 15, 18, 16)) +
  geom_text(aes(label = Replicate), nudge_y = 3) +
  labs(title="O. taurus RNA-seq Samples - Salmon",x = "PC1 (22.5%)", y = "PC2 (11.6%)") +
  theme_bw()

ggsave(Ot_PCA, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_PCA_plot_salmon.pdf",
       width = 6, height = 5)


### Fore tibia PCA ####
cols <- c( "F" = "#C1272D", "M" = "#2166AC") 
Ot_L_info <- filter(Ot_sample_table, Trait == "L")
Ot_L_counts <- Ot_counts %>% select(contains("L")) 

dds_Ot_L_dg <- DESeqDataSetFromMatrix(countData = round(Ot_L_counts),
                                      colData = Ot_L_info,
                                      design = ~ 1)
rlog_dds_Ot_L <- rlog(dds_Ot_L_dg) 
rlog_counts_Ot_L <- as.data.frame(assay(rlog_dds_Ot_L))
rlog_countsT_Ot_L <- as.data.frame(t(rlog_counts_Ot_L))

PCA_Ot_L <- prcomp(rlog_countsT_Ot_L, center = TRUE)
PCs_Ot_L <- as.data.frame(PCA_Ot_L$x)
PCs_Ot_L$Sample <- row.names(PCs_Ot_L)
PCs_Ot_L <- inner_join(PCs_Ot_L, Ot_L_info, by='Sample') 
summary(PCA_Ot_L)

PCA_Ot_L <- ggplot(data = PCs_Ot_L, 
                   aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus RNAseq Fore Tibia",x = "PC1 (25.7%)", y = "PC2 (20.0%)") +
  theme_bw()

ggsave(PCA_Ot_L, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/PCA_Ot_L.pdf",
       width = 4.30, height = 3.77)

### Posterior head PCA ####
Ot_PH_info <- filter(Ot_sample_table, Trait == "PH")
Ot_PH_counts <- Ot_counts %>% select(contains("PH")) 

dds_Ot_PH_dg <- DESeqDataSetFromMatrix(countData = round(Ot_PH_counts),
                                       colData = Ot_PH_info,
                                       design = ~ 1)
rlog_dds_Ot_PH <- rlog(dds_Ot_PH_dg) 
rlog_counts_Ot_PH <- as.data.frame(assay(rlog_dds_Ot_PH))
rlog_countsT_Ot_PH <- as.data.frame(t(rlog_counts_Ot_PH))

PCA_Ot_PH <- prcomp(rlog_countsT_Ot_PH, center = TRUE)
PCs_Ot_PH <- as.data.frame(PCA_Ot_PH$x)
PCs_Ot_PH$Sample <- row.names(PCs_Ot_PH)
PCs_Ot_PH <- inner_join(PCs_Ot_PH, Ot_PH_info, by='Sample') 
summary(PCA_Ot_PH)

PCA_Ot_PH <- ggplot(data = PCs_Ot_PH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus RNAseq PH",x = "PC1 (27.9%)", y = "PC2 (12.4%)") +
  theme_bw()

ggsave(PCA_Ot_PH, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/PCA_Ot_PH.pdf",
       width = 4.30, height = 3.77)

### Anterior head PCA ####
Ot_AH_info <- filter(Ot_sample_table, Trait == "AH")
Ot_AH_counts <- Ot_counts %>% select(contains("AH")) 

dds_Ot_AH_dg <- DESeqDataSetFromMatrix(countData = round(Ot_AH_counts),
                                       colData = Ot_AH_info,
                                       design = ~ 1)
rlog_dds_Ot_AH <- rlog(dds_Ot_AH_dg) 
rlog_counts_Ot_AH <- as.data.frame(assay(rlog_dds_Ot_AH))
rlog_countsT_Ot_AH <- as.data.frame(t(rlog_counts_Ot_AH))

PCA_Ot_AH <- prcomp(rlog_countsT_Ot_AH, center = TRUE)
PCs_Ot_AH <- as.data.frame(PCA_Ot_AH$x)
PCs_Ot_AH$Sample <- row.names(PCs_Ot_AH)
PCs_Ot_AH <- inner_join(PCs_Ot_AH, Ot_AH_info, by='Sample') 
summary(PCA_Ot_AH)

PCA_Ot_AH <- ggplot(data = PCs_Ot_AH, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus RNAseq AH",x = "PC1 (31.9%)", y = "PC2 (13.6%)") +
  theme_bw()

ggsave(PCA_Ot_AH, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/PCA_Ot_AH.pdf",
       width = 4.30, height = 3.77)


### Genitalia PCA ####
Ot_G_info <- filter(Ot_sample_table, Trait == "G")
Ot_G_counts <- Ot_counts %>% select(contains("G")) 

dds_Ot_G_dg <- DESeqDataSetFromMatrix(countData = round(Ot_G_counts),
                                      colData = Ot_G_info,
                                      design = ~ 1)
rlog_dds_Ot_G <- rlog(dds_Ot_G_dg) 
rlog_counts_Ot_G <- as.data.frame(assay(rlog_dds_Ot_G))
rlog_countsT_Ot_G <- as.data.frame(t(rlog_counts_Ot_G))

PCA_Ot_G <- prcomp(rlog_countsT_Ot_G, center = TRUE)
PCs_Ot_G <- as.data.frame(PCA_Ot_G$x)
PCs_Ot_G$Sample <- row.names(PCs_Ot_G)
PCs_Ot_G <- inner_join(PCs_Ot_G, Ot_G_info, by='Sample') 
summary(PCA_Ot_G)

PCA_Ot_G <- ggplot(data = PCs_Ot_G, 
                   aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus RNAseq Genitalia",x = "PC1 (42.6%)", y = "PC2 (14.0%)") +
  theme_bw()

ggsave(PCA_Ot_G, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/PCA_Ot_G.pdf",
       width = 4.30, height = 3.77)

### Elytra PCA ####
Ot_E_info <- filter(Ot_sample_table, Trait == "E")
Ot_E_counts <- Ot_counts %>% select(contains("-E")) 

dds_Ot_E_dg <- DESeqDataSetFromMatrix(countData = round(Ot_E_counts),
                                      colData = Ot_E_info,
                                      design = ~ 1)
rlog_dds_Ot_E <- rlog(dds_Ot_E_dg) 
rlog_counts_Ot_E <- as.data.frame(assay(rlog_dds_Ot_E))
rlog_countsT_Ot_E <- as.data.frame(t(rlog_counts_Ot_E))

PCA_Ot_E <- prcomp(rlog_countsT_Ot_E, center = TRUE)
PCs_Ot_E <- as.data.frame(PCA_Ot_E$x)
PCs_Ot_E$Sample <- row.names(PCs_Ot_E)
PCs_Ot_E <- inner_join(PCs_Ot_E, Ot_E_info, by='Sample') 
summary(PCA_Ot_E)

PCA_Ot_E <- ggplot(data = PCs_Ot_E, 
                   aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus RNAseq Elytra",x = "PC1 (34.8%)", y = "PC2 (16.0%)") +
  theme_bw()

ggsave(PCA_Ot_E, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/PCA_Ot_E.pdf",
       width = 4.30, height = 3.77)

##### Step 3: analyze differential gene expression #####
### generate DESeq dataset with grouping variable design ####
Ot_dds <- DESeqDataSetFromTximport(Ot_salmon_tx, Ot_sample_table, design = ~ Sex + Trait)
Ot_dds$group <- factor(paste0(Ot_dds$Sex, Ot_dds$Trait))
design(Ot_dds) <- ~ group

# filter out genes with very low read count across most samples
Ot_dds <- Ot_dds[rowSums(counts(Ot_dds)) > 5, ]

Ot_dds <- DESeq(Ot_dds)
resultsNames(Ot_dds)

# save normalized counts output
Ot_norm_counts <- counts(Ot_dds, normalized=TRUE)
colnames(Ot_norm_counts) <- Ot_sample_table$Sample

### perform LFC shrink on results ####
Ot_dds_shrink <- DESeq(Ot_dds, betaPrior=TRUE)
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","MPH","FPH"))) %>% filter(padj <= 0.1)) 
# 1730 
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","MAH","FAH"))) %>% filter(padj <= 0.1)) 
# 1954
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","MG","FG"))) %>% filter(padj <= 0.1)) 
# 1740
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","ML","FL"))) %>% filter(padj <= 0.1)) 
# 790 
tally(as.data.frame(results(Ot_dds_shrink, contrast=c("group","ME","FE"))) %>% filter(padj <= 0.1)) 
# 1204

### save results of each comparison as a dataframe #####
# posterior head
Ot_PH_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","MPH","FPH")))
Ot_PH_MvF_gene <- rownames_to_column(Ot_PH_MvF_gene)
names(Ot_PH_MvF_gene)[1] <- "gene"
Ot_PH_MvF_gene <- Ot_PH_MvF_gene %>% 
  mutate(PH_upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                           ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                                  "ns")))

Ot_PH_MvF_deg <- Ot_PH_MvF_gene %>% filter(padj <= 0.1)
nrow(Ot_PH_MvF_deg)
# anterior head
Ot_AH_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","MAH","FAH")))
Ot_AH_MvF_gene <- rownames_to_column(Ot_AH_MvF_gene)
names(Ot_AH_MvF_gene)[1] <- "gene"
Ot_AH_MvF_gene <- Ot_AH_MvF_gene %>% 
  mutate(AH_upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                           ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                                  "ns")))

Ot_AH_MvF_deg <- Ot_AH_MvF_gene %>% filter(padj <= 0.1)

# elytra
Ot_E_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","ME","FE")))
Ot_E_MvF_gene <- rownames_to_column(Ot_E_MvF_gene)
names(Ot_E_MvF_gene)[1] <- "gene"
Ot_E_MvF_gene <- Ot_E_MvF_gene %>% 
  mutate(E_upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                                 "ns")))

Ot_E_MvF_deg <- Ot_E_MvF_gene %>% filter(padj <= 0.1)

# fore tibia
Ot_L_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","ML","FL")))
Ot_L_MvF_gene <- rownames_to_column(Ot_L_MvF_gene)
names(Ot_L_MvF_gene)[1] <- "gene"
Ot_L_MvF_gene <- Ot_L_MvF_gene %>% 
  mutate(L_upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                                 "ns")))

Ot_L_MvF_deg <- Ot_L_MvF_gene %>% filter(padj <= 0.1)

# genitalia
Ot_G_MvF_gene <- as.data.frame(results(Ot_dds_shrink, contrast=c("group","MG","FG")))
Ot_G_MvF_gene <- rownames_to_column(Ot_G_MvF_gene)
names(Ot_G_MvF_gene)[1] <- "gene"
Ot_G_MvF_gene <- Ot_G_MvF_gene %>% 
  mutate(G_upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                                 "ns")))

Ot_G_MvF_deg <- Ot_G_MvF_gene %>% filter(padj <= 0.1)

##### Step 4: plotting differential gene expression results #####
### gene expression volcano plots ######
### posterior head 
# annotate DEG groups 
Ot_PH_MvF_gene <- Ot_PH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))

Ot_PH_MvF_gene %>%
  dplyr::count(upreg)

cols <- c("F" = "#C1272D", "M" = "#2166AC", "ns" = "grey") 
sizes <- c("F" = 2, "M" = 2, "ns" = 1) 
alphas <- c("F" = 0.7, "M" = 0.7, "ns" = 1)

Ot_PH_vol_plot <- Ot_PH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = upreg,    
             size = upreg,
             alpha = upreg)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) + # 4 male, 1 female outliers log padj > 15
  xlim(-5, 5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Posterior Head")
Ot_PH_vol_plot

ggsave(Ot_PH_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_PH_vol_plot.pdf",
       width = 4, height = 4)

### anterior head 
Ot_AH_MvF_gene <- Ot_AH_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Ot_AH_MvF_gene %>%
  dplyr::count(upreg)

Ot_AH_vol_plot <- Ot_AH_MvF_gene %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = upreg,    
             size = upreg,
             alpha = upreg)) + 
  geom_point(shape = 21, 
             colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) + 
  ylim(0,12) +
  xlim(-5, 5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Anterior Head")
Ot_AH_vol_plot

ggsave(Ot_AH_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_AH_vol_plot.pdf",
       width = 4, height = 4)

### genitalia
Ot_G_MvF_gene <- Ot_G_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))

Ot_G_MvF_gene %>%
  dplyr::count(upreg)

Ot_G_vol_plot <- Ot_G_MvF_gene %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = upreg,    
             size = upreg,
             alpha = upreg)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) +
  xlim(-5, 5) + # two female outliers lfc > 5, one female outlier >12 padj
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Genitalia")
Ot_G_vol_plot

ggsave(Ot_G_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_G_vol_plot.pdf",
       width = 4, height = 4)

### fore tibia
Ot_L_MvF_gene <- Ot_L_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))
Ot_L_MvF_gene %>%
  dplyr::count(upreg)

Ot_L_vol_plot <- Ot_L_MvF_gene %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = upreg,    
             size = upreg,
             alpha = upreg)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) + # one outlier female lfc > 12
  xlim(-5, 5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Fore Tibia")
Ot_L_vol_plot

ggsave(Ot_L_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_L_vol_plot.pdf",
       width = 4, height = 4)

### elytra 
Ot_E_MvF_gene <- Ot_E_MvF_gene %>% 
  mutate(upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                        ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                               "ns")))

Ot_E_MvF_gene %>%
  dplyr::count(upreg)

Ot_E_vol_plot <- Ot_E_MvF_gene %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = upreg,    
             size = upreg,
             alpha = upreg)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) +
  xlim(-5, 5) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(title="Elytra")
Ot_E_vol_plot

ggsave(Ot_E_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_E_vol_plot.pdf",
       width = 4, height = 4)

### clustered heatmaps of DEGs #####
# posterior head 
# subset count matrix by DEGs and trait for heatmap 
Ot_PH_norm_counts <- select(as.data.frame(Ot_norm_counts), contains("PH"))
Ot_PH_norm_counts <- rownames_to_column(Ot_PH_norm_counts)
names(Ot_PH_norm_counts)[1] <- "gene"
Ot_PH_deg <- inner_join(Ot_PH_MvF_deg, Ot_PH_norm_counts, by="gene")

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(Ot_PH_deg[,9:20]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_PH_deg[,9:20], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_PH_deg_heatmap <- heatmap.2(as.matrix(Ot_PH_deg[,9:20]), 
                               Rowv=as.dendrogram(clustRows), 
                               Colv=as.dendrogram(clustColumns),
                               RowSideColors=module.color,
                               col=myheatcolors, scale='row', labRow=NA,
                               density.info="none", trace="none",)

# anterior head
Ot_AH_norm_counts <- select(as.data.frame(Ot_norm_counts), contains("AH"))
Ot_AH_norm_counts <- rownames_to_column(Ot_AH_norm_counts)
names(Ot_AH_norm_counts)[1] <- "gene"
Ot_AH_deg <- inner_join(Ot_AH_MvF_deg, Ot_AH_norm_counts, by="gene")

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(Ot_AH_deg[,9:20]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_AH_deg[,9:20], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_AH_deg_heatmap <- heatmap.2(as.matrix(Ot_AH_deg[,9:20]), 
                               Rowv=as.dendrogram(clustRows), 
                               Colv=as.dendrogram(clustColumns),
                               RowSideColors=module.color,
                               col=myheatcolors, scale='row', labRow=NA,
                               density.info="none", trace="none",)

# genitalia 
Ot_G_norm_counts <- select(as.data.frame(Ot_norm_counts), contains("G"))
Ot_G_norm_counts <- rownames_to_column(Ot_G_norm_counts)
names(Ot_G_norm_counts)[1] <- "gene"
Ot_G_deg <- inner_join(Ot_G_MvF_deg, Ot_G_norm_counts, by="gene")

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(Ot_G_deg[,9:20]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_G_deg[,9:20], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_G_deg_heatmap <- heatmap.2(as.matrix(Ot_G_deg[,9:20]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

Ot_G_deg_heatmap
# fore tibia
Ot_L_norm_counts <- select(as.data.frame(Ot_norm_counts), contains("L"))
Ot_L_norm_counts <- rownames_to_column(Ot_L_norm_counts)
names(Ot_L_norm_counts)[1] <- "gene"
Ot_L_deg <- inner_join(Ot_L_MvF_deg, Ot_L_norm_counts, by="gene")

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(Ot_L_deg[,9:20]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_L_deg[,9:20], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_L_deg_heatmap <- heatmap.2(as.matrix(Ot_L_deg[,9:20]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

# elytra
Ot_E_norm_counts <- select(as.data.frame(Ot_norm_counts), contains("-E"))
Ot_E_norm_counts <- rownames_to_column(Ot_E_norm_counts)
names(Ot_E_norm_counts)[1] <- "gene"
Ot_E_deg <- inner_join(Ot_E_MvF_deg, Ot_E_norm_counts, by="gene")

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(Ot_E_deg[,9:19]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_E_deg[,9:19], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_E_deg_heatmap <- heatmap.2(as.matrix(Ot_E_deg[,9:19]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

### bar chart - total number of DEGs ####
sex <- c('F', 'M','F', 'M','F', 'M','F', 'M','F', 'M') 
trait <- c('PH', 'PH','G','G','L','L','AH','AH','E','E')
upreg.genes <- c(1014,716,1185,555,435,355,1088,866,696,508)

Ot_deg_numbers <- data.frame(sex,trait,upreg.genes)
head(Ot_deg_numbers)
Ot_deg_numbers$trait = factor(Ot_deg_numbers$trait, levels = c("E","L","G",'AH',"PH"), ordered = TRUE)

Ot_bar_deg <-ggplot(data=Ot_deg_numbers, aes(x=trait, y=upreg.genes, fill=sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylim(0, 1800) +
  geom_text(aes(label=upreg.genes), hjust = -0.3, size=3.5, position = position_dodge(0.9))+ # outside bars
  scale_fill_manual(values=c("#C1272D","#2166AC"))+
  labs(title="O. taurus Sex-biased genes", y= "number of significantly upregulated genes")+
  theme_classic()+
  coord_flip()
Ot_bar_deg

ggsave(Ot_bar_deg, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/figures/Ot_bar_deg.pdf",
       width = 6, height = 3)
##### Step 5: GO term enrichment ###########
### generate gene to GO mapping by joining table of O. taurus UniProt BLAST hits to UniProt GO table ######
uniprot_hits <- read_delim("/Users/ericanadolski/Documents/Sexual_dimorphism_project/GO_enrich/uniprot-blast-hits2.txt", col_names = c("gene","Entry.Name"))
# uniprot GO database
uniprot_go <- read.delim("/Users/ericanadolski/Documents/Sexual_dimorphism_project/GO_enrich/uniprot_reviewed_db.tsv")
go_terms_full <- left_join(uniprot_hits, uniprot_go, by="Entry.Name") 
go_terms <- go_terms_full %>% dplyr::select(gene, Entry.Name, Gene.Ontology.IDs)
head(go_terms)
write.table(go_terms, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/GO_enrich/Ot_go_terms.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
# reformat to  geneID2GO in BBEdit

### read in custom Otau GO annotations in geneID2GO format ####
OtgeneID2GO <- readMappings(file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/GO_enrich/geneID2GO.txt")
str(head(OtgeneID2GO))

# provide predefined list of genes of interest as factor (0,1) for each trait and sex
# for each trait run separate GO analyses for female-upregulated genes and male-upregulated genes

### posterior head female upreg #####
OtPHF_upreg <- Ot_PH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))
head(OtPHF_upreg)
OtPHF_geneList <- as.factor(OtPHF_upreg$F.upreg)
names(OtPHF_geneList) <- Ot_F_upreg$gene
length(OtPHF_geneList) # total number of genes
summary(OtPHF_geneList) # gives # of genes of interest

# build GOdata object
OtPHF_GOdata <- new("topGOdata",
                    description = "GO analysis of Otau PH Female upreg genes",
                    ontology = "BP", 
                    allGenes = OtPHF_geneList, 
                    annot = annFUN.gene2GO, 
                    gene2GO = OtgeneID2GO, 
                    nodeSize = 5)

# run test for significance using the weight01 algorithm (default) with fisher 
OtPHF_weight_fisher_result=runTest(OtPHF_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
OtPHF_allGO=usedGO(OtPHF_GOdata)
OtPHF_all_res=GenTable(OtPHF_GOdata, weightFisher=OtPHF_weight_fisher_result, orderBy='weightFisher', topNodes=length(OtPHF_allGO))
OtPHF_all_res$weightFisher <- as.numeric(OtPHF_all_res$weightFisher)

#get list of significant GO terms
OtPHF_results.table.p= OtPHF_all_res[which(OtPHF_all_res$weightFisher<=0.001),]
dim(OtPHF_results.table.p) # 7

### posterior head male upreg #########
Ot_M_upreg <- Ot_PH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "1",
                                 "0")))
head(Ot_M_upreg)
PH_M_geneList <- as.factor(Ot_M_upreg$M.upreg)
names(PH_M_geneList) <- Ot_M_upreg$gene

PH_M_GOdata <- new("topGOdata",
                   description = "GO analysis of Otau PH Male upreg genes",
                   ontology = "BP", 
                   allGenes = PH_M_geneList, 
                   annot = annFUN.gene2GO, 
                   gene2GO = OtgeneID2GO, 
                   nodeSize = 5)

# run test for significance
PH_M_weight_fisher_result=runTest(PH_M_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
PH_M_allGO=usedGO(PH_M_GOdata)
PH_M_all_res=GenTable(PH_M_GOdata, weightFisher=PH_M_weight_fisher_result, orderBy='weightFisher', topNodes=length(PH_M_allGO))
PH_M_all_res$weightFisher <- as.numeric(PH_M_all_res$weightFisher)

#get list of significant GO 
PH_M_results.table.p= PH_M_all_res[which(PH_M_all_res$weightFisher<=0.001),]
dim(PH_M_results.table.p) # 12 

### anterior head male upreg #########
Ot_AH_M_upreg <- Ot_AH_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "1",
                                 "0")))
AH_M_geneList <- as.factor(Ot_AH_M_upreg$M.upreg)
names(AH_M_geneList) <- Ot_AH_M_upreg$gene

# build GOdata object
AH_M_GOdata <- new("topGOdata",
                   description = "GO analysis of Otau AH Male upreg genes",
                   ontology = "BP", 
                   allGenes = AH_M_geneList, 
                   annot = annFUN.gene2GO, 
                   gene2GO = OtgeneID2GO, 
                   nodeSize = 5)

# run test for significance
AH_M_weight_fisher_result=runTest(AH_M_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
AH_M_allGO=usedGO(AH_M_GOdata)
AH_M_all_res=GenTable(AH_M_GOdata, weightFisher=AH_M_weight_fisher_result, orderBy='weightFisher', topNodes=length(AH_M_allGO))
AH_M_all_res$weightFisher <- as.numeric(AH_M_all_res$weightFisher)

#get list of significant GO 
AH_M_results.table.p= AH_M_all_res[which(AH_M_all_res$weightFisher<=0.001),]
dim(AH_M_results.table.p) # 22 

### anterior female upreg #########
Ot_AH_F_upreg <- Ot_AH_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))
AH_F_geneList <- as.factor(Ot_AH_F_upreg$F.upreg)
names(AH_F_geneList) <- Ot_AH_F_upreg$gene
# build GOdata object
AH_F_GOdata <- new("topGOdata",
                   description = "GO analysis of Otau AH Female upreg genes",
                   ontology = "BP", 
                   allGenes = AH_F_geneList, 
                   annot = annFUN.gene2GO, 
                   gene2GO = OtgeneID2GO, 
                   nodeSize = 5)

# run test for significance using the weight01 algorithm (default) with fisher 
AH_F_weight_fisher_result=runTest(AH_F_GOdata, algorithm='weight01', statistic='fisher') 

# generate a table of results
AH_F_allGO=usedGO(AH_F_GOdata)
AH_F_all_res=GenTable(AH_F_GOdata, weightFisher=AH_F_weight_fisher_result, orderBy='weightFisher', topNodes=length(AH_F_allGO))
AH_F_all_res$weightFisher <- as.numeric(AH_F_all_res$weightFisher)

#get list of significant GO 
AH_F_results.table.p= AH_F_all_res[which(AH_F_all_res$weightFisher<=0.001),]
dim(AH_F_results.table.p) # 14 

### genitalia male upreg #########
Ot_G_M_upreg <- Ot_G_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "1",
                                 "0")))
G_M_geneList <- as.factor(Ot_G_M_upreg$M.upreg)
names(G_M_geneList) <- Ot_G_M_upreg$gene

# build GOdata object
G_M_GOdata <- new("topGOdata",
                  description = "GO analysis of Otau G Male upreg genes",
                  ontology = "BP", 
                  allGenes = G_M_geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = OtgeneID2GO, 
                  nodeSize = 5)

# run test for significance
G_M_weight_fisher_result=runTest(G_M_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
G_M_allGO=usedGO(G_M_GOdata)
G_M_all_res=GenTable(G_M_GOdata, weightFisher=G_M_weight_fisher_result, orderBy='weightFisher', topNodes=length(G_M_allGO))
G_M_all_res$weightFisher <- as.numeric(G_M_all_res$weightFisher)

#get list of significant GO 
G_M_results.table.p= G_M_all_res[which(G_M_all_res$weightFisher<=0.001),]
dim(G_M_results.table.p) # 16 

### genitalia female upreg #########
Ot_G_F_upreg <- Ot_G_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))
G_F_geneList <- as.factor(Ot_G_F_upreg$F.upreg)
names(G_F_geneList) <- Ot_G_F_upreg$gene

# build GOdata object
G_F_GOdata <- new("topGOdata",
                  description = "GO analysis of Otau G Female upreg genes",
                  ontology = "BP", 
                  allGenes = G_F_geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = OtgeneID2GO, 
                  nodeSize = 5)

# run test for significance 
G_F_weight_fisher_result=runTest(G_F_GOdata, algorithm='weight01', statistic='fisher') 

# generate table of results
G_F_allGO=usedGO(G_F_GOdata)
G_F_all_res=GenTable(G_F_GOdata, weightFisher=G_F_weight_fisher_result, orderBy='weightFisher', topNodes=length(G_F_allGO))
G_F_all_res$weightFisher <- as.numeric(G_F_all_res$weightFisher)

#get list of significant GO 
G_F_results.table.p= G_F_all_res[which(G_F_all_res$weightFisher<=0.001),]
dim(G_F_results.table.p) # 15 

### fore tibia male upreg #########
Ot_L_M_upreg <- Ot_L_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "1",
                                 "0")))
L_M_geneList <- as.factor(Ot_L_M_upreg$M.upreg)
names(L_M_geneList) <- Ot_L_M_upreg$gene

# build GOdata object
L_M_GOdata <- new("topGOdata",
                  description = "GO analysis of Otau L Male upreg genes",
                  ontology = "BP", 
                  allGenes = L_M_geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = OtgeneID2GO, 
                  nodeSize = 5)

# run test for significance 
L_M_weight_fisher_result=runTest(L_M_GOdata, algorithm='weight01', statistic='fisher') 

# generate a table of results
L_M_allGO=usedGO(L_M_GOdata)
L_M_all_res=GenTable(L_M_GOdata, weightFisher=L_M_weight_fisher_result, orderBy='weightFisher', topNodes=length(L_M_allGO))
L_M_all_res$weightFisher <- as.numeric(L_M_all_res$weightFisher)

#get list of significant GO 
L_M_results.table.p= L_M_all_res[which(L_M_all_res$weightFisher<=0.001),]
dim(L_M_results.table.p) # 16 

### fore tibia female upreg #########
Ot_L_F_upreg <- Ot_L_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))
L_F_geneList <- as.factor(Ot_L_F_upreg$F.upreg)
names(L_F_geneList) <- Ot_L_F_upreg$gene

# build GOdata object
L_F_GOdata <- new("topGOdata",
                  description = "GO analysis of Otau L Female upreg genes",
                  ontology = "BP", 
                  allGenes = L_F_geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = OtgeneID2GO, 
                  nodeSize = 5)

# run test for significance 
L_F_weight_fisher_result=runTest(L_F_GOdata, algorithm='weight01', statistic='fisher') 

# generate a table of results
L_F_allGO=usedGO(L_F_GOdata)
L_F_all_res=GenTable(L_F_GOdata, weightFisher=L_F_weight_fisher_result, orderBy='weightFisher', topNodes=length(L_F_allGO))
L_F_all_res$weightFisher <- as.numeric(L_F_all_res$weightFisher)

#get list of significant GO 
L_F_results.table.p= L_F_all_res[which(L_F_all_res$weightFisher<=0.001),]
dim(L_F_results.table.p)# 1 

### elytra male upreg #########
Ot_E_M_upreg <- Ot_E_MvF_gene %>% 
  mutate(M.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "0",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "1",
                                 "0")))
E_M_geneList <- as.factor(Ot_E_M_upreg$M.upreg)
names(E_M_geneList) <- Ot_E_M_upreg$gene

# build GOdata object
E_M_GOdata <- new("topGOdata",
                  description = "GO analysis of Otau E Male upreg genes",
                  ontology = "BP", 
                  allGenes = E_M_geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = OtgeneID2GO, 
                  nodeSize = 5)

# run test for significance
E_M_weight_fisher_result=runTest(E_M_GOdata, algorithm='weight01', statistic='fisher') 

# generate a table of results
E_M_allGO=usedGO(E_M_GOdata)
E_M_all_res=GenTable(E_M_GOdata, weightFisher=E_M_weight_fisher_result, orderBy='weightFisher', topNodes=length(E_M_allGO))
E_M_all_res$weightFisher <- as.numeric(E_M_all_res$weightFisher)

#get list of significant GO 
E_M_results.table.p= E_M_all_res[which(E_M_all_res$weightFisher<=0.001),]
dim(E_M_results.table.p) # 21 @ 0.001

### elytra female upreg ########
Ot_E_F_upreg <- Ot_E_MvF_gene %>% 
  mutate(F.upreg = ifelse(padj <= 0.1 & log2FoldChange < 0, "1",
                          ifelse(padj <= 0.1 & log2FoldChange > 0, "0",
                                 "0")))
E_F_geneList <- as.factor(Ot_E_F_upreg$F.upreg)
names(E_F_geneList) <- Ot_E_F_upreg$gene

# build GOdata object
E_F_GOdata <- new("topGOdata",
                  description = "GO analysis of Otau E Female upreg genes",
                  ontology = "BP", 
                  allGenes = E_F_geneList, 
                  annot = annFUN.gene2GO, 
                  gene2GO = OtgeneID2GO, 
                  nodeSize = 5)

# run test for significance
E_F_weight_fisher_result=runTest(E_F_GOdata, algorithm='weight01', statistic='fisher') 

# generate a table of results
E_F_allGO=usedGO(E_F_GOdata)
E_F_all_res=GenTable(E_F_GOdata, weightFisher=E_F_weight_fisher_result, orderBy='weightFisher', topNodes=length(E_F_allGO))
E_F_all_res$weightFisher <- as.numeric(E_F_all_res$weightFisher)

#get list of significant GO 
E_F_results.table.p= E_F_all_res[which(E_F_all_res$weightFisher<=0.001),]
dim(E_F_results.table.p) # 8 

############################# STEP 8 - EXTRA DATA PLOTS
##### Step 6: get gene annotations via UniProt best hits #######
### Otau3 blast best hits to Otau2 
Otau3_prot_hits_clean <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau3_prot_hits_clean.txt")
dim(Otau3_prot_hits_clean) # 242025 rows
Otau3_prot_hits_clean$eval <- as.numeric(Otau3_prot_hits_clean$eval)
Otau3_prot_hits_best <- Otau3_prot_hits_clean %>% 
  group_by(OT3_ID) %>% 
  top_n(-1,eval) %>% 
  dplyr::slice(which.max(pident))
dim(Otau3_prot_hits_best) # 16055 rows

Otau2_prot_names <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau2_prot_names.txt")
Otau3_prot_anno <- left_join(Otau3_prot_hits_best, Otau2_prot_names)
colnames(Otau3_prot_anno)<- c("gene", "OT2_ID", "pident", "eval", "OT2_description")

### UniProt annotations
# go terms dataframe has Protein.names annotation
uniprot_anno <- go_terms_full %>% dplyr::select(gene, Protein.names, Organism)

#####Step 7: annotate and export trait-specific DEG tables ####
deg_Ot_PH_MvF_anno <- join_all(list(Ot_PH_MvF_deg, uniprot_anno, Otau3_prot_anno), by='gene', type='left')
write.table(deg_Ot_PH_MvF_anno, 
            file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/Ot_PH_MvF_deg.txt", sep = "\t", quote = FALSE, row.names = FALSE)

deg_Ot_AH_MvF_anno <- join_all(list(Ot_AH_MvF_deg, Otau3_prot_anno, uniprot_anno), by='gene', type='left')
write.table(deg_Ot_AH_MvF_anno, 
            file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/Ot_AH_MvF_deg.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
deg_Ot_E_MvF_anno <- join_all(list(Ot_E_MvF_deg, Otau3_prot_anno, uniprot_anno), by='gene', type='left')
write.table(deg_Ot_E_MvF_anno, 
            file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/Ot_E_MvF_deg.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
deg_Ot_L_MvF_anno <- join_all(list(Ot_L_MvF_deg, Otau3_prot_anno, uniprot_anno), by='gene', type='left')
write.table(deg_Ot_L_MvF_anno, 
            file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/Ot_L_MvF_deg.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
deg_Ot_G_MvF_anno <- join_all(list(Ot_G_MvF_deg, Otau3_prot_anno, uniprot_anno), by='gene', type='left')
write.table(deg_Ot_G_MvF_anno, 
            file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/RNAseq/Ot_G_MvF_deg.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

########## Step 8: plotting doublesex expression levels #####
### plotting counts per million using'cpm' function from EdgeR
# make the DGEList
Ot_eds_r <- DGEList(Ot_counts)
# calculate TMM normalization factors
Ot_eds_r <- calcNormFactors(Ot_eds_r)
#get the normalized counts
Ot_cpm <- cpm(Ot_eds_r, log=FALSE) # matrix output

# subset counts for dsx isoforms
Ot_dsx <- as.data.frame(t(Ot_cpm[c("jg10020.t1", "jg10020.t2", "jg10020.t3", "jg10020.t4", "jg10020.t5", "jg10020.t6", "jg10020.t7"),]))
Ot_dsx$Sex_Trait <- Ot_sample_table$Sex_Trait

# convert to long data format
Ot_dsx <- rename(Ot_dsx, "F1"="jg10020.t1", "F2"="jg10020.t2", "F3"="jg10020.t3", "F4"="jg10020.t4", "F5"="jg10020.t5", "M"="jg10020.t6", "C"="jg10020.t7")
Ot_dsx_long <- Ot_dsx %>% rownames_to_column("sample") %>% pivot_longer(cols=c("F1", "F2", "F3", "F4", "F5", "M", "C"),names_to='isoform',values_to='counts')

dsx_plot <- ggplot(Ot_dsx_long) +
  geom_boxplot(aes(x=Sex_Trait, y=counts, color=isoform)) +
  scale_x_discrete(limits=c("M_G","F_G","M_PH","F_PH","M_AH","F_AH","M_L","F_L","M_E","F_E")) +
  labs(title="O.tau dsx isoform expression across sample types",x="Sample Type", y = "CPM") +
  theme_classic() + 
  scale_color_manual(values=c("#999999", "#FF808F","#CD9600", "#7CAE00", "#00BE67","#00BFC4","#668FFF"))

ggsave(dsx_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/figures/dsx_expression.pdf",
       width = 5, height = 4)

### Step 9: plotting ventral veinless expression levels #####
#### EdgeR 'cpm' function to get counts per million
Ot_vvl <- as.data.frame(Ot_cpm["jg3443.t1",])
Ot_vvl$Sex_Trait <- Ot_sample_table$Sex_Trait
colnames(Ot_vvl) <- c("jg3443.t1","Sex_Trait")

vvl_plot <- ggplot(Ot_vvl, aes(x=Sex_Trait, y=jg3443.t1)) +
  geom_boxplot() +
  scale_x_discrete(limits=c("F_AH","M_AH","F_PH","M_PH","F_L","M_L","F_G","M_G","F_E","M_E")) +
  labs(title="O.tau ventral veinless mRNA read counts",x="Sample Type", y = "CPM") + theme_classic()

ggsave(vvl_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/figures/vvl_expression.pdf",
       width = 5, height = 4)

# II. ATACseq data analyses ############################################################################
###### Step 1: import ATACseq read counts from bedtools multicov #####
Ot_counts_ATAC <- read.delim("/Users/ericanadolski/GitHub/Beetle-sexual-dimorphism/peak_counts/Ot_counts_05.txt", header=FALSE)
colnames(Ot_counts_ATAC) <- c("chr", "start", "end", "peak", "Ot_F1_AH", "Ot_F1_E", "Ot_F1_G", "Ot_F1_L","Ot_F1_PH", "Ot_F2_AH", "Ot_F2_E", "Ot_F2_G", "Ot_F2_L", "Ot_F2_PH", "Ot_F3_AH", "Ot_F3_E", "Ot_F3_G", "Ot_F3_L", "Ot_F3_PH", "Ot_F4_AH", "Ot_F4_E", "Ot_F4_G", "Ot_F4_L", "Ot_F4_PH", "Ot_F5_AH", "Ot_F5_E", "Ot_F5_G", "Ot_F5_L", "Ot_F5_PH", "Ot_M1_AH", "Ot_M1_E", "Ot_M1_G", "Ot_M1_L", "Ot_M1_PH", "Ot_M2_AH", "Ot_M2_E", "Ot_M2_G", "Ot_M2_L", "Ot_M2_PH", "Ot_M3_AH", "Ot_M3_E", "Ot_M3_G", "Ot_M3_L", "Ot_M3_PH", "Ot_M4_AH", "Ot_M4_E", "Ot_M4_G", "Ot_M4_L", "Ot_M4_PH", "Ot_M5_AH","Ot_M5_E", "Ot_M5_G", "Ot_M5_L", "Ot_M5_PH")

# sample info table
Ot_info_ATAC <- read.delim("/Users/ericanadolski/GitHub/Beetle-sexual-dimorphism/Ot_sample_info_table.txt", row.names = NULL)
Ot_info_ATAC <- Ot_info_ATAC %>% mutate(Sex_Trait = paste0(Sex, "_", Trait),
                              Species_Sex = paste0(Species, "_", Sex),
                              Species_Trait = paste0(Species, "_", Trait),) 
###### Step 2: normalize X chromosome counts and filter out low read counts #####
# write function to multiply by 2
double <- function(observed) {
  result <- (observed * 2)
}

# subset out row of X chromosome
Ot_counts_X <- Ot_counts_ATAC %>%
  filter(chr == "Scaffold8")

# subset female columns
Ot_counts_X_F <- Ot_counts_X[,1:29]

# subset male columns and multiply read counts by 2 to normalize to female's two copies of the X
Ot_counts_X_M_2 <- Ot_counts_X[,30:54] %>%
  mutate_all(list(double))

# add back normalized male to female columns 
Ot_counts_X_norm <- cbind(Ot_counts_X_M_2, Ot_counts_X_F)

# then append X chr rows back
Ot_counts_norm <- Ot_counts_ATAC %>%
  filter(!chr == "Scaffold8")
Ot_counts_norm <- rbind(Ot_counts_norm, Ot_counts_X_norm)

# filter out OCRs with extremely low read counts 
cpm_Ot_norm <- cpm(Ot_counts_norm[c(5:54)]) # 175884
countcheck_Ot_norm <- cpm_Ot_norm > 3
keep_Ot_norm <- which(rowSums(countcheck_Ot_norm) >= 5)
Ot_filtered_counts_norm <- Ot_counts_norm[keep_Ot_norm,] #88405 OCRs

###### Step 3: data visualization ######
### principal component analysis all traits ######
dds_Ot_ATAC <- DESeqDataSetFromMatrix(countData = Ot_filtered_counts_norm[5:54],
                                 colData = Ot_info_ATAC,
                                 design = ~ 1)
rlog_dds_Ot_ATAC <- rlog(dds_Ot_ATAC)
rlog_counts_Ot_ATAC <- as.data.frame(assay(rlog_dds_Ot_ATAC))
rlog_countsT_Ot_ATAC <- as.data.frame(t(rlog_counts_Ot_ATAC))

PCA_Ot_ATAC <- prcomp(rlog_countsT_Ot_ATAC, center = TRUE)
PCs_Ot_ATAC <- as.data.frame(PCA_Ot_ATAC$x)
PCs_Ot_ATAC$Sample <- row.names(PCs_Ot_ATAC)
PCs_Ot_ATAC <- inner_join(PCs_Ot_ATAC, Ot_info_ATAC, by='Sample') 
summary(PCA_Ot_ATAC)

Ot_PCA_ATAC <- ggplot(data = PCs_Ot_ATAC, 
                    aes(x = PC1, y = PC2, color = Sex, shape=Trait)) + 
  geom_point(size = 3) +
  scale_fill_manual(values=c("#C1272D", "#2166AC")) + 
  scale_shape_manual(values = c(17, 3, 15, 18, 16)) +
  geom_text(aes(label = Replicate), nudge_y = 10) +
  labs(title="O. taurus ATAC-seq OCRs",x = "PC1 (28.3%)", y = "PC2 (12.2%)") +
  theme_bw()

ggsave(Ot_PCA_ATAC, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/Ot_PCA_plot_ATAC.pdf",
       width = 6, height = 5)

### PCA genitalia ####
G_info_ATAC <- filter(Ot_info_ATAC, Trait == "G")
G_filtered_counts_ATAC <- Ot_filtered_counts_norm %>% select(contains("G")) 

dds_G_ATAC <- DESeqDataSetFromMatrix(countData = G_filtered_counts_ATAC,
                                     colData = G_info_ATAC,
                                     design = ~ 1)
rlog_dds_G_ATAC <- rlog(dds_G_ATAC) 
rlog_counts_G_ATAC <- as.data.frame(assay(rlog_dds_G_ATAC))
rlog_countsT_G_ATAC <- as.data.frame(t(rlog_counts_G_ATAC))

PCA_G_ATAC <- prcomp(rlog_countsT_G_ATAC, center = TRUE)
PCs_G_ATAC <- as.data.frame(PCA_G_ATAC$x)
PCs_G_ATAC$Sample <- row.names(PCs_G_ATAC)
PCs_G_ATAC <- inner_join(PCs_G_ATAC, Ot_info_ATAC, by='Sample') 
summary(PCA_G_ATAC)

cols <- c( "F" = "#C1272D", "M" = "#2166AC") 

PCA_G_ATAC <- ggplot(data = PCs_G_ATAC, 
                     aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Genitalia",x = "PC1 (40.8%)", y = "PC2 (19.3%)") +
  theme_bw()

ggsave(PCA_G_ATAC, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/PCA_G_ATAC.pdf", width = 4.30, height = 3.77)

### PCA posterior head ####
PH_info_ATAC <- filter(Ot_info_ATAC, Trait == "PH")
PH_filtered_counts_ATAC <- Ot_filtered_counts_norm %>% select(contains("PH")) 

dds_PH_ATAC <- DESeqDataSetFromMatrix(countData = PH_filtered_counts_ATAC,
                                    colData = PH_info_ATAC,
                                    design = ~ 1)
rlog_dds_PH_ATAC <- rlog(dds_PH_ATAC) 
rlog_counts_PH_ATAC <- as.data.frame(assay(rlog_dds_PH_ATAC))
rlog_countsT_PH_ATAC <- as.data.frame(t(rlog_counts_PH_ATAC))

PCA_PH_ATAC <- prcomp(rlog_countsT_PH_ATAC, center = TRUE)
PCs_PH_ATAC <- as.data.frame(PCA_PH_ATAC$x)
PCs_PH_ATAC$Sample <- row.names(PCs_PH_ATAC)
PCs_PH_ATAC <- inner_join(PCs_PH_ATAC, PH_info_ATAC, by='Sample') 
summary(PCA_PH_ATAC)

PCA_PH_ATAC <- ggplot(data = PCs_PH_ATAC, 
                    aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Posterior Head",x = "PC1 (40.8%)", y = "PC2 (19.3%)") +
  theme_bw()

ggsave(PCA_PH_ATAC, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/PCA_PH_ATAC.pdf", width = 4.30, height = 3.77)

### PCA anterior head ####
AH_info_ATAC <- filter(Ot_info_ATAC, Trait == "AH")
AH_filtered_counts_ATAC <- Ot_filtered_counts_norm %>% select(contains("AH")) 

dds_AH_ATAC <- DESeqDataSetFromMatrix(countData = AH_filtered_counts_ATAC,
                                      colData = AH_info_ATAC,
                                      design = ~ 1)
rlog_dds_AH_ATAC <- rlog(dds_AH_ATAC) 
rlog_counts_AH_ATAC <- as.data.frame(assay(rlog_dds_AH_ATAC))
rlog_countsT_AH_ATAC <- as.data.frame(t(rlog_counts_AH_ATAC))

PCA_AH_ATAC <- prcomp(rlog_countsT_AH_ATAC, center = TRUE)
PCs_AH_ATAC <- as.data.frame(PCA_AH_ATAC$x)
PCs_AH_ATAC$Sample <- row.names(PCs_AH_ATAC)
PCs_AH_ATAC <- inner_join(PCs_AH_ATAC, AH_info_ATAC, by='Sample') 
summary(PCA_AH_ATAC)

PCA_AH_ATAC <- ggplot(data = PCs_AH_ATAC, 
                      aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Anterior Head",x = "PC1 (40.8%)", y = "PC2 (19.3%)") +
  theme_bw()

ggsave(PCA_AH_ATAC, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/PCA_AH_ATAC.pdf", width = 4.30, height = 3.77)

### PCA fore tibia ####
L_info_ATAC <- filter(Ot_info_ATAC, Trait == "L")
L_filtered_counts_ATAC <- Ot_filtered_counts_norm %>% select(contains("L")) 

dds_L_ATAC <- DESeqDataSetFromMatrix(countData = L_filtered_counts_ATAC,
                                     colData = L_info_ATAC,
                                     design = ~ 1)
rlog_dds_L_ATAC <- rlog(dds_L_ATAC) 
rlog_counts_L_ATAC <- as.data.frame(assay(rlog_dds_L_ATAC))
rlog_countsT_L_ATAC <- as.data.frame(t(rlog_counts_L_ATAC))

PCA_L_ATAC <- prcomp(rlog_countsT_L_ATAC, center = TRUE)
PCs_L_ATAC <- as.data.frame(PCA_L_ATAC$x)
PCs_L_ATAC$Sample <- row.names(PCs_L_ATAC)
PCs_L_ATAC <- inner_join(PCs_L_ATAC, L_info_ATAC, by='Sample') 
summary(PCA_L_ATAC)

PCA_L_ATAC <- ggplot(data = PCs_L_ATAC, 
                     aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Fore Tibia",x = "PC1 (40.8%)", y = "PC2 (19.3%)") +
  theme_bw()

ggsave(PCA_L_ATAC, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/PCA_L_ATAC.pdf",
       width = 4.30, height = 3.77)

### PCA posterior head ####
E_info_ATAC <- filter(Ot_info_ATAC, Trait == "E")
E_filtered_counts_ATAC <- Ot_filtered_counts_norm %>% select(contains("_E")) 

dds_E_ATAC <- DESeqDataSetFromMatrix(countData = E_filtered_counts_ATAC,
                                     colData = E_info_ATAC,
                                     design = ~ 1)
rlog_dds_E_ATAC <- rlog(dds_E_ATAC) 
rlog_counts_E_ATAC <- as.data.frame(assay(rlog_dds_E_ATAC))
rlog_countsT_E_ATAC <- as.data.frame(t(rlog_counts_E_ATAC))

PCA_E_ATAC <- prcomp(rlog_countsT_E_ATAC, center = TRUE)
PCs_E_ATAC <- as.data.frame(PCA_E_ATAC$x)
PCs_E_ATAC$Sample <- row.names(PCs_E_ATAC)
PCs_E_ATAC <- inner_join(PCs_E_ATAC, E_info_ATAC, by='Sample') 
summary(PCA_E_ATAC)

PCA_E_ATAC <- ggplot(data = PCs_E_ATAC, 
                     aes(x = PC1, y = PC2, color=Sex)) + 
  geom_point(size = 6, shape=16) +
  scale_color_manual(values = cols) +
  labs(title="O.taurus ATAC-seq Elytra",x = "PC1 (40.8%)", y = "PC2 (19.3%)") +
  theme_bw()

ggsave(PCA_E_ATAC, file = "/Users/ericanadolski/Documents/SexuaE_dimorphism_project/ATACseq/figures/PCA_E_ATAC.pdf",
       width = 4.30, height = 3.77)

### heatmap all traits  ######
Ot_sampleDists_ATAC <- dist(rlog_countsT_Ot_ATAC)
Ot_sampleDistMatrix_ATAC <- as.matrix(Ot_sampleDists_ATAC) 
rownames(Ot_sampleDistMatrix_ATAC) <- paste(Ot_info_ATAC$Sample) 
colnames(Ot_sampleDistMatrix_ATAC) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(Ot_sampleDistMatrix_ATAC,
         clustering_distance_rows=Ot_sampleDists_ATAC,
         clustering_distance_cols=Ot_sampleDists_ATAC,
         col=colors)

##### Step 3: differential chromatin accessibility analysis ##########
### differential accessibility between the sexes ####
### genitalia 
dds_G_ATAC <- DESeqDataSetFromMatrix(countData = G_filtered_counts_ATAC,
                                     colData = G_info_ATAC,
                                     design = ~ Sex)
dds_G_ATAC <- DESeq(dds_G_ATAC)
Ot_G_DA_resLFC <- lfcShrink(dds_G_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(Ot_G_DA_resLFC, independentFiltering=FALSE) %>% 
        filter(padj <= 0.1)) # 1078

### posterior head
dds_PH_ATAC <- DESeqDataSetFromMatrix(countData = PH_filtered_counts_ATAC,
                                      colData = PH_info_ATAC,
                                      design = ~ Sex)
dds_PH_ATAC <- DESeq(dds_PH_ATAC)
Ot_PH_DA_resLFC <- lfcShrink(dds_PH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(Ot_PH_DA_resLFC, independentFiltering=FALSE) %>% 
        filter(padj <= 0.1)) # 3150

### anterior head 
dds_AH_ATAC <- DESeqDataSetFromMatrix(countData = AH_filtered_counts_ATAC,
                                      colData = AH_info_ATAC,
                                      design = ~ Sex)
dds_AH_ATAC <- DESeq(dds_AH_ATAC)
Ot_AH_DA_resLFC <- lfcShrink(dds_AH_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(Ot_AH_DA_resLFC, independentFiltering=FALSE) %>% 
        filter(padj <= 0.1)) # 28

### fore tibia
dds_L_ATAC <- DESeqDataSetFromMatrix(countData = L_filtered_counts_ATAC,
                                     colData = L_info_ATAC,
                                     design = ~ Sex)
dds_L_ATAC <- DESeq(dds_L_ATAC)
Ot_L_DA_resLFC <- lfcShrink(dds_L_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(Ot_L_DA_resLFC, independentFiltering=FALSE) %>% 
        filter(padj <= 0.1)) # 234

### elytra
dds_E_ATAC <- DESeqDataSetFromMatrix(countData = E_filtered_counts_ATAC,
                                     colData = E_info_ATAC,
                                     design = ~ Sex)
dds_E_ATAC <- DESeq(dds_E_ATAC)
Ot_E_DA_resLFC <- lfcShrink(dds_E_ATAC, coef="Sex_M_vs_F", type="normal")
tally(as.data.frame(Ot_E_DA_resLFC, independentFiltering=FALSE) %>% 
        filter(padj <= 0.1)) # 69

### save pairwise comparisons into dataframes #####
# genitalia 
OtG_sex_res_OCR <- as.data.frame(Ot_G_DA_resLFC, independentFiltering=FALSE)
OtG_sex_res_OCR <- cbind(Ot_filtered_counts_norm[1:4], OtG_sex_res_OCR)
OtG_sex_res_OCR_sig <- OtG_sex_res_OCR %>% filter(padj <= 0.1)

#posterior head
OtPH_sex_res_OCR <- as.data.frame(Ot_PH_DA_resLFC, independentFiltering=FALSE)
OtPH_sex_res_OCR <- cbind(Ot_filtered_counts_norm[1:4], OtPH_sex_res_OCR)
OtPH_sex_res_OCR_sig <- OtPH_sex_res_OCR %>% filter(padj <= 0.1)

# anterior head 
OtAH_sex_res_OCR<- as.data.frame(Ot_AH_DA_resLFC, independentFiltering=FALSE)
OtAH_sex_res_OCR <- cbind(Ot_filtered_counts_norm[1:4], OtAH_sex_res_OCR)
OtAH_sex_res_OCR_sig <- OtAH_sex_res_OCR %>% filter(padj <= 0.1)

# fore tibia
OtL_sex_res_OCR<- as.data.frame(Ot_L_DA_resLFC, independentFiltering=FALSE)
OtL_sex_res_OCR <- cbind(Ot_filtered_counts_norm[1:4], OtL_sex_res_OCR)
OtL_sex_res_OCR_sig <- OtL_sex_res_OCR %>% filter(padj <= 0.1)

# elytra
OtE_sex_res_OCR<- as.data.frame(Ot_E_DA_resLFC, independentFiltering=FALSE)
OtE_sex_res_OCR <- cbind(Ot_filtered_counts_norm[1:4], OtE_sex_res_OCR)
OtE_sex_res_OCR_sig <- OtE_sex_res_OCR %>% filter(padj <= 0.1)

### volcano plots of sex-responsive peaks ####
cols <- c("F" = "#C1272D", "M" = "#2166AC", "ns" = "grey") 
sizes <- c("F" = 2, "M" = 2, "ns" = 1) 
alphas <- c("F" = 0.5, "M" = 0.5, "ns" = 1)

### genitalia 
# annotate DA peaks as more open in females or males
OtG_sex_res_OCR <- OtG_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
OtG_sex_res_OCR %>% dplyr::count(DA)

Ot_G_DA_vol_plot <- OtG_sex_res_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ylim(0,12) + 
  xlim(-2, 2) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Genitalia")

ggsave(Ot_G_DA_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/Ot_G_DA_vol_plot.pdf",
       width = 4, height = 4)

### posterior head
OtPH_sex_res_OCR <- OtPH_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
OtPH_sex_res_OCR %>% dplyr::count(DA)

Ot_PH_DA_vol_plot <- OtPH_sex_res_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + scale_size_manual(values = sizes) + scale_alpha_manual(values = alphas) + 
  ylim(0,12) + # 2 F, 1 M outliers lfc > 3
  xlim(-2, 2) + theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Posterior Head")
Ot_PH_DA_vol_plot

ggsave(Ot_PH_DA_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/Ot_PH_DA_vol_plot.pdf",
       width = 4, height = 4)

### anterior head 
OtAH_sex_res_OCR <- OtAH_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
OtAH_sex_res_OCR %>% dplyr::count(DA)

Ot_AH_DA_vol_plot <- OtAH_sex_res_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + scale_size_manual(values = sizes) + scale_alpha_manual(values = alphas) + 
  ylim(0,12) + 
  xlim(-2, 2) + theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5), 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Anterior Head")
Ot_AH_DA_vol_plot

ggsave(Ot_AH_DA_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/Ot_AH_DA_vol_plot.pdf",
       width = 4, height = 4)

### fore tibia
OtL_sex_res_OCR <- OtL_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
OtL_sex_res_OCR %>% dplyr::count(DA)

Ot_L_DA_vol_plot <- OtL_sex_res_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + scale_size_manual(values = sizes) + scale_alpha_manual(values = alphas) + 
  ylim(0,12) + # 3 M outliers padj > 12
  xlim(-2, 2) + theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Fore Tibia")
Ot_L_DA_vol_plot

ggsave(Ot_L_DA_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/Ot_L_DA_vol_plot.pdf",
       width = 4, height = 4)

### elytra
OtE_sex_res_OCR <- OtE_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))
OtE_sex_res_OCR %>% dplyr::count(DA)

Ot_E_DA_vol_plot <- OtE_sex_res_OCR %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = DA, size = DA, alpha = DA)) + 
  geom_point(shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = cols) + scale_size_manual(values = sizes) + scale_alpha_manual(values = alphas) + 
  ylim(0,12) + # 2 M outliers log padj > 12
  xlim(-2, 2) + theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="Elytra")
Ot_E_DA_vol_plot

ggsave(Ot_E_DA_vol_plot, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/Ot_E_DA_vol_plot.pdf",
       width = 4, height = 4)

### clustered heatmap of sex-responsive peaks ####
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))

# genitalia 
# subset count matrix by DEGs and trait for heatmap 
G_filtered_counts_ATAC <- cbind(Ot_filtered_counts_norm[1:4], G_filtered_counts_ATAC)
Ot_G_DA_ocr <- inner_join(OtG_sex_res_OCR_sig, G_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_G_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_G_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_G_DA_heatmap <- heatmap.2(as.matrix(Ot_G_DA_ocr[,14:23]), 
                             Rowv=as.dendrogram(clustRows), 
                             Colv=as.dendrogram(clustColumns),
                             RowSideColors=module.color,
                             col=myheatcolors, scale='row', labRow=NA,
                             density.info="none", trace="none",)

### posterior head 
PH_filtered_counts_ATAC <- cbind(Ot_filtered_counts_norm[1:4], PH_filtered_counts_ATAC)
Ot_PH_DA_ocr <- inner_join(OtPH_sex_res_OCR_sig, PH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_PH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_PH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_PH_DA_heatmap <- heatmap.2(as.matrix(Ot_PH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

### anterior head
AH_filtered_counts_ATAC <- cbind(Ot_filtered_counts_norm[1:4], AH_filtered_counts_ATAC)
Ot_AH_DA_ocr <- inner_join(OtAH_sex_res_OCR_sig, AH_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_AH_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_AH_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_AH_DA_heatmap <- heatmap.2(as.matrix(Ot_AH_DA_ocr[,14:23]), 
                              Rowv=as.dendrogram(clustRows), 
                              Colv=as.dendrogram(clustColumns),
                              RowSideColors=module.color,
                              col=myheatcolors, scale='row', labRow=NA,
                              density.info="none", trace="none",)

### fore tibiae
L_filtered_counts_ATAC <- cbind(Ot_filtered_counts_norm[1:4], L_filtered_counts_ATAC)
Ot_L_DA_ocr <- inner_join(OtL_sex_res_OCR_sig, L_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_L_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_L_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_L_DA_heatmap <- heatmap.2(as.matrix(Ot_L_DA_ocr[,14:23]), 
                             Rowv=as.dendrogram(clustRows), 
                             Colv=as.dendrogram(clustColumns),
                             RowSideColors=module.color,
                             col=myheatcolors, scale='row', labRow=NA,
                             density.info="none", trace="none",)

### elytra
E_filtered_counts_ATAC <- cbind(Ot_filtered_counts_norm[1:4], E_filtered_counts_ATAC)
Ot_E_DA_ocr <- inner_join(OtE_sex_res_OCR_sig, E_filtered_counts_ATAC, by="peak")

clustRows <- hclust(as.dist(1-cor(t(Ot_E_DA_ocr[,14:23]), method="pearson")), method="complete")
clustColumns <- hclust(as.dist(1-cor(Ot_E_DA_ocr[,14:23], method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
Ot_E_DA_heatmap <- heatmap.2(as.matrix(Ot_E_DA_ocr[,14:23]), 
                             Rowv=as.dendrogram(clustRows), 
                             Colv=as.dendrogram(clustColumns),
                             RowSideColors=module.color,
                             col=myheatcolors, scale='row', labRow=NA,
                             density.info="none", trace="none",)

### differential accessibility across traits ####
dds_Ot<- DESeqDataSetFromMatrix(countData = Ot_filtered_counts_norm[5:54], 
                                colData = Ot_info, design = ~Trait)
dds_Ot<- DESeq(dds_Ot)

# PH vs AH
Ot_AHPH_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","AH","PH"), type="normal")
Ot_AHPH_shrunk <- as.data.frame(Ot_AHPH_DA_resLFC, independentFiltering=FALSE)
Ot_AHPH_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_AHPH_shrunk)
Ot_AHPH_sig <- Ot_AHPH_shrunk %>% filter(padj <= 0.1)

# PH vs L
Ot_LPH_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","L","PH"), type="normal")
Ot_LPH_shrunk <- as.data.frame(Ot_LPH_DA_resLFC, independentFiltering=FALSE)
Ot_LPH_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_LPH_shrunk)
Ot_LPH_sig <- Ot_LPH_shrunk %>% filter(padj <= 0.1)

# PH vs E
Ot_EPH_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","E","PH"), type="normal")
Ot_EPH_shrunk <- as.data.frame(Ot_EPH_DA_resLFC, independentFiltering=FALSE)
Ot_EPH_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_EPH_shrunk)
Ot_EPH_sig <- Ot_EPH_shrunk %>% filter(padj <= 0.1)

# PH vs G
Ot_GPH_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","G","PH"), type="normal")
Ot_GPH_shrunk <- as.data.frame(Ot_GPH_DA_resLFC, independentFiltering=FALSE)
Ot_GPH_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_GPH_shrunk)
Ot_GPH_sig <- Ot_GPH_shrunk %>% filter(padj <= 0.1)

# AH vs G
Ot_GAH_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","G","AH"), type="normal")
Ot_GAH_shrunk <- as.data.frame(Ot_GAH_DA_resLFC, independentFiltering=FALSE)
Ot_GAH_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_GAH_shrunk)
Ot_GAH_sig <- Ot_GAH_shrunk %>% filter(padj <= 0.1)

# AH vs E
Ot_EAH_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","E","AH"), type="normal")
Ot_EAH_shrunk <- as.data.frame(Ot_EAH_DA_resLFC, independentFiltering=FALSE)
Ot_EAH_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_EAH_shrunk)
Ot_EAH_sig <- Ot_EAH_shrunk %>% filter(padj <= 0.1)

# AH vs L
Ot_LAH_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","L","AH"), type="normal")
Ot_LAH_shrunk <- as.data.frame(Ot_LAH_DA_resLFC, independentFiltering=FALSE)
Ot_LAH_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_LAH_shrunk)
Ot_LAH_sig <- Ot_LAH_shrunk %>% filter(padj <= 0.1)

# E vs L
Ot_LE_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","L","E"), type="normal")
Ot_LE_shrunk <- as.data.frame(Ot_LE_DA_resLFC, independentFiltering=FALSE)
Ot_LE_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_LE_shrunk)
Ot_LE_sig <- Ot_LE_shrunk %>% filter(padj <= 0.1)

# E vs G
Ot_GE_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","G","E"), type="normal")
Ot_GE_shrunk <- as.data.frame(Ot_GE_DA_resLFC, independentFiltering=FALSE)
Ot_GE_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_GE_shrunk)
Ot_GE_sig <- Ot_GE_shrunk %>% filter(padj <= 0.1)

# L vs G
Ot_GL_DA_resLFC <- lfcShrink(dds_Ot, contrast=c("Trait","G","L"), type="normal")
Ot_GL_shrunk <- as.data.frame(Ot_GL_DA_resLFC, independentFiltering=FALSE)
Ot_GL_shrunk <- cbind(Ot_filtered_counts_norm[1:4], Ot_GL_shrunk)
Ot_GL_sig <- Ot_GL_shrunk %>% filter(padj <= 0.1)

##### Step 4: annotate peaks with information on nearest genes ######
#### read in genomic transcript coordinates ###
Ot_transcript_coords <- read.delim("/Users/ericanadolski/Documents/Genomes/Otau3/Otau_transcript_coords.txt")
# correct start and end positions for the negative strand transcripts
Ot_transcript_starts <- Ot_transcript_coords %>% 
  mutate(start_correct = ifelse(orientation == "+", start, end),
         end_correct = ifelse(orientation == "+",end, start)) %>% 
  select(gene, chr, start_gene = start_correct, orientation)

### anterior head ####
# denote which condition showed higher accessibility
OtAH_sex_res_OCR <- OtAH_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))

OtAH_OCR_sig_location  <- OtAH_sex_res_OCR %>%
  filter(padj <= 0.1) %>% mutate(mid = (end-start)/2 + start) %>% select(peak, chr, DA, start, mid, end)

# combine peak list with nearby genes
OtAH_DA_peak_gene_map <- 
  left_join(OtAH_OCR_sig_location, Ot_transcript_starts, by = "chr", relationship = "many-to-many") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% # peak should be within 25kb of gene
  group_by(peak) %>% top_n(-2, abs(peakGeneDist)) #selecting 2 closest genes

# combine peak-gene map with protein annotations
OtAH_DA_peak_gene_map_anno <- left_join(OtAH_DA_peak_gene_map, Otau3_prot_anno, by = "gene")
OtAH_F_peak_map_anno <- OtAH_DA_peak_gene_map_anno %>% filter(DA == "F")
OtAH_M_peak_map_anno <- OtAH_DA_peak_gene_map_anno %>% filter(DA == "M")

### posterior head ####
OtPH_sex_res_OCR <- OtPH_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))

OtPH_OCR_sig_location  <- OtPH_sex_res_OCR %>%
  filter(padj <= 0.1) %>% mutate(mid = (end-start)/2 + start) %>% select(peak, chr, DA, start, mid, end)

# combine peak list with nearby genes
OtPH_DA_peak_gene_map <- 
  left_join(OtPH_OCR_sig_location, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% 
  group_by(peak) %>% top_n(-2, abs(peakGeneDist)) 

# add protein annotations
OtPH_DA_peak_gene_map_anno <- left_join(OtPH_DA_peak_gene_map, Otau3_prot_anno, by = "gene")
OtPH_F_peak_map_anno <- OtPH_DA_peak_gene_map_anno %>% filter(DA == "F")
OtPH_M_peak_map_anno <- OtPH_DA_peak_gene_map_anno %>% filter(DA == "M")

### genitalia ####
OtG_sex_res_OCR <- OtG_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))

OtG_OCR_sig_location  <- OtG_sex_res_OCR %>%
  filter(padj <= 0.1) %>% mutate(mid = (end-start)/2 + start) %>% select(peak, chr, DA, start, mid, end)

# map peaks to nearby genes
OtG_DA_peak_gene_map <- 
  left_join(OtG_OCR_sig_location, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% 
  group_by(peak) %>% top_n(-2, abs(peakGeneDist)) 

# add protein annotations
OtG_DA_peak_gene_map_anno <- left_join(OtG_DA_peak_gene_map, Otau3_prot_anno, by = "gene")
OtG_F_peak_map_anno <- OtG_DA_peak_gene_map_anno %>% filter(DA == "F")
OtG_M_peak_map_anno <- OtG_DA_peak_gene_map_anno %>% filter(DA == "M")

### fore tibia ####
OtL_sex_res_OCR <- OtL_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))

OtL_OCR_sig_location  <- OtL_sex_res_OCR %>%
  filter(padj <= 0.1) %>% mutate(mid = (end-start)/2 + start) %>% select(peak, chr, DA, start, mid, end)

# map peaks to nearby genes
OtL_DA_peak_gene_map <- 
  left_join(OtL_OCR_sig_location, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% 
  group_by(peak) %>% top_n(-2, abs(peakGeneDist)) 

# add protein annotations
OtL_DA_peak_gene_map_anno <- left_join(OtL_DA_peak_gene_map, Otau3_prot_anno, by = "gene")
OtL_F_peak_map_anno <- OtL_DA_peak_gene_map_anno %>% filter(DA == "F")
OtL_M_peak_map_anno <- OtL_DA_peak_gene_map_anno %>% filter(DA == "M")

### elytra ####
OtE_sex_res_OCR <- OtE_sex_res_OCR %>% 
  mutate(DA = ifelse(padj <= 0.1 & log2FoldChange < 0, "F",
                     ifelse(padj <= 0.1 & log2FoldChange > 0, "M",
                            "ns")))

OtE_OCR_sig_location  <- OtE_sex_res_OCR %>%
  filter(padj <= 0.1) %>% mutate(mid = (end-start)/2 + start) %>% select(peak, chr, DA, start, mid, end)

# map peaks to nearby genes
OtE_DA_peak_gene_map <- 
  left_join(OtE_OCR_sig_location, Ot_transcript_starts, by = "chr") %>% 
  mutate(peakGeneDist = mid-start_gene) %>% 
  filter(abs(peakGeneDist) <= 25000) %>% 
  group_by(peak) %>% top_n(-2, abs(peakGeneDist)) 

# add protein annotations
OtE_DA_peak_gene_map_anno <- left_join(OtE_DA_peak_gene_map, Otau3_prot_anno, by = "gene")
OtE_F_peak_map_anno <- OtE_DA_peak_gene_map_anno %>% filter(DA == "F")
OtE_M_peak_map_anno <- OtE_DA_peak_gene_map_anno %>% filter(DA == "M")

##### Step 5: plot bar chart of sex-responsive OCRs across traits ####
sex <- c('F', 'M','F', 'M','F', 'M','F', 'M','F', 'M') 
trait <- c('PH', 'PH','G','G','L','L','AH','AH','E','E')
DA_peaks <- c(1237,1913,236,842,58,176,7,21,22,47)

Ot_DA_numbers <- data.frame(sex,trait,DA_peaks)
head(Ot_DA_numbers)
Ot_DA_numbers$trait = factor(Ot_DA_numbers$trait, levels = c("E","L","G",'AH',"PH"), ordered = TRUE)

# Otaurus
Ot_bar_DA <-ggplot(data=Ot_DA_numbers, aes(x=trait, y=DA_peaks, fill=sex)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  ylim(0, 2000) +
  geom_text(aes(label=DA_peaks), hjust = -0.3, size=3.5, position = position_dodge(0.9))+ # outside bars
  scale_fill_manual(values=c("#C1272D","#2166AC"))+
  labs(title="O. taurus Sex-responsive Chromatin", y= "number of significantly more accessible OCRs")+
  theme_classic()+
  coord_flip()
Ot_bar_DA

ggsave(Ot_bar_DA, file = "/Users/ericanadolski/Documents/Sexual_dimorphism_project/ATACseq/figures/Ot_bar_DA.pdf",
       width = 6, height = 3)

##### STEP 7: assess peaks near doublesex ######
## dsx start & end coordinates = Scaffold 7, 7503131, 7570866
# create list of all peaks near dsx (within 25kb up or downstream)
near_dsx <- Ot_counts_5 %>% 
  select(1:4) %>% 
  filter(chr == "Scaffold7") %>% 
  filter(start > 7478131 & start < 7595866)

# filter by list of trait_responsive peaks
trait_res_list <- read.delim("/Users/ericanadolski/GitHub/Beetle-sexual-dimorphism/peak-sets/Otau-sig/trait_res_sorted.txt", col.names = c("peak"))
trait_res_near_dsx <- inner_join(trait_res_list, near_dsx, by = "peak")
nrow(trait_res_near_dsx)

# filter by list of sex_responsive peaks
sex_res_list <- read.delim("/Users/ericanadolski/GitHub/Beetle-sexual-dimorphism/peak-sets/Otau-sig/sex_res_sorted.txt", col.names = c("peak"))
sex_res_near_dsx <- inner_join(sex_res_list, near_dsx, by = "peak")
nrow(sex_res_near_dsx)

##### STEP 8: EXPORT significant peak tables ######
### export sex-responsive peaks sets ####
write.table(OtPH_sex_res_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtPH_sex_res_peaks.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtAH_sex_res_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtAH_sex_res_peaks.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtL_sex_res_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtL_sex_res_peaks.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtG_sex_res_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtG_sex_res_peaks.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtE_sex_res_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtE_sex_res_peaks.txt", sep = "\t", quote = FALSE, row.names = FALSE,)

### export sex-responsive peak sets annotated with closest genes ####
write.table(OtAH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtAH_F_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtAH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtAH_M_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtE_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtE_F_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtE_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtE_M_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtG_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtG_F_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtG_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtG_M_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtL_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtL_F_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtL_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtL_M_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtPH_F_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtPH_F_peak_map_anno.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(OtPH_M_peak_map_anno, 
            file = "./GitHub/Beetle-sexual-dimorphism/DAP-tables-annotated/OtPH_M_peak_map_anno.txt",sep = "\t", quote = FALSE, row.names = FALSE)

# export sex-specific OCR lists - bed file for motif enrichments
OtPH_peaks_M_up <- filter(OtPH_sex_res_sig_anno, DA =="M")
OtPH_peaks_F_up <- filter(OtPH_sex_res_sig_anno, DA =="F")
OtAH_peaks_M_up <- filter(OtAH_sex_res_sig_anno, DA =="M")
OtAH_peaks_F_up <- filter(OtAH_sex_res_sig_anno, DA =="F")
OtL_peaks_M_up <- filter(OtL_sex_res_sig_anno, DA =="M")
OtL_peaks_F_up <- filter(OtL_sex_res_sig_anno, DA =="F")
OtG_peaks_M_up <- filter(OtG_sex_res_sig_anno, DA =="M")
OtG_peaks_F_up <- filter(OtG_sex_res_sig_anno, DA =="F")
OtE_peaks_M_up <- filter(OtE_sex_res_sig_anno, DA =="M")
OtE_peaks_F_up <- filter(OtE_sex_res_sig_anno, DA =="F")

write.table(OtPH_peaks_M_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtPH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtPH_peaks_F_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtPH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtAH_peaks_M_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtAH_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtAH_peaks_F_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtAH_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtL_peaks_M_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtL_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtL_peaks_F_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtL_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtG_peaks_M_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtG_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtG_peaks_F_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtG_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtE_peaks_M_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtE_peaks_M_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(OtE_peaks_F_up, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/OtE_peaks_F_up.txt", sep = "\t", quote = FALSE, row.names = FALSE,)

### export all trait-responsive peak sets ####
write.table(Ot_AHPH_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_AHPH_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_EPH_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_EPH_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_LPH_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_LPH_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_GPH_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_GPH_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_EAH_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_EAH_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_GAH_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_GAH_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_LAH_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_LAH_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_LE_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_LE_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_GE_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_GE_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)
write.table(Ot_GL_sig, file = "./GitHub/Beetle-sexual-dimorphism/peak-sets/Ot_GL_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE,)