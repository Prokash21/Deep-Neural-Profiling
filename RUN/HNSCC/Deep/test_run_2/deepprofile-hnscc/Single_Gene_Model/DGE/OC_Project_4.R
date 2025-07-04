
BiocManager::install("DESeq2")
install.packages("pheatmap")
install.packages("umap")
BiocManager::install("WGCNA")



# Load libraries
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(umap)
library(ggplot2)
library(WGCNA)

# Load Data
#
# Expression Data
# count_data_x_ev.csv

# load the count data
count_data <- read.csv("GSE290057_QC_clean.csv", header = TRUE, row.names = 1)
head(count_data)

# Metadata
# meta_data_x_ev.csv

# load the sample info
sample_info <- read.csv("GSE290057_meta_data.csv", header = TRUE, row.names = 1)
head(sample_info)
# Ensure the number of samples match
if (ncol(count_data) != nrow(sample_info)) {
  stop("Number of samples in count_data and sample_info do not match!")
}

# Convert the Col name of user to "Treatment"
colnames(sample_info)[colnames(sample_info) == colnames(sample_info)] <- "Treatment"

# Convert Condition to factor
sample_info$Treatment <- factor(sample_info$Treatment)

##############################################################################################################

##############################################################################################################

# Quality Control: detect outlier genes
library(WGCNA)
gsg <- goodSamplesGenes(t(count_data)) # check the data format
summary(gsg)

# Remove outlier genes
data <- count_data[gsg$goodGenes == TRUE,]

################################################# PCA ########################################################

#***** Make A Loop *****#

# Check for NA or Infinite values
summary(data)
is.na(data)

is.infinite(data)
# Replace NA and Infinite values with zero

data[is.na(data)] <- 0
data[is.infinite(data)] <- 0

# Verify no NA or Infinite values remain
summary(data)

#***** Make A Loop *****#


# Remove non-numeric columns for PCA
# Remove the gene ID column
data_numeric <- data[, sapply(data, is.numeric)]

# data_numeric <- data[,1:12]

# Perform PCA
pca <- prcomp(t(data_numeric))

# View the PCA results
summary(pca)

# Prepare PCA data for plotting
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, digits = 2)

# Merge PCA data with metadata
pca.dat <- cbind(pca.dat, sample_info)

library(ggplot2)

ggplot(pca.dat, aes(PC1, PC2, color = Treatment)) +
  geom_point(size = 3) +  # Adjust point size to match t-SNE
  geom_text(aes(label = rownames(pca.dat)), hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) +  # Bold, smaller labels
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %')) +
  theme_minimal() +
  
  # Adjust the legend
  theme(legend.position = "bottom",                            # Move the legend to the bottom
        legend.title = element_text(size = 10),                # Legend title size
        legend.text = element_text(size = 8),                  # Legend text size
        legend.key.size = unit(0.5, "cm")) +                   # Reduce legend key size
  guides(color = guide_legend(nrow = 2, byrow = TRUE))          # Organize legend into 2 rows
################################################# UMAP #####################################################


# Load libraries
library(umap)

# Load Data
# data <- read.csv("count_data.csv")
# Remove non-numeric columns for Umap
# Remove the gene ID column
# data_numeric <- data[, sapply(data, is.numeric)]

# load the metadata
# sample_info <- read.csv("meta_data.csv", header =TRUE,row.names = 1)

# Set random seed for reproducibility
set.seed(123)
#to do
n_samples <- nrow(data_numeric)  # Number of samples in your dataset
min_dist <- 1 / sqrt(n_samples)
n_neighbors <- round(sqrt(n_samples))

# Perform UMAP dimensionality reduction
umap_result <- umap(t(data_numeric), n_neighbors = 5, min_dist = 0.5)

# Extract UMAP coordinates and combine with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  sample_info
)


# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = Treatment)) +
  geom_point(size = 3) +  # Adjust point size
  geom_text(aes(label = rownames(sample_info)), hjust = .5, vjust = 1.5, size = 3, show.legend = FALSE) +  # Adjust label size and position
  labs(title = "UMAP",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  
  # Adjust the legend
  theme(legend.position = "bottom",                            # Move the legend to the bottom
        legend.title = element_text(size = 10),                # Legend title size
        legend.text = element_text(size = 8),                  # Legend text size
        legend.key.size = unit(0.5, "cm")) +                   # Reduce legend key size
  guides(color = guide_legend(nrow = 2, byrow = TRUE))          # Organize legend into 2 rows


################################################# t-SNE ###############################################

library(Rtsne)

# Load Data
# data <- read.csv("count_data.csv")
# Remove non-numeric columns for t-SNE
# Remove the gene ID column
# data_numeric <- data[, sapply(data, is.numeric)]

# load the metadata
# sample_info <- read.csv("meta_data.csv", header =TRUE,row.names = 1)

set.seed(123)

# Perform t-SNE
num_samples <- length(rownames(sample_info))
perplexity=num_samples*.26
# to do # date: updated on 6th Aug**
tsne_result <- Rtsne(t(data_numeric), dims=2, perplexity=perplexity, verbose=TRUE, max_iter=500)
# Create a data frame for plotting
tsne_data <- data.frame(
  X = tsne_result$Y[,1],
  Y = tsne_result$Y[,2],
  sample_info
)


# Plot the t-SNE results using ggplot2 with smaller sample name labels
ggplot(tsne_data, aes(x = X, y = Y, color = Treatment)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(sample_info)), hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) + # Reduced text size
  theme_minimal() +
  ggtitle("t-SNE Plot") +
  xlab("t-SNE 1") +
  ylab("t-SNE 2") +
  
  # Adjust the legend
  theme(legend.position = "bottom",                            # Move the legend to the bottom
        legend.title = element_text(size = 10),                # Legend title size
        legend.text = element_text(size = 8),                  # Legend text size
        legend.key.size = unit(0.5, "cm")) +                   # Reduce legend key size
  guides(color = guide_legend(nrow = 2, byrow = TRUE))          # Organize legend into 2 rows


##################################################### k means ##################################################
# Normalize the data (optional but recommended)
data_numeric_scaled <- t(scale(data_numeric))

num_conditions <- as.numeric(length(unique(sample_info$Treatment)))

# Combine numeric data and treatment data
combined_data <- cbind(data_numeric_scaled, Treatment_num = sample_info$Treatment_num)

# Set seed for reproducibility
set.seed(123)

# Perform K-means clustering using num_conditions as the number of clusters
num_clusters <- num_conditions  # Use num_conditions here
kmeans_result <- kmeans(combined_data, centers = num_clusters)

# Create a data frame for plotting
kmeans_data_plot <- data.frame(
  X = combined_data[, 1],  # Assuming you want to use the first feature for plotting
  Y = combined_data[, 2],  # Assuming you want to use the second feature for plotting
  Cluster = as.factor(kmeans_result$cluster),
  Treatment = sample_info$Treatment
)

# Plot K-means clustering results using ggplot2
ggplot(kmeans_data_plot, aes(x = X, y = Y, color = Treatment)) +
  geom_point(size = 3) +  # Adjust point size
  geom_text(aes(label = rownames(kmeans_data_plot)), hjust = 0.5, vjust = 1.5, size = 3, show.legend = FALSE) +  # Add labels
  labs(title = "K-means Clustering",  # Set title
       x = "Feature 1",               # X-axis label
       y = "Feature 2",               # Y-axis label
       color = "Treatment") +          # Legend title
  theme_minimal() +                   # Minimal theme
  
  # Adjust the legend
  theme(legend.position = "bottom",                            # Move the legend to the bottom
        legend.title = element_text(size = 10),                # Legend title size
        legend.text = element_text(size = 8),                  # Legend text size
        legend.key.size = unit(0.5, "cm")) +                   # Reduce legend key size
  guides(color = guide_legend(nrow = 2, byrow = TRUE))          # Organize legend into 2 rows


########################################### Phylogenetic Tree ###############################################
set.seed(123)
install.packages("ape")
library(ape)


# Detect outlier samples using hierarchical clustering
htree <- hclust(dist(t(data_numeric)), method = "average")

phylo_tree <- as.phylo(htree)


# Assuming 'htree' and 'metadata' are correctly defined

# Convert labels to factors or unique numeric indices for colors
group_colors <- as.factor(sample_info$Treatment)
tip_colors <- as.numeric(group_colors)


# Plot the dendrogram with colored labels

plot.phylo(phylo_tree,
           type = "phylogram",
           tip.color = tip_colors,
           cex = 0.8,
           main = "Hierarchical Clustering Dendrogram"
)

# Add a legend for the treatment groups
legend("topleft",
       legend = levels(group_colors),
       col = 1:length(levels(group_colors)),
       pch = 19,
       cex = 0.8,
       bty = "n"
)


# Load Data
# data <- read.csv("count_data.csv")
# Remove non-numeric columns for htree
# Remove the gene ID column
# data_numeric <- data[, sapply(data, is.numeric)]

# load the metadata
# sample_info <- read.csv("meta_data.csv", header =TRUE,row.names = 1)

# Detect outlier samples using hierarchical clustering
tree <- dist(t(data_numeric))
htree <- hclust(tree, method = "average")

# Convert labels to factors or unique numeric indices for colors
label_colors <- as.numeric(factor(sample_info$Treatment))

# Plot the dendrogram with colored labels
plot(htree, labels = rownames(sample_info), main = "Hierarchical Clustering Dendrogram", sub = "", xlab = "", ylab = "")

##**********************#

##**********************##

########################################## Normalization | NODE 3 ###################################################


# Convert non-integer values to integers in count data
data <- round(data)
head(data)

# Create a new count data object
new_data <- as.matrix(data)
head(new_data)

# Display dimensions for verification
cat("Dimensions of data:", dim(data), "\n")
cat("Dimensions of new_data:", dim(new_data), "\n")
cat("Dimensions of sample_info:", dim(sample_info), "\n")

# Start Normalization

library(DESeq2)


# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_data, colData = sample_info, design = ~ Treatment)

# Set the factor levels for the Treatment column based on unique values
condition <- unique(sample_info$Treatment)
dds$Treatment <- factor(dds$Treatment, levels = condition)

# Input all factor from metadata # TO do

# Filter genes with low counts (less than 75% of sample number)
threshold <- round(dim(sample_info)[1] * 0.70)
keep <- rowSums(counts(dds)) >= threshold
dds <- dds[keep,]

# Perform DESeq2 analysis
dds <- DESeq(dds)

# save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"Normalized_Count_Data.csv")


##################################################################################################################################################
# Log2 transformation for count data
count_matrix <- counts(dds) + 1  # Adding 1 to avoid log(0)
log2_count_matrix <- log2(count_matrix)

boxplot(log2_count_matrix, outline = FALSE, main = "Boxplot of Log2-transformed Count Data",
        cex.main = 0.9,  # Make title size smaller
        ylab = "Log2-transformed Counts",
        cex.axis = 0.7, las = 1, font.axis = 1, xaxt = "n")  # Set font.axis to 1 for normal

# Rotate x-axis labels to 45 degrees
text(x = 1:ncol(log2_count_matrix), 
     y = par("usr")[3] - 0.5,  # Adjust y position as needed
     labels = colnames(log2_count_matrix),
     srt = 45, adj = 1, xpd = TRUE, cex = 0.7)  # Removed font for normal

# Log2 transformation for normalized count data
normalized_counts <- counts(dds, normalized = TRUE)
log2_normalized_counts <- log2(normalized_counts + 1)  # Adding 1 to avoid log(0)

boxplot(log2_normalized_counts, outline = FALSE,
        main = "Boxplot of Log2-transformed Normalized Count Data",
        cex.main = 0.9,  # Make title size smaller
        ylab = "Log2-transformed Counts",
        cex.axis = 0.7, las = 1, font.axis = 1, xaxt = "n")  # Set font.axis to 1 for normal

# Rotate x-axis labels to 45 degrees
text(x = 1:ncol(log2_normalized_counts), 
     y = par("usr")[3] - 0.5,  # Adjust y position as needed
     labels = colnames(log2_normalized_counts),
     srt = 45, adj = 1, xpd = TRUE, cex = 0.7)  # Removed font for normal

################################################### DEGs, LFC, FDR | NODE 5 ###################################################

# Find out P values, Log Fold Change(LFC) Values, False Discovery Rate (FDR) Values 

condition <- as.data.frame(condition)

print(condition)

# condition
# 1                    mock
# 2 MPXV_clade_IIa_infected
# 3 MPXV_clade_IIb_infected
# 4   MPXV_clade_I_infected

# Ensure the reference level is a character string
ref_level <- as.character(condition[2,]) #User selection

# set the reference/ control for the treatment factor
dds$Treatment <- relevel(dds$Treatment, ref = ref_level)# User input ref = "....."

# Perform DESeq2 analysis
dds <- DESeq(dds)

# Identify available coefficient names
coeff_names <- as.data.frame(resultsNames(dds))

# Print the coefficient names
print(coeff_names)

# 1                      Intercept
# 2     Treatment_eADMSC_vs_eBMMSC
# 3      Treatment_BMMSC_vs_eBMMSC
# 4      Treatment_ADMSC_vs_eBMMSC
# 5       Treatment_MCF7_vs_eBMMSC
# 6       Treatment_HeLa_vs_eBMMSC
# 7   Treatment_MDAMB231_vs_eBMMSC
# 8        Treatment_TM6_vs_eBMMSC
# 9      Treatment_eHeLa_vs_eBMMSC
# 10     Treatment_eMCF7_vs_eBMMSC
# 11 Treatment_eMDAMB231_vs_eBMMSC
# 12      Treatment_A549_vs_eBMMSC
# 13     Treatment_H1975_vs_eBMMSC
# 14     Treatment_eA549_vs_eBMMSC
# 15    Treatment_eH1975_vs_eBMMSC

# 9, 10, 11, 14, 15
# For User: Select 1 or 2 or 3 or 4
# BiocManager::install("apeglm")
X <- coeff_names[2,]


resLFC <- lfcShrink(dds, coef =X  , type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)
is.na(resLFC)
# Assuming resLFC is your data frame
#resLFC_filtered <- resLFC[!apply(is.na(resLFC), 1, any), ]
resLFC <- resLFC[!apply(is.na(resLFC), 1, any), ]
resLFC <- as.data.frame(resLFC)

# Researchers are often interested in minimizing the number of false discoveries. 
# Only Keep the significant genes padj (FDR) values is less than 0.05
# resLFC_p_cut <- resLFC[resLFC$padj < 0.05,]

##************** remove the NA rows ************##

# create histogram plot of p-values
hist(resLFC$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "Frequencies of padj-values")

summary(resLFC)


# Calculate the number of genes with padj < 0.05
num_significant_genes <- sum(resLFC$padj < 0.05)

# Create histogram plot of p-values
hist(resLFC$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "Frequencies of padj-values")

# Add text to indicate the number of significant genes
mtext(paste("Number of genes with padj < 0.05:", num_significant_genes), side = 3, line = 0.5, adj = 0.5)

# to do # date: updated on 6th Aug**
write.csv(resLFC, file = paste0('resLFC_', X,'.csv'), row.names = TRUE)

################################################################################


# Upregulated genes
Upregulated <- resLFC[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05, ]
Upregulated_padj <- Upregulated[order(Upregulated$padj ),]
write.csv(Upregulated_padj, file = paste0('Upregulated_padj_', X,'.csv'), row.names = TRUE)

# Downregulated genes
Downregulated <- resLFC[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05, ]
Downregulated_padj <- Downregulated[order(Downregulated$padj),]
write.csv(Downregulated_padj, file = paste0('Downregulated_padj_', X,'.csv'), row.names = TRUE )

################################################### Volcano Plot | NODE 6 ###########################################################

################################################################################

library(dplyr)
down_num <- nrow(Downregulated)
up_num <- nrow(Upregulated)
# Create a Volcano Plot

library(ggplot2)

# User will select specific gene names to label and highlight with borders
# genes_to_highlight <- c("URS00004A2461_Y1", "URS00004A2461_Y3", "URS00004A2461_Y4","URS00004A2461_Y5")  
# genes_to_highlight <- c("URS000007CA8F_Y1", "URS000014AF0D_Y1", "URS00002273D1_Y1", 
#                         "URS000048D319_Y1", "URS000048F2DE_Y1", "URS0000103047_Y1", 
#                         "URS00003AF995_Y1", "URS0000598D4C_Y1", "URS00002017A0_Y1", 
#                         "URS000016432E_Y3", "URS000017CEAC_Y3", "URS00001B6FFF_Y3", 
#                         "URS0000244EDB_Y3", "URS000026BC5F_Y3", "URS0000278337_Y3", 
#                         "URS0000294A73_Y3", "URS000029BC51_Y3", "URS000029DE41_Y3", 
#                         "URS000033C665_Y3", "URS00003D1FB8_Y3", "URS00005CF03F_Y3", 
#                         "URS00003FA741_Y3", "URS000018E9D7_Y4", "URS0000007D24_Y4", 
#                         "URS0000188F7D_Y4", "URS00004A2461_Y5")  
# genes_to_highlight <- c("URS0000103047_Y1", "URS00005CF03F_Y3", "URS0000007D24_Y4","URS00004A2461_Y5")  

################################################################################
genes_to_highlight <- c("SLC9A7", "RAP1GAP2")  

# Adding gene names as a column
resLFC$gene <- row.names(resLFC)

volcano_plot <- ggplot(resLFC, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange > 1.0 & -log10(padj) > 1.3, "Upregulated",
                                ifelse(log2FoldChange < -1.0 & -log10(padj) > 1.3, "Downregulated", "Not Significant"))),
             size = 2.5, alpha = 0.5) +
  
  # Add horizontal dashed line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  
  # Add vertical dashed lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  # Customize plot labels and add the header
  labs(
    title = paste0("Volcano Plot of ", X, " (Up: ", up_num, ", Down: ", down_num, ")"),
    x = "Log2 Fold Change",
    y = "-log10(padj)",
    color = "Regulation"
  ) +
  
  # Customize color palette for regulation categories
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  
  # Use a minimal theme for the plot
  theme_minimal() +
  
  # Add borders around specific genes and label them
  geom_point(data = subset(resLFC, gene %in% genes_to_highlight),
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, size = 4, stroke = 1, color = "black", fill = NA) +  # shape = 21 for a circle with a border
  geom_text(data = subset(resLFC, gene %in% genes_to_highlight),
            aes(label = gene),
            vjust = -0.5, hjust = 0.5, size = 3, color = "black")

# Save the plot immediately after creating it
file_name <- paste0("Volcano_Plot_", X, ".pdf")
ggsave(filename = file_name, plot = volcano_plot, width = 8, height = 6)

# Print confirmation message
print(paste("Plot saved as:", file_name))





















































































































#
#
#
#
# Create a Volcano Plot
ggplot(resLFC,aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange > 1.0 & -log10(padj) > 1.3, "Upregulated",
                                ifelse(log2FoldChange < -1.0 & -log10(padj) > 1.3, "Downregulated", "Not Significant"))),
             size = 2.5, alpha = 0.5) +
  
  # Add horizontal dashed line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  
  # Add vertical dashed lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  # Customize plot labels and add the header
  labs(
    title = paste0("Volcano Plot of ", X, " (Up: ", up_num, ", Down: ", down_num, ")"),# Add the header here
    x = "Log2 Fold Change",
    y = "-log10(padj)",
    color = "Regulation"
  ) +
  
  # Customize color palette for regulation categories
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  
  # Use a minimal theme for the plot
  theme_minimal()

################################################################################
# Create a Volcano Plot
volcano_plot <- ggplot(resLFC, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange > 1.0 & -log10(padj) > 1.3, "Upregulated",
                                ifelse(log2FoldChange < -1.0 & -log10(padj) > 1.3, "Downregulated", "Not Significant"))),
             size = 2.5, alpha = 0.5) +
  
  # Add horizontal dashed line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  
  # Add vertical dashed lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  # Customize plot labels
  labs(
    title = paste0("Volcano Plot of ", X, " (Up: ", up_num, ", Down: ", down_num, ")"),
    x = "Log2 Fold Change",
    y = "-log10(padj)",
    color = "Regulation"
  ) +
  
  # Customize color palette for regulation categories
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  
  # Use a minimal theme for the plot
  theme_minimal()

# Save the plot immediately after creating it
file_name <- paste0("Volcano_Plot_", X, ".pdf")
ggsave(filename = file_name, plot = volcano_plot, width = 8, height = 6)

# Print confirmation message
print(paste("Plot saved as:", file_name))

###################################################################################################################