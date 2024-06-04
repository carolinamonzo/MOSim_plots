## Kendall correlations for 3C

# We do log1p then the mean between cells of the same celltype so we have one 
# column per celltype. We compare the profiles the true associations with 
# Kendall (TRUE). We remove all the genes from the clusters from the 
# expression dataframe and run Kendall on the other combinations (FALSE).

devtools::load_all("~/workspace/1_conesalab/MOSim")
suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
  library(WGCNA)
  library(Seurat)
  library(RColorConesa)
  library(acorde)
  library(corrr)
})
library(pbmcMultiome.SeuratData)

set.seed(000)

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells8clus8000_scMOSim_2groups_022.rds")

## Plotting big sample
#sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_few_2602.rds")

# create cell-to-cell-type ID table
ct <- tibble::tibble(Cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

rna <- as.data.frame(acorde::scale_isoforms(sim$Group_1$Rep_1$`sim_scRNA-seq`@counts, 
                              isoform_col = NULL))

rna[is.na(rna)] <- 0

atac <- as.data.frame(acorde::scale_isoforms(sim$Group_1$Rep_1$`sim_scATAC-seq`@counts, 
                              isoform_col = NULL))

atac[is.na(atac)] <- 0

rownames(rna) <- rna$transcript
rna$transcript <- NULL
rownames(atac) <- atac$transcript
atac$transcript <- NULL

asocG2 <- sim$AssociationMatrices$AssociationMatrix_Group_2[sim$AssociationMatrices$AssociationMatrix_Group_2$Peak_cluster %in% c(1:11),]
asocG2 <- asocG2[asocG2$Gene_cluster %in% c(1:11),]

calculate_mean_per_list_df <- function(df, named_lists) {
  means <- list()
  for (name in names(named_lists)) {
    columns <- named_lists[[name]]
    means[[name]] <- rowMeans(df[, columns, drop = FALSE])
  }
  # Combine the list of means into a dataframe
  means_df <- do.call(cbind, means)
  # Add column names
  colnames(means_df) <- names(named_lists)
  return(means_df)
}

cell_types <- sim$cellTypes

means_celltype_rna <- calculate_mean_per_list_df(rna, cell_types)
means_celltype_atac <- calculate_mean_per_list_df(atac, cell_types)

## Remove flat profiles. Genes or peaks where the rowsum of absolute values is less than 0.1
row_sums <- rowSums(abs(means_celltype_rna))
means_celltype_rna <- means_celltype_rna[row_sums >= 0.1, ]

row_sums <- rowSums(abs(means_celltype_atac))
means_celltype_atac <- means_celltype_atac[row_sums >= 0.1, ]

genes <- means_celltype_rna[match(asocG2$Gene_ID, rownames(means_celltype_rna)),]
peaks <- means_celltype_atac[match(asocG2$Peak_ID, rownames(means_celltype_atac)),]

# Drop rows of NA
genes <- genes[rowSums(is.na(genes)) != ncol(genes), ]
peaks <- peaks[rowSums(is.na(peaks)) != ncol(peaks), ]

# Function to compute Kendall correlation between rows
compute_kendall <- function(gene_row, peak_row) {
  return(cor(gene_row, peak_row, method = "kendall"))
}

# Initialize a matrix to store the correlation values
correlation_matrix <- matrix(nrow = nrow(genes), ncol = nrow(peaks))
rownames(correlation_matrix) <- rownames(genes)
colnames(correlation_matrix) <- rownames(peaks)

# Compute Kendall correlation for each gene and peak
for (i in 1:nrow(genes)) {
  for (j in 1:nrow(peaks)) {
    correlation_matrix[i, j] <- compute_kendall(as.numeric(genes[i, ]), as.numeric(peaks[j, ]))
  }
}

# Convert the correlation matrix to a data frame
correlation_df <- as.data.frame(correlation_matrix)

# Remove rows and columns where all are NA
correlation_df <- correlation_df[rowSums(is.na(correlation_df)) != ncol(correlation_df),]
correlation_df <- correlation_df[, colSums(is.na(correlation_df)) != nrow(correlation_df)]

# Wide to long
correlation_df$Gene_ID <- row.names(correlation_df)

long_trues <- correlation_df %>% pivot_longer(!Gene_ID, names_to = "Peak_ID", values_to = "cors")

# Bring the cluster names
long_trues <- merge(long_trues, asocG2[, c("Gene_ID", "Gene_cluster")], by = c("Gene_ID"))
long_trues <- merge(long_trues, asocG2[, c("Peak_ID", "Peak_cluster")], by = c("Peak_ID"))

# Keep only those that belong to the same cluster or opposite.
trues1 <- long_trues[long_trues$Gene_cluster==long_trues$Peak_cluster, ]

