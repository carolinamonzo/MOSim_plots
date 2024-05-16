devtools::load_all("~/workspace/1_conesalab/MOSim")
suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
  library(WGCNA)
  library(Seurat)
  library(RColorConesa)
  library(acorde)
})
library(pbmcMultiome.SeuratData)

set.seed(000)

#sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells11clus8800_scMOSim_2groups_000.rds")
sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_few_2602.rds")

# Normalize using Seurat's method
#dat <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
#dat <- Seurat::SCTransform(dat, verbose = FALSE, return.only.var.genes = FALSE)

# Create cell-tocelltype ID table
rna <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
rna <- Seurat::SCTransform(rna, verbose = FALSE, return.only.var.genes = FALSE)

atac <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scATAC-seq`@counts)
atac <- Seurat::SCTransform(atac, verbose = FALSE, return.only.var.genes = FALSE)

dat <- rbind(as.data.frame(rna[["SCT"]]$data), as.data.frame(atac[["SCT"]]$data))

ct <- tibble::tibble(cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

coexpr.log_corr_percentil <- stats::cor(percentile_expr(dat, ct, percentile_no = 10))

row.names(coexpr.log_corr_percentil) <- row.names(dat)
colnames(coexpr.log_corr_percentil) <- row.names(dat)

saveRDS(object = coexpr.log_corr_percentil, file = '~/workspace/1_conesalab/test_scMOSim/paper_plots/percentil_celltypes_few.rds')


## For Kmeans we will need this, where columns are genes and rows are the percentils, 10 per celltype
percentil <- acorde::percentile_expr(dat, ct, percentile_no = 9)
colnames(percentil) <- row.names(dat)









rna <- log1p(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)

atac <- log1p(sim$Group_2$Rep_1$`sim_scATAC-seq`@counts)

dat <- rbind(as.data.frame(rna), as.data.frame(atac))

ct <- tibble::tibble(cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

coexpr.log_corr_percentil <- acorde::percentile_cor(dat, ct, percentile_no = 9)

row.names(coexpr.log_corr_percentil) <- row.names(dat)
colnames(coexpr.log_corr_percentil) <- row.names(dat)

asocG2 <- sim$AssociationMatrices$AssociationMatrix_Group_2[sim$AssociationMatrices$AssociationMatrix_Group_2$Peak_cluster %in% c(1:11),]
asocG2 <- asocG2[asocG2$Gene_cluster %in% c(1:11),]

# Subset rows as genes and columns as peaks.

coexpr.log_corr_percentil <- coexpr.log_corr_percentil[row.names(coexpr.log_corr_percentil) %in% row.names(rna),]
coexpr.log_corr_percentil <- as.data.frame(coexpr.log_corr_percentil[, colnames(coexpr.log_corr_percentil) %in% row.names(atac)])
coexpr.log_corr_percentil$Gene_ID <- row.names(coexpr.log_corr_percentil)

# Subset the genes and peaks from the association matrix
trues <- coexpr.log_corr_percentil[row.names(coexpr.log_corr_percentil) %in% asocG2$Gene_ID, ]
trues <- coexpr.log_corr_percentil[,colnames(coexpr.log_corr_percentil) %in% asocG2$Peak_ID, ]





long <- as.data.frame(coexpr.log_corr_percentil) %>% pivot_longer(!Gene_ID, names_to = "Peak_ID", values_to = "corperc")

# This is super weird, only 32 rows are common between dataframes, it should be all 3000
trues <- merge(long, asocG2[, c("Gene_ID", "Peak_ID")], by = c("Gene_ID", "Peak_ID"))




