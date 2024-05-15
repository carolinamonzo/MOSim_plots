devtools::load_all("/home/cmonzo/workspace/MOSim")
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

sim <- readRDS("~/workspace/mosim_paper/sim_6cells11clus8800_scMOSim_2groups_000.rds")

# Normalize using Seurat's method
#dat <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
#dat <- Seurat::SCTransform(dat, verbose = FALSE, return.only.var.genes = FALSE)

# Create cell-tocelltype ID table
ct <- tibble::tibble(Cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

rna <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
rna <- Seurat::SCTransform(rna, verbose = FALSE, return.only.var.genes = FALSE)

atac <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scATAC-seq`@counts)
atac <- Seurat::SCTransform(atac, verbose = FALSE, return.only.var.genes = FALSE)

dat <- rbind(as.data.frame(rna[["SCT"]]$data), as.data.frame(atac[["SCT"]]$data))

ct <- tibble::tibble(cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

coexpr.log_corr_percentil <- acorde::percentile_cor(dat, ct, percentile_no = 10)

row.names(coexpr.log_corr_percentil) <- row.names(dat)
colnames(coexpr.log_corr_percentil) <- row.names(dat)

saveRDS(object = coexpr.log_corr_percentil, file = '~/workspace/mosim_paper/percentil_celltypes_many.rds')

