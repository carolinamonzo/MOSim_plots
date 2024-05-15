## correlations

setwd("~/workspace/1_conesalab/MOSim")
devtools::load_all()
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

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells11clus8800_scMOSim_2groups_000.rds")

numcells <- "few_RNA_G1R1"

# Normalize using Seurat's method
dat <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
dat <- Seurat::SCTransform(dat, verbose = FALSE, return.only.var.genes = TRUE)

datATAC <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scATAC-seq`@counts)
datATAC <- Seurat::SCTransform(datATAC, verbose = FALSE, return.only.var.genes = TRUE)

settings <- scOmicSettings(sim)

# Create cell-tocelltype ID table
ct <- tibble::tibble(Cell = colnames(dat[["SCT"]]$data),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

# Subset the association matrix to only variable genes
associationMatrix <- as.data.frame(settings$AssociationMatrix_Group_2[settings$AssociationMatrix_Group_2$Gene_ID %in% rownames(dat[['SCT']]$data),])
# And variable features
associationMatrix <- as.data.frame(associationMatrix[associationMatrix$Peak_ID %in% rownames(datATAC[['SCT']]$data),])

# Check out the variable associated features

testRNA1 <- dat[["SCT"]]$data[rownames(dat[["SCT"]]$data) %in% c("MOB3C"), ]



df_plot <- data.frame(x = 1:length(testRNA1), y = testRNA1)
ggplot(df_plot, aes(x = x, y = y)) + geom_line()



