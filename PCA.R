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

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_2groups_000.rds")

numcells <- "many_ATAC_G2R1"

# Normalize using Seurat's method
dat <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scATAC-seq`@counts)
dat <- Seurat::SCTransform(dat, verbose = FALSE, return.only.var.genes = FALSE)

# Create cell-tocelltype ID table
ct <- tibble::tibble(Cell = colnames(dat[["SCT"]]$data),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

library(factoextra)


pca_sc <- prcomp(t(dat[["SCT"]]$data))
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + geom_point() + theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

ggsave(width = 3, height = 2, filename = paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/", numcells, "PCA_PC1_PC2.pdf"))

# Since the third component is also important to visualize it
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC2, y = PC3, colour = cell_type)) + geom_point() + theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

ggsave(width = 3, height = 2, filename = paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/", numcells, "PCA_PC2_PC3.pdf"))

fviz_eig(pca_sc)
ggsave(width = 3, height = 2, filename = paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/", numcells, "PCA_scree.pdf"))


