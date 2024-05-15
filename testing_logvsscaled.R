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

#sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells11clus8800_scMOSim_2groups_000.rds")

## Plotting big sample
sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_few_2602.rds")

results <- scOmicResults(sim)
settings <- scOmicSettings(sim)


## Let's try Angeles' visualization.

st <- "G2R1"

# create cell-to-cell-type ID table
ct <- tibble::tibble(Cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))


#### Doing the acorde scaling

coexpr.scaled <- acorde::scale_isoforms(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts, 
                                        isoform_col = NULL)

coexpr.scaled[is.na(coexpr.scaled)] <- 0

# Check some profiles: 
df <- zoo(coexpr.scaled[, 2:5])
plot(df)

df <- zoo(coexpr.scaled[, c(2,12, 25,37)])
plot(df)

# Make a correlation of all cells
coexpr.scaled_corr <- coexpr.scaled %>% select(-transcript) %>% cor(method = "spearman")
# Visualise the correlations between the first 5 samples
coexpr.scaled_corr[1:5, 1:5]

library(corrr)
# While the strongest correlations are within cell-type, some celltypes are similar
# To others, mostly because they are all lymphocytes
rplot(coexpr.scaled_corr) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Make PCA
library(factoextra)


pca_sc <- prcomp(t(coexpr.scaled %>% select(-transcript)))
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + geom_point() + theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

#ggsave(width = 3, height = 2, filename = paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/", numcells, "PCA_PC1_PC2.pdf"))

# Since the third component is also important to visualize it
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC2, y = PC3, colour = cell_type)) + geom_point() + theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

#ggsave(width = 3, height = 2, filename = paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/", numcells, "PCA_PC2_PC3.pdf"))

fviz_eig(pca_sc)
#ggsave(width = 3, height = 2, filename = paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/", numcells, "PCA_scree.pdf"))




##### Doing the stuff on log1 data


# PCA on the log1(1p) ori data
pca_raw<- prcomp(log(t(as.data.frame(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)) + 1))
pca_raw$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + geom_point() + theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

dat <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
dat <- Seurat::SCTransform(dat, verbose = FALSE, return.only.var.genes = FALSE)

pca_sc <- prcomp(t(dat[["SCT"]]$data))
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + geom_point() + theme_classic(base_size = 10) +
  scale_color_conesa(palette = "main")

library(factoextra)
fviz_eig(pca_raw)
fviz_eig(pca_sc)


# MAKE HEATMAPS OF CORRELATION BETWEEN CELLS
coexpr.log_corr <- dat[["SCT"]]$data %>% cor(method = "spearman")
# Visualise the correlations between the first 5 cells
coexpr.log_corr[1:5, 1:5]

library(corrr)
# While the strongest correlations are within cell-type, some celltypes are similar
# To others, mostly because they are all lymphocytes
rplot(coexpr.log_corr) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

coexpr.log_corr_kendall <- dat[["SCT"]]$data %>% cor(method = "kendall")
# Visualise the correlations between the first 5 cells
coexpr.log_corr_kendall[1:5, 1:5]
# While the strongest correlations are within cell-type, some celltypes are similar
# To others, mostly because they are all lymphocytes
rplot(coexpr.log_corr_kendall) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))