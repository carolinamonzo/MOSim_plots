setwd("~/workspace/1_conesalab/MOSim")
devtools::load_all()
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

######## ANGELES SCALING ##########
# Make a correlation of all cells
coexpr.scaled_corr <- coexpr.scaled %>% select(-transcript) %>% cor(method = "spearman")
# Visualise the correlations between the first 5 samples
coexpr.scaled_corr[1:5, 1:5]


# While the strongest correlations are within cell-type, some celltypes are similar
# To others, mostly because they are all lymphocytes
rplot(coexpr.scaled_corr) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Make a correlation of all cells
coexpr.scaled_corr <- coexpr.scaled %>% select(-transcript) %>% cor(method = "kendall")
# Visualise the correlations between the first 5 samples
coexpr.scaled_corr[1:5, 1:5]

# While the strongest correlations are within cell-type, some celltypes are similar
# To others, mostly because they are all lymphocytes
rplot(coexpr.scaled_corr) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######## LOG1 ##########

dat <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
dat <- Seurat::SCTransform(dat, verbose = FALSE, return.only.var.genes = FALSE)

coexpr.log_corr <- dat[["SCT"]]$data %>% cor(method = "spearman")
# Visualise the correlations between the first 5 cells
coexpr.log_corr[1:5, 1:5]


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