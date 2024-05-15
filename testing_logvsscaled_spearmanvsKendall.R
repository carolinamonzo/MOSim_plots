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

## PERCENTIL CORRELATION
percentile_expr <- function(data, id_table, percentile_no = 10, isoform_col = NULL){
  
  # handle rownames and data type
  if(is.null(isoform_col) == TRUE){
    data <- data %>% as.data.frame %>% rownames_to_column("transcript")
  }
  
  # split cell IDs by cell type labels
  cells_split <- split(id_table$cell, id_table$cell_type)
  
  # check that percentile_no is between 4 and 100
  if(percentile_no < 4 | percentile_no > 100){
    warning("percentile_no is not between 4 (quantiles) and 100 (percentiles).")
  }
  
  # generate step for probabilities in quantile function
  step <- 1/percentile_no
  # calculate percentiles by cell type for each transcript
  percentile_list <- map(cells_split, ~(select(data, all_of(.)) %>%
                                          apply(1, quantile, seq(0, 1, step)) %>% as.data.frame))
  percentiles <- bind_rows(percentile_list)
  colnames(percentiles) <- data %>% select(isoform_col) %>% unlist
  
  return(percentiles)
}

ct <- tibble::tibble(cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))
coexpr.log_corr_percentil <- stats::cor(percentile_expr(coexpr.log_corr, ct, percentile_no = 10))
coexpr.log_corr_percentil[1:5, 1:5]
# While the strongest correlations are within cell-type, some celltypes are similar
# To others, mostly because they are all lymphocytes
rplot(coexpr.log_corr_percentil) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
