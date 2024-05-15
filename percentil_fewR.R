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
ct <- tibble::tibble(Cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

rna <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
rna <- Seurat::SCTransform(rna, verbose = FALSE, return.only.var.genes = FALSE)

atac <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scATAC-seq`@counts)
atac <- Seurat::SCTransform(atac, verbose = FALSE, return.only.var.genes = FALSE)

dat <- rbind(as.data.frame(rna[["SCT"]]$data), as.data.frame(atac[["SCT"]]$data))

ct <- tibble::tibble(cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

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

coexpr.log_corr_percentil <- stats::cor(percentile_expr(dat, ct, percentile_no = 10))

row.names(coexpr.log_corr_percentil) <- row.names(dat)
colnames(coexpr.log_corr_percentil) <- row.names(dat)

saveRDS(object = coexpr.log_corr_percentil, file = '~/workspace/1_conesalab/test_scMOSim/paper_plots/percentil_celltypes_few.rds')



