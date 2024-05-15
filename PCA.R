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

numcells <- "RNA_G2R1"

# Normalize using Seurat's method
dat <- Seurat::CreateSeuratObject(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)
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



## Now K-means starting from the log2 normalized
asocG2 <- sim$AssociationMatrices$AssociationMatrix_Group_2[sim$AssociationMatrices$AssociationMatrix_Group_2$Gene_cluster %in% c(1:11),]
data <- as.data.frame(dat[["SCT"]]$data)
data <- data[rownames(data) %in% asocG2$Gene_ID,]

coexpr.scaled.celltype <- data

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

coexpr.scaled.celltype <- calculate_mean_per_list_df(coexpr.scaled.celltype, cell_types)

max_itr <-  100
n_clust  <-  11  ## number of cluster 

## Remove NAs

coexpr.scaled.celltype[is.na(coexpr.scaled.celltype)] <- 0

kmeans_out  <- kmeans(as.matrix(coexpr.scaled.celltype),n_clust,iter.max = max_itr)
# With cluster info from the kmeans

data_with_cust_info <- as.data.frame(coexpr.scaled.celltype) %>% 
  mutate(clust = paste("clust_", kmeans_out$cluster,sep = ""))

mean_expression <- data_with_cust_info %>%
  group_by(clust) %>%
  summarise_all(mean)

## visualise  each cluster 

df <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 


data_with_cust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(data_with_cust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = df, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 4)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/2A_", numcells, "_kmeans.pdf"))


## Now let's check the cluster assignment from kmeans compared to the forced cluster
# in the association matrix

data <- as.data.frame(dat[["SCT"]]$data)
data <- calculate_mean_per_list_df(data, cell_types)
data <- as.data.frame(data[rownames(data) %in% asocG2$Gene_ID,])

subAsoc <- asocG2[, c("Gene_ID", "Gene_cluster")]
subAsoc <- subAsoc %>% mutate(ori_clust = paste0("ori_clust_", asocG2$Gene_cluster))
subAsoc <- subAsoc[, c("Gene_ID", "ori_clust")]

data$Gene_ID <- row.names(data)
data_with_oricust_info <- merge(data, subAsoc, by = "Gene_ID")
data_with_oricust_info$Gene_ID <- NULL

mean_expression <- data_with_oricust_info %>%
  group_by(ori_clust) %>%
  summarise_all(mean)

## visualise  each cluster 

df <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 

data_with_oricust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(data_with_oricust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = df, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~ori_clust, ncol = 4)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/original_2A_", numcells, "_kmeans.pdf"))
