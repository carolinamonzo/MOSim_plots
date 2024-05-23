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
  library(cluster)
})
library(pbmcMultiome.SeuratData)

set.seed(000)

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells11clus8800_scMOSim_2groups_000.rds")

ct <- tibble::tibble(Cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

rna <- log1p(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)

atac <- log1p(sim$Group_2$Rep_1$`sim_scATAC-seq`@counts)

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

asocG2 <- sim$AssociationMatrices$AssociationMatrix_Group_2[sim$AssociationMatrices$AssociationMatrix_Group_2$Gene_cluster %in% c(1:11),]
means_celltype_rna <- as.data.frame(means_celltype_rna[rownames(means_celltype_rna) %in% asocG2$Gene_ID,])

# Cluster our 11 clusters of interest
rna_gower_dist <- cluster::daisy(means_celltype_rna, metric = "gower")
rna_pam <- pam(rna_gower_dist, diss = TRUE, k = 11)

pam_cluster <- as.data.frame(rna_pam$clustering)
colnames(pam_cluster) <- c("pam_cluster")

pam_cluster <- as.data.frame(pam_cluster) %>% 
  mutate(clust = paste("clust_", pam_cluster,sep = ""))
pam_cluster$pam_cluster <- NULL

rna_with_clust_info <- merge(means_celltype_rna, pam_cluster, by = 0)

rownames(rna_with_clust_info) <- rna_with_clust_info$Row.names
rna_with_clust_info$Row.names <- NULL


mean_expression <- rna_with_clust_info %>%
  group_by(clust) %>%
  summarise_all(mean)

## visualise  each cluster 

df <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 


rna_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(rna_with_clust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = df, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 4)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/3A_gowerPAM.pdf"))

## Plot this nicer for the supplementary
# sil_width <- c(NA)
# 
# for(i in 2:15){
#   
#   pam_fit <- pam(rna_gower_dist,
#                  diss = TRUE,
#                  k = i)
#   
#   sil_width[i] <- pam_fit$silinfo$avg.width
#   
# }
# 
# plot(1:15, sil_width,
#      xlab = "Number of clusters",
#      ylab = "Silhouette Width")
# lines(1:15, sil_width)



####### Doing the same but scaling with acorde, looks very good
rna <- as.data.frame(acorde::scale_isoforms(sim$Group_1$Rep_1$`sim_scRNA-seq`@counts, 
                              isoform_col = NULL))

rna[is.na(rna)] <- 0

cell_types <- sim$cellTypes
rownames(rna) <- rna$transcript
rna$transcript <- NULL

means_celltype_rna <- calculate_mean_per_list_df(rna, cell_types)

asocG2 <- sim$AssociationMatrices$AssociationMatrix_Group_2[sim$AssociationMatrices$AssociationMatrix_Group_2$Gene_cluster %in% c(1:11),]
means_celltype_rna <- as.data.frame(means_celltype_rna[rownames(means_celltype_rna) %in% asocG2$Gene_ID,])

# Cluster our 11 clusters of interest
rna_gower_dist <- cluster::daisy(means_celltype_rna, metric = "gower")
rna_pam <- pam(rna_gower_dist, diss = TRUE, k = 11)

pam_cluster <- as.data.frame(rna_pam$clustering)
colnames(pam_cluster) <- c("pam_cluster")

pam_cluster <- as.data.frame(pam_cluster) %>% 
  mutate(clust = paste("clust_", pam_cluster,sep = ""))
pam_cluster$pam_cluster <- NULL

rna_with_clust_info <- merge(means_celltype_rna, pam_cluster, by = 0)

rownames(rna_with_clust_info) <- rna_with_clust_info$Row.names
rna_with_clust_info$Row.names <- NULL


mean_expression <- rna_with_clust_info %>%
  group_by(clust) %>%
  summarise_all(mean)

## visualise  each cluster 

df <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 


rna_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(ncol(rna_with_clust_info))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = df, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 4)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/3A_G1_gowerPAM_acordeScaled.pdf"))


## Checking the original clustering
means_celltype_rna_ori <- calculate_mean_per_list_df(rna, cell_types)

asocG2 <- sim$AssociationMatrices$AssociationMatrix_Group_2[sim$AssociationMatrices$AssociationMatrix_Group_2$Gene_cluster %in% c(1:11),]
means_celltype_rna_ori <- as.data.frame(means_celltype_rna_ori[rownames(means_celltype_rna_ori) %in% asocG2$Gene_ID,])

subAsoc <- asocG2[, c("Gene_ID", "Gene_cluster")]
subAsoc <- subAsoc %>% mutate(ori_clust = paste0("ori_clust_", asocG2$Gene_cluster))
subAsoc <- subAsoc[, c("Gene_ID", "ori_clust")]

means_celltype_rna_ori$Gene_ID <- row.names(means_celltype_rna_ori)
means_celltype_rna_ori <- merge(means_celltype_rna_ori, subAsoc, by = "Gene_ID")
means_celltype_rna_ori$Gene_ID <- NULL

mean_expression <- means_celltype_rna_ori %>%
  group_by(ori_clust) %>%
  summarise_all(mean)

## visualise  each cluster 

df <- mean_expression %>% gather(key = "variable", value = "value", -c(1)) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) 

means_celltype_rna_ori %>% 
  gather(key = "variable" , value = "value", -c(ncol(means_celltype_rna_ori))) %>%  ###the last column is the cluster
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num)) +   
  geom_line(alpha = 1, color = "grey") + 
  geom_line(data = df, aes(x = variable, y= value, group = row_num), linewidth = 1, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~ori_clust, ncol = 4)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/3A_oriCluster_acordeScaling.pdf"))
