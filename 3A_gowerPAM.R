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
library(cowplot)

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

set.seed(000)

### read the simulated dataset

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells11clus8800_scMOSim_2groups_000.rds")

settings <- scOmicSettings(sim, TF = TRUE)
cell_types <- sim$cellTypes

### SCALE USING ACORDE

rna <- acorde::scale_isoforms(sim$Group_1$Rep_1$`sim_scRNA-seq`@counts, 
                                            isoform_col = NULL)

rna[is.na(rna)] <- 0

ct <- tibble::tibble(Cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

### PLOT ORIGINAL CLUSTERS

## Extract clusters list to plot originals
clusters_list <- settings$AssociationMatrix_Group_1[settings$AssociationMatrix_Group_1$Gene_cluster %in% c(1:11), c("Gene_ID", "Gene_cluster")]
clusters_list <- split(clusters_list$Gene_ID, clusters_list$Gene_cluster)


# compute average-by-cell type cluster patterns
cluster_patterns <- map(clusters_list,
                        ~acorde::calculate_cluster_profile(rna,
                                                           isoform_ids = .,
                                                           id_table = ct,
                                                           isoform_col = "transcript"))


theme_set(theme_cowplot())

pattern_plots <- map(cluster_patterns,
                     plot_cluster_profile,
                     ct_labels = c("CD16 Mono", "NK", "Treg", "cDC", "CD4 TEM", "Memory B"))

plot_grid(plotlist = pattern_plots, 
          labels = NULL, 
          ncol = 4)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/Supp_OriClusters_acordeScaled.pdf"))

####### Clustering genes with acorde scaling, manhattan and kmedoids


rna <- as.data.frame(rna)
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


##### Plot scree of the medoids

sil_width <- c(NA)
 
for(i in 2:15){
  
  pam_fit <- pam(rna_gower_dist,
                 diss = TRUE,
                 k = i)
  
  sil_width[i] <- pam_fit$silinfo$avg.width
  
}

sil_width <- data.frame(
  Position = 1:length(sil_width), 
  Value = sil_width)

ggplot(sil_width, aes(x = Position, y = Value)) + geom_point() + geom_line() +
  labs(x = "Number of medoids", y = "Pam average fit") + theme_classic()
ggsave(width = 3, height = 2, filename = "~/workspace/1_conesalab/test_scMOSim/paper_plots/Kmedoids_scree.pdf")


