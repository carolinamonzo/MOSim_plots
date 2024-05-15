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

set.seed(000)
# ## Plots for paper

#sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_few_2602.rds")

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_6cells11clus8800_scMOSim_2groups_000.rds")

numcells <- "many_RNA_G2R1"
results <- scOmicResults(sim)
settings <- scOmicSettings(sim, TF = TRUE)

## SC plots for paper
# Figure 2A representation of k-means clusters for scRNAseq data forced to 
# follow a co-expression pattern.

# Extract group one
data <- as.data.frame(results$Group_2$Rep_1$`sim_scRNA-seq`@counts)
# Keep only genes in co-expression clusters
asocG1 <- settings$AssociationMatrix_Group_2[settings$AssociationMatrix_Group_2$Gene_cluster %in% c(1:9),]
data <- data[asocG1$Gene_ID,]
# Scale data and make means per celltype
coexpr.scaled <- acorde::scale_isoforms(data, isoform_col = NULL)
#coexpr.scaled[is.na(coexpr.scaled)] <- 0

# Remove the transcripts names so we can work on the matrix
coexpr.scaled.celltype <- coexpr.scaled %>% select(-transcript)

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


### kmeans
max_itr <-  50
n_clust  <-  9  ## number of cluster 
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
  geom_line(data = df, aes(x = variable, y= value, group = row_num), size = 1, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 3)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/many_2A_", numcells, "_kmeans.pdf"))


#####################


## Now let's check the cluster assignment from kmeans compared to the forced cluster
# in the association matrix

data_with_oricust_info <- as.data.frame(coexpr.scaled.celltype) %>% 
  mutate(ori_clust = paste0("ori_clust_", asocG1$Gene_cluster))

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
  geom_line(data = df, aes(x = variable, y= value, group = row_num), size = 1, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~ori_clust, ncol = 3)
ggsave(paste0("~/workspace/1_conesalab/test_scMOSim/paper_plots/original_many_2A_", numcells, "_kmeans.pdf"))



#Extract the association matrix count the number of times the cluster is the same as it should be
#and the number of times it is not.
