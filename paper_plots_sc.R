suppressPackageStartupMessages({
  library(tidyverse)
  library(scater)
  library(WGCNA)
  library(Seurat)
  library(RColorConesa)
  library(acorde)
})

set.seed(123)
# ## Plots for paper

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_few_2602.rds")
results <- scOmicResults(sim)
settings <- scOmicSettings(sim)

## SC plots for paper
# Figure 2A representation of k-means clusters for scRNAseq data forced to 
# follow a co-expression pattern.

# Extract group one
data <- as.data.frame(results$Group_1$Rep_1$`sim_scRNA-seq`@counts)
# Keep only genes in co-expression clusters
asocG1 <- settings$AssociationMatrix_Group_1[settings$AssociationMatrix_Group_1$Gene_cluster %in% c(1:6),]
data <- data[asocG1$Gene_ID,]
# Scale data and make means per celltype
coexpr.scaled <- acorde::scale_isoforms(data, isoform_col = NULL)
coexpr.scaled[is.na(coexpr.scaled)] <- 0

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

coexpr.scaled.celltype <- calculate_mean_per_list_df(coexpr.scaled.celltype, cell_types)

#rownames(coexpr.scaled.celltype) <- coexpr.scaled$transcript

### kmeans
max_itr <-  50
n_clust  <-  8  ## number of cluster 
set.seed(123) ## reproduce the cluster 
kmeans_out  <- kmeans(coexpr.scaled.celltype,n_clust,iter.max = max_itr)
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
  geom_line(data = df, aes(x = variable, y= value, group = row_num), size = 1.5, color = "red3") +
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 45 , vjust = 0.4)) +
  facet_wrap(~clust, ncol = 3)
ggsave("~/workspace/1_conesalab/test_scMOSim/paper_plots/2A_kmeans.pdf")

## Now let's check the cluster assignment from kmeans compared to the forced cluster
# in the association matrix

# First we need to make Angeles' visualization



