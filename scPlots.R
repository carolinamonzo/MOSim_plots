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

sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_few_2602.rds")

## Plotting big sample
#sim <- readRDS("~/workspace/1_conesalab/test_scMOSim/paper_plots/sim_scMOSim_2groups_000.rds")

results <- scOmicResults(sim)
settings <- scOmicSettings(sim)


## Let's try Angeles' visualization.

st <- "G3R1"

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

# Check covariation between cells
coexpr.scaled %>% ggplot(aes(Treg_cell1, Treg_cell2)) + geom_point() + geom_abline(colour = "brown")

# create cell-to-cell-type ID table
ct <- tibble::tibble(Cell = colnames(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts),
                     cell_type = rep(names(sim$cellTypes), times = lengths(sim$cellTypes)))

#Make PCA
pca_sc <- prcomp(t(coexpr.scaled %>% select(-transcript)))
pca_sc$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + geom_point()

# PCA on the unscaled data === This looks good

pca_raw<- prcomp(log(t(as.data.frame(sim$Group_2$Rep_1$`sim_scRNA-seq`@counts)) + 1))
pca_raw$x %>% as_tibble(rownames = "Cell") %>% full_join(ct, by = "Cell") %>% 
  ggplot(aes(x = PC1, y = PC2, colour = cell_type)) + geom_point()

library(factoextra)

#Check how much variability each PC explains
fviz_eig(pca_raw)


# Wide to long for plotting
coexpr.long <- coexpr.scaled %>% pivot_longer(cols = `CD16 Mono_cell1`:`Memory B_cell371`,
                                              names_to = "Cell",
                                              values_to = "cts")

# Join the celltype info
coexpr.long <- full_join(coexpr.long, ct, by = "Cell")

# Plot distribution
coexpr.long %>% ggplot(aes(cts, colour = cell_type)) + geom_freqpoly(binwidth = 1)

# Extract DE genes to visualize expression trends
candidate_genes <- settings$AssociationMatrix_Group_2 %>%
  filter(Gene_DE %in% c("Up", "Down")) %>%
  pull(Gene_ID) %>% unique()

# Filter our long table to retain genes of interest and summarize counts per celltype
coexpr.mean <- coexpr.long %>% filter(transcript %in% candidate_genes) %>%
  group_by(transcript, cell_type) %>% summarise(mean_cts = mean(cts),
                                                nrep = n()) %>% ungroup()

# Plot the expression pattern
coexpr.mean %>% ggplot(aes(cell_type, mean_cts)) +
  geom_line(aes(group = transcript), alpha = 0.3)

##### Do the same with only up genes
# Extract DE genes to visualize expression trends
UPcandidate_genes <- settings$AssociationMatrix_Group_2 %>%
  filter(Gene_DE %in% c("Up")) %>%
  pull(Gene_ID) %>% unique()

# Filter our long table to retain genes of interest and summarize counts per celltype
UPcoexpr.mean <- coexpr.long %>% filter(transcript %in% UPcandidate_genes) %>%
  group_by(transcript, cell_type) %>% summarise(mean_cts = mean(cts),
                                                nrep = n()) %>% ungroup()

# Plot the expression pattern
UPcoexpr.mean %>% ggplot(aes(cell_type, mean_cts)) +
  geom_line(aes(group = transcript), alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "brown", linetype = "dashed")

##### Do the same with only down genes
# Extract DE genes to visualize expression trends
DOWNcandidate_genes <- settings$AssociationMatrix_Group_2 %>%
  filter(Gene_DE %in% c("Down")) %>%
  pull(Gene_ID) %>% unique()

# Filter our long table to retain genes of interest and summarize counts per celltype
DOWNcoexpr.mean <- coexpr.long %>% filter(transcript %in% DOWNcandidate_genes) %>%
  group_by(transcript, cell_type) %>% summarise(mean_cts = mean(cts),
                                                nrep = n()) %>% ungroup()

# Plot the expression pattern
DOWNcoexpr.mean %>% ggplot(aes(cell_type, mean_cts)) +
  geom_line(aes(group = transcript), alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "brown", linetype = "dashed")

## Working on clustering
hclust_matrix <- coexpr.scaled %>% select(-transcript) %>% as.matrix()
rownames(hclust_matrix) <- coexpr.scaled$transcript

hclust_matrix <- hclust_matrix[candidate_genes, ]
hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

hclust_matrix[is.na(hclust_matrix)] <- 0

gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 10, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

## Clustering using hierarchical clustering of genes.

gene_cluster <- cutree(gene_hclust, k = 8) %>% 
  # turn the named vector into a tibble
  enframe()

colnames(gene_cluster) <- c("transcript", "cluster")

head(gene_cluster)

trans_cts_cluster <- coexpr.mean %>% 
  inner_join(gene_cluster, by = "transcript")

head(trans_cts_cluster)

trans_cts_cluster %>% 
  ggplot(aes(cell_type, mean_cts)) +
  geom_line(aes(group = transcript), alpha = 0.3) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +
  facet_grid(cols = vars(cluster))


library(ComplexHeatmap)
#Heatmap(hclust_matrix, show_row_names = FALSE)

## Extract the assigned pattern (remove rows with regulator or anything not gene):
pat <- settings$AssociationMatrix_Group_2$Gene_cluster[
  !is.na(settings$AssociationMatrix_Group_2$Gene_cluster)]

# Remove genes with no expression pattern
pat <- pat[ !pat == 0]


## Here is where we have the problem, pat is a list

# compute average-by-cell type cluster patterns
cluster_patterns <- map( pat,
                        ~acorde::calculate_cluster_profile(coexpr.scaled,
                                                           isoform_ids = .,
                                                           id_table = ct,
                                                           isoform_col = "transcript"))


# plot patterns
library(cowplot)

theme_set(theme_cowplot())


pattern_plots <- map(cluster_patterns,
                     acorde::plot_cluster_profile,
                     ct_labels = c("CD4_TEM", "cDC", "Memory_B", "Treg"))

plot_grid(plotlist = pattern_plots, 
          labels = NULL, 
          ncol = 4)
