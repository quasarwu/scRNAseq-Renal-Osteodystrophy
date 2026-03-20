library(Seurat)
library(harmony)
library(ggplot2)

# 1. Merge objects and prepare layers for Seurat v5
sc_merged <- merge(sc_list_clean[[1]], 
                   y = sc_list_clean[-1], 
                   add.cell.ids = names(sc_list_clean), 
                   project = "ROD_Project")
sc_merged <- JoinLayers(sc_merged)

# 2. Variable features, scaling, PCA and Harmony integration
sc_merged <- sc_merged %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunHarmony(group.by.vars = "orig.ident", plot_convergence = FALSE)

# 3. Quantitative calculation of optimal PC dimension
pct <- sc_merged[["pca"]]@stdev / sum(sc_merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

cond1 <- which(cumu > 90 & pct < 5)[1]
cond2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
valid_pcs <- c(cond1, cond2)


# 4. Visualize and save the Elbow Plot with the threshold line
p_elbow <- ElbowPlot(sc_merged)$data %>% 
  ggplot() +
  geom_point(aes(x = dims, y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw() + 
  labs(title = paste0("Elbow plot: Quantitative approach (Optimal PC = ", pc.use, ")"))

# 5. Downstream clustering using the calculated pc.use
sc_merged <- sc_merged %>%
  RunUMAP(reduction = "harmony", dims = 1:pc.use, verbose = FALSE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:pc.use, verbose = FALSE) %>%
  FindClusters(resolution = 0.6, verbose = FALSE)

# 6. Final visualization
p_umap <- DimPlot(sc_merged, reduction = "umap", label = TRUE, repel = TRUE)
print(p_umap)
ggsave("UMAP_Clusters.pdf", plot = p_umap, width = 7, height = 6)

# Save the final processed object
saveRDS(sc_merged, "03_sc_merged_clustered.rds")
