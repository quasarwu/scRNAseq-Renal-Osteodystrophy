library(Seurat)
library(harmony)
library(ggplot2)
library(clustree)

# Merge the list of Seurat objects into a single object
sc_merged <- merge(sc_list_clean[[1]], 
                   y = sc_list_clean[2:length(sc_list_clean)], 
                   add.cell.ids = names(sc_list_clean), 
                   project = "ROD_Project")

# Required for Seurat v5 to perform proper downstream scaling and PCA
sc_merged <- JoinLayers(sc_merged)

sc_merged <- FindVariableFeatures(sc_merged, selection.method = "vst", nfeatures = 2500, verbose = FALSE)
hvgs <- VariableFeatures(sc_merged)

sc_merged <- sc_merged %>% 
  ScaleData(verbose = FALSE) %>%  
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident", plot_convergence = FALSE)

pct <- sc_merged[["pca"]]@stdev / sum(sc_merged[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

cond1 <- which(cumu > 90 & pct < 5)[1]
cond2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
valid_pcs <- c(cond1, cond2)
valid_pcs <- valid_pcs[!is.na(valid_pcs)]

if (length(valid_pcs) > 0) {
  pc.use <- min(valid_pcs)
} else {
  pc.use <- 20
  warning("Strict threshold not met, defaulting to PC = 20")
}

cat("Optimal PC dimension calculated:", pc.use, "\n")

p_elbow <- ElbowPlot(sc_merged)$data %>% 
  ggplot() +
  geom_point(aes(x = dims, y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw() + 
  labs(title = paste0("Elbow plot: Quantitative approach (Optimal PC = ", pc.use, ")"))
print(p_elbow)

sc_merged <- sc_merged %>%
  RunUMAP(reduction = "harmony", dims = 1:pc.use, verbose = FALSE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:pc.use, verbose = FALSE)

sc_merged <- FindClusters(object = sc_merged, resolution = seq(0.1, 1.2, 0.1), verbose = FALSE)

# clustree(sc_merged@meta.data, prefix = "RNA_snn_res.")

Idents(sc_merged) <- "RNA_snn_res.0.6"

DimPlot(sc_merged, reduction = "umap", label = TRUE)
