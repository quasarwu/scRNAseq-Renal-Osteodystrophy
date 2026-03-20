library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# 1. Subset Macrophages and prepare features
Macrophage <- subset(sc_merged, celltype == "Macrophages")
print(table(Macrophage$group))

scale_genes <- rownames(Macrophage)

# 2. Feature selection, scaling, PCA, and Integration
Macrophage <- Macrophage %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500, verbose = FALSE) %>% 
  ScaleData(features = scale_genes, verbose = FALSE) %>% 
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident", plot_convergence = FALSE)

# 3. Dynamic PC dimension calculation (with safety fallback)
pct <- Macrophage[["pca"]]@stdev / sum(Macrophage[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

cond1 <- which(cumu > 90 & pct < 5)[1]
cond2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
valid_pcs <- c(cond1, cond2)
valid_pcs <- valid_pcs[!is.na(valid_pcs)]

pc.use <- ifelse(length(valid_pcs) > 0, min(valid_pcs), 15)
cat("Optimal PC dimension for Macrophages:", pc.use, "\n")

# 4. Dimensionality reduction and Clustering
Macrophage <- Macrophage %>%
  RunUMAP(reduction = "harmony", dims = 1:pc.use, verbose = FALSE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:pc.use, verbose = FALSE) %>% 
  FindClusters(resolution = 0.2, verbose = FALSE)

# 5. Marker identification and verification
mac_markers <- FindAllMarkers(Macrophage, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 6. Assign functional subtypes using a dictionary
mac_annotations <- c(
  "0" = "Transitional Mac",
  "1" = "Transitional Mac",
  "2" = "Inflammatory Mac",
  "3" = "Inflammatory Mac",
  "6" = "Inflammatory Mac",
  "4" = "Lipid-Associated Mac",
  "5" = "Repair Mac",
  "7" = "Repair Mac"
)

Macrophage <- RenameIdents(Macrophage, mac_annotations)
Macrophage$stim <- Idents(Macrophage)

# 7. Final visualization and export
print(table(Macrophage$stim))

p_mac_umap <- DimPlot(Macrophage, reduction = "umap", label = TRUE, repel = TRUE) + 
  ggtitle("Macrophage Subpopulations")

saveRDS(Macrophage, "05_Macrophage_Subclustered.rds")
