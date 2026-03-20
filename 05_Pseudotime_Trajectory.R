library(monocle)
library(Seurat)
library(ggplot2)

# 1. Safely extract sparse count matrix (Crucial for memory management in Seurat V5)
seurat_obj <- Macrophage
expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

p_data <- seurat_obj@meta.data
f_data <- data.frame(gene_short_name = rownames(seurat_obj), row.names = rownames(seurat_obj))

pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)

# 2. Initialize Monocle CellDataSet
cds <- newCellDataSet(
  expr_matrix,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 3. Select ordering genes using Seurat's FindAllMarkers
marker_genes <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0.5)
write.csv(marker_genes, file = "Monocle_Input_Markers.csv", row.names = FALSE)

expressed_genes <- unique(marker_genes$gene)
cds <- cds[expressed_genes, ]

# 4. Perform differential gene test to find significant trajectory genes
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~stim", cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))


# 5. Dimensionality reduction and trajectory construction
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

saveRDS(cds, file = "06_Macrophage_Monocle2.rds")

plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE)

plot_cell_trajectory(cds, color_by = "stim", size = 1, show_backbone = TRUE) +

plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) 
