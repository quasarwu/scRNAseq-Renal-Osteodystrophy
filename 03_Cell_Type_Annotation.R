# 1. Identify all markers for each cluster
cell_markers <- FindAllMarkers(sc_merged, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(cell_markers, file = "CellMarkers_merge.csv", row.names = FALSE)

# 2. Define canonical markers
markers <- c(
  "Cd3d", "Cd3e", "Ccl5",
  "Cd19", "Cd79a", "Ms4a1",
  "Sdc1", "Prdm1", "Jchain",
  "Alpl", "Sp7", "Bglap", "Ibsp",
  "Hba-a1", "Hbb-bs", "Hbb-bt",
  "Mcpt8", "Prss34", "Cpa3",
  "Cd34", "Adgrg1", "Prom1",
  "Siglech", "Irf8", "Bst2",
  "Elane", "Mpo", "Prtn3",
  "Adgre1", "Cd68", "Apoe",
  "Ly6c2", "Csf1r", "Ccr2",
  "Acp5", "Ctsk", "Calcr", "Atp6v0d2",
  "S100a8", "S100a9", "Ly6g", "Cxcr2"
)

# 3. Visualize markers and current clusters
Idents(sc_merged) <- "RNA_snn_res.0.6"
DotPlot(sc_merged, features = unique(markers), scale = TRUE) + RotatedAxis()
DimPlot(sc_merged, label = TRUE, reduction = "umap")

# 4. Annotate clusters using a mapped dictionary
cluster_annotations <- c(
  "0" = "Neutrophils", "1" = "Neutrophils", "2" = "Neutrophils", "4" = "Neutrophils", "13" = "Neutrophils", "14" = "Neutrophils",
  "3" = "B cells", "7" = "B cells", "12" = "B cells", "21" = "B cells",
  "5" = "Macrophages",
  "6" = "T cells", "19" = "T cells", "26" = "T cells",
  "8" = "NMPs", "11" = "NMPs", "25" = "NMPs",
  "9" = "Monocytes",
  "10" = "Erythroblasts", "15" = "Erythroblasts", "17" = "Erythroblasts", "24" = "Erythroblasts",
  "16" = "HSC",
  "18" = "Osteoblasts",
  "20" = "pDCs",
  "22" = "Basophils",
  "23" = "Plasma Cells"
)

sc_merged <- RenameIdents(sc_merged, cluster_annotations)
sc_merged$celltype <- Idents(sc_merged)

# 5. Review annotations and save the finalized object
print(table(Idents(sc_merged)))
saveRDS(sc_merged, "04_sc_merged_annotated.rds")
