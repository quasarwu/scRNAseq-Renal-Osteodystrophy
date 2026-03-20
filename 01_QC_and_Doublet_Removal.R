samples <- c("con1", "con2", "con3", "ROD1", "ROD2", "ROD3")
sample_names <- c("CON1", "CON2", "CON3", "ROD1", "ROD2", "ROD3")
groups <- c(rep("Con", 3), rep("ROD", 3))

sc_list <- list()
for (i in seq_along(samples)) {
  counts <- Read10X(data.dir = paste0(samples[i], "/")) 
  sc_list[[sample_names[i]]] <- CreateSeuratObject(
    counts = counts, 
    project = sample_names[i], 
    min.cells = 3, 
    min.features = 200
  )
  sc_list[[sample_names[i]]]$group <- groups[i]
}

run_qc_and_doublets <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  
  genes_in_data <- rownames(obj)
  
  obj <- subset(obj, subset = nFeature_RNA > 250 & percent.mt < 10 & 
                  nCount_RNA > 500 & nFeature_RNA < 6000 & 
                  nCount_RNA < 50000)
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  mt_genes <- grep("^mt-", genes_in_data, value = TRUE)
  ribo_genes <- grep("^Rp[sl]", genes_in_data, value = TRUE)
  exclude_genes <- unique(c(globin_features, mt_genes, ribo_genes))
  
  VariableFeatures(obj) <- setdiff(VariableFeatures(obj), exclude_genes)
  
  obj <- ScaleData(obj, verbose = FALSE) %>% 
    RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:20, verbose = FALSE) %>% 
    FindNeighbors(dims = 1:20, verbose = FALSE) %>% 
    FindClusters(resolution = 0.5, verbose = FALSE)
  
  sweep.res.list <- paramSweep(obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk_best <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  doublet_rate <- (ncol(obj) / 1000) * 0.008
  nExp_poi <- round(doublet_rate * ncol(obj))
  homotypic.prop <- modelHomotypic(obj$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  obj <- doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = pk_best, 
                       nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  obj <- subset(obj, cells = rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Singlet"])
  
  obj@meta.data <- obj@meta.data[, !grepl("^DF\\.classifications|^pANN|^RNA_snn|^seurat_clusters", colnames(obj@meta.data))]
  
  cat("Sample", unique(obj$orig.ident), "QC & Doublet removal complete. Remaining cells:", ncol(obj), "\n")
  return(obj)
}

sc_list_clean <- lapply(sc_list, run_qc_and_doublets)

saveRDS(sc_list_clean, "02_sc_list_clean.rds")
