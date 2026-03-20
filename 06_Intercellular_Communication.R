library(Seurat)
library(CellChat)
library(patchwork)

# 1. Efficiently subset target cell populations
target_cells <- c("Transitional Mac", "Inflammatory Mac", "Lipid-Associated Mac", 
                  "Repair Mac", "Monocytes", "Osteoblasts")
sc_cellchat <- subset(sc_merged, celltype %in% target_cells)

# Prepare normalized data for CellChat
sc_cellchat <- JoinLayers(sc_cellchat)
sc_cellchat <- NormalizeData(sc_cellchat)
sc_cellchat$labels <- sc_cellchat$celltype

# 2. Define a streamlined function for the CellChat pipeline
run_cellchat_pipeline <- function(seurat_obj, group_name) {
  cat("Processing group:", group_name, "\n")
  
  # Subset by group condition
  group_data <- subset(seurat_obj, subset = group == group_name)
  
  # Initialize CellChat object
  data_input <- GetAssayData(group_data[["RNA"]], layer = "data")
  meta <- group_data@meta.data
  
  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.mouse
  
  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  # Compute communication probabilities
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Compute network centrality directly
  cellchat <- netAnalysis_computeCentrality(cellchat)
  
  # Extract and save dataframes
  df_net <- subsetCommunication(cellchat)
  df_pathway <- subsetCommunication(cellchat, slot.name = "netP")
  
  write.csv(df_net, file = paste0("CellChat_Net_", group_name, ".csv"), row.names = FALSE)
  write.csv(df_pathway, file = paste0("CellChat_Pathway_", group_name, ".csv"), row.names = FALSE)
  
  # Save individual group object
  saveRDS(cellchat, file = paste0("07_CellChat_", group_name, ".rds"))
  
  return(cellchat)
}

# 3. Execute the pipeline for both groups
cc_Con <- run_cellchat_pipeline(sc_cellchat, "Con")
cc_ROD  <- run_cellchat_pipeline(sc_cellchat, "ROD")

# 4. Merge the datasets for comparative analysis
cc_list <- list(Con = cc_Con, ROD = cc_ROD)
cellchat_merged <- mergeCellChat(cc_list, add.names = names(cc_list), cell.prefix = TRUE)

# Save the final comparative object
saveRDS(cellchat_merged, file = "08_CellChat_Merged_Comparison.rds")
