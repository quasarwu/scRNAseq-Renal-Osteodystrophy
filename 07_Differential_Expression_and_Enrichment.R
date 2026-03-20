library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)

# 1. Subset and Prepare Data
Inflammatory <- subset(Macrophage, stim == "Inflammatory Mac")
Idents(Inflammatory) <- "group"

# 2. Differential Expression Analysis
MYDEG <- FindMarkers(Inflammatory, 
                     ident.1 = "ROD", 
                     ident.2 = "Con", 
                     verbose = FALSE, 
                     test.use = "wilcox", 
                     min.pct = 0.1)

MYDEG$Gene <- rownames(MYDEG)
MYDEG <- subset(MYDEG, pct.1 > 0 & pct.2 > 0)

logFC_cutoff <- 1
p_cutoff <- 0.05

MYDEG$change <- "NOT"
MYDEG$change[MYDEG$p_val < p_cutoff & MYDEG$avg_log2FC > logFC_cutoff] <- "UP"
MYDEG$change[MYDEG$p_val < p_cutoff & MYDEG$avg_log2FC < -logFC_cutoff] <- "DOWN"

write.csv(MYDEG, file = 'MYDEG_filtered.csv', row.names = TRUE)

# 3. Minimalist Volcano Plot
min_nonzero_p <- min(MYDEG$p_val[MYDEG$p_val > 0])
MYDEG$p_val[MYDEG$p_val == 0] <- min_nonzero_p * 0.1
MYDEG$neg_log10_pvalue <- -log10(MYDEG$p_val)


p_volcano <- ggplot(MYDEG, aes(x = avg_log2FC, y = neg_log10_pvalue, color = change)) +
  geom_point(alpha = 0.8, size = 2, stroke = 0) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "grey50") +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)")

ggsave("Volcano_Plot_Minimal.pdf", p_volcano, width = 6, height = 5)

# 4. Enrichment Analysis Preparation
deg_genes_up <- MYDEG$Gene[MYDEG$change == "UP"]
deg_entrez <- bitr(deg_genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# 5. KEGG Analysis & Minimalist Plot
kegg_enrich <- enrichKEGG(gene = deg_entrez$ENTREZID, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH")
kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
write.csv(as.data.frame(kegg_enrich), "KEGG_results.csv", row.names = FALSE)

kegg_selected <- c(
  "Oxidative phosphorylation - Mus musculus (house mouse)",
  "HIF-1 signaling pathway - Mus musculus (house mouse)", 
  "Carbon metabolism - Mus musculus (house mouse)", 
  "Glycolysis / Gluconeogenesis - Mus musculus (house mouse)", 
  "Pentose phosphate pathway - Mus musculus (house mouse)",
  "Chemokine signaling pathway - Mus musculus (house mouse)"
)

kegg_df <- as.data.frame(kegg_enrich) %>%
  filter(Description %in% kegg_selected) %>%
  mutate(log10P = -log10(pvalue))

p_kegg <- ggplot(kegg_df, aes(x = reorder(Description, log10P), y = log10P))
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "-Log10(P-value)")

ggsave("KEGG_Selected_Minimal.pdf", plot = p_kegg, height = 4, width = 8)

# 6. GO Analysis & Minimalist Plot
go_enrich <- enrichGO(gene = deg_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "all", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
write.csv(as.data.frame(go_enrich), file = 'GO_results.csv', row.names = FALSE)

go_selected <- c(
  "cellular response to oxygen levels", 
  "cellular response to decreased oxygen levels",
  "cellular respiration", 
  "ATP metabolic process",
  "ATP biosynthetic process",
  "response to oxidative stress"
)

go_df <- as.data.frame(go_enrich) %>%
  filter(Description %in% go_selected) %>%
  mutate(log10P = -log10(pvalue))

p_go <- ggplot(go_df, aes(x = reorder(Description, log10P), y = log10P))
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "-Log10(P-value)")

ggsave("GO_Selected_Minimal.pdf", plot = p_go, height = 4, width = 8)
