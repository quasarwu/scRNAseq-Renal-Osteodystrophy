# scRNAseq-Renal-Osteodystrophy
Customized R pipelines for single-cell RNA sequencing analysis, characterizing macrophage heterogeneity and metabolic reprogramming in renal osteodystrophy.

## Pipeline Structure
The analysis pipeline is organized into sequentially numbered R scripts to ensure complete computational reproducibility:

* `01_QC_and_Doublet_Removal.R`: Raw data loading, stringent quality control metrics, and doublet identification using *DoubletFinder*.
* `02_Integration_and_Clustering.R`: Data normalization, exclusion of dissociation-induced stress genes, batch effect correction via *Harmony*, and unsupervised clustering.
* `03_Cell_Type_Annotation.R`: Marker gene identification and assignment of distinct cell lineages.
* `04_Macrophage_Subclustering.R`: High-resolution sub-clustering of the macrophage population and functional enrichment analysis.
* `05_Pseudotime_Trajectory.R`: Cell state transition modeling using *Monocle*.
* `06_Intercellular_Communication.R`: Ligand-receptor network inference using *CellChat*.
* `07_Differential_Expression_and_Enrichment.R`: Statistical identification of DEGs and functional pathway enrichment (GO/KEGG).

## System Requirements and Dependencies
The analytical scripts were developed and executed in **R (version [4.3.1])**. 

**Core R packages required:**
* `Seurat` (v[5.0.1])
* `Harmony` (v[1.0.3])
* `DoubletFinder` (v[2.0.3])
* `clusterProfiler` (v[4.8.1])
* `CellChat` (v[1.6.1])
* `Monocle` (v[2.28.0] or Monocle3)
* `AUCell` / `GSVA` (for pathway scoring)
