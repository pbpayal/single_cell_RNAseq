Raw fastq files from single cell
sequencing were processed using the CellRanger (10X Genomics Cell
Ranger 7.0.1) pipeline. The reads were aligned to the
CellRanger_GRCh38_2020-A (https://support.10xgenomics.com/
single-cell-gene-expression/software/release-notes/build#GRCh38
_2020A) human reference genome. The filtered feature-barcode
matrices produced from the 10X pipeline were used for further
downstream clustering and analysis. Unsupervised cell clustering
was performed by Seurat (4.1.1) (9) in R. (4.1.0). For each sample,
the filtered feature-barcode matrix produced from the CellRanger
pipeline was read and converted to Seurat object using “Read10X”
and “CreateSeuratObject” functions, respectively. Criteria were used
to identify gel bead-in-emulsions (GEMs) likely containing mRNAs
derived only from a single cell, where GEMs were retained if more
than 200 genes and less than 2500 genes were detected. Genes were
retained for downstream analysis if they were detected in minimum
of three cells. For each sample, the ‘‘NormalizeData’’ was performed
to normalize for gene expression followed by ‘‘FindVariableGene’’
function to identify a subset of genes, in this case limited to 2000
genes, that exhibit high cell-to-cell variation. The 8 samples in
controland experiment groups were integrated using the
‘‘FindIntergrationAnchors’’ and ‘‘IntegrateData’’ functions with the dimension parameter set to 20. Next, the integrated dataset was
scaled, and PCA was performed on this dataset. The first 20 principal
components (PCs) were used to perform UMAP to place similar cells
together in low-dimensional space. Then the shared nearest-neighbor
graph (SNN) was constructed using the “FindNeighbors” function,
using the first 20 PCs. The “FindClusters” function, which
implements a graph-based Louvain algorithm was applied to
identify 11 distinct clusters. The resolution was set to 0.2 for this
dataset. The clusters were annotated with the top markers using
“FindAllMarkers” function. Classic cell markers were selected as
reference after extensive literature and database search. (example -
Human Protein Atlas (10) and PanglaoDB (11)). The following genes
were used as markers for cell type identity - MS4A1(B cells), CD8B,
CD8A (T cells), LYZ (Dendritic cells), GNLY (NK cells), IL7R,
BCL11B, MAL (T memory cells/T cells), SLC25A37 (Erythroid-like
and erythroid precursor cells), GZMB, PTGDS, PPP1R14B
(Plasmacytoid dendritic cells), HBB, HBA1, HBA2 (Red blood
cells) and PPBP (Platelets). The differentially expressed genes
between the cells from control and the children with MtD were
calculated by “FindMarkers”.

Paper DOI - https://doi.org/10.3389/fimmu.2023.1142634
