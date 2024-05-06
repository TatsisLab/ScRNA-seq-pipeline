## Assessing the similarity of cell types between datasets by MetaNeighbor

**1**. Install MetaNeighbor[^1]
```bash
## installation
devtools::install_git('https://github.com/gillislab/MetaNeighbor')
library(MetaNeighbor)
library(SingleCellExperiment)
library(Seurat)
```

**2**. Assess the similarity of cell types
```bash
## read input from Seurat object
leaf<-readRDS("leaf.Rds")
bud<-readRDS("bud.Rds")
## transform Seurat object to sce object
leaf_sce<-as.SingleCellExperiment(leaf)
bud_sce<-as.SingleCellExperiment(bud)
## merge two objects to a list object
sce<-list(leaf_sce=leaf_sce,bud_sce=bud_sce)
sce$leaf_sce$"cell type" <- sce$leaf_sce$seurat_clusters
sce$bud_sce$"cell type" <- sce$bud_sce$seurat_clusters
fused_data<-mergeSCE(sce)
## identify variable genes
global_hvgs <- variableGenes(dat=fused_data,exp_labels=fused_data$study_id)
## generate aurocs
aurocs <- MetaNeighborUS(var_genes = global_hvgs,dat = fused_data,study_id = fused_d
ata$study_id,cell_type = fused_data$"cell type",fast_version=TRUE)
## generate cell clusters above 0.9 score
top_hits = topHits(aurocs,dat = fused_data,study_id = fused_data$study_id,cell_type
= fused_data$"cell type",threshold = 0.9)
```
[^1]: Characterizing the replicability of cell types defined by single cell RNA-sequencing data using MetaNeighbor. Crow, M., Paul, A., Ballouz, S. et al. Nat Commun 9, 884 (2018). https://doi.org/10.1038/s41467-018-03282-0
