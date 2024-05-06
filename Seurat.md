## Processing of single-cell RNA-seq data with Seurat

**1**. Load the matrix from Cellranger or alevin [^1],
```bash
## load the R package
library(Seurat)
## set the dir and load the martrix
dir="~/res/cellranger/outs/filtered_feature_bc_matrix"
list.files(dir)
counts<-Read10X(data.dir = dir)
## filter the cells with less than 500 genes expressed
sce<-CreateSeuratObject(counts,min.cells = 1,min.features = 500) 
```
**2**. Filter and normalize the data
```bash
## filter the cells with more than 3000 genes expressed
sce<-subset(sce,subset= nFeature_RNA<3000)
## normalize the data
sce<-NormalizeData(sce,normalization.method = "LogNormalize",scale.factor = 10000)
## find variable features
sce<-FindVariableFeatures(sce,selection.method = "vst",nfeatures = 2000)
## scale the data
all.genes<-rownames(sce)
sce<-ScaleData(sce,features=all.genes)
## optional: scTransform to replace the classic nomalization method
BiocManager::install("glmGamPoi")
sce <- SCTransform(sce, method = "glmGamPoi", verbose = FALSE)
```

**3**. PCA dimensionality reduction
```bash
sce_PCA<-RunPCA(sce,features=VariableFeatures(object=sce))
## determine the dimention
pct <- sce_PCA [["pca"]]@stdev / sum( sce_PCA [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
elbow<-ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + geom_text() + 
geom_vline(xintercept = 90, color = "grey") + geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +  
theme_bw()
pdf("elbow.pdf")
print(elbow)
dev.off()
## check the elbow plot to determine the dimension
```
**4**. Determine the resolution
```bash
sce_findNei<-FindNeighbors(sce_PCA,dims=1:16)
## find best resolution
library(clustree)
sce_findclu<-FindClusters(sce_findNei,resolution=c(seq(0.1,1,.1)))
clustree<-clustree(sce_findclu@meta.data, prefix = "RNA_snn_res.")
sce_findClu_025<-FindClusters(sce_findNei,resolution=0.25)
## UMAP reduction optional: you can also use t-SNE
sce_umap<-RunUMAP(sce_findClu_025,dims=1:16)
## check the cell numbers for each cluster
table(Idents(sce_umap))
```
**5**. Remove the potential doublets
```bash
## function to find doublets
Find_doublet <- function(data){sweep.res.list <- paramSweep_v3(data, PCs = 1:16, sct = FALSE)
                               sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
                               bcmvn <- find.pK(sweep.stats)
                               nExp_poi <- round(0.05*ncol(data))
                               p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
                               data <- doubletFinder_v3(data, PCs = 1:16, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
                               colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
                               return(data)
                               }
## use the reduction data 
sce_doublet<-Find_doublet(sce_umap)
## check the ration of singlets and doublets
doublet_pdf<-DimPlot(sce_doublet,reduction="umap",group.by="doublet_info",cols=c("black","gold"))
## remove the doublets
sce_single<-subset(sce_doublet,subset=doublet_info=="Singlet")
```
**6**. Find the most variable genes as markers
```bash
##find marker genes
sce.markers <- FindAllMarkers(sce_single, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(sce.markers,"marker.csv",col.names=T,row.names=F,sep=",", )
## subset a spercific cluster like hyper
sce_hyper<-subset(sce_single,idents=11)
```
**7**. Visulization
```bash
## classic Dimplot
umap<-DimPlot(sce_umap, reduction = "umap",label=T)
## customized Dimplot plot
library(RColorBrewer)
my_cols <-c('9'='#884caa','0'='#5f2b8f','10'='#b06ec5','1'='#ffb8ff','11'='#FF4848','3'='#7ed6a8','6'='#209658','8'='#ffba00','4'='#ff7b00','2'='#73cede','7'='#0566ad','5'='#00a8e1')
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
umap<-DimPlot(leaf_single, reduction = "umap",label=F,cols=alpha(my_cols2,0.8),pt.size=0.5)

## classic Featureplot
gene<-FeaturePlot(sce_umap, features = c("geneID"))
## custoimzed Featureplot
library(scCustomize)
gene<-FeaturePlot_scCustom(seurat_object=sce_umap, features = "geneID",colors_use=viridis_plasma_dark_high,pt.size=2,alpha_exp = 0.8,alpha_na_exp= 0.8)
## Genes density plot
gene_density<- Plot_Density_Joint_Only(seurat_object = sce_umap, features = c("GeneID1","GeneID2","GeneID3"),viridis_palette = "plasma")

## classic Dot plot
dot<-DotPlot(sce_umap,features=c("GeneID"))+RotatedAxis()
## Customized Dot plot
library(scCustomize)
dot<-DotPlot_scCustom(seurat_object=sce_umap,features=c("GeneID"),colors_use=viridis_plasma_dark_high)+RotatedAxis()
```

### References

[^1]: Alevin-fry unlocks rapid, accurate and memory-frugal quantification of single-cell RNA-seq data. He, D., Zakeri, M., Sarkar, H., Soneson, C., Srivastava, A., and Patro, R. Nat Methods 19, 316–322 (2022). https://doi.org/10.1038/s41592-022-01408-3
