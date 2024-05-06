## Customized annotation of cell-types based on ScType

**1**. Retrieve a marker gene list from single-cell database, like **[PlantscRNAdb, ](http://ibi.zju.edu.cn/plantscrnadb/index.php)**[^1]
**[scPlantDB, ](https://biobigdata.nju.edu.cn/scplantdb)**[^2]
**[PsctH, ](http://jinlab.hzau.edu.cn/PsctH/)**[^3]
**[PCMDB.](https://www.tobaccodb.org/pcmdb/homePage)**[^4]
```bash
## marker gene list from PlantscRNAdb
$ head ara_markerdegree_v3.txt
Arabidopsis thaliana	Leaf	Spongy mesophyll cell	AT1G01140	Marker #2
Arabidopsis thaliana	Leaf	Spongy mesophyll cell	AT1G01320	Marker #2
Arabidopsis thaliana	Leaf	Spongy mesophyll cell	AT1G01470	Marker #2
Arabidopsis thaliana	Leaf	Spongy mesophyll cell	AT1G02660	Marker #2
Arabidopsis thaliana	Leaf	Spongy mesophyll cell	AT1G06680	Marker #2
...
## simplify the list to two columns (geneID, cell type)
$ awk 'BEGIN{FS="\t";OFS="\t"}{print $3,$4}' ara_markerdegree_v3.txt > ara_markerdegree_v3_2c.txt
$ head ara_markerdegree_v3_2c.txt
AT1G01120	Lower epidermal cell
AT1G01120	Upper epidermal cell
AT1G01120	Upper epidermal cell
AT1G01120	Upper epidermal cell
AT1G01120	Upper epidermal cell
...
## calculate the syntenic relationship between your genome with Arabidopsis thaliana
$ head At_Hp.syn
AT1G01020	Hper_g105835
AT1G01030	Hper_g105836
AT1G01040	Hper_g105837
AT1G01050	Hper_g105838
AT1G01060	Hper_g105840
...
## organize a ScTypeDB format file according to https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypeDB_full.xlsx in Excel.
```
**2**. Annotate using ScType
```bash
## start R
R
## load the module
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("~/res/sctype/script/optimized/gene_sets_prepare.R")
source("~/res/sctype/script/optimized/sctype_score_.R")
## load the customized database file
db_="~/res/sctype/ScTypeDB.xlsx"
## set the tissue type
tissue="Leaf"
gs_list=gene_sets_prepare(db_,tissue)
es.max = sctype_score(scRNAseqData = sce_single[["RNA"]]@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, NULL)
cL_resutls = do.call("rbind", lapply(unique(sce_single@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sce_single@meta.data[sce_single@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sce_single@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
## check scores
print(sctype_scores[,1:3])
## add customclassif attributes
sce_single@meta.data$customclassif = ""for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  sce_single@meta.data$customclassif[sce_single@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
dim_label<-DimPlot(sce_single, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')

```





### References

[^1]: PlantscRNAdb: A database for plant single-cell RNA analysis. Chen H, Yin X, Guo L, Yao J, Ding Y, Xu X, Liu L, Zhu QH, Chu Q, Fan L. Mol Plant. 2021 Jun 7;14(6):855-857. doi: 10.1016/j.molp.2021.05.002. Epub 2021 May 4. PMID: 33962062.
[^2]: ScPlantDB: a comprehensive database for exploring cell types and markers of plant cell atlasesZhaohui He, Yuting Luo, Xinkai Zhou, Tao Zhu, Yangming Lan, Dijun Chen, Nucleic Acids Research, Volume 52, Issue D1, 5 January 2024, Pages D1629–D1638, https://doi.org/10.1093/nar/gkad706
[^3]: Plant Single Cell Transcriptome Hub (PsctH): an integrated online tool to explore the plant single-cell transcriptome landscapeZhongping. Xu et.al. Plant Biotechnology Journal (2021), DOI: 10.1111/pbi.13725
[^4]: PCMDB: a curated and comprehensive resource of plant cell markers. Jin J, Lu P, Xu Y, Tao J, Li Z, Wang S, Yu S, Wang C, Xie X, Gao J, Chen Q, Wang L, Pu W, Cao P. Nucleic Acids Res. 2022 Jan 7;50(D1):D1448-D1455. doi: 10.1093/nar/gkab949. 
[^5]: Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Ianevski, A., Giri, A.K. & Aittokallio, T. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w
