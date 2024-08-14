# ScRNA-seq pipeline

This is a pipeline that handles the single-cell RNA-seq data of **St. John's wort** (*Hypericum perforatum*)[^1]🌿. 

### Working Background
Single-cell RNA sequencing  has emerged as a revolutionary technology in genomics, enabling the analysis of the genomic and transcriptomic landscapes at an unprecedented resolution. This technology is particularly popular in human and mouse studies due to its ability to unravel complex cellular heterogeneities and molecular dynamics that are crucial for understanding various biological processes and diseases.

In recent years, the application of single-cell sequencing has extended to the field of plant biology, enhancing our understanding of plant development, physiology, and interaction with the environment. Initially, this technology was predominantly applied to model plants, such as *Arabidopsis thaliana*, where comprehensive genomic resources and easier cell isolation techniques facilitated early and extensive studies.

However, the expansion of single-cell sequencing to non-model plants, especially medicinal plants, faces several challenges. First, the baseline biological research on these plants is relatively underdeveloped. Many non-model plants have complex genomes and unique physiological properties that complicate cell isolation and preparation for sequencing. Furthermore, after sequencing, the annotation of genomic data in non-model plants presents additional hurdles due to the lack of experimentally proved marker genes and limited genetic information. This gap significantly hinders the potential for detailed molecular insights that could drive pharmacological innovations and conservation efforts.

### Annotation of cell types
Annotation of cell types in a plant single-cell atlas is usually based on single or a set of marker genes. However, such markers are absent for *H. perforatum*. 

There is, however, an extensive database of genes from mainstream model plants used as cell-type markers. Identification of functional orthologs between *A. thaliana* and *H. perforatum* by sequence alignment and homology are prone to false positives. Instead, synteny is one of the most reliable criteria for establishing the functional relationships between orthologous genes or genomic regions across different species.

To annotate the cell clusters with an unbiased method, we performed a fully automated annotation pipeline comprising: 
- Compiling an extensive marker-gene dataset for *A. thaliana* from database like **[PlantscRNAdb](http://ibi.zju.edu.cn/plantscrnadb/index.php)**[^2]
* Correlating the orthologous marker genes between the two species and generated a marker-gene dataset for *H. perforatum* through syntenic analysis between *A. thaliana* and *H. perforatum* by **[MCscan](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version))**[^3]
+ Performing the **[ScType](https://github.com/IanevskiAleksandr/sc-type)**[^4] platform, an automated cell-type annotation based on a pool of markers
![Scheme](/Scheme_annotation.png)
> Scheme of rapid and automated cell-type annotation in non-model plants. Adopted from ScType.[^4]
### Dependencies

**1**. From the treatment of single-cell RNA-seq raw data to plots like UMAP figures were most done in R

```bash
## It's better to create a specific enviroment with conda
conda create -n hyper -c conda-forge r-base=4.1.2 r-rgeos
## activate the conda enviroment
conda activate hyper
## start R console
R
```

**2**. Install the R packages

```bash
## packages for processing single-cell raw data
install.packages('Seurat')
install.packages("clustree")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
## packages for plotting
install.packages('ggplot2')
install.packages('ggsci')
Install.packages('RColorBrewer')
install.packages('scCustomize')
```

**3**. Install the MCscan python version

```bash
pip install jcvi
## download and install the dependencies of jcvi
## http://last.cbrc.jp/
export PATH="$PATH:/yourdir/last/src"
## https://github.com/tanghaibao/jcvi-bin/blob/master/bin/scip
export PATH="$PATH:/yourdir/Scip"
## more details please see https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)#dependencies
```



> [!TIP]
> If you had any questions during processing please post issues or send me the Email (wusong@cemps.ac.cn)(songwu2024@163.com).

> [!NOTE]
> If you found this flow useful please cite our work of hyperforin biosynthesis 😃[^1] and all applied softwares.[^2][^3][^4][^5][^6]  

### References
[^1]: Single cell RNA sequencing facilitates the elucidation of the complete biosynthesis of the antidepressant hyperforin in St. John's wort. Song Wu, Ana Luisa Malaco Morotti, Jun Yang, Ertao Wang, Evangelos C. Tatsis Mol Plant 2024 Aug 12:S1674-2052(24)00258-2; doi: 10.1016/j.molp.2024.08.003
[^2]: PlantscRNAdb: A database for plant single-cell RNA analysis. Chen H, Yin X, Guo L, Yao J, Ding Y, Xu X, Liu L, Zhu QH, Chu Q, Fan L. Mol Plant. 2021 Jun 7;14(6):855-857. doi: 10.1016/j.molp.2021.05.002. Epub 2021 May 4. PMID: 33962062.
[^3]: Synteny and collinearity in plant genomes. Tang H, Bowers JE, Wang X, Ming R, Alam M, Paterson AH. Science. 2008 Apr 25;320(5875):486-8. doi: 10.1126/science.1153917. PMID: 18436778.
[^4]: Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Ianevski, A., Giri, A.K. & Aittokallio, T. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w
[^5]: Clustering trees: a visualization for evaluating clusterings at multiple resolutions. Zappia L, Oshlack A. Gigascience. 2018;7. DOI:gigascience/giy083.
[^6]: Characterizing the replicability of cell types defined by single cell RNA-sequencing data using MetaNeighbor. Crow, M., Paul, A., Ballouz, S. et al. Nat Commun 9, 884 (2018). https://doi.org/10.1038/s41467-018-03282-0

