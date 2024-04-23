# **ScRNA-seq pipeline**

This is a pipeline that handles the single-cell RNA-seq data of **St. John's wort** (*Hypericum perforatum*)[^1]🌿. 

### Working Background

### Annotation of cell types
Annotation of cell types in a plant single-cell atlas is usually based on single or a set of marker genes. However, such markers are absent for *H. perforatum*. 

There is, however, an extensive database of genes from mainstream model plants used as cell-type markers. Identification of functional orthologs between *A. thaliana* and *H. perforatum* by sequence alignment and homology are prone to false positives. Instead, synteny is one of the most reliable criteria for establishing the functional relationships between orthologous genes or genomic regions across different species.

To annotate the cell clusters with an unbiased method, we performed a fully automated annotation pipeline comprising: 
- Compiling an extensive marker-gene dataset for *A. thaliana* from database like **[PlantscRNAdb](http://ibi.zju.edu.cn/plantscrnadb/index.php)**[^2]
* Correlating the orthologous marker genes between the two species and generated a marker-gene dataset for *H. perforatum* through syntenic analysis between *A. thaliana* and *H. perforatum* by **[MCScan](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version))**[^3]
+ Performing the **[ScType](https://github.com/IanevskiAleksandr/sc-type)**[^4] platform, an automated cell-type annotation based on a pool of markers

### Dependencies

说明如何安装和设置你的项目。例如：

```bash
git clone https://github.com/yourusername/yourprojectname.git
cd yourprojectname
pip install -r requirements.txt
```
### References
[^1]: Single-cell RNA-seq based elucidation of the antidepressant hyperforin biosynthesis de novo in St. John’s wort. Song Wu, Ana Luisa Malaco Morotti, Jun Yang, Ertao Wang, Evangelos C. Tatsis bioRxiv 2024.01.24.577018; doi: https://doi.org/10.1101/2024.01.24.577018
[^2]: PlantscRNAdb: A database for plant single-cell RNA analysis. Chen H, Yin X, Guo L, Yao J, Ding Y, Xu X, Liu L, Zhu QH, Chu Q, Fan L. Mol Plant. 2021 Jun 7;14(6):855-857. doi: 10.1016/j.molp.2021.05.002. Epub 2021 May 4. PMID: 33962062.
[^3]: Synteny and collinearity in plant genomes. Tang H, Bowers JE, Wang X, Ming R, Alam M, Paterson AH. Science. 2008 Apr 25;320(5875):486-8. doi: 10.1126/science.1153917. PMID: 18436778.
[^4]: Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data. Ianevski, A., Giri, A.K. & Aittokallio, T. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w
> [!NOTE]
> If you found the mistakes of the scripts or had the ERROR report during processing please post issues or send me the Email (wusong@cemps.ac.cn).
> If you found this flow useful please cite our work of hyperforin biosynthesis.
