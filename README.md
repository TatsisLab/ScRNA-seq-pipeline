# **ScRNA-seq pipeline**

This is a pipeline that handles the single-cell RNA-seq data of **St. John's wort** (*Hypericum perforatum*)🌿. 

## Working Background

## Annotation of cell types
Annotation of cell types in a plant single-cell atlas is usually based on single or a set of marker genes. However, such markers are absent for *H. perforatum*. There is, however, an extensive database of genes from mainstream model plants used as cell-type markers. Identification of functional orthologs between *A. thaliana* and *H. perforatum* by sequence alignment and homology are prone to false positives. Instead, synteny is one of the most reliable criteria for establishing the functional
relationships between orthologous genes or genomic regions across different species. To annotate the cell clusters with an unbiased method, we performed a fully automated annotation pipeline comprising: 
- Compiling an extensive marker-gene dataset for *A. thaliana* from database like **[PlantscRNAdb](http://ibi.zju.edu.cn/plantscrnadb/index.php)**
* Correlating the orthologous marker genes between the two species and generated a marker-gene dataset for *H. perforatum* through syntenic analysis between *A. thaliana* and *H. perforatum* by **[MCScan](https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version))**
+ Performing the **[ScType](https://github.com/IanevskiAleksandr/sc-type)** platform, an automated cell-type annotation based on a pool of markers

## 安装

说明如何安装和设置你的项目。例如：

```bash
git clone https://github.com/yourusername/yourprojectname.git
cd yourprojectname
pip install -r requirements.txt
