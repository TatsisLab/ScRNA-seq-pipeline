# **ScRNA-seq pipeline**

This is a pipeline that handles the single-cell RNA-seq data of St. John's wort (*Hypericum perforatum*)🌿. 

## Working Background

Annotation of cell types in a plant single-cell atlas is based on single or a set of marker genes.
201 However, such markers are absent for H. perforatum. There is, however, an extensive database of
genes from mainstream model plants used as cell-type markers59-61 202 . Identification of functional
203 orthologs between A. thaliana and H. perforatum by sequence alignment and homology are prone
204 to false positives. Instead, synteny is one of the most reliable criteria for establishing the functional
relationships between orthologous genes or genomic regions across different species62 205 . To annotate
206 the cell clusters with an unbiased method, we performed a fully automated annotation pipeline
comprising: (1) compiling an extensive marker-gene dataset for A. thaliana59 207 ; (2) through syntenic
208 analysis between A. thaliana and H. perforatum, we correlated the orthologous marker genes
209 between the two species and generated a marker-gene dataset for H. perforatum; and (3) using the
ScType63 210 platform, an automated cell-type annotation was performed based on a pool of markers

## 安装

说明如何安装和设置你的项目。例如：

```bash
git clone https://github.com/yourusername/yourprojectname.git
cd yourprojectname
pip install -r requirements.txt
