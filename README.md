# ClusterMatch aligns single-cell RNA-seqencing data at the multi-scale cluster level via stable matching
## 1.	Introduction to ClusterMatch
### 1.1	Overview of ClusterMatch
Unsupervised clustering of single-cell RNA sequencing (scRNA-seq) data holds the promise of known and novel cell type characterization in various biological and clinical contexts. However, multiple sources of variability in the data poses challenges to deal with high-dimensional and high-noise characteristics, batch effects, and intrinsic multi-scale clustering resolutions. Here, we present ClusterMatch, a stable match optimization model to account the statistical uncertainty by aligning scRNA-seq data at the cluster level. In one hand, ClusterMatch leverages the mutual correspondence by canonical correlation analysis (CCA) and multi-scale Louvain clustering algorithms to identify cluster with optimized resolutions. In the other hand it utilizes stable matching framework to align scRNA-seq data in the latent space while maintaining interpretability with overlapped marker gene set. Through extensive experiments, we demonstrate the efficacy of ClusterMatch in data integration, cell type annotation, and cross-species/timepoint alignment scenarios. Our results show ClusterMatch's ability to utilize both global and local information of scRNA-seq data, sets the appropriate resolution of multi-scale clustering, and offers interpretability by utilizing marker genes.
![Figure1](https://github.com/BateerCoder/ClusterMatch/assets/150581842/ec50430c-58f1-4498-8da7-c7a246f4dae5)
### 1.2 Installing R package
```
library(devtools)
devtools::install_github("BateerCoder/ClusterMatch")
```
## 2.	Main functions of ClusterMatch
### 2.1 Load data
```
dendritic_batch1 <- read.csv(".../data/dendritic/batch1.csv", row.names = 1)
dendritic_batch2 <- read.csv(".../data/dendritic/batch2.csv", row.names = 1)
dendritic_celltype <- read.csv(".../data/dendritic/celltype.csv")
human_dLGN <- read.csv(".../data/dLGN/human_dLGN.csv", row.names = 1)
macaque_dLGN <- read.csv(".../data/dLGN/macaque_dLGN.csv", row.names = 1)
human_celltype <- read.csv(".../data/dLGN/human_celltype.csv")
rownames(human_celltype) <- human_celltype$cells
macaque_celltype <- read.csv(".../data/dLGN/macaque_celltype.csv")
rownames(macaque_celltype) <- macaque_celltype$cells
```
The format of the input file is as follows
1. The row names: gene symbols.
2. The column names: cell IDs.
3. Other place: the expression values (counts or TPM) for a gene in a cell.

As an example, we will use two batches of human dendritic cell datasets in sections 2.2, 2.3, and 2.4. Batch 1 consists of 288 cells with three annotated cell types: CD141 DC, double-negative, and plasmacytoid DC. Batch 2 also contains 288 cells with three annotated cell types: CD1C DC, double-negative, and plasmacytoid DC.

We will also use human and macaque dorsal lateral geniculate nucleus (dLGN) datasets as an example in section 2.5. The human dLGN dataset consists of 952 cells, while the macaque dLGN dataset consists of 1,723 cells. The annotations have been refined through classification and manual review, resulting in four cell types: koniocellular (K), magnocellular and parvocellular projection neurons (MP), and two GABAergic cells (GABA1 and GABA2).
### 2.2 Find optimal resolutions
```
dendritic_res <- ClusterMatch_resolution(dendritic_batch1, dendritic_batch2, ref.norm = FALSE, que.norm = FALSE)
```
The optimal resolutions of batch1 and batch2
```
dendritic_res$D1_ref_res
[1] 1.1
dendritic_res$D2_que_res
[1] 1.4
```
### 2.3	Match scRNA-seq data at the cluster level
```
dendritic_matching <- ClusterMatch_matching(dendritic_batch1, dendritic_batch2, ref.res = dendritic_res$D1_ref_res, que.res = dendritic_res$D2_que_res, ref.norm = FALSE, que.norm = FALSE, random_PCC = 1.3)
```
Matching matrix
```
dendritic_matching$Matching_matrix
     D1_1 D1_0 D1_2
D2_1    0    2    0
D2_2    0    0    2
D2_0    0    0    0
```
### 2.4	Integrate scRNA-seq datasets
```
dendritic_integration <- ClusterMatch_integration(dendritic_batch1, dendritic_batch2, ref.res = dendritic_res$D1_ref_res,
                                                  que.res = dendritic_res$D2_que_res, ref.norm = FALSE, que.norm = FALSE, random_PCC = 1.3, distance_diff = 3, distance_same = 1)
```
UMAP visualization colored by the annotated batches

```
umap_df <- ClusterMatch_UMAP(embedding = dendritic_integration$cell_embedding, celltype = dendritic_celltype)

batch_colour=c("#E64540","#3F81BB")
celltype_colour=c("#E64136","#5F78A3","#EDA6C3","#96C561")
library(ggplot2)
ggplot(umap_df,aes(X1,X2,color=batch)) + 
  scale_color_manual(values = batch_colour)+
  geom_point() + theme_bw() +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size = 20)) +
  labs(x="UMAP_1",y="UMAP_2",
       title = "ClusterMatch")
```
<img width="725" alt="截屏2023-12-07 08 51 54" src="https://github.com/BateerCoder/ClusterMatch/assets/150581842/c128fb10-1349-4b2d-8ac6-208ac2e3d39c">

UMAP visualization colored by the annotated cell types
```
ggplot(umap_df,aes(X1,X2,color=label)) + 
  scale_color_manual(values = celltype_colour)+
  geom_point() + theme_bw() +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),text = element_text(size = 20)) +
  labs(x="UMAP_1",y="UMAP_2",
       title = "ClusterMatch")
```
<img width="725" alt="截屏2023-12-07 08 51 36" src="https://github.com/BateerCoder/ClusterMatch/assets/150581842/9c54dc3e-99a1-4159-8495-4fcb5db58ab5">

### 2.5	Annotate the query data based on the reference data labels
```
dLGN_annotation <- ClusterMatch_annotation(human_dLGN, macaque_dLGN, human_celltype, que.res=2, ref.norm = FALSE, que.norm = FALSE)
```
```
table(dLGN_annotation$D2_query$predicted_celltypes==macaque_celltype$celltype)

FALSE  TRUE 
   31  1693
```


