# Cluster/Cell similarity with `jaccard`, `hclust` and [metaNeighbor](https://github.com/maggiecrow/MetaNeighbor)
- **Brief:** 分组数据间各个cluster的相似性
- **Fature** 优化输出可视化sankeyplot加颜色，heatmap.pdf的大小关联分群数量
- **Log:**
  - v1.1.4
    - 250905 sankey_order设置太麻烦了，修改为非必填项，无输入则可自动构建order; 修改pdf输出大小的判断，会根据亚群数目设置pdf宽度，环状热图经常报错已经取消，有需要可以在仓库的PLOT找到plot.ipynb作为个性化可视化
    - 250827 跟新Description
    - 250815(v1.1.4) 增添了`only_metaNeighbor`的判断默认只运行`metaNeighbor`，检查metaNeighbor运行的矩阵对象为`counts`，解决了可视化热图颜色块不显示的问题
- **Tradition:** `cell_similarity` `jm`

---
# Input
- **Variable**
  - `rds` 包含分组信息和cluster信息的rds对象，只跑metaNeighbor只需要rds包含`RNA@counts`,提取`counts`作为singlecellexperiment对象
  - `prefix` 每个`rds`运行之后输出文件的前缀，顺序对应
  - `sankey_order` 将每个`rds`里面`batch_key`用'|'连接用于设定sankey plot的顺序
  - `batch_key` 分组信息键，默认每个`rds`文件的`batch_key`一致
  - `cluster_key` 分群信息键，默认每个`rds`文件的`cluster_key`一致
  - `mem_similarity` 运行similarity的资源(GB)
- **Example** [download](https://github.com/ydgenomics/WDL/blob/main/Similarity/v1.1.4/Similarity_v1.1.4.csv)

| EntityID | rds | prefix | sankey_order | batch_key | cluster_key | mem_similarity |
|-|-|-|-|-|-|-|
| test_peanut | /Files/yangdong/wdl/SCP/Dataget/W202508040017201/01_dataget/peanut/peanut_merge.rds | peanut] | H1314|H2014 | biosample | leiden_res_0.50 | 32 |

---
# Output
主要看metaNeighbor结果，hclust和jaccard的结果仅供参考
- **Frame**
```shell
tree /data/input/Files/yangdong/wdl/SCP/Similarity/W202508060107189
/data/input/Files/yangdong/wdl/SCP/Similarity/W202508060107189
├── input.json
└── peanut
    ├── peanut_metaNeighbor.csv
    ├── peanut_metaNeighbor.pdf
    ├── peanut_metaNeighbor_tophits.csv
    ├── peanut_sanky.html
    ├── peanut_sct.rds
    └── seq.txt

2 directories, 7 files
```
- **Next**
  - Anno 细胞注释
  - Integration_scIB 整合去批次和评测
- **Interpretation**
  - 热图可视化相似性AUROC
  - sankey图默认仅展示相似性大于0.8的条目`slimit=0.80`
  - .csv文件为矩阵文件，可供个性化分析可视化

---
# Detail
- **Overview** 提取rds的RNA@counts为sce对象，找sce对象的批次的hvgs做相似性计算
- **Software**
  - **metaNeighbor**: The output is an AUROC matrix, where each value represents the similarity between cell types across batches or datasets. Higher AUROC values indicate greater similarity and more consistent annotation between datasets.
- **Script**
  - metaNeighbor.R
  - sanky_plot.py
- **Image**
  - **metaNeighbor--08** 供`git clone`维护流程用
  - **metaNeighbor--07**
- **test**
`as.SingleCellExperiment()` 是 Seurat 提供的方法，用于将 Seurat 对象转换为 SingleCellExperiment 对象（Bioconductor 生态常用的单细胞数据结构）。转换后的 SingleCellExperiment 对象会保留：表达矩阵（默认是 RNA assay 的 data 槽，即 log-normalized 数据;细胞元数据（colData，对应 Seurat 的 meta.data）; 基因元数据（rowData，对应 Seurat 的 feature metadata）；降维结果（如 PCA、UMAP，存储在 reducedDims 中）
```R
# 方法 1：直接提取 Seurat的counts 并构建 SCE
sce <- SingleCellExperiment(
  assays = list(counts = GetAssayData(sdata, slot = "counts")),
  colData = sdata@meta.data
)
# 方法 2：使用 Seurat::as.SingleCellExperiment() 并指定 slot
sce <- as.SingleCellExperiment(sdata, assay = "RNA", slot = "counts")
```
```R
#使用MetaNeighbor计算每个批次中细胞类型之间的相关性
library(MetaNeighbor)
Aurocs_matrix = MetaNeighborUS(var_genes = global_hvgs, 
                               dat = cca.results.sce, 
                               study_id = cca.results.sce$batch, 
                               cell_type = cca.results.sce$celltype, 
                               fast_version = T)
```

---
# Reference & Citation
- [Article](https://drive.google.com/file/d/17fc8paCD7v6RA6GfwRVNcqle85cx52my/view?usp=drive_link)
- [[R包] MetaNeighbor 第一期 评估不同数据集中细胞类型注释的一致性](https://mp.weixin.qq.com/s/cb9DWJm8zNc1J9wEUNTUVg)
- [常被提起的Jaccard指数是什么？怎么在单细胞中运用和实现Jaccard相似性比较？](https://mp.weixin.qq.com/s/-6iM2phNUh2Qo0wbN0Azpw)
- [练习R：hclust()函数层次聚类分析](https://mp.weixin.qq.com/s/-AvRPX7DG5fzyAVmv8wg7Q)


---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [WDL/Similarity](https://github.com/ydgenomics/WDL/tree/main/Similarity)