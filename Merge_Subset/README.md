# `Merge_Subset` is used to merge mutiple seurat objects(including subseted data)
- **Brief:** 将多个Seurat对象(.rds)整合为一个对象同时满足取子集(subset)后再整合，整合后重新标准化、降维、聚类和可视化
- **Fature** 
- **Log:** 
  - v1.0.0 
    - 250827 更新Description
    - 250818 第一次提交
- **Tradition:** Integration

---
# Input
- **Variable**
  - `input_rds` Array[Array[File]] 输入rds文件。分大list和小list，小list放要整合在一起的元素，设计大list考虑到只做subset而不做整合的需求
  - `cluster_key` Array[Array[String]] 取子集参考的键名，与上面的rds文件一一对应
  - `cluster_value` Array[Array[String]] 取子集参考键名下的值，与上面一一对应
  - `plot_keys ` Array [String] 期待可视化键名，与上面的大list一一对应，小list里面可视化多个键用,连接； 例如"sample,leiden_res_0.50"
  - `r_value` Array [Float] 降维的resolution，与上面的大list一一对应
  - `prefix` Array [String] 整合后文件输出的前缀
  - `mem_merge_subset` Int 内存资源(GB)
- **Example** [download](https://github.com/ydgenomics/WDL/blob/main/Merge_Subset/v1.0.0/Merge_Subset_v1.0.0.csv)

| EntityID | input_rds1 | input_rds2 | input_rds3 | input_rds4 | input_rds5 | cluster_key1 | cluster_key2 | cluster_key3 | cluster_key4 | cluster_key5 | cluster_value1 | cluster_value2 | cluster_value3 | cluster_value4 | cluster_value5 | plot_keys | r_value | prefix | mem_merge_subset |
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
| cotton | /Files/yangdong/cotton/NB2025081415270346626070/cotton_C1.hr.rds | /Files/yangdong/cotton/NB2025081415270346626070/cotton_D3.hr.rds | /Files/yangdong/cotton/NB2025081415270346626070/cotton_E1.hr.rds | /Files/yangdong/cotton/NB2025081415270346626070/cotton_G3.hr.rds | /Files/yangdong/cotton/NB2025081415270346626070/cotton_K2.hr.rds | leiden_res_0.50 | leiden_res_0.50 | leiden_res_0.50 | leiden_res_0.50 | leiden_res_0.50 | 5 | 3 | 2 | 5 | 5 | sample\|leiden_res_0.50 | 0.5 | cotton_fibre_0.5 | 8 |

---
# Output
- **Frame**
```shell
tree /data/input/Files/yangdong/wdl/SCP/Merge_Subset/W202508270022225
/data/input/Files/yangdong/wdl/SCP/Merge_Subset/W202508270022225
├── cotton_fibre_0.5
│   ├── cotton_C1.hr.rds.rds
│   ├── cotton_C1.hr.rds_umap.pdf
│   ├── cotton_D3.hr.rds.rds
│   ├── cotton_D3.hr.rds_umap.pdf
│   ├── cotton_E1.hr.rds.rds
│   ├── cotton_E1.hr.rds_umap.pdf
│   ├── cotton_fibre_0.5.rds
│   ├── cotton_fibre_0.5_umap.pdf
│   ├── cotton_G3.hr.rds.rds
│   ├── cotton_G3.hr.rds_umap.pdf
│   ├── cotton_K2.hr.rds.rds
│   └── cotton_K2.hr.rds_umap.pdf
└── input.json

2 directories, 13 files
```
- **Next**
  - Track-pseudotime 取子集的数据做拟时序分析
  - Similarity 整合后的数据查看分组间的cluster相似性
- **Interpretation**
  - `umap.pdf` 包含整合后降维可视化umap图和单独对象的umap图
  - `.rds` 包含整合后的.rds和单独对象标准化流程的.rds

---
# Detail
- **Pipeline**
  1. 读取rds文件地址(以|连接)
  2. 读取第一个元素，`cluster_key`必须存在于`meta.data`，判断`cluster_value`第一个元素是不是`all`，是则全部读取而不取子集，取子集按`cluster_value`元素来，所以必须都存在于`meta.data[[cluster_key]]`值中。
  3. 确定好对象是否取子集后就运行Seurat标准流程(标准化、降维、聚类和可视化)
  4. 判断rds文件地址数量是否大于1，大于1则做Merge，Merge过程中同第二步的处理
  5. 将单独降维得到的分辨率的值重命名为`*_0`，为整个对象降维留出键，最后可视化整个对象
- **Software**
  - Seurat
- **Script**
  - merge_subset.R
- **Image**
  - metaNeighbor--08 
- **Test**
```R
seurat_pipeline <- function(seu, r_value, plot_keys, prefix){
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, nfeatures = 3000)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, features = VariableFeatures(object = seu), verbose = FALSE)
    seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
    seu <- FindClusters(seu, resolution = as.numeric(r_value)) # RNA_snn_res.0.5
    seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
    pdf(paste0(prefix, "_umap.pdf"), width = 10, height = 8)
    for (plot_key in plot_keys){
        p1 <- DimPlot(seu, reduction = "umap", group.by = plot_key, shuffle = TRUE, label = TRUE)
        print(p1)
    }
    dev.off()
    saveRDS(seu, paste0(prefix,".rds"))
    return(seu)
}
```

---
# Reference & Citation
> [单细胞亚群取子集后的细分亚群再命名的两个难题](https://mp.weixin.qq.com/s/IfLm9MOpm6RKjCsXy_1AWw)
> [提取单细胞亚群进行后续再分析](https://mp.weixin.qq.com/s/dtsrP6AZx4_fHiVxeHxfZQ)

---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [WDL/Merge_Subset](https://github.com/ydgenomics/WDL/tree/main/Merge_Subset)