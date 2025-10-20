# Using [SoupX](https://github.com/constantAmateur/SoupX) and [scrublet](https://github.com/swolock/scrublet) do QC
- **Brief:** SoupX去除环境污染，scrublet评估双胞并去除预测为双胞的细胞，质控后用scanpy做标准化、降维、聚类、Marker和可视化
- **Fature:** 使用更新更适合数据的质控软件
- **Log:**
  - 1.2.2
    - 0827 更新Description
    - 0814 soupx后的数据会与两两不一样，增加添加splice和unsplice的元素不一样的判断(0814); doublet判断的图像缺失；缺少splice文件的问题；figure config的问题，已经注释掉(0814)
    - 1.2.2  250806 为`merge`增加了长度判断，对于无对照组即分组小于2的输入则不做merge
    - 250805
      - scrublet的结果存储在`.obs['predicted_doublet']` & `.obs['doublet_score']`去除`predicted_doublet`中为`True`的细胞，可视化中若`predicted_doublet`的umap为空则说明全是`False`，scanpy可视化umap至少两个及其以上变量数目。
      - 当缺少splice和unsplice矩阵或后续分析不涉及RNA Velocity时，可以不输入`SpliceMatrix`和`UnspliceMatrix`，为保证流程运行会将`FilterMatrix`/`SoupX结果`同时作为`SpliceMatrix`和`UnspliceMatrix`的输入。
      - `.X`里面是标准化后的矩阵，如果转rds应该先`adata.X=adata.layers['counts'].copy()`保证`.X`为原始数据。
      - 最新的脚本代码仓库[Dataget/v1.2.2](https://github.com/ydgenomics/WDL/tree/main/Dataget/v1.2.2)
  - 1.2.1
    - 修改为Dataget流程，支持多分组数据一次投递，新增子任务`merge`做.h5ad转.rds并merge做Seurat的标准化，该对象可用于后面做**Similarity**分析
  - 1.2.0
    - 0606 修改因`CreateSeuratObject()`自动更改基因名中'_'为'-'的问题，将task封装为函数即`run_*`；另外在流程部署上取消了脚本封装，避免多次保存环境和公布流程引起的维护问题，样本间concat也是取并集，尽量保存多的特征信息；另外在三个矩阵合并时基因取并集，至于细胞感觉也取并集，保留更多信息。
    - 20250516 统一了输出的marker基因csv包含的列`gene,cluster,p_val_adj,avg_log2FC`，便于下游分析；另外对多个resolution的marker基因的pdf和csv进行了保存`0.5, 0.8, 1.0`
    - 20250507 修改了三个矩阵存在细胞数不同的情况(Soupx处理后的矩阵)--取交集，修改了可视化pct_counts_mt的判断
    - 20250429 修改了三个矩阵整合为取基因的交集，另外为scrublet_estimate_doublecell.py运行添加了` > log.txt 2>&1`，用于保存运行过程信息
    - 20250417 优化了三个矩阵得到一个对象的基因选择，都以FilterMatrix为基准
    - 20250414 1.引入了splice和unsplice矩阵到anndata对象的layers中，有利于后面的RNA velocity分析; 2.将sample名作为后缀加到细胞名后面，保证了每个样本的细胞名不重复; 3.根据SoupX的默认参数maxrho为0.2，根据样本实际情况调整; 4.放弃了原先的大目录检索，之前的不利于流程维护，更加推荐大家使用表格投递任务
  - 1.0.0
    - 20250305 修复了无线粒体基因和有线粒体基因数据在QC质控的判断
- **Tradition:** dataget_scRNAseq


---
# Input
- **Variable:**
  - `RawMatrix` Array[Array[File]] 原始矩阵，同一分组下的数据在小list里面，大list代表不同分组数据
  - `FilterMatrix` Array[Array[File]] 质控过的矩阵，同上一一对应
  - `SpliceMatrix` Array[Array[File]] 剪切矩阵，同上一一对应，如无可以用FilterMatrix作为输入来保证流程正常运行，结果在.layers[splice]
  - `UnspliceMatrix` Array[Array[File]] 非剪切矩阵，同上一一对应，如无可以用FilterMatrix作为输入来保证流程正常运行,结果在.layers[unsplice]
  - `sample_value` Array[Array[String]] 自定义的样本信息，同上一一对应
  - `biosample_value` Array [String] 自定义的分组信息，对应大list
  - `species` String 自定义名称作为输出文件后缀
  - `mem_soupx` Int soupx的资源(GB)
  - `mem_scrublet` Int scrublet的资源(GB)
  - `mem_merge` Int merge的资源(GB)
  - `mitogenes_txt` File 线粒体基因列表
  - `mito_threshold` Int 线粒体基因最高占比，默认为5


- **csv** [download](https://github.com/ydgenomics/WDL/blob/main/Dataget/v1.2.2/Dataget_v1.2.2.csv) 
- **Example** 

| EntityID | RawMatrix1 | RawMatrix2 | FilterMatrix1 | FilterMatrix2 | SpliceMatrix1 | SpliceMatrix2 | UnspliceMatrix1 | UnspliceMatrix2 | sample_value1 | sample_value2 | biosample_value | species |
|-|-|-|-|-|-|-|-|-|-|-|-|-|
| test_peanut | /Files/husasa/anther-241218/H1314/V3RNA24120200006-2/HS-V3RNA24120200006/output/raw_matrix | /Files/husasa/anther-241218/H1314/V3RNA24120200007/HS-V3RNA24120200007/output/raw_matrix | /Files/husasa/anther-241218/H1314/V3RNA24120200006-2/HS-V3RNA24120200006/output/filter_matrix | /Files/husasa/anther-241218/H1314/V3RNA24120200007/HS-V3RNA24120200007/output/filter_matrix | /Files/husasa/anther-241218/H1314/V3RNA24120200006-2/HS-V3RNA24120200006/output/attachment/splice_matrix | /Files/husasa/anther-241218/H1314/V3RNA24120200007/HS-V3RNA24120200007/output/attachment/splice_matrix | /Files/husasa/anther-241218/H1314/V3RNA24120200006-2/HS-V3RNA24120200006/output/attachment/RNAvelocity_matrix | /Files/husasa/anther-241218/H1314/V3RNA24120200007/HS-V3RNA24120200007/output/attachment/RNAvelocity_matrix | V3RNA24120200006 | V3RNA24120200007 | H1314 | peanut |
| test_peanut | /Files/husasa/anther-241218/H2014/V3RNA24120200008_2/HS-V3RNA24120200008/output/raw_matrix | /Files/husasa/anther-241218/H2014/V3RNA24120200009/HS-V3RNA24120200009/output/raw_matrix | /Files/husasa/anther-241218/H2014/V3RNA24120200008_2/HS-V3RNA24120200008/output/filter_matrix | /Files/husasa/anther-241218/H2014/V3RNA24120200009/HS-V3RNA24120200009/output/filter_matrix | /Files/husasa/anther-241218/H2014/V3RNA24120200008_2/HS-V3RNA24120200008/output/attachment/splice_matrix | /Files/husasa/anther-241218/H2014/V3RNA24120200009/HS-V3RNA24120200009/output/attachment/splice_matrix | /Files/husasa/anther-241218/H2014/V3RNA24120200008_2/HS-V3RNA24120200008/output/attachment/RNAvelocity_matrix | /Files/husasa/anther-241218/H2014/V3RNA24120200009/HS-V3RNA24120200009/output/attachment/RNAvelocity_matrix | V3RNA24120200008 | V3RNA24120200009 | H2014 |  |


---
# Output
- **Frame**
```shell
tree /data/input/Files/yangdong/wdl/SCP/Dataget/W202508040017201
/data/input/Files/yangdong/wdl/SCP/Dataget/W202508040017201
├── 01_dataget
│   ├── H1314_dataget
│   │   ├── figures
│   │   │   ├── dotplot_leiden_res_0.50_marker.pdf
│   │   │   ├── dotplot_leiden_res_0.80_marker.pdf
│   │   │   ├── dotplot_leiden_res_1.00_marker.pdf
│   │   │   ├── pca_potentially_undesired_features.pdf
│   │   │   ├── umap_batch.pdf
│   │   │   ├── umap_leiden_clus.pdf
│   │   │   └── umap_quality.pdf
│   │   ├── H1314.h5ad
│   │   ├── marker_csv
│   │   │   ├── leiden_res_0.50.markers.csv
│   │   │   ├── leiden_res_0.80.markers.csv
│   │   │   └── leiden_res_1.00.markers.csv
│   │   ├── qc.pdf
│   │   └── summary.txt
│   ├── H1314_soupx_dataget
│   │   ├── figures
│   │   │   ├── dotplot_leiden_res_0.50_marker.pdf
│   │   │   ├── dotplot_leiden_res_0.80_marker.pdf
│   │   │   ├── dotplot_leiden_res_1.00_marker.pdf
│   │   │   ├── pca_potentially_undesired_features.pdf
│   │   │   ├── umap_batch.pdf
│   │   │   ├── umap_leiden_clus.pdf
│   │   │   └── umap_quality.pdf
│   │   ├── H1314_soupx.h5ad
│   │   ├── marker_csv
│   │   │   ├── leiden_res_0.50.markers.csv
│   │   │   ├── leiden_res_0.80.markers.csv
│   │   │   └── leiden_res_1.00.markers.csv
│   │   ├── qc.pdf
│   │   ├── summary.txt
│   │   ├── V3RNA24120200006_rho.pdf
│   │   ├── V3RNA24120200006_soupx_rho.txt
│   │   ├── V3RNA24120200007_rho.pdf
│   │   └── V3RNA24120200007_soupx_rho.txt
│   ├── H2014_dataget
│   │   ├── figures
│   │   │   ├── dotplot_leiden_res_0.50_marker.pdf
│   │   │   ├── dotplot_leiden_res_0.80_marker.pdf
│   │   │   ├── dotplot_leiden_res_1.00_marker.pdf
│   │   │   ├── pca_potentially_undesired_features.pdf
│   │   │   ├── umap_batch.pdf
│   │   │   ├── umap_leiden_clus.pdf
│   │   │   └── umap_quality.pdf
│   │   ├── H2014.h5ad
│   │   ├── marker_csv
│   │   │   ├── leiden_res_0.50.markers.csv
│   │   │   ├── leiden_res_0.80.markers.csv
│   │   │   └── leiden_res_1.00.markers.csv
│   │   ├── qc.pdf
│   │   └── summary.txt
│   ├── H2014_soupx_dataget
│   │   ├── figures
│   │   │   ├── dotplot_leiden_res_0.50_marker.pdf
│   │   │   ├── dotplot_leiden_res_0.80_marker.pdf
│   │   │   ├── dotplot_leiden_res_1.00_marker.pdf
│   │   │   ├── pca_potentially_undesired_features.pdf
│   │   │   ├── umap_batch.pdf
│   │   │   ├── umap_leiden_clus.pdf
│   │   │   └── umap_quality.pdf
│   │   ├── H2014_soupx.h5ad
│   │   ├── marker_csv
│   │   │   ├── leiden_res_0.50.markers.csv
│   │   │   ├── leiden_res_0.80.markers.csv
│   │   │   └── leiden_res_1.00.markers.csv
│   │   ├── qc.pdf
│   │   ├── summary.txt
│   │   ├── V3RNA24120200008_rho.pdf
│   │   ├── V3RNA24120200008_soupx_rho.txt
│   │   ├── V3RNA24120200009_rho.pdf
│   │   └── V3RNA24120200009_soupx_rho.txt
│   ├── peanut
│   │   ├── H1314.hr.rds
│   │   ├── H2014.hr.rds
│   │   ├── peanut_merge.pdf
│   │   ├── peanut_merge.rds
│   │   ├── saved_layers.txt
│   │   └── saved_paths.txt
│   └── peanut_soupx
│       ├── H1314_soupx.hr.rds
│       ├── H2014_soupx.hr.rds
│       ├── peanut_soupx_merge.pdf
│       ├── peanut_soupx_merge.rds
│       ├── saved_layers.txt
│       └── saved_paths.txt
└── input.json

16 directories, 73 files
```

- **Next**
  - Anno 细胞注释
  - Similarity 分组数据cluster间的相似性
  - anno_sctype 以细胞类型对应的高质量marker基因(.csv)为基础用sctype做注释加新键sctype
  - anno_singler
  - Integration_scIB


- **Interpretation**
  - 每个样本污染值评估(.pdf & .txt)
  - scanpy做去双胞后标准化流程(标准化、降维、聚类、找分群marker基因和可视化)，以dataget结尾的目录是未做SoupX处理的，soupx_dataget结尾的目录是做了SoupX再做的scrublet [Interpretation of results](https://mp.weixin.qq.com/s/xsxtCRFCi-y_3unfOkT-kQ)
  - 其它目录是子任务merge做的h5ad转rds后的整合(Seurat), rds都做过标准化，可以直接用于后续的anno系列流程


---
# Workflow
- **Pipeline:** 可选的用SoupX去除单细胞实验存在的环境污染(实验室环境RNA或细胞提前破碎导致的环境RNA等)，单细胞分离时可能存在双胞捕获为一个细胞，通过scrublet评估每个细胞为双胞的得分，最后基于质控后的数据集使用scanpy做降维聚类和可视化。多分组数据转rds之后做merge，便于后续的注释和Similarity分析
  - 路线1：评估环境污染后去污，再去除双胞，质控后做降维聚类可视化；
  - 路线2：只做去除双胞（考虑到去污效果差异和过处理），质控后做降维聚类可视化；

- **Software:**
  - SoupX：是一个R包，去除背景 RNA 污染——SoupX 利用空液滴（empty droplets）中的游离 RNA 和聚类信息来对表达量进行矫正，从而达到去噪效果。一个液滴捕获的数据是细胞内源 mRNA UMI 总和 + 游离 mRNA 的 UMI 总和 [demo](https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html) [SoupX tutorial](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html)
  - scrublet 是一个用于单细胞 RNA 测序（scRNA-seq）数据中检测双细胞（doublets）的 Python 工具。双细胞是指在实验过程中，两个或多个细胞被错误地封装在同一个液滴中，导致测序结果中出现混合的转录组信号。scrublet 通过模拟双细胞并使用 k-最近邻分类器来计算每个细胞的双细胞得分（doublet score）[demo](https://github.com/swolock/scrublet/blob/master/examples)

- **Script:**
  - scrublet_estimate.py
  - soupx.R
  - merge_subset.R

- **Image:**
  - soupx-R--03, SoupX-R--02
  - scrublet-py--05, scrublet-py--04
  - sceasy-schard-10, sceasy-schard--02

```shell
conda create -n r r-base=4.4 -y
conda activate r
conda install conda-forge::r-seurat -y
conda install conda-forge::r-soupx -y
conda install bioconda::bioconductor-decontx -y
conda install bioconda::presto -y
conda install bioconda::bioconductor-dropletutils -y
conda install conda-forge::r-optparse -y
conda install bioconda::r-sceasy -y
conda install conda-forge::r-reticulate -y
conda install conda-forge::r-devtools -y
# devtools::install_github("cellgeni/schard")
conda create -n py python=3.12 -y
conda activate py
conda install conda-forge::scanpy -y
conda install bioconda::scrublet -y
conda install conda-forge::leidenalg -y
```


---
# Reference & Citation
> **Sincerely thank all teachers and researchers who provide open source resources**
> 1. [SoupX——去除RNA污染](https://mp.weixin.qq.com/s/7g9Zo6IPqTafSjKCeAFNIQ)
> 2. [使用DecontX预测和去除单细胞转录组的环境游离RNA污染](https://mp.weixin.qq.com/s/ndt9Fsgg5dNxIOh9m7j9Bw)
> 3. [是否细胞周期矫正，去除双细胞和环境RNA污染——单细胞入门到进阶(初级篇2）](https://mp.weixin.qq.com/s/HgTVwfDfE4lzBXJKihlknA)
> 4. [*生信钱同学*·全代码干货奉上——多样本多方案去除单细胞环境RNA污染——这次把这个聊清楚](https://mp.weixin.qq.com/s/1eJq3u-aKpQaL9CM7bV94g)
> 4. [单细胞去噪工具一览](https://mp.weixin.qq.com/s/78RC4qH_Kw_eb-rql_QGjg)
![scCDC对各个去污工具的评估](../png/scCDC_ability.png)
> 5. [还在纠结双细胞质控方法吗！一文说清楚](https://mp.weixin.qq.com/s/64hB2cj-NwojuZbdiyEGzg)
![doublecell](../png/doublecell_ability.png)


---
# Coder
  - **Coder:** yangdong(yangdong@genomics.cn)
  - **Github:** [ydgenomics](https://github.com/ydgenomics)
  - **Prospect:** Do interesting and competitive works, open source, make progress and collaboration!
  - **Repository:** [WDL/Dataget]()