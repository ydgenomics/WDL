# Using multi-methods do Integration and Eliminate batch effects (Integration)
- **Brief:** 先对h5ad对象做整合加新键biosample，然后以biosample作为批次键值做去批次，最后可视化感兴趣的键
- **Fature**
- **Log:**
  - v1.0.3 250812 取消可视化`sample`，python的umap添加`legend_loc='on data'`，是图更易阅读
  - 250724 优化了任务运行条件判断`whether_sct`
- **Related:** *scIB* *Integration_scIB* *multi_samples_scRNAseq_integration*


---
# Input
- **Variable**
  - h5ad(Array[File])：待整合.h5ad文件，其中layers['counts']为原始数据 || layers['counts']不存在时.X为原始数据
  - biosample_value(Array[String]): 整合后对象'biosample'列对应的值
  - species(String): 物种名作为输出文件的prefix
  - whether_sct(String): 是否使用SCTransform.CCA和SCTransform.harmony，"yes"即使用，"no"则不适用
  - other1_key(String): 展示在umap图上信息的键
  - other2_key(String): 展示在umap图上信息的键
  - resolution(Float): 分辨率控制分群数量resolution
  - mem_concat(Int): concat的资源
  - mem_scvi(Int): scvi的资源
  - mem_integration2(Int): integration2的资源
  - mem_scdatacg(Int): scdatacg的资源
  - mem_integration4(Int): integration4的资源，如果whether_sct为"yes"需要提供更多的资源
  - mem_dealplus(Int): dealplus的资源

- **Table submit**
[.csv download]()


---
# Output
- **Interpretation**
  - `.pdf` 可视化umap
  - `.rds/.h5ad` 整合去批次输出的rds/h5ad文件


- **Downstream analysis**
  - scIB
  - DEA

---
# Workflow
- **Overview**: 实验差异引起的非生物学差异(批次效应)干扰下游分析，不同整合去批次方法在多样的数据中表现效果各异，做多种方法的整合去批次供选择。六种整合方法scVI, harmony, rliger.inmf, bbknnr, SCTransform.CCA和SCTransform.harmony，其中前4种默认都要运行。concat后的对象做scanpy的降维聚类得到unintegrated.h5ad，可以对比整合去批次和未整合去批次结果。
  1. 将多个.h5ad文件的原始文件`layers['counts']`做`concat`拿到一个scanpy对象(anndata)
  2. 直接运行harmonypy和scVI整合
  3. `scdatacg`将anndata转化为seurat对象
  4. 拿到转化后的.rds文件直接运行rliger,bbknnR,SCTransform.CCA,SCTransform.harmony
  5. `dealplus` 可视化感兴趣的键

- **Software**
  - [rliger](https://github.com/welch-lab/liger)
  - [harmony](https://github.com/immunogenomics/harmony)
  - [bbknnR](https://github.com/ycli1995/bbknnR)
  - [CCA]()
  - [scVI](https://github.com/scverse/scvi-tools)

- **Image**
  - Integration

---
# Reference
> 


---
# Coder info
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration.
- **Repository:** [WDL/Integration]()

---
