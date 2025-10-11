# **Convert** maybe easy
- **Brief:** 可能比较轻松的rds和h5ad转换流程(支持多矩阵转换)
- **Fature** 
- **Log:** 
  - v1.0.1
    - 251011 规避`CreateAssayObject`修改基因名'_'为'-'的问题
    - 250827 完善Description
    - 250819 第一次提交


---
# Input
- **Variable**
  - `input_file` 输入的.rds或.h5ad文件，后缀必须为.rds或.h5ad
  - `layers` 要转换的矩阵，如果转所有矩阵则输入`all`，如果只保留.layers['counts']或RNA@counts则输入`RNA`
  - `mem_scdatacg` scdatacg需要的内存(GB)

**Note:** h5ad转rds时，会默认先用layers['counts']赋值给.X后再转为RNA$counts，保证了输出的h5ad为原始文件, 如果没有layers['counts']才直接.X转

- **csv** [download](https://github.com/ydgenomics/WDL/blob/main/Convert/v1.0.1/Convert_v1.0.1.csv)
- **Example** 

| EntityID | input_file | layers | mem_scdatacg |
|-|-|-|-|
| cotton | /Files/yangdong/wdl/SCP/Merge_Subset/W202508270022225/cotton_fibre_0.5/cotton_C1.hr.rds.rds | RNA | 4 |
| cotton | /Files/yangdong/wdl/SCP/Merge_Subset/W202508270022225/cotton_fibre_0.5/cotton_D3.hr.rds.rds |   |   |
| cotton | /Files/yangdong/wdl/SCP/Merge_Subset/W202508270022225/cotton_fibre_0.5/cotton_E1.hr.rds.rds |   |   |
| cotton | /Files/yangdong/wdl/SCP/Merge_Subset/W202508270022225/cotton_fibre_0.5/cotton_G3.hr.rds.rds |   |   |
| cotton | /Files/yangdong/wdl/SCP/Merge_Subset/W202508270022225/cotton_fibre_0.5/cotton_K2.hr.rds.rds |   |   |


---
# Output
- **Frame**
```shell
tree /data/input/Files/yangdong/wdl/SCP/Convert/W202508270039982
/data/input/Files/yangdong/wdl/SCP/Convert/W202508270039982
├── cotton_C1.hr.rds.rds
├── cotton_C1.hr.rds.rh.h5ad
├── cotton_D3.hr.rds.rds
├── cotton_D3.hr.rds.rh.h5ad
├── cotton_E1.hr.rds.rds
├── cotton_E1.hr.rds.rh.h5ad
├── cotton_G3.hr.rds.rds
├── cotton_G3.hr.rds.rh.h5ad
├── cotton_K2.hr.rds.rds
├── cotton_K2.hr.rds.rh.h5ad
└── input.json

1 directory, 11 files
```

- **Interpretation**
  - `.rh.h5ad` rds转h5ad的输出文件
  - `.hr.rds` h5ad转rds的输出文件


---
# Detail
- **Pipeline**
  - **rds2h5ad** 先根据Assay进行拆分为多个rds文件，然后对每个Assay的counts转h5ad，最后将多个h5ad按layers合并，RNA@counts默认为.X。可以通过指定`layers`来控制需要转的矩阵，`all`即转全部
  - **h5ad2rds** 先根据layers进行拆分为多个h5ad文件，然后对每个layers的h5ad转rds, 最后将多个rds按Assay合并，RNA@counts默认为.layers['counts']。同上可以进行矩阵指定。*值得注意的是h5ad转rds会补齐counts和data矩阵，此时data矩阵其实也是原始矩阵==counts，如果后续分析需要标准化的矩阵，请做处理后使用*。
- **Software**
  - sceasy
  - schard
- **Script**
  - deal_layers_ydgenomics.py
  - convert_rdsAh5ad.R
  - convert_rdsAh5ad2.R
  - convert.sh
- **Image**
  - sceasy-schard--10 sceasy-schard--02

---
# Reference & Citation
> [100w个细胞的单细胞h5ad对象转换为seurat：python版本](https://mp.weixin.qq.com/s/DxgZ-p9b51mzGASpxyhPNw)
> [认识单细胞分析中的各种数据结构](https://mp.weixin.qq.com/s/D7a2RlbJttaOlUNDAU5P3Q)

---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [WDL/Convert](https://github.com/ydgenomics/WDL/tree/main/Convert)