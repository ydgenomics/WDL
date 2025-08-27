# DEAs: Revealing the difference genes of multiple clusters(DEA-Seurat/DEA-memento/DEA-FindMarkers)
- **Brief:**
  - DEA-Seurat 基于Seurat找各个cluster的marker基因并做可视化
  - DEA-memento 考虑到多样性(方差)变化做两分组数据的差异分析 
  - DEA-FindMarkers	基于Seurat做两分组数据的差异分析
- **Fature** 提供更好的可视化方案
- **Log:**
  - 250827 第一次完善Description

---
# Input
- **Variable**
  - DEA-Seurat
    - `rds`
    - `prefix`
    - `cluster_key`
    - `only_pos`
    - `mem_markers`
    - `batch_key`
  - DEA-memento
    - `input_h5ad` 
    - `group_key` 
    - `control_value`
    - `sample_key` 
    - `top_number` 
    - `mem_memento` 
  - DEA-FindMarkers
    - `rds` 
    - `prefix`
    - `ident_1` 
    - `ident_2`
    - `cluster_key`
    - `mem_findmarkers`




- **Example** 

DEA-Seurat [download]()

| EntityID | rds | prefix | cluster_key | only_pos | mem_markers | batch_key |
|-|-|-|-|-|-|-|
| test_peanut | /Files/yangdong/wdl/SCP/Dataget/W202508040017201/01_dataget/peanut/peanut_merge.rds | peanut | leiden_res_0.50 | yes | 16 | sample |

DEA-memento [download]()

| EntityID | input_h5ad | group_key | control_value | sample_key | top_number | mem_memento |
|-|-|-|-|-|-|-|
| yd_test | /Files/yangdong/wdl/SCP/Integration/W202508120010920/03_integration/peanut.h5ad | biosample | H1314 | sample | 2 | 16 |

DEA-FindMarkers [download]()

| EntityID | rds | prefix | ident_1 | ident_2 | cluster_key | mem_findmarkers |
|-|-|-|-|-|-|-|
| yd_test | /Files/yangdong/wdl/SCP/Merge_Subset/W202508190026068/peanut_0/peanut_0.rds | peanut_0 | H1314 | H2014 | biosample | 8 |


---
# Output
- **Frame**

DEA-Seurat
```shell
tree /data/input/Files/yangdong/wdl/SCP/DEA-Seurat/W202508190044072
/data/input/Files/yangdong/wdl/SCP/DEA-Seurat/W202508190044072
├── input.json
└── peanut_markers
    ├── allmarkers_peanut_merge.rds.csv
    ├── conserved_markers_peanut_merge.rds.csv
    ├── volcano_allmarkers_peanut_merge.rds.csv_0.01.csv
    ├── volcano_allmarkers_peanut_merge.rds.csv_0.01.pdf
    ├── volcano_conserved_markers_peanut_merge.rds.csv_0.01.csv
    └── volcano_conserved_markers_peanut_merge.rds.csv_0.01.pdf

2 directories, 7 files
```
DEA-memento
```shell
tree /data/input/Files/yangdong/wdl/SCP/DEA-memento/W202508190034508
/data/input/Files/yangdong/wdl/SCP/DEA-memento/W202508190034508
├── input.json
└── memento
    ├── differential_expression_replicate_peanut.pdf
    ├── result_1d_replicate_peanut.txt
    └── topgenes_boxplot
        ├── output_arahy.Tifrunner.gnm2.ann2.Ah16g512800_boxplot.pdf
        ├── output_arahy.Tifrunner.gnm2.ann2.Ah17g300100_boxplot.pdf
        ├── output_arahy.Tifrunner.gnm2.ann2.Ah19g485300_boxplot.pdf
        └── output_STRG.4445_boxplot.pdf

3 directories, 7 files
```
DEA-FindMarkers
```shell
tree /data/input/Files/yangdong/wdl/SCP/DEA-FindMarkers/W202508190044611
/data/input/Files/yangdong/wdl/SCP/DEA-FindMarkers/W202508190044611
├── input.json
└── peanut_0_findmarkers
    ├── markers_peanut_0.rds_H1314_vs_H2014.csv
    └── volcano_plot_H1314_vs_H2014.pdf

2 directories, 3 files
```

- **Next**
  - Enrich 富集分析


- **Interpretation**



---
# Detail
- **Pipeline**

- **Software**

- **Script**

- **Image**

---
# Reference & Citation
> 


---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [WDL/DEAs](https://github.com/ydgenomics/WDL/tree/main/DEAs) [DEAs](https://github.com/ydgenomics/DEAs)