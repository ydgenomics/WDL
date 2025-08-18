# Method of moments framework for differential expression analysis by [memento]()
- **Brief:** 除了解释平均基因表达的差异还挖掘表达在每个细胞多样性的差异
- **Fature** 理解结果的各个指标，选择最优的阈值和可视化方案
- **Log:** 
  - v1.0.2 250816
- **Tradition:** DEA-memento

---
# Input
- **Variable**

---
# Output
- **Interpretation**
  - `.txt` 潜在差异基因的各项指标，包括未校正的p值和平均表达系数和多样性系数
  - `.pdf` 散点图，注释前几个低pval的，黑色圆点的是小于p_threshold的基因
  - `boxplot` 箱线图可视化前几个基因在不同处理下不同样本的分布
  - `2d text` 不建议运行，计算前几个基因与其它基因的相关性

---
# Detail
- **Overview**
关注表达和多样性都变高的基因，按de_pval做小到大取前十做可视化

- **Software**
  1.超几何分布恢复基因表达的真实分布(依赖于测序的`capture_rate`)；2.考虑到协变量(sample)的差异平均表达和差异多样性(方差); 3.使用的p值是uncorrected
`de_coef` and `de_pval` refer to the **mean effect size** and **uncorrected differential mean p-value** respectively.
`ve_coef` and `dv_pval` refer to the **variability effect size** and **uncorrected differential variability p-value** respectively.

- **Image**
```bash
pip install memento-de
```

- **Plot**
  [Nature 复现 | 多分组差异基因展示 log2FC-FC 散点图](https://mp.weixin.qq.com/s/UJRE4cX7zjQAT0kj9XEZDw)

---
# Reference & Citation
> 


---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [Scripts/enrich_scRNAseq](https://github.com/ydgenomics/Scripts/tree/main/enrich_scRNAseq)
---
