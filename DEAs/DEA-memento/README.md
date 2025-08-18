# Method of moments framework for differential expression analysis by [memento](https://github.com/yelabucsf/scrna-parameter-estimation)
- **Brief:** 除了解释平均基因表达的差异还挖掘表达在每个细胞多样性的差异
- **Fature** 理解结果的各个指标，选择最优的阈值和可视化方案
- **Log:** 
  - v1.0.2 250816
- **Related:** DEA-FindMarkers

---
# Input
- **Variable**
  - 

---
# Output
- **Interpretation** Differential Variability Genes (DV) 的含义：差异变异性基因（Differential Variability Genes, DV）是指在不同条件下，基因表达的变异性存在显著差异的基因。具体来说，这些基因在某一条件下的表达变异性可能显著高于或低于其他条件下的变异性。这可能意味着这些基因在不同条件下的表达调控机制存在差异，或者受到不同因素的影响，导致其表达的稳定性不同。例如，在某种疾病状态下，某些基因的表达变异性可能增加，这可能与疾病的发生发展或异质性有关。
  - `.txt` 潜在差异基因的各项指标，包括未校正的p值和平均表达系数和多样性系数
  - `.pdf` 散点图，注释前几个`top_number`低pval的，黑色圆点的是小于`P_threshold=0.05`的基因(蓝色是最低de_pval的基因, 红色是最低dv_pval的基因)
  - `boxplot` 箱线图可视化前几个基因在不同处理下不同样本的分布
  - `2d text` 不建议运行，计算前几个基因与其它基因的相关性

---
# Detail
- **Overview**
关注表达和多样性都变高的基因，按de_pval做小到大取前十做可视化

- **Software**
  - `memento` is a Python package for estimating the mean, variability, and gene correlation from scRNA-seq data as well as contructing a framework for hypothesis testing of differences in these parameters between groups of cells. Method-of-moments estimators are used for parameter estimation, and efficient resampling is used to construct confidence intervals and establish statistical significance.
  1.超几何分布恢复基因表达的真实分布(依赖于测序的`capture_rate`)；2.考虑到协变量(sample)的差异平均表达和差异多样性(方差); 3.使用的p值是uncorrected
`de_coef` and `de_pval` refer to the **mean effect size** and **uncorrected differential mean p-value** respectively.
`ve_coef` and `dv_pval` refer to the **variability effect size** and **uncorrected differential variability p-value** respectively.

- **Image**
memento--12 /usr/bin/python
```bash
pip install memento-de
```

- **Plot**
  [Nature 复现 | 多分组差异基因展示 log2FC-FC 散点图](https://mp.weixin.qq.com/s/UJRE4cX7zjQAT0kj9XEZDw)

- **Script**
```shell
input_h5ad='/data/work/Single-Cell-Pipeline/DEA-memento/peanut.h5ad'
group_key="biosample"
control_value='H1314'
sample_key="sample"
capture_rate=0.07
min_perc=0.7
pval_threshold=0.05
n_cpu=2
top_number=1
perform_2d_test='yes'
/usr/bin/python /data/work/Single-Cell-Pipeline/DEA-memento/dea_memento.py \
--input_h5ad $input_h5ad --group_key $group_key --control_value $control_value --sample_key $sample_key \
--capture_rate $capture_rate --min_perc $min_perc --pval_threshold $pval_threshold \
--n_cpu $n_cpu --top_number $top_number --perform_2d_test $perform_2d_test
```

---
# Reference & Citation
> [Cell || Resource || 2024 || Memento: 用于单细胞 RNA 测序数据差异表达分析的矩量框架方法](https://mp.weixin.qq.com/s/WGl51Y8DiIQHybq4cM_bRA)
> [差异基因找的不好？Cell刚发的这个单细胞差异统计的方法，可以用到咱们自己的数据上](https://mp.weixin.qq.com/s/5kwfiapfBKVSQmo5Zuh6EA)
> [【综述】Nature Methods | 干货！一文读懂单细胞转录组分析的现状和问题！](https://mp.weixin.qq.com/s/pF95r2p0KK1LSW-l5YFrsw)
> [单细胞转录组差异分析的8大痛点](https://mp.weixin.qq.com/s/4VThQEdclByz6jaJ3EQOdQ)
> [单细胞下游分析 | 基础知识 | ③差异分析](https://mp.weixin.qq.com/s/2sQ2cbEQv2VWmKpXCTmTrQ)
> [特异性基因筛选方法比较](https://mp.weixin.qq.com/s/8g4lSm8az7dyo__RZALL1w)

---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [DEAs/DEA-memento]()


---