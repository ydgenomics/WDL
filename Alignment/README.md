# Alignment: Find the best pairs with one2one
- **Brief:**蛋白质的序列比对并筛选出最优的一对一关系
- **Fature**
- **Log:** 
  - v1.0.0
    - 251014 新增`type`区分输入的是nucleotide还是protein，使用blastn做核酸序列比对
    - 250828 提交Description

---
# Input
- **Variable**
  - `fasta1` 待比对的蛋白质序列fasta文件
  - `fasta2` 待比对的蛋白质序列fasta文件
  - `method` 使用的序列比对软件可选 "diamond" or "blastp", 默认比对速度更快的`diamond`
  - `mem_alignment` 比对使用的内存GB
- **Example** [download](https://github.com/ydgenomics/WDL/blob/main/Alignment/v1.0.0/Alignment_v1.0.0.csv)

| EntityID | fasta1 | fasta2 | method | mem_alignment |
|-|-|-|-|-|
| test | /Files/yangdong/wdl/blast/Ghirsutum_gene_peptide_trimmed.fasta | /Files/yangdong/wdl/blast/TM-1_V2.1.gene.pep.fa | diamond | 8 |

---
# Output
- **Frame**
```shell
reciprocal_best.txt
```
- **Next**
  - Anno-singler 数据映射
- **Interpretation**
  - `reciprocal_best.txt` 最佳的一对一比对关系

---
# Detail
- **Pipeline** 蛋白质序列单独建库，然后分别比对到数据库，拿到两个多对多的关系，通过筛选确定最佳的一对一关系
- **Software**
  - diamond
  - blastp
- **Script**
  - diamond_blastp.sh
- **Image**
  - Alignment--01, Alignment

```shell
source /opt/software/miniconda3/bin/activate
conda create -n alignment -y
conda activate alignment
conda install bioconda::blast -y
conda install bioconda::diamond -y
```

---
# Reference & Citation
> [同源基因比对：顶刊都在用的跨物种分析](https://mp.weixin.qq.com/s/jv2Z8NVWZwzeVjmm5c9NNg)
> [史上最详细的blast安装附视频](https://mp.weixin.qq.com/s/rEBqjN-fGOp_loTmyEuMJA)
> [Clustal Omega—广泛使用的多序列比对工具](https://mp.weixin.qq.com/s/f9pEFWJJoNCqlFEfd77aOA)
> [DIAMOND:快又准的蛋白序列比对软件](https://mp.weixin.qq.com/s/5UhthY9PHfN7zxZbJdZaJA)


---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [WDL/Alignment](https://github.com/ydgenomics/WDL/tree/main/Alignment)

