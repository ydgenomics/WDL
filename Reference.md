# A repository to display source of public data in DCS cloud
- **Brief:** 公共资源非常重要，一方面可以了解研究进展，另一方面是支撑进一步研究的宝贵资源。将这些资源整理并便于使用是有意义的，清楚的文件介绍和地址说明，避免重复下载，提升效率
- **Fature:** 以项目驱动，整理对应研究方向的资源
- **Log:**
  - 1.0.0 250809 第一次总结

---
# Database

---
# 拟南芥(Arabidopsis thaliana) [Tair](https://www.arabidopsis.org/)
  - **Protein**: Tiar下载蛋白质序列并对序列处理之后用eggnog-mapper做基因富集

|来源|云平台地址|035个性分析地址|
|-|-|-|
|[Tair](https://www.arabidopsis.org/download/list?dir=Proteins%2FTAIR10_protein_lists)|/Files/yangdong/Reference/Arabidopsis_thaliana/|/data/users/yangdong/yangdong_8632f88957bb4c4f85daf59edaf6b059/online/Reference/Arabidopsis_thaliana/TAIR10_pep_unique.fa|
|[eggnog-mapper](http://eggnog-mapper.embl.de/)|||


<details>
  <summary>蛋白质序列名处理脚本</summary>

```shell
# 蛋白质序列名是转录本名，去掉后缀拿到基因名
# 去掉序列名中.后面的部分包括.且只保留唯一序列名
awk '/^>/ {
    sub(/\.[0-9]+.*$/, "", $0)   # 去掉 .数字及之后内容
    if (!seen[$0]++) print;      # 第一次出现才打印
    next
}
{ print }' /data/work/Reference/TAIR10_pep_20101214.fa > TAIR10_pep_unique.fa
```
</details>

  - **scRNA**: scplantdb下载拟南芥Inflorescence的rds文件(SRP320285)

|ID|组织|细胞/基因|细胞类型|来源|云平台地址|035个性分析地址|
|-|-|-|-|-|-|-|
|SRP320285|Inflorescence|15662/78894|Celltype:"Shoot system epidermis","Vascular cambium","G2/M phase","Flower meristem","Mesophyll","Unknow","Cortex","S phase","Xylem"|[scplantdb](https://biobigdata.nju.edu.cn/scplantdb/dataset)|/Files/yangdong/Reference/Arabidopsis_thaliana/SRP320285.rds|/data/users/yangdong/yangdong_8632f88957bb4c4f85daf59edaf6b059/online/Reference/Arabidopsis_thaliana/SRP320285.rds| 

---
# 水稻(Oryza sativa) 

**GO**

> [我查完了，真的没有了 | 2025年中最全汇总之水稻常用数据库](https://mp.weixin.qq.com/s/ypXhQgVWu9-vYeLrJ8s34g)

---
# 玉米(Zea mays)(Maize)


> [育种工具｜玉米研究常用数据库](https://mp.weixin.qq.com/s/lkqxc27mEdnBgvb-GlyVSQ)

---
# 棉花(Zea mays)(Cotton)


> [中棉所生物信息中心与多家单位联合发布棉花多组学数据库](https://mp.weixin.qq.com/s/5fjQTzFnlPh1oNsrTyD5OQ)

---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [WDL/Reference.md](https://github.com/ydgenomics/WDL/blob/main/Reference.md)
---