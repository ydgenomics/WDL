# **Convert** maybe easy
- **Brief:** 可能比较轻松的rds和h5ad转换流程(支持多矩阵转换)
- **Fature**
- **Log:** 
  - v1.0.1
    - 250819 第一次提交


---
# Input
- **Variable**
  - `input_file` 输入的.rds或.h5ad文件，后缀必须为.rds或.h5ad
  - `tools` 选择工具`single_convert`或`multi_convert`，默认使用`multi_convert`
  - `layers` 要转的矩阵，如果转所有矩阵则输入`all`，如果只保留.layers['counts']或RNA@counts则输入`RNA`

h5ad转rds时，如果选的multi_convert,layers='RNA',会默认先用layers['counts']赋值给.X后再转为RNA$counts，保证了输出的h5ad为原始文件, 如果没有layers['counts']才直接.X转；而如果选的single_convert,默认直接从.X转，.X为非原始矩阵的话会导致RNA$counts也为非原始矩阵

---
# Output
- **Interpretation**
  - `.rh.h5ad` rds转h5ad的输出文件
  - `.hr.rds` h5ad转rds的输出文件


---
# Detail
- **Overview**

- **Software**
  - sceasy
  - schard

- **Image**
  - sceasy-schard--10 sceasy-schard--02

---
# Reference & Citation
> 


---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [Convert]()


---
