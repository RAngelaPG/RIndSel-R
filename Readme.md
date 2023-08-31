## RIndSel: An R Package for calculated some different indexes

The RIndSel Package ([Pacheco & et. Al, 2016](https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10854)) calculated some different index like: Base Linear Phenotypic
Selection Index (BLPSI), Linear Phenotypic Selection Index (LPSI), Eigen Selection Index
Method (ESIM), Predetermined Proportional Gain Linear Phenotypic Selection Index (PPGLPSI),
Restricted Linear Phenotypic Selection Index (RLPSI), Restricted Eigen Selection
Index Method (RESIM), Linear Marker Selection Index (LMSI), Molecular Eigen Selection
Index Method (MESIM) and Linear Genomic Selection Index (LGSI). In this repository we maintain the latest
version beta version. The latest stable release can be downloaded from [CRAN](https://cran.r-project.org/web/packages/BGLR/index.html).


#### Citation

Please cite [Pacheco & et. Al, 2016](https://data.cimmyt.org/dataset.xhtml?persistentId=hdl:11529/10854).


#### Installation

**From CRAN (stable release)**.

```R
  install.packages(pkg='RIndSel',repos='https://cran.r-project.org/')
```

**From GitHub (development version, added features)**.


```R
   install.packages(pkg='devtools',repos='https://cran.r-project.org/')  #1# install devtools
   library(devtools)                                                     #2# load the library
   install_git('https://github.com/RAngelaPG/RIndSel-R')                 #3# install BGLR from GitHub
```
## Indexes
----------------------------------------------------------------

**Examples function**
----------------------------------------------------------------

  - [1. Linear Phenotypic Selection Index](https://github.com/RAngelaPG/RIndSel-R/blob/master/inst/md/LPSI.md)
