# R code to reproduce the analysis of the paper

A benchmarking study of individual somatic variant callers and voting-based ensembles for whole-exome sequencing


### dependencies

* R version >= 4.2.1
* library(circlize)
* library(ComplexHeatmap)
* library(ggpubr)
* library(stringr)
* library(e1071)
* library(ComplexHeatmap)
* library(ggrepel)
* library(future)

### Install

R docker with all the dependencies

```
sudo singularity build Rsom.simg docker://ngsom/bioconductor-base
```

### Data



### Run

```
singularity shell Rsom.simg
```
