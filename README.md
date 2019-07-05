# RiboPlotR
### Introduction
RiboPlotR package for Ribo-plot

```R

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
BiocManager::install("GenomicAlignments")
BiocManager::install("rtracklayer") 
BiocManager::install("Rsamtools") 

library(devtools)  #install devtools if necessary
install_github("hsinyenwu/RiboPlotR")
```
### Basic functions

### Examples
```R
library(devtools)
install_github("hsinyenwu/RiboPlotR")

library(RiboPlotR)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)

ath <- system.file("extdata", "TAIR10.29_part.gtf", package = "RiboPlotR", mustWork = TRUE)
RRNA <- system.file("extdata", "Root_test.bam", package = "RiboPlotR", mustWork = TRUE)
SRNA <- system.file("extdata", "Shoot_test.bam", package = "RiboPlotR", mustWork = TRUE)
RRibo <- system.file("extdata", "riboRoot.bed", package = "RiboPlotR", mustWork = TRUE)
SRibo <- system.file("extdata", "riboShoot.bed", package = "RiboPlotR", mustWork = TRUE)

gene.structure(annotation=ath, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
uorf.structure(uorf_annotation=ath, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
rna_bam.ribo(rna1=RRNA,rna2=SRNA,ribo1=RRibo,ribo2=SRibo,
             RNAlab1="Shoot_RNA",
             RNAlab2="Root_RNA",
             Ribolab1="Shoot_Ribo",
             Ribolab2="Root_Ribo",
             RNAseqBamPaired="single")

PLOTc2("AT4G21910")
PLOTc("AT3G02470",uORF = "AT3G02468",NAME=" SAMDC")

```


### Citation
