# RiboPlotR
### Introduction
RiboPlotR package for Ribo-plot

Install required packages.
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
BiocManager::install("GenomicAlignments")
BiocManager::install("rtracklayer") 
BiocManager::install("Rsamtools") 

#Install RiboPlotR
library(devtools)
install_github("hsinyenwu/RiboPlotR")
```

### Examples
```R
# Load RiboPlotR and essential packages
library(RiboPlotR)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)

# Load datasets
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "RiboPlotR", mustWork = TRUE) #Annotation
ugtf <- system.file("extdata", "AT3G02468.gtf", package = "RiboPlotR", mustWork = TRUE) #uORF annotation
RRNA <- system.file("extdata", "Root_test_PE.bam", package = "RiboPlotR", mustWork = TRUE) #Root RNA-seq data
SRNA <- system.file("extdata", "Shoot_test_PE.bam", package = "RiboPlotR", mustWork = TRUE) #Shoot RNA-seq data
RRibo <- system.file("extdata", "riboRoot.bed", package = "RiboPlotR", mustWork = TRUE) #Root Ribo-seq data
SRibo <- system.file("extdata", "riboShoot.bed", package = "RiboPlotR", mustWork = TRUE) #Shoot Ribo-seq data

# Run gene.structure function to load gtf
gene.structure(annotation=agtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

# Run uorf.structure to load uORF gtf
uorf.structure(uorf_annotation=ugtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

# Run rna_bam.ribo to load root and shoot RNA-seq and Ribo-seq data sets
# Here root is the first dataset and shoot is the second dataset 
rna_bam.ribo(Ribo1=RRibo,
             RNAseqBam1=RRNA,
             RNAlab1="RNA count",
             Ribolab1="Ribo count",
             S_NAME1="Root",
             Ribo2=SRibo,
             RNAseqBam2=SRNA,
             RNAlab2="RNA count",
             Ribolab2="Ribo count",
             S_NAME2="Shoot",
             RNAseqBamPaired="paired")

#Plot AT4G21910 
PLOTc2("AT4G21910") #default using first isoform
```
![AT4G21910](https://github.com/hsinyenwu/RiboPlotR/blob/master/image/AT4G21910.png)  

```R
AT4G21910.png
PLOTc2("AT4G21910",isoform=2)

#Plot Root data (PLOTc uses the first RNA-seq and Ribo-seq dataset by default. Here the first dataset is the Root dataset.) 
PLOTc("AT3G02470",uORF = "AT3G02468",NAME=" SAMDC")
```

![SAMDC_uORF](https://github.com/hsinyenwu/RiboPlotR/blob/master/image/SAMDC_uORF.png)

```R
#Plot Shoot data (Here is an example how to plot the second dataset using PLOTc)
PLOTc("AT3G02470",uORF="AT3G02468",NAME=" SAMDC",RNAbam1 = RNAseqBam2, ribo1 = Ribo2, SAMPLE1 = "Shoot")

#Plot both dataset wiht PLOTC2
PLOTc2("AT3G02470",uORF = "AT3G02468",NAME=" SAMDC",isoform=3)
```

### Citation
Visualizing the periodic Ribo-seq reads with RiboPlotR
https://www.biorxiv.org/content/10.1101/694646v1
### Session Info
```R
sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils    
[7] datasets  methods   base     

other attached packages:
 [1] GenomicAlignments_1.20.1    Rsamtools_2.0.0            
 [3] Biostrings_2.52.0           XVector_0.24.0             
 [5] SummarizedExperiment_1.14.0 DelayedArray_0.10.0        
 [7] BiocParallel_1.18.0         matrixStats_0.54.0         
 [9] GenomicFeatures_1.36.3      AnnotationDbi_1.46.0       
[11] Biobase_2.44.0              GenomicRanges_1.36.0       
[13] GenomeInfoDb_1.20.0         IRanges_2.18.1             
[15] S4Vectors_0.22.0            BiocGenerics_0.30.0        
[17] RiboPlotR_0.1.0             usethis_1.5.1              
[19] devtools_2.0.2             

loaded via a namespace (and not attached):
 [1] progress_1.2.2         remotes_2.1.0         
 [3] lattice_0.20-38        rtracklayer_1.44.0    
 [5] yaml_2.2.0             blob_1.1.1            
 [7] XML_3.98-1.20          rlang_0.4.0           
 [9] pkgbuild_1.0.3         glue_1.3.1            
[11] withr_2.1.2            DBI_1.0.0             
[13] bit64_0.9-7            sessioninfo_1.1.1     
[15] GenomeInfoDbData_1.2.1 stringr_1.4.0         
[17] zlibbioc_1.30.0        memoise_1.1.0         
[19] callr_3.3.0            biomaRt_2.40.1        
[21] ps_1.3.0               curl_3.3              
[23] Rcpp_1.0.1             backports_1.1.4       
[25] desc_1.2.0             pkgload_1.0.2         
[27] fs_1.3.1               bit_1.1-14            
[29] hms_0.4.2              digest_0.6.20         
[31] stringi_1.4.3          processx_3.4.0        
[33] grid_3.6.0             rprojroot_1.3-2       
[35] cli_1.1.0              tools_3.6.0           
[37] bitops_1.0-6           magrittr_1.5          
[39] RCurl_1.95-4.12        RSQLite_2.1.1         
[41] crayon_1.3.4           pkgconfig_2.0.2       
[43] Matrix_1.2-17          prettyunits_1.0.2     
[45] assertthat_0.2.1       httr_1.4.0            
[47] rstudioapi_0.10        R6_2.4.0              
[49] compiler_3.6.0      
```
