Repitope: Epitope immunogenicity prediction from interresidue contact potential profiling
===============================================

The 'Repitope' package provides a structured framework of quantitative prediction of the immunogenicity of MHC-I-loaded peptides through interresidue contact potential profiling.
Note: Repitope is currently under construction.

Installation
------------------------
Install the latest version from [GitHub](https://github.com/masato-ogishi/Repitope) as follows:
``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/plotUtility")
devtools::install_github("masato-ogishi/Repitope")
```
-   You might be prompted to install some packages before installling Repitope. Follow the message(s).

Usage
------------------
0. Working environment setup
``` r
options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
library(tidyverse)
library(Repitope)
```
1. Dataset
-   The following datasets are compiled.
``` r
MHCI_Human
MHCI_Rodents
MHCI_Primates
MHCII_Human
MHCII_Rodents
MHCII_Primates
```
-   Also, the compilation process itself is defined as a single function.
``` r
# Epitope datasets [MHC class I]
MHCI_Human <- Epitope_Import(
  system.file("IEDB_Assay_MHCI_Human.csv.gz", package="Repitope"),
  OtherFileNames=list(
    system.file("Calis1.csv", package="Repitope"),    ## Supplementary dataset from Calis et al., 2013.
    system.file("Calis2.csv", package="Repitope"),    ## Supplementary dataset from Calis et al., 2013.
    system.file("Chowell.csv", package="Repitope"),   ## Supplementary dataset from Chowell et al., 2015.
    system.file("EPIMHC.csv", package="Repitope"),    ## ClassI, Human, Annotated with T-cell activity
    system.file("HCV.csv", package="Repitope"),       ## LANL HCV epitope dataset.
    system.file("HIV.csv", package="Repitope"),       ## LANL HIV epitope dataset. ("best-defined")
    system.file("IMMA2.csv", package="Repitope"),     ## POPISK paper (Tung et al., 2011.), http://140.113.239.45/POPISK/download.php
    system.file("MHCBN.csv", package="Repitope"),     ## ClassI, Human, Annotated with T-cell activity
    system.file("TANTIGEN.csv", package="Repitope")   ## TANTIGEN T cell epitope dataset; entries annotated by in vitro or in vivo experiments are retained (but not MS experiments)
  ),
  peptideLengthSet=8:11
)
## 1873/21162 (8.85%) peptides have contradicting annotations.
## 6957/21162 (32.9%) peptides are immunogenic in at least one of the observations.
# Epitope datasets [MHC class II]
MHCII_Human <- Epitope_Import(
  system.file("IEDB_Assay_MHCII_Human.csv.gz", package="Repitope"),
  peptideLengthSet=11:30
)
## 4505/31693 (14.2%) peptides have contradicting annotations.
## 16642/31693 (52.5%) peptides are immunogenic in at least one of the observations.
```
-   Create TCR fragment library from the internally stored public TCR clonotype set.
``` r
fragLibDT <- CPP_FragmentLibrary(TCRSet_Public, fragLenSet=3:11, maxFragDepth=100000, seedSet=1:5)
fst::write_fst(fragLibDT, "./Path/To/Your/Directory/FragmentLibrary_TCRSet_Public.fst", compress=0)
```

Reference
------------------------
Ogishi, M and Yotsuyanagi, H. (2018) "The landscapes of T cell epitope immunogenicity in sequence space." bioRxiv. doi: https://doi.org/10.1101/155317
