Repitope: Epitope immunogenicity prediction through in silico TCR-peptide contact potential profiling
===============================================

The 'Repitope' package provides a structured framework of quantitative prediction of immunogenicity for peptides bound to MHC class I and class II through TCR-peptide contact potential profiling (CPP).

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
2. Features
-   Features can be calculated as follows.
``` r
# Features [MHC I]
FeatureDFList_MHCI <- Features(
  peptideSet=unique(c(MHCI_Human$Peptide, MHCI_Rodents$Peptide, MHCI_Primates$Peptide)),
  fragLib="./Path/To/Your/Directory/FragmentLibrary_TCRSet_Public.fst",
  aaIndexIDSet="all",
  fragLenSet=3:8,
  fragDepthSet=10000,
  fragLibTypeSet="Weighted",
  seedSet=1:5,                                   ## must be the same seed set for the fragment library
  coreN=parallel::detectCores()                  ## parallelization
  tmpDir="./Path/To/Your/Temporary/Directory/"   ## where intermediate files are stored
)
saveFeatureDFList(FeatureDFList_MHCI, "./Path/To/Your/Directory/MHCI/FeatureDF_")
# Features [MHC II]
FeatureDFList_MHCII <- Features(
  peptideSet=unique(c(MHCII_Human$Peptide, MHCII_Rodents$Peptide, MHCII_Primates$Peptide)),
  fragLib="./Path/To/Your/Directory/FragmentLibrary_TCRSet_Public.fst",
  aaIndexIDSet="all",
  fragLenSet=3:11,
  fragDepthSet=10000,
  fragLibTypeSet="Weighted",
  seedSet=1:5,                                   ## must be the same seed set for the fragment library
  coreN=parallel::detectCores()                  ## parallelization
  tmpDir="./Path/To/Your/Temporary/Directory/"   ## where intermediate files are stored
)
saveFeatureDFList(FeatureDFList_MHCII, "./Path/To/Your/Directory/MHCII/FeatureDF_")
```
-   Feature selection can be performed as follows.
``` r
# Feature selection [MHC I]
dt_feature_MHCI <- fst::read_fst("./Path/To/Your/Directory/MHCI/FeatureDF_Weighted.10000.fst", as.data.table=T)
minFeatureSet_MHCI_Human <- Features_MinimumFeatures(
  featureDFList=list(dt_feature_MHCI[Peptide%in%MHCI_Human$Peptide,]),
  metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)][,Cluster:=.I],
  seedSet=1:5,
  corThreshold=0.75,
  featureNSet=100,
  criteria="intersect",
  returnImpDFList=T
)
saveRDS(minFeatureSet_MHCI_Human, "./Path/To/Your/Directory/MHCI/MinimumFeatureSet_MHCI_Human.rds")
# Feature selection [MHC II]
dt_feature_MHCII <- fst::read_fst("./Path/To/Your/Directory/MHCII/FeatureDF_Weighted.10000.fst", as.data.table=T)
minFeatureSet_MHCII_Human <- Features_MinimumFeatures(
  featureDFList=list(dt_feature_MHCII[Peptide%in%MHCII_Human$Peptide,]),
  metadataDF=MHCII_Human[,.(Peptide,Immunogenicity)][,Cluster:=.I],
  seedSet=1:5,
  corThreshold=0.75,
  featureNSet=100,
  criteria="intersect",
  returnImpDFList=T
)
saveRDS(minFeatureSet_MHCII_Human, "./Path/To/Your/Directory/MHCII/MinimumFeatureSet_MHCII_Human.rds")
```
3. Immunogenicity scores
-   Probability estimates from multiple extremely randomized trees (ERTs) are summrized into a single numerical scale.
-   Predictions can be performed as follows.
``` r
# Parameters & datasets
featureDF_MHCI <- fst::read_fst("./Path/To/Your/Directory/MHCI/FeatureDF_Weighted.10000.fst", as.data.table=T)
featureDF_MHCI_Human <- featureDF_MHCI[Peptide%in%MHCI_Human$Peptide, ]
featureDF_MHCII <- fst::read_fst("D:/Research/Immunogenicity/MHCII/FeatureDF_Weighted.10000.fst", as.data.table=T)
featureDF_MHCII_Human <- featureDF_MHCII[Peptide%in%MHCII_Human$Peptide, ]
# MHCI
scoreDF_MHCI_Human <- Immunogenicity_Score(
  featureDF=featureDF_MHCI_Human,
  metadataDF=MHCI_Human[,.(Peptide, Immunogenicity)],
  featureSet=MHCI_Human_MinimumFeatureSet,
  seedSet=1:5
)
readr::write_csv(scoreDF_MHCI_Human, "./Path/To/Your/Directory/MHCI/ScoreDF_MHCI_Human.csv")
# MHCII
scoreDF_MHCII_Human <- Immunogenicity_Score(
  featureDF=featureDF_MHCII_Human,
  metadataDF=MHCII_Human[,.(Peptide, Immunogenicity)],
  featureSet=MHCII_Human_MinimumFeatureSet,
  seedSet=1:5
)
readr::write_csv(scoreDF_MHCII_Human, "./Path/To/Your/Directory/MHCII/ScoreDF_MHCII_Human.csv")
```

Reference
------------------------
Ogishi, M and Yotsuyanagi, H. (2018) "The landscape of T cell epitope immunogenicity in sequence space." bioRxiv. doi: https://doi.org/10.1101/155317
