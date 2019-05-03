Repitope: Epitope immunogenicity prediction through in silico TCR-peptide contact potential profiling
===============================================

The 'Repitope' package provides a structured framework of quantitative prediction of immunogenicity and escape potential for a given set of peptides presented onto MHC class I and class II molecules by (approximately) simulating the TCR-peptide intermolecular interactions in silico.  
Note: This is the original Repitope package repository in R. For Python scripts, see [Repitope-Python](https://github.com/masato-ogishi/Repitope-Python).

Concept
------------------------
The underlying concept of the Repitope framework is simple. We assume that there are some hidden, non-linear patterns in the sequences of epitopes with high population-wide immunogenicity (i.e., prone to trigger T cell activation in multiple individuals), and that human TCR sequences in turn have some patterns optimized for recognition of such immunogenic epitopes. Because TCR-epitope interaction is governed by the physicochemical principles like other protein-protein interactions, we assume that amino acid pairwise contact potential matrices, derived from studies on protein folding, can be utilized to describe the necessary features for TCR-epitope interaction. In this scenario, we could utilize TCR sequences as "baits" to probe highly immunogenic epitopes. Utilizing minimally pre-selected features (32 and 27 features for MHC-I and MHC-II epitopes, respectively), we conducted machine learning and computed a probabilistic estimate of immunogenicity (which we termed "immunogenicity score") for each epitope. Please note that our pipeline utilizes a large set of TCR sequences as a whole, and therefore cannot specifically identify TCR clones recognizing specific epitopes. In other words, the immunogenicity score is intended for predicting a population-level, not an individual-level, probability of effective T cell response.  
The concept of escape potential is also simple. Even if one epitope has a high immunogenicity by itself, single amino acid mutations (substitutions and indels) could drastically alter the immunogenicity. Such epitopes are considerd "unstable" in terms of T cell recognition (i.e., prone to escaping), and therefore may not be appropriate for therapuetic use. Although we focus only on sequence alteration here, post-translational modification (e.g., glycosylation) could have the same impact in terms of T cell recognition. We compute probabilistic estimates of imunogenicity for both an original epitope and all of its single-AA variants, and measure the maximal difference of immunogencity between the original and the variants. In short, the smaller escape potential could be a nindicator for more stable immune response.

Installation
------------------------
Install the latest version as follows:
``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/plotUtility")
devtools::install_github("masato-ogishi/Repitope")
```
-  You may be prompted to install some packages before installing Repitope. Follow the messages.  
-  This package depends on some packages in the [*Bioconductor*](https://www.bioconductor.org/) (e.g., BioStrings) that cannot be automatically installed through devtools::install_github. Install packages manually as required.  
-  You need an appropriate rJava setting beforehand.

Usage
------------------
0. Working environment
``` r
options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
library(tidyverse)
library(data.table)
library(Repitope)
```
1. Datasets
-  The following datasets are included in the package.
``` r
# A summary table for the peptide sequences with T cell assay annotations
MHCI_Human
MHCI_Rodents
MHCI_Primates
MHCII_Human
MHCII_Rodents
MHCII_Primates

# A compiled set of TCR CDR3b sequences commonly identified in multiple TCR datasets (see the Reference for details)
TCRSet_Public

# Minimal feature sets
MHCI_Human_MinimumFeatureSet    ## identified using MHCI_Human
MHCII_Human_MinimumFeatureSet   ## identified using MHCII_Human
```
-  The compilation process of epitope datasets itself is provided as a single function so that users can compile their own epitope datasets from IEDB and other sources.
``` r
# Epitope datasets [MHC-I]
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

# Epitope datasets [MHC-II]
MHCII_Human <- Epitope_Import(
  system.file("IEDB_Assay_MHCII_Human.csv.gz", package="Repitope"),
  peptideLengthSet=11:30
)
## 4505/31693 (14.2%) peptides have contradicting annotations.
## 16642/31693 (52.5%) peptides are immunogenic in at least one of the observations.
```
-  Create TCR fragment library from the internally stored public TCR clonotype set. Note: you need to use the same set of random seeds throughout the pipeline.
``` r
fragLibDT <- CPP_FragmentLibrary(TCRSet_Public, fragLenSet=3:11, maxFragDepth=100000, seedSet=1:5)
fst::write_fst(fragLibDT, "./Path/To/Your/Directory/FragmentLibrary.fst", compress=0)
```
2. Features
-  Features can be calculated as follows. Note: This computation is time-consuming and resource-intensive. Computation can be resumed if temporary files are stored in the temporary directory provided.  
-  You may find the pre-computed feature dataframes for the example datasets (i.e., human, rodents, and primates) in [Mendeley Data](https://data.mendeley.com/datasets/2hp96k6m2c/3). Note: these data must be downloaded with a good internet connection. Otherwise the downloaded files could be broken, and subsequent codes would not work properly. The expected file sizes are 1.3 and 2.6GB for MHC-I and MHC-II, respectively.
``` r
# Features [MHC-I]
featureDFList_MHCI <- Features(
  peptideSet=unique(c(MHCI_Human$Peptide, MHCI_Rodents$Peptide, MHCI_Primates$Peptide)),
  fragLib="./Path/To/Your/Directory/FragmentLibrary.fst",
  aaIndexIDSet="all",
  fragLenSet=3:8,
  fragDepthSet=10000,
  fragLibTypeSet="Weighted",
  seedSet=1:5,                                   ## must be the same random seeds used for preparing the fragment library
  coreN=parallel::detectCores(logical=F),        ## parallelization
  tmpDir="./Path/To/Your/Temporary/Directory/"   ## where intermediate files are stored
)
saveFeatureDFList(featureDFList_MHCI, "./Path/To/Your/Directory/MHCI/FeatureDF_")

# Features [MHC-II]
featureDFList_MHCII <- Features(
  peptideSet=unique(c(MHCII_Human$Peptide, MHCII_Rodents$Peptide, MHCII_Primates$Peptide)),
  fragLib="./Path/To/Your/Directory/FragmentLibrary.fst",
  aaIndexIDSet="all",
  fragLenSet=3:11,
  fragDepthSet=10000,
  fragLibTypeSet="Weighted",
  seedSet=1:5,                                   ## must be the same seed set for the fragment library
  coreN=parallel::detectCores(logical=F),        ## parallelization
  tmpDir="./Path/To/Your/Temporary/Directory/"   ## where intermediate files are stored
)
saveFeatureDFList(featureDFList_MHCII, "./Path/To/Your/Directory/MHCII/FeatureDF_")
```
-  Feature selection can be performed as follows.
``` r
# Feature selection [MHC-I]
featureDF_MHCI <- fst::read_fst("./Path/To/Your/Directory/MHCI/FeatureDF_Weighted.10000.fst", as.data.table=T)
minFeatureSet_MHCI_Human <- Features_MinimumFeatures(
  featureDFList=list(featureDF_MHCI[Peptide%in%MHCI_Human$Peptide,]),
  metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)][,Cluster:=.I],
  seedSet=1:5,
  corThreshold=0.75,
  featureNSet=100,
  criteria="intersect",
  returnImpDFList=T
)
saveRDS(minFeatureSet_MHCI_Human, "./Path/To/Your/Directory/MHCI/MinimumFeatureSet_MHCI_Human.rds")

# Feature selection [MHC-II]
featureDF_MHCII <- fst::read_fst("./Path/To/Your/Directory/MHCII/FeatureDF_Weighted.10000.fst", as.data.table=T)
minFeatureSet_MHCII_Human <- Features_MinimumFeatures(
  featureDFList=list(featureDF_MHCII[Peptide%in%MHCII_Human$Peptide,]),
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
-  Probability estimates from multiple extremely randomized trees (ERTs) are summrized into a single numerical scale, termed "immunogenicity score".  
-  Prediction can be performed as follows.
``` r
# Datasets
featureDF_MHCI <- fst::read_fst("./Path/To/Your/Directory/MHCI/FeatureDF_Weighted.10000.fst", as.data.table=T)
featureDF_MHCII <- fst::read_fst(""./Path/To/Your/Directory/MHCII/FeatureDF_Weighted.10000.fst", as.data.table=T)

# Probability estimation [MHC-I]
scoreDF_MHCI_Human <- Immunogenicity_Score(
  featureDF=featureDF_MHCI[Peptide%in%MHCI_Human$Peptide,],
  metadataDF=MHCI_Human[,.(Peptide, Immunogenicity)],
  featureSet=MHCI_Human_MinimumFeatureSet,
  seedSet=1:5
)
readr::write_csv(scoreDF_MHCI_Human, "./Path/To/Your/Directory/MHCI/ScoreDF_MHCI_Human.csv")

# Probability estimation [MHC-II]
scoreDF_MHCII_Human <- Immunogenicity_Score(
  featureDF=featureDF_MHCII[Peptide%in%MHCII_Human$Peptide,],
  metadataDF=MHCII_Human[,.(Peptide, Immunogenicity)],
  featureSet=MHCII_Human_MinimumFeatureSet,
  seedSet=1:5
)
readr::write_csv(scoreDF_MHCII_Human, "./Path/To/Your/Directory/MHCII/ScoreDF_MHCII_Human.csv")
```

4. Epitope prioritization with immunogenicity scores and escape potentials
-  A straitforward wrapper function to compute the two metrics, immunogenicity score and escape potential, for a given set of peptides is provided.  
-  Prediction can be performed as follows.
``` r
# Datasets
fragLibDT <- fst::read_fst("./Path/To/Your/Directory/FragmentLibrary.fst", as.data.table=T)
featureDF_MHCI <- fst::read_fst("./Path/To/Your/Directory/MHCI/FeatureDF_Weighted.10000.fst", as.data.table=T)
featureDF_MHCII <- fst::read_fst("./Path/To/Your/Directory/MHCII/FeatureDF_Weighted.10000.fst", as.data.table=T)

# Prioritization [MHC-I]
res_MHCI <- EpitopePrioritization(
  featureDF=featureDF_MHCI[Peptide%in%MHCI_Human$Peptide,], 
  metadataDF=MHCI_Human[,.(Peptide,Immunogenicity)],
  peptideSet=peptideSet_of_interest,
  fragLib=fragLibDT,
  aaIndexIDSet="all",
  fragLenSet=3:8,
  fragDepthSet=10000,
  fragLibTypeSet="Weighted",
  featureSet=MHCI_Human_MinimumFeatureSet,
  seedSet=1:5,
  coreN=parallel::detectCores(logical=F),        ## parallelization
  outDir="./Output"  ## Intermediate and final output files will be stored under this directory
)

# Prioritization [MHC-II]
res_MHCII <- EpitopePrioritization(
  featureDF=featureDF_MHCII[Peptide%in%MHCII_Human$Peptide,], 
  metadataDF=MHCII_Human[,.(Peptide,Immunogenicity)],
  peptideSet=peptideSet_of_interest,
  fragLib=fragLibDT,
  aaIndexIDSet="all",
  fragLenSet=3:11,
  fragDepthSet=10000,
  fragLibTypeSet="Weighted",
  featureSet=MHCII_Human_MinimumFeatureSet,
  seedSet=1:5,
  coreN=parallel::detectCores(logical=F),
  outDir="./Output"  ## Intermediate and final output files will be stored under this directory
)
```

Reference
------------------------
Ogishi, M and Yotsuyanagi, H. (2019) "Quantitative prediction of the landscape of T cell epitope immunogenicity in sequence space." *Frontiers in Immunology*. https://www.frontiersin.org/articles/10.3389/fimmu.2019.00827/full
