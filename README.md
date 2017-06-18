Repitope: Epitope prediction via repertoire-wide TCR-peptide contact profiles
===============================================

The 'Repitope' package provides some easy-to-use functions for conducting epitope prediction analysis via repertoire-wide TCR-peptide contact profiles.

Installation
------------------------

-   Install the latest version from [GitHub](https://github.com/masato-ogishi/Repitope) as follows:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/Repitope")
```

-   When you encounter errors, error might be resolved after installing and updating dependent packages.

Loading
------------------

``` r
# Loading the Repitope package
library(Repitope)

# Datasets included
data(Chowell)
data(HCV)
data(HIV)
data(WellEstablishedNeoepitopes)
data(Stronen_Best)
data(Stronen_All)
data(Calis)
data(Rizvi)
data(vanAllen)
data(TCGA)
```

Usage
-----------------------------------
``` r
# Calculating rTPCPs for a given set of input peptides
# :::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(12345)
code <- code
# :::::::::::::::::::::::::::::::::::::::::::::::::::

# Calculating mrTPCPs for a given set of input peptides[9-mers]
# :::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(12345)
code <- code
# :::::::::::::::::::::::::::::::::::::::::::::::::::

# Machine learning
# :::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(12345)
code <- code
# :::::::::::::::::::::::::::::::::::::::::::::::::::

# Calculating rTPCPs/mrTPCPs with a given TCR repertoire
# :::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(12345)
code <- code
# :::::::::::::::::::::::::::::::::::::::::::::::::::

# Predicting immunogenicity using the pre-trained support vector machine classifier
# :::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(12345)
code <- code
# :::::::::::::::::::::::::::::::::::::::::::::::::::

```
