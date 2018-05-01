Repitope: Epitope immunogenicity prediction from interresidue contact potential profiling
===============================================

The 'Repitope' package provides a structured framework of quantitative prediction of the immunogenicity of MHC-I-loaded peptides through interresidue contact potential profiling.
Note: Repitope is currently under construction.

Installation
------------------------

-   Install the latest version from [GitHub](https://github.com/masato-ogishi/Repitope) as follows:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/plotUtility")
devtools::install_github("masato-ogishi/Repitope")
```

-   You might be prompted to install some packages before installling Repitope. Follow the message(s).

Loading
------------------

``` r
options(java.parameters="-Xmx60G")  ## allow the JAVA session to use larger memory space
library(tidyverse)
library(Repitope)
```

Reference
------------------------

Ogishi, M and Yotsuyanagi, H. (2018) "Immunogenicity landscape of cytotoxic T lymphocyte epitopes in sequence space." (Manuscript in preparation)
