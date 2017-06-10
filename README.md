Repitope: Epitope prediction via repertoire-wide TCR-peptide contact profiles
===============================================

The 'Repitope' package provides some easy-to-use functions for conducting epitope prediction analysis via repertoire-wide TCR-peptide contact profiles.

Installation and loading
------------------------

-   Install the latest version from [GitHub](https://github.com/masato-ogishi/Repitope) as follows:

``` r
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
```

Distribution
------------

``` r
library(Repitope)
# Example code 1
# :::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(1234)
wdata = data.frame(
   sex = factor(rep(c("F", "M"), each=200)),
   weight = c(rnorm(200, 55), rnorm(200, 58)))
head(wdata, 4)
#>   sex   weight
#> 1   F 53.79293
#> 2   F 55.27743
#> 3   F 56.08444
#> 4   F 52.65430

# Example code 2
# :::::::::::::::::::::::::::::::::::::::::::::::::::
# Change outline and fill colors by groups ("sex")
# Use custom palette
ggdensity(wdata, x = "weight",
   add = "mean", rug = TRUE,
   color = "sex", fill = "sex",
   palette = c("#00AFBB", "#E7B800"))
```
