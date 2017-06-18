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
# Setting up the working environment
set.seed(12345)
library(tidyverse)
library(ggsci)
library(caret)

# Enable parallel computing
library(doParallel)
cl <- makePSOCKcluster(7)
registerDoParallel(cl)

# The dataset originally constructed by Chowell et al.
data(Chowell)

# Calculating rTPCPs for a given set of input peptides
df.rTPCP.short <- rTPCP(Chowell$Peptide[1:10]) ## Choose the first ten peptide to reduce calculation burden

# Calculating mrTPCPs for a given set of input peptides[9-mers]
df.chowell <- Chowell %>% 
  dplyr::filter(nchar(Peptide)==9) %>%
  dplyr::distinct(Peptide, .keep_all=T) %>% 
  DescTools::Sort(ord="Peptide", decreasing=F)
df.mrTPCP <- dplyr::bind_cols(
  dplyr::select(df.chowell, Immunogenicity), 
  mrTPCP(df.chowell$Peptide, winds = 4:5)
) %>%
dplyr::mutate(Immunogenicity=factor(Immunogenicity, levels=c("Positive", "Negative")))

# Machine learning
trainIndex <- createDataPartition(df.mrTPCP$"Immunogenicity", p=2/3, list=F, times=1)
df.train <- df.mrTPCP[trainIndex,]
df.test <- df.mrTPCP[-trainIndex,]
table(df.train$Immunogenicity)
table(df.test$Immunogenicity)
ml.options <- trainControl(
  ## 10-fold CV
  method="repeatedcv",
  number=10,
  ## repeated ten times
  repeats=10,
  savePredictions="final", 
  verboseIter=T, 
  classProbs=T
)
ml.model.preProc <- caret::preProcess(df.train, method = c("center", "scale"))
ml.model <- caret::train(Immunogenicity~., data=predict(ml.model.preProc, newdata=df.train), method="svmPoly", trControl=ml.options)
ml.pred <- predict(ml.model, predict(ml.model.preProc, newdata=df.test))
caret::confusionMatrix(ml.pred, df.test$Immunogenicity)

# Predicting immunogenicity using the pre-trained support vector machine classifier
# :::::::::::::::::::::::::::::::::::::::::::::::::::
ml.pred <- predict_PreTrained(df.mrTPCP, variableType="mrTPCP")
caret::confusionMatrix(ml.pred, df.mrTPCP$Immunogenicity)
# :::::::::::::::::::::::::::::::::::::::::::::::::::

```
