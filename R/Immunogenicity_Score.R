#' Immunogenicity score.
#'
#' \code{Immunogenicity_Score} calculates immunogenicity scores.
#'
#' @param featureDFFileNames_train A character vector of the fst file names of feature dataframes for model training.
#' @param featureDFFileNames_predict A character vector of the fst file names of feature dataframes for predicting immunogenicity scores. If \code{NULL}, model training data was splitted into training and prediction subdatasets.
#' @param metadata_train A dataframe containig "Peptide" and "Immunogenicity" columns for model training.
#' @param featureSet A set of features for model training.
#' @param seedSet A set of random seeds for bootstrapping.
#' @param maxJavaMemory The upper limit of memory for Java virtual machine in megabytes.
#' @param coreN The number of cores to be used for parallelization. Set \code{NULL} to skip parallelization.
#' @importFrom dplyr %>%
#' @param dplyr select
#' @param data.table :=
#' @param data.table rbindlist
#' @param data.table as.data.table
#' @param data.table setcolorder
#' @param parallel detectCores
#' @param pbapply pblapply
#' @param caret preProcess
#' @param extraTrees extraTrees
#' @param fst read_fst
#' @param BBmisc chunk
#' @export
#' @rdname Immunogenicity_Score
#' @name Immunogenicity_Score
Immunogenicity_Score <- function(
  featureDFFileNames_train,
  featureDFFileNames_predict=NULL,
  metadata_train,
  featureSet="all",
  seedSet=1:10,
  maxJavaMemory="6G",
  coreN=parallel::detectCores()
){
  # Working function
  Immunogenicity_Score_Single <- function(
    featureDFFileNames_train,
    featureDFFileNames_predict=NULL,
    metadata_train,
    featureSet="all",
    seed=12345,
    maxJavaMemory="6G",
    coreN=parallel::detectCores()
  ){
    set.seed(seed)
    options(java.parameters=paste0("-Xmx", maxJavaMemory))

    # Data for model training
    leng <- length(featureDFFileNames_train)
    dt_train_list <- lapply(featureDFFileNames_train, fst::read_fst, as.data.table=T)
    if(identical(featureSet, "all")){
      featureSet <- sort(unique(unlist(
        lapply(dt_train_list, function(dt){setdiff(colnames(dt), "Peptide")})
      )))
    }
    dt_train_list <- lapply(dt_train_list, function(dt){dplyr::select(dt, c("Peptide", featureSet))})

    # Data for predicting immunogenicity
    if(is.null(featureDFFileNames_predict)){
      peptideSet <- sort(unique(unlist(
        lapply(dt_train_list, function(dt){dt$"Peptide"})
      )))
      peptideSetList <- BBmisc::chunk(peptideSet, n.chunks=leng, shuffle=T)
      dt_pred_list <- lapply(1:leng, function(i){
        dt <- dt_train_list[[i]]
        dt <- dt[Peptide %in% peptideSetList[[i]],,]
        dt <- dt[,DataSource:=featureDFFileNames_train[i]]
        return(dt)
      })
      dt_train_list <- lapply(1:leng, function(i){
        dt <- dt_train_list[[i]]
        dt[!Peptide %in% peptideSetList[[i]],,]
      })
    }else{
      dt_pred <- data.table::rbindlist(lapply(featureDFFileNames_predict, function(f){
        dt <- fst::read_fst(f, as.data.table=T)
        dt <- dplyr::select(dt, c("Peptide", featureSet))
        dt <- dt[,DataSource:=f]
        return(dt)
      }))
    }

    predDTList <- pbapply::pblapply(
      1:leng,
      function(i){
        set.seed(i)

        ## Preprocessing
        dt_train <- dt_train_list[[i]]
        pp_train <- caret::preProcess(dplyr::select(dt_train, -Peptide), method=c("center", "scale"))
        dt_train <- predict(pp_train, dt_train)

        ## Model training
        trgt <- merge(dt_train, metadata_train, by="Peptide", sort=F)$"Immunogenicity"
        tab <- as.numeric(table(trgt))
        w <- 1/tab[trgt]
        mat <- as.matrix(dplyr::select(dt_train, -Peptide))
        et <- extraTrees::extraTrees(
          x=mat, y=trgt,
          mtry=35, numRandomCuts=2, weights=w,
          numThreads=ifelse(is.null(coreN), 1, coreN)
        )

        ## Prediction
        if(is.null(featureDFFileNames_predict)) dt_pred <- dt_pred_list[[i]]
        mat <- as.matrix(predict(pp_train, dplyr::select(dt_pred, -Peptide, -DataSource)))
        pred <- cbind(
          dplyr::select(dt_pred, Peptide, DataSource),
          data.frame(PredictedImmunogenicity=predict(et, mat, probability=F)),
          as.data.frame(predict(et, mat, probability=T))
        )
        pred <- data.table::as.data.table(pred)

        ## Output
        gc();gc()
        return(pred)
      }
    )

    set.seed(seed)
    predDT <- data.table::rbindlist(predDTList)
    if(is.null(featureDFFileNames_predict)){
      predDT <- merge(predDT, metadata_train, by="Peptide", sort=F)
    }
    predDT[,Seed:=seed]
    return(predDT)
  }

  # Batch workflow
  dt <- data.table::rbindlist(lapply(1:length(seedSet), function(i){
    s <- seedSet[i]
    cat(i, "/", length(seedSet), " | Seed = ", s, "\n", sep="")
    dt <- Immunogenicity_Score_Single(
      featureDFFileNames_train,
      featureDFFileNames_predict,
      metadata_train,
      featureSet,
      seed=s,
      maxJavaMemory,
      coreN
    )
    gc();gc()
    return(dt)
  }))
  if(is.null(featureDFFileNames_predict)){
    data.table::setcolorder(dt, c("Peptide","Immunogenicity","PredictedImmunogenicity","Positive","Negative","Seed","DataSource"))
  }else{
    data.table::setcolorder(dt, c("Peptide","PredictedImmunogenicity","Positive","Negative","Seed","DataSource"))
  }
  return(dt)
}
