#' Immunogenicity score.
#'
#' @param featureDFFileNames_train A character vector of the fst file names of feature dataframes for model training.
#' @param featureDFFileNames_predict A character vector of the fst file names of feature dataframes for predicting immunogenicity scores. If the length equals to that of \code{featureDFFileNames_train}, each dataframe will be subject to the prediction by the classifier trained from the corresponding training dataframe. Otherwise, they will be subject to the prediction in an all-by-all manner. If \code{NULL}, model training dataframe will be splitted into training and prediction dataframes.
#' @param metadata_train A dataframe containig "Peptide" and "Immunogenicity" columns for model training.
#' @param featureSet A set of features for model training.
#' @param seedSet A set of random seeds for bootstrapping. Can be a single value when \code{featureDFFileNames_predict} is provided.
#' @param coreN The number of cores to be used for parallelization. Set \code{NULL} to skip parallelization.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom data.table :=
#' @importFrom data.table rbindlist
#' @importFrom data.table as.data.table
#' @importFrom data.table setcolorder
#' @importFrom parallel detectCores
#' @importFrom pbapply pblapply
#' @importFrom caret preProcess
#' @importFrom extraTrees extraTrees
#' @importFrom fst read_fst
#' @importFrom BBmisc chunk
#' @export
#' @rdname Immunogenicity_Score
#' @name Immunogenicity_Score
Immunogenicity_Score <- function(
  featureDFFileNames_train,
  featureDFFileNames_predict=NULL,
  metadata_train,
  featureSet="all",
  seedSet=1:10,
  coreN=parallel::detectCores()
){
  # Working functions
  dt_lists <- function(
    featureDFFileNames_train,
    featureDFFileNames_predict=NULL,
    featureSet="all",
    seed=12345
  ){
    set.seed(seed)

    ## Data for model training
    leng <- length(featureDFFileNames_train)
    dt_train_list <- lapply(featureDFFileNames_train, fst::read_fst, as.data.table=T)
    if(identical(featureSet, "all")){
      featureSet <- sort(unique(unlist(
        lapply(dt_train_list, function(dt){setdiff(colnames(dt), c("Peptide", "Immunogenicity", "DataType", "Cluster"))})
      )))
    }
    dt_train_list <- lapply(1:leng, function(i){
      dt_train_list[[i]] %>%
        dplyr::select(c("Peptide", featureSet)) %>%
        dplyr::mutate(DataSource_Train=featureDFFileNames_train[i],
                      Seed=NA) %>%
        data.table::as.data.table()
    })

    ## Data for predicting immunogenicity
    if(is.null(featureDFFileNames_predict)){
      peptideSet <- sort(unique(unlist(
        lapply(dt_train_list, function(dt){dt$"Peptide"})
      )))
      peptideSetList <- BBmisc::chunk(peptideSet, n.chunks=leng, shuffle=T)  ## requires a random seed
      dt_pred_list <- lapply(1:leng, function(i){
        dt_train_list[[i]] %>%
          dplyr::filter(Peptide %in% peptideSetList[[i]]) %>%
          dplyr::mutate(DataSource_Pred=featureDFFileNames_train[i],
                        Seed=seed) %>%
          data.table::as.data.table()
      })
      dt_train_list <- lapply(1:leng, function(i){
        dt_train_list[[i]] %>%
          dplyr::filter(!Peptide %in% peptideSetList[[i]]) %>%
          dplyr::mutate(DataSource_Train=featureDFFileNames_train[i],
                        Seed=seed) %>%
          data.table::as.data.table()
      })
    }else{
      dt_pred_list <- lapply(featureDFFileNames_predict, function(f){
        fst::read_fst(f) %>%
          dplyr::select(c("Peptide", featureSet)) %>%
          dplyr::mutate(DataSource_Pred=f,
                        Seed=NA) %>%
          data.table::as.data.table()
      })
    }
    gc();gc()
    list(dt_train_list, dt_pred_list)
  }

  immscore <- function(dt_train, dt_pred){
    ## Preprocessing
    pp_train <- caret::preProcess(dplyr::select(dt_train, featureSet), method=c("center", "scale"))
    dt_train <- predict(pp_train, dt_train)

    ## Model training
    trgt <- merge(dt_train, metadata_train, by="Peptide", sort=F)$"Immunogenicity"
    tab <- as.numeric(table(trgt))
    w <- 1/tab[trgt]
    mat <- as.matrix(dplyr::select(dt_train, featureSet))
    et <- extraTrees::extraTrees(
      x=mat, y=trgt,
      mtry=35, numRandomCuts=2, weights=w,
      numThreads=ifelse(is.null(coreN), 1, coreN)
    )

    ## Prediction
    mat <- as.matrix(predict(pp_train, dplyr::select(dt_pred, featureSet)))
    predDT <- cbind(
      dplyr::select(dt_pred, Peptide, DataSource_Pred, Seed),
      data.frame(PredictedImmunogenicity=predict(et, mat, probability=F)),
      as.data.frame(predict(et, mat, probability=T))
    ) %>%
      dplyr::mutate(DataSource_Train=unique(dt_train$DataSource_Train)) %>%
      data.table::as.data.table()
    data.table::setcolorder(predDT, c("Peptide","PredictedImmunogenicity","Positive","Negative","DataSource_Train","DataSource_Pred","Seed"))

    ## Output
    gc();gc()
    return(predDT)
  }

  immscore_batch <- function(dt_train_list, dt_pred_list){
    if(length(dt_train_list)==length(dt_pred_list)){
      predDTList <- pbapply::pblapply(
        1:length(dt_train_list),
        function(i){immscore(dt_train_list[[i]], dt_pred_list[[i]])}
      )
    }else{
      comb <- expand.grid(1:length(dt_train_list), 1:length(dt_pred_list))
      predDTList <- pbapply::pblapply(
        1:nrow(comb),
        function(i){immscore(dt_train_list[[comb[i, 1]]], dt_pred_list[[comb[i, 2]]])}
      )
    }
    predDT <- data.table::rbindlist(predDTList)
    if(is.null(featureDFFileNames_predict)){
      predDT <- merge(predDT, metadata_train, by="Peptide", sort=F)
      data.table::setcolorder(predDT, c("Peptide","Immunogenicity","PredictedImmunogenicity","Positive","Negative","DataSource_Train","DataSource_Pred","Seed"))
    }
    return(predDT)
  }

  # Batch analysis
  immScoreDT <- data.table::rbindlist(lapply(1:length(seedSet), function(i){
    s <- seedSet[i]
    cat(i, "/", length(seedSet), " | Seed = ", s, "\n", sep="")
    dts <- dt_lists(
      featureDFFileNames_train,
      featureDFFileNames_predict,
      featureSet,
      seed=s
    )
    immScoreDT <- immscore_batch(dts[[1]], dts[[2]])
    rm(list=setdiff(ls(), "immScoreDT"))
    gc();gc()
    return(immScoreDT)
  }))

  # Output
  rm(list=setdiff(ls(), "immScoreDT"))
  gc();gc()
  return(immScoreDT)
}
