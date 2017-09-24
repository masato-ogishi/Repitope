#' Machine learning for immunogenicity prediction
#' 
#' \code{immunogenicityPrediction_MachineLearning.DataFormat} generates machine-learning-ready datasets.\cr
#' \code{immunogenicityPrediction_MachineLearning.Single} performs machine learning with a single algorithm and a single random seed.\cr
#' \code{immunogenicityPrediction_MachineLearning} performs machine learning with multiple combinations of algorithms and random seeds.
#' 
#' @param featureDF A feature dataframe. A "Peptide" column is required.
#' @param featureDFSet A set of feature dataframes.
#' @param metadataDF A metadata dataframe labeled as "Peptide", "Immunogenicity", "MHC", and "Cluster". No extra metadata is allowed (mistaken as features).
#' @param formattedDataList A formatted list of dataframes returned by \code{immunogenicityPrediction_MachineLearning.DataFormat}.
#' @param trainControlOptions A trainControl object.
#' @param mlAlgorithm A string indicating the caret machine learning algorithm.
#' @param mlAlgorithmSet A set of strings indicating the caret machine learning algorithms.
#' @param mlAlgorithmLabel A label for the caret machine learning algorithm.
#' @param mlAlgorithmLabelSet A set of labels for the caret machine learning algorithms.
#' @param seed A random seed.
#' @param seedSet A set of random seeds.
#' @param outDir A directory for exporting the results.
#' @importFrom dplyr %>%
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr bind_rows
#' @importFrom tidyr drop_na
#' @importFrom DescTools Sort
#' @importFrom caret createDataPartition
#' @importFrom caret preProcess
#' @importFrom caret train
#' @importFrom caret upSample
#' @importFrom caret downSample
#' @importFrom caret confusionMatrix
#' @importFrom classifierplots classifierplots_folder
#' @export
#' @rdname immunogenicityPrediction_MachineLearning
#' @name immunogenicityPrediction_MachineLearning
immunogenicityPrediction_MachineLearning.DataFormat <- function(featureDF, metadataDF, seed=12345){
  # Combine metadata
  df <- dplyr::left_join(metadataDF, featureDF, by="Peptide")
  
  # Randomly choose one epitope per each cluster
  set.seed(seed)
  ML_Data <- df %>%
    (function(d){d[sample.int(nrow(d)),]}) %>%
    dplyr::distinct(Cluster, .keep_all=T) %>%
    dplyr::distinct(Peptide, .keep_all=T) %>%
    DescTools::Sort(ord=c("Peptide", "Cluster")) %>%
    dplyr::mutate(Immunogenicity=factor(Immunogenicity, levels=c("Positive", "Negative")))

  # Random data splitting
  set.seed(seed)
  trainID <- caret::createDataPartition(ML_Data$"Immunogenicity", p=7/10, list=F)
  ML_Train <- ML_Data[trainID,] ## Proportion == 0.7
  ML_Test <- ML_Data[-trainID,]
  testID <- caret::createDataPartition(ML_Test$"Immunogenicity", p=2/3, list=F)
  ML_Test <- ML_Test[testID,]   ## Proportion == 0.2
  ML_Valid <- ML_Test[-testID,] ## Proportion == 0.1

  # Pre-processing
  ML_Train_Features <- dplyr::select(ML_Train, -Peptide, -Immunogenicity, -MHC, -Cluster)
  ML_Train_Features <- predict(caret::preProcess(ML_Train_Features, method=c("nzv", "corr")), ML_Train_Features)
  featureSet <- colnames(ML_Train_Features)
  pp_Train <- caret::preProcess(ML_Train_Features, method=c("center", "scale"))
  ML_Train <- predict(pp_Train, dplyr::select(ML_Train, Peptide, Immunogenicity, MHC, Cluster, dplyr::one_of(featureSet)))
  ML_Test <- predict(pp_Train, dplyr::select(ML_Test, Peptide, Immunogenicity, MHC, Cluster, dplyr::one_of(featureSet)))
  ML_Valid <- predict(pp_Train, dplyr::select(ML_Valid, Peptide, Immunogenicity, MHC, Cluster, dplyr::one_of(featureSet)))

  # Output
  list("ML_Data"=ML_Data, "ML_Train"=ML_Train, "ML_Test"=ML_Test, "ML_Valid"=ML_Valid)
}
#' @export
#' @rdname immunogenicityPrediction_MachineLearning
#' @name immunogenicityPrediction_MachineLearning
immunogenicityPrediction_MachineLearning.Single <- function(
  formattedDataList, trainControlOptions,
  mlAlgorithm, mlAlgorithmLabel, seed=12345,
  outDir="./classifiers/"
){
  # Train
  set.seed(seed)
  sink(tempfile())
  model <- try(caret::train(Immunogenicity~.,
                            dplyr::select(formattedDataList$"ML_Train", -Peptide, -MHC, -Cluster),
                            method=mlAlgorithm,
                            trControl=trainControlOptions),
               silent=T)
  sink()
  if(class(model)=="try-error"){
    print("Model training process encountered an unexpected error!")
    return("Error")
  }
  cm <- caret::confusionMatrix(predict(model, formattedDataList$"ML_Test"), formattedDataList$"ML_Test"$"Immunogenicity")

  # Save
  out <- paste0(outDir, "/", mlAlgorithmLabel, "_Seed", seed)
  dir.create(out, showWarnings=F, recursive=T)
  dir.create(paste0(out, "/Classifierplots"), showWarnings=F, recursive=T)
  labels <- as.character(formattedDataList$"ML_Test"$"Immunogenicity")
  labels[labels=="Positive"] <- 1
  labels[labels=="Negative"] <- 0
  labels <- as.numeric(labels)
  classifierplots::classifierplots_folder(
    labels, predict(model, formattedDataList$"ML_Test", type="prob")[[1]],
    paste0(out, "/Classifierplots")
  )
  saveRDS(model, paste0(out, "/Classifier.rds"))
  sink(paste0(out, "/ConfusionMatrix.txt"))
  print(cm)
  sink()
  gc();gc()

  # Output
  list("ML_Model"=model, "ML_ConfMat"=cm)
}
#' @export
#' @rdname immunogenicityPrediction_MachineLearning
#' @name immunogenicityPrediction_MachineLearning
immunogenicityPrediction_MachineLearning <- function(
  featureDFSet, metadataDF, trainControlOptions,
  mlAlgorithmSet, mlAlgorithmLabelSet, seedSet=1:5,
  outDir="./classifiers/"
){
  cmSummaryList <- as.list(numeric(length(seedSet)))
  for(i in 1:length(seedSet)){
    set.seed(seedSet[[i]])
    data <- immunogenicityPrediction_MachineLearning.DataFormat(featureDFSet[[i]], metadataDF, seedSet[[i]])
    dir.create(paste0(outDir, "/"), showWarnings=F, recursive=T)
    saveRDS(data, paste0(outDir, "/Data_Seed", seedSet[[i]], ".rds"))
    mlList <- mapply(
      function(x, y){
        immunogenicityPrediction_MachineLearning.Single(
          data, trainControlOptions,
          mlAlgorithm=x, mlAlgorithmLabel=y, seed=seedSet[[i]], outDir
        )
      }, mlAlgorithmSet, mlAlgorithmLabelSet, SIMPLIFY=F)
    errors <- mlList=="Error"
    mlList <- mlList[!errors]
    cmList <- lapply(mlList, function(ml){ml$"ML_ConfMat"})
    modelList <- lapply(mlList, function(ml){ml$"ML_Model"})
    cmSummary <- as.data.frame(sapply(cmList, function(cm){c(cm$"overall", cm$"byClass")}))
    colnames(cmSummary) <- mlAlgorithmLabelSet[!errors]
    cmSummary[["Metric"]] <- rownames(cmSummary)
    cmSummary[["RandomSeed.ML"]] <- seedSet[[i]]
    cmSummary <- cmSummary %>% 
      tidyr::gather(Algorithm, Value, -Metric, -RandomSeed.ML) %>% 
      dplyr::select(RandomSeed.ML, Algorithm, Metric, Value) %>%
      tidyr::spread(Algorithm, Value)
    cmSummaryList[[i]] <- cmSummary
  }
  cmSummary <- dplyr::bind_rows(cmSummaryList)
  write.csv(cmSummary, paste0(outDir, "/SummaryTable.csv"), row.names=F)
  list("ML_ModelList"=modelList, "ML_ConfMat_SummaryTable"=cmSummary)
}
