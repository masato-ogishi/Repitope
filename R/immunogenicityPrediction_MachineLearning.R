#' Machine learning for immunogenicity prediction
#' 
#' \code{immunogenicityPrediction_MachineLearning.DataFormat} generates machine-learning-ready datasets.\cr
#' \code{immunogenicityPrediction_MachineLearning.Single} performs machine learning with a single algorithm and a single random seed.\cr
#' \code{immunogenicityPrediction_MachineLearning} performs machine learning with multiple combinations of algorithms and random seeds.
#' 
#' @param featureDF A feature dataframe. A "Peptide" column is required.
#' @param featureDFSet A list of feature dataframes. The names will be used as parameter sets by default.
#' @param metadataDF A metadata dataframe labeled as "Peptide", "Immunogenicity", and "Cluster". No extra metadata is allowed. Feature dataframes are downsized based on the "Peptide" column in this metadata. In other words, peptides not contained in the metadata dataframe are omitted.
#' @param formattedDataList A formatted list of dataframes returned by \code{immunogenicityPrediction_MachineLearning.DataFormat}.
#' @param trainControlOptions A trainControl object.
#' @param weighted Logical. Whether class imbalance should be adjusted by class-specific weighting.
#' @param mlAlgorithm A string indicating the caret machine learning algorithm.
#' @param mlAlgorithmSet A set of strings indicating the caret machine learning algorithms.
#' @param mlAlgorithmLabel A label for the caret machine learning algorithm.
#' @param mlAlgorithmLabelSet A set of labels for the caret machine learning algorithms.
#' @param seed A random seed.
#' @param parameterSet A set of parameters to generate features. Ex. "global-local.3.10.1" (AlignmentMethod.FragLen.Depth.Seed). A random seed value is reused for machine learning.
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
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pbapply
#' @importFrom pbapply pblapply
#' @export
#' @rdname immunogenicityPrediction_MachineLearning
#' @name immunogenicityPrediction_MachineLearning
immunogenicityPrediction_MachineLearning.DataFormat <- function(featureDF, metadataDF, seed=12345){
  # Combine metadata
  df <- dplyr::left_join(metadataDF, featureDF, by="Peptide") %>%
    tidyr::drop_na() %>%
    dplyr::mutate(Immunogenicity=factor(Immunogenicity, levels=c("Positive", "Negative")))
  
  # Randomly choose one epitope per each cluster [Minor class is prioritized]
  set.seed(seed)
  outcome_table <- table(df$"Immunogenicity")
  if(outcome_table["Positive"]>=outcome_table["Negative"]){ 
    outcome_minor_major <- c("Negative", "Positive")
  }else{ 
    outcome_minor_major <- c("Positive", "Negative")
  }
  ML_Data <- df %>%
    dplyr::mutate(Immunogenicity=factor(Immunogenicity, levels=outcome_minor_major)) %>%
    (function(d){d[sample.int(nrow(d)),]}) %>%
    dplyr::distinct(Immunogenicity, Cluster, .keep_all=T) %>%
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
  ML_Train_Features <- dplyr::select(ML_Train, -Peptide, -Immunogenicity, -Cluster)
  ML_Train_Features <- predict(caret::preProcess(ML_Train_Features, method=c("nzv", "corr")), ML_Train_Features)
  featureSet <- colnames(ML_Train_Features)
  pp_Train <- caret::preProcess(ML_Train_Features, method=c("center", "scale"))
  
  ML_Data <- predict(pp_Train, dplyr::select(ML_Data, Peptide, Immunogenicity, Cluster, dplyr::one_of(featureSet)))
  ML_Train <- predict(pp_Train, dplyr::select(ML_Train, Peptide, Immunogenicity, Cluster, dplyr::one_of(featureSet)))
  ML_Test <- predict(pp_Train, dplyr::select(ML_Test, Peptide, Immunogenicity, Cluster, dplyr::one_of(featureSet)))
  ML_Valid <- predict(pp_Train, dplyr::select(ML_Valid, Peptide, Immunogenicity, Cluster, dplyr::one_of(featureSet)))

  # Output
  list("ML_Data"=ML_Data, "ML_Train"=ML_Train, "ML_Test"=ML_Test, "ML_Valid"=ML_Valid)
}
#' @export
#' @rdname immunogenicityPrediction_MachineLearning
#' @name immunogenicityPrediction_MachineLearning
immunogenicityPrediction_MachineLearning.TrainAndEvaluate <- function(
  formattedDataList, trainControlOptions, weighted=F,
  mlAlgorithm, mlAlgorithmLabel, 
  seed=12345,
  outDir="./classifiers/"
){
  # (Optional) Weighting for class imbalance
  if(weighted==T){
    outcome_table <- table(formattedDataList$"ML_Train"$"Immunogenicity")
    weights <- formattedDataList$"ML_Train"$"Immunogenicity"
    levels(weights) <- c(as.numeric(outcome_table[names(outcome_table)[2]]/nrow(formattedDataList$"ML_Train")), 
                         as.numeric(outcome_table[names(outcome_table)[1]]/nrow(formattedDataList$"ML_Train")))
    weights <- as.numeric(as.vector(weights))
    trainControlOptions$weights <- weights
  }
  
  # Train
  set.seed(seed)
  sink(tempfile())
  model <- try(caret::train(
    Immunogenicity~.,
    dplyr::select(formattedDataList$"ML_Train", -Peptide, -Cluster),
    method=mlAlgorithm,
    trControl=trainControlOptions,
    metric="Kappa"
  ), silent=T)
  sink()
  if(class(model)=="try-error"){
    print("Model training process encountered an unexpected error!")
    return("Error")
  }
  
  # Evaluate & Save
  out <- paste0(outDir, "/", mlAlgorithmLabel)
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
  cm <- caret::confusionMatrix(predict(model, formattedDataList$"ML_Test"), formattedDataList$"ML_Test"$"Immunogenicity")
  sink(paste0(out, "/ConfusionMatrix.txt"))
  print(cm)
  sink()

  # Output
  gc();gc()
  return(cm)
}
#' @export
#' @rdname immunogenicityPrediction_MachineLearning
#' @name immunogenicityPrediction_MachineLearning
immunogenicityPrediction_MachineLearning.Single <- function(
  featureDF, metadataDF, trainControlOptions, weighted=F,
  mlAlgorithmSet, mlAlgorithmLabelSet, 
  seed=12345,
  outDir="./classifiers/"
){
  dir.create(paste0(outDir, "/"), showWarnings=F, recursive=T)
  set.seed(seed)
  mlData <- immunogenicityPrediction_MachineLearning.DataFormat(featureDF, metadataDF, seed)
  saveRDS(mlData, paste0(outDir, "/FormattedData.rds"))
  cmList <- mapply(
    function(x, y){
      immunogenicityPrediction_MachineLearning.TrainAndEvaluate(
        mlData, trainControlOptions, weighted,
        mlAlgorithm=x, mlAlgorithmLabel=y, 
        seed=seed, outDir=outDir
      )
    }, mlAlgorithmSet, mlAlgorithmLabelSet, SIMPLIFY=F)
  errors <- cmList=="Error"
  cmList <- cmList[!errors]
  cmSummary <- as.data.frame(sapply(cmList, function(cm){c(cm$"overall", cm$"byClass")}))
  colnames(cmSummary) <- mlAlgorithmLabelSet[!errors]
  cmSummary[["Metric"]] <- rownames(cmSummary)
  cmSummary[["Seed"]] <- seed
  cmSummary <- cmSummary %>% 
    tidyr::gather(Algorithm, Value, -Metric, -Seed) %>% 
    dplyr::select(Seed, Algorithm, Metric, Value) %>%
    tidyr::spread(Algorithm, Value)
  return(cmSummary)
}
#' @export
#' @rdname immunogenicityPrediction_MachineLearning
#' @name immunogenicityPrediction_MachineLearning
immunogenicityPrediction_MachineLearning <- function(
  featureDFSet, metadataDF, trainControlOptions, weighted=F,
  mlAlgorithmSet, mlAlgorithmLabelSet, 
  parameterSet=names(featureDFSet),
  outDir="./classifiers/"
){
  time.start <- proc.time()
  
  # Reconstructing parameter sets
  parameterDF <- as.data.frame(stringr::str_split(parameterSet, stringr::fixed("."), simplify=T))
  colnames(parameterDF) <- c("Alignment", "FragLen", "Depth", "Seed")
  parameterDF$"Alignment" <- as.character(parameterDF$"Alignment")
  parameterDF$"FragLen" <- as.character(parameterDF$"FragLen")
  parameterDF$"Depth" <- as.character(parameterDF$"Depth")
  parameterDF$"Seed" <- as.numeric(as.character(parameterDF$"Seed"))
  
  # Non-parallelized machine learning
  message(paste0("Number of parameter combinations = ", nrow(parameterDF)))
  message(paste0("Non-parallelized machine learning was started. (Memory occupied = ", memory.size(), "[Mb])"))
  cmSummary <- dplyr::bind_rows(
    lapply(1:nrow(parameterDF), 
           function(i){
             immunogenicityPrediction_MachineLearning.Single(
               featureDFSet[[i]], metadataDF, trainControlOptions, weighted,
               mlAlgorithmSet, mlAlgorithmLabelSet,
               seed=parameterDF$"Seed"[[i]],
               outDir=file.path(outDir, paste0(c("", "FranLen", "Depth", "Seed"), unlist(parameterDF[i,]), collapse="_"))
             )
           }
          )
  )
  write.csv(cmSummary, paste0(outDir, "/SummaryTable.csv"), row.names=F)
  message(paste0("Non-parallelized machine learning was finished. (Memory occupied = ", memory.size(), "[Mb])"))
  gc();gc();
  
  # Output
  time.end <- proc.time()
  message(paste0("Overall time required = ", (time.end-time.start)[3], "[sec]"))
  return(cmSummary)
}
