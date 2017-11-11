#' Learning curve analysis for the caret machine learning pipeline
#' 
#' @param trainingData.preProc A pre-processed dataframe for model training.
#' @param testingData.preProc A pre-processed dataframe for model testing.
#' @param outcomeLabelName A string for the outcome label name.
#' @param mlAlgorithm A string indicating the caret machine learning algorithm.
#' @param trainControlOptions A trainControl object.
#' @param weighted Logical. Whether class imbalance should be adjusted by class-specific weighting.
#' @param proportionSet A numeric vector indicating the set of proportions for learning curve analysis.
#' @param seed A random seed.
#' @param ... Arguments passed to \code{caret::train}
#' @importFrom caret train
#' @importFrom caret postResample
#' @importFrom ggsci scale_color_npg
#' @import ggplot2
#' @export
#' @rdname learningCurve
#' @name learningCurve
learningCurve <- function(trainingData.preProc, testingData.preProc, 
                          outcomeLabelName, mlAlgorithm="xgbLinear", trainControlOptions, weighted=F, 
                          proportionSet=(1:10)/10, seed=12345) {
  set.seed(seed)
  proportionSet <- sort(unique(proportionSet))
  resampled <- vector(mode="list", length=length(proportionSet))
  tested <- vector(mode="list", length=length(proportionSet))
  apparent <- vector(mode="list", length=length(proportionSet))
  for(i in seq(along=proportionSet)) {
    datasize <- floor(nrow(trainingData.preProc)*proportionSet[i])
    cat("Training for ", round(proportionSet[i]*100, 1), "% (n = ", datasize, ")\n", sep="")
    if(proportionSet[i]<1){
      in_mod <- sample(1:nrow(trainingData.preProc), datasize)
    }else{
      in_mod <- 1:nrow(trainingData.preProc)
    }
    trainingData.sub <- trainingData.preProc[in_mod,]
    if(weighted==T){
      outcome_table <- table(trainingData.sub[[outcomeLabelName]])
      weights <- trainingData.sub[[outcomeLabelName]]
      levels(weights) <- c(as.numeric(outcome_table[names(outcome_table)[2]]/nrow(trainingData.sub)), 
                           as.numeric(outcome_table[names(outcome_table)[1]]/nrow(trainingData.sub)))
      weights <- as.numeric(as.vector(weights))
      trainControlOptions$weights <- weights
    }
    form <- as.formula(paste0(outcomeLabelName, "~."))
    mod <- caret::train(form, trainingData.sub, method=mlAlgorithm, trControl=trainControlOptions, metric="Kappa")
    if(i==1) perf_names <- mod$perfNames
    
    # Resampling
    resampled[[i]] <- merge(mod$resample, mod$bestTune)
    resampled[[i]]$Training_Size <- length(in_mod)
    
    # Accuracy in training set
    app_perf <- caret::postResample(
      predict(mod, newdata=trainingData.sub, type="raw"), 
      trainingData.sub[[outcomeLabelName]]
    )
    app_perf <- as.data.frame(t(app_perf))
    app_perf$Training_Size <- length(in_mod)    
    apparent[[i]] <- app_perf
    
    # Accuracy in testing set
    test_perf <- caret::postResample(
      predict(mod, newdata=testingData.preProc, type="raw"), 
      testingData.preProc[[outcomeLabelName]]
    )
    test_perf <- as.data.frame(t(test_perf))
    test_perf$Training_Size <- length(in_mod)
    tested[[i]] <- test_perf
    
    cat("-Resampling: ", mean(resampled[[i]]$Accuracy), "\n",
        "-Training: ", mean(apparent[[i]]$Accuracy), "\n", 
        "-Testing: ", mean(tested[[i]]$Accuracy), "\n",
        sep="")
  }
  
  # Summarize
  resampled <- do.call("rbind", resampled)
  resampled <- resampled[, c(perf_names, "Training_Size")]
  resampled$Type <- "Resampling"
  apparent <- do.call("rbind", apparent)
  apparent <- apparent[, c(perf_names, "Training_Size")]
  apparent$Type <- "Training"
  out <- rbind(resampled, apparent)
  tested <- do.call("rbind", tested)
  tested <- tested[, c(perf_names, "Training_Size")]
  tested$Type <- "Testing"
  out <- rbind(out, tested)
  out[["Type"]] <- factor(out[["Type"]], levels=c("Training","Testing","Resampling"))
  print(out)
  
  # Visualize
  lc.plot <- ggplot(out, aes_string(x="Training_Size", y="Kappa", color="Type")) + 
    geom_smooth(method=loess, span=.8) + 
    ggsci::scale_color_npg() +
    theme_Publication()
  
  return(list("learningCurve.DF"=out, "learningCurve.Plot"=lc.plot))
}

## https://github.com/tobigithub/caret-machine-learning/blob/master/caret-parallel/learning-curve-plots-caret-parallel.R
