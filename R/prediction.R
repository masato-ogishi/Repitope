#' @title Predict immunogenicity using pre-trained classifiers
#' @description Predict immunogenicity from the rTPCP/mrTPCP matrix calculated from user-provided peptide sequences using pre-trained classifiers
#' @param variableDataFrame A dataframe containing either rTPCP or mrTPCP variables. Currently, only MIYS990106-derived variables with window sizes of 4 and 5 are supported. In case of mrTPCP, currently only nonapeptides are supported.
#' @param variableType Either "rTPCP" or "mrTPCP".
#' @importFrom dplyr bind_cols
#' @importFrom stats predict
#' @importFrom caret predict.train
#' @export
predict_PreTrained <- function(variableDataFrame, variableType){
  if(variableType=="rTPCP"){
    model.preProc <- rTPCP_MIYS990106_svmPoly_34512_PreProcessModel
    model <- rTPCP_MIYS990106_svmPoly_34512
  }
  if(variableType=="mrTPCP"){
    model.preProc <- mrTPCP_MIYS990106_svmPoly_12345_PreProcessModel
    model <- mrTPCP_MIYS990106_svmPoly_12345
  }
  df.pred <- predict(model.preProc, newdata=variableDataFrame)
  df.pred <- dplyr::bind_cols(
    predict(model, newdata=df.pred, type="raw"),
    predict(model, newdata=df.pred, type="prob")
  )
  return(df.pred)
}
