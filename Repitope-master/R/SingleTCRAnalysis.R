#' Single TCR-peptide contact potential profiling.
#'
#' @param weightStringSet A set of contact weights. E.g. "1|2|3|2|1"
#' @param peptideSet A set of peptide sequences.
#' @param tcrSet A set of corresponding TCR sequences.
#' @param aaIndexIDSet A set of AAIndex IDs indicating the AACP scales to be used. Set "all" to shortcut the selection of all available AACP scales.
#' @param fragLenSet A set of sliding window sizes. Must be between 3 and 8.
#' @param seed A random seed.
#' @param coreN The number of cores to be used for parallelization.
#' @param dt_deltaG_sCPP A data.table containing the "DeltaG" and sCPP feature columns.
#' @param dt_deltaG_sCPP_outlier Optional. If provided, it would be used as outliers.
#' @param dt_deltaG_sCPP_new Optional. If provided, TCR affnities would be predicted for this dataset. Otherwise, predicted TCR affinities for the original dataset would be returned.
#' @param corVarName A variable name used for univariate analysis.
#' @param corVarNames A set of variable names used for multivariate regression.
#' @param dt_univ An output data.table of univariate analysis.
#' @param sig A siginificance cutoff for univariate feature selection.
#' @param bestOnly Logical. If True, only the best sCPP feature per AAIndex will be returned. If False, all sCPP features derived from selected AAIndices will be returned.
#' @export
#' @rdname SingleTCRAnalysis
#' @name contactFootprintDensityPlot
contactFootprintDensityPlot <- function(weightStringSet){
  bs <- lapply(strsplit(weightStringSet, "|", fixed=T), as.numeric)
  bs <- bs[!is.na(bs)]
  bs <- lapply(bs, function(w){
    p <- 1:length(w)-weighted.mean(1:length(w), w=w)
    mapply(function(p, w){rep(p, w)}, p, w, SIMPLIFY=F)
  })
  bs <- unlist(bs)
  bs <- bs[!is.na(bs)]
  propFormat <- function(numericalVector){
    paste0(formatC(numericalVector, format="f", digits=1), "%")
  }
  proportion.1 <- sum(bs>=-1 & bs<=1)/length(bs)*100
  proportion.2 <- sum(bs>=-2 & bs<=2)/length(bs)*100 - proportion.1
  proportion.3 <- 100 - proportion.1 - proportion.2
  bd <- density(bs)
  bd <- data.table::data.table("x"=bd$x, "y"=bd$y)
  cl <- RColorBrewer::brewer.pal(3, "Accent")
  ggplot(NULL) +
    geom_area(data=bd, aes(x=x, y=y), fill="grey50") +
    geom_area(data=bd[x>=-2 & x<=2, ], aes(x=x, y=y), fill=cl[1]) +
    geom_area(data=bd[x>=-1 & x<=1, ], aes(x=x, y=y), fill=cl[3]) +
    geom_vline(xintercept=-2:2, linetype=3, size=1, color="black") +
    geom_line(data=bd, aes(x=x, y=y), color="black") +
    scale_x_continuous(name="Aligned contact sites", breaks=seq(-10, 10, by=1)) +
    ylab("density") +
    ggpubr::theme_pubr(base_size=16) +
    annotate("text", label=propFormat(proportion.1), color=cl[3], x=3, hjust=0, y=0.2, size=5) +
    annotate("text", label=propFormat(proportion.2), color=cl[1], x=3, hjust=0, y=0.15, size=5) +
    annotate("text", label=propFormat(proportion.3), color="grey50", x=3, hjust=0, y=0.10, size=5)
}

#' @export
#' @rdname SingleTCRAnalysis
#' @name Features_sCPP
Features_sCPP <- function(
  peptideSet, tcrSet,
  aaIndexIDSet=c("BETM990101inv","KESO980102inv","KOLA930101","MICC010101","SIMK990104","VENM980101","VENM980101inv"),
  fragLenSet=3:8,
  seed=12345,
  coreN=parallel::detectCores(logical=F)
){
  if(length(peptideSet)!=length(tcrSet)){
    message("The number of peptides and that of TCRs must be the same!")
    return(NULL)
  }

  AACPMatrixList <- CPP_AACPMatrix()
  paramDT <- data.table::CJ(
    "FragLen"=fragLenSet,
    "AAIndexID"=aaIndexIDSet
  )

  cl <- parallel::makeCluster(coreN, type="SOCK")
  snow::clusterSetupRNGstream(cl, seed=rep(seed, 6))
  doSNOW::registerDoSNOW(cl)
  sink(tempfile())
  pb <- pbapply::timerProgressBar(max=length(peptideSet), style=1)
  sink()
  opts <- list(progress=function(n) {
    pbapply::setTimerProgressBar(pb, n)
  })
  dt_scpp <- foreach::foreach(i=1:length(peptideSet), .inorder=F, .packages=c("tidyverse","data.table","Repitope"), .options.snow=opts)%dopar%{
    pept <- peptideSet[i]
    tcr <- tcrSet[i]
    res <- lapply(1:nrow(paramDT), function(j){
      fragLen <- paramDT$FragLen[j]
      aaIndexID <- paramDT$AAIndexID[j]
      aacpMat <- AACPMatrixList[[aaIndexID]]
      fragSet <- sequenceSlidingWindow(tcr, w=fragLen)
      fragSet <- c(fragSet, stringi::stri_reverse(fragSet))
      al <- Biostrings::pairwiseAlignment(
        subject=pept, pattern=fragSet,
        substitutionMatrix=aacpMat,
        type="global-local", gapOpening=100, gapExtension=100,
        scoreOnly=T
      )
      statSet <- c("max","min","mean","sd","skew","kurtosis")
      stat <- psych::describe(al, interp=F, skew=T, type=3, ranges=T, IQR=F)[statSet]
      as.numeric(stat)
    }) %>%
      data.table::as.data.table() %>%
      data.table::transpose()
    colnames(res) <- c("Max","Min","Mean","SD","Skew","Kurt")
    res[,Peptide:=pept][,CDR3B:=tcr]
    res <- cbind(paramDT, res)
    return(res)
  } %>%
    data.table::rbindlist()
  data.table::setorder(dt_scpp, Peptide, CDR3B, FragLen, AAIndexID)
  close(pb)
  parallel::stopCluster(cl)
  gc();gc()

  col_id <- c("Peptide", "CDR3B", "FragLen", "AAIndexID")
  col_val <- setdiff(colnames(dt_scpp), col_id)
  dt_scpp <- data.table::melt.data.table(dt_scpp, id=col_id, measure=col_val, variable.name="Stat", value.name="Value")
  dt_scpp[,Feature:=paste0("sCPP_",AAIndexID,"_",Stat,"_",FragLen)]
  dt_scpp <- data.table::dcast.data.table(dt_scpp, Peptide+CDR3B~Feature, value.var="Value", fun=mean)
  return(dt_scpp)
}

#' @export
#' @rdname SingleTCRAnalysis
#' @name univariateCorrelationAnalysis
univariateCorrelationAnalysis <- function(
  dt_deltaG_sCPP,
  coreN=parallel::detectCores(logical=F)
){
  cl <- parallel::makeCluster(coreN, type="SOCK")
  featureSet <- grep("sCPP_", colnames(dt_deltaG_sCPP), value=T)
  dt_deltaG_sCPP_univ <- pbapply::pblapply(featureSet, function(col){
    lm.model <- lm(as.formula(paste0("DeltaG~", col)), data=dt_deltaG_sCPP)
    lm.summary <- summary(lm.model)
    adj.r.sq <- lm.summary$adj.r.squared
    f.stat <- lm.summary$fstatistic
    p <- pf(f.stat[1], f.stat[2], f.stat[3], lower.tail=F)
    data.table::data.table("Feature"=col, "AdjRSq"=adj.r.sq, "PValue"=p)
  }, cl=cl) %>%
    data.table::rbindlist() %>%
    data.table::setorder(PValue)
  dt_deltaG_sCPP_univ[,AAIndexID:=sapply(stringr::str_split(Feature, "_"), dplyr::nth, 2)]
  parallel::stopCluster(cl)
  gc();gc()
  return(dt_deltaG_sCPP_univ)
}

#' @export
#' @rdname SingleTCRAnalysis
#' @name univariateCorrelationPlot
univariateCorrelationPlot <- function(
  dt_deltaG_sCPP,
  dt_deltaG_sCPP_outlier=NULL,
  corVarName="sCPP_BETM990101inv_Skew_4"
){
  summary.model <- summary(lm(as.formula(paste0("DeltaG~", corVarName)), data=dt_deltaG_sCPP))
  adj.r.sq <- summary.model$adj.r.squared
  f.stat <- summary.model$fstatistic
  p <- pf(f.stat[1], f.stat[2], f.stat[3], lower.tail=F)
  pl <- ggpubr::ggscatter(
    dt_deltaG_sCPP, x="DeltaG", y=corVarName, xlab="Delta G (kcal/mol)",
    color="black", shape=21, size=3,
    add="reg.line",
    add.params=list(color="blue", fill="lightgray"),
    conf.int=T
  )
  if(!is.null(dt_deltaG_sCPP_outlier)){
    pl <- pl + geom_point(data=dt_deltaG_sCPP_outlier,
                          aes_string(x="DeltaG", y=corVarName),
                          shape=21, fill="firebrick1", size=3)
  }
  pl <- pl + ggpubr::theme_pubr(base_size=16)
  x_range <- range(dt_deltaG_sCPP$DeltaG)
  y_range <- range(dt_deltaG_sCPP[[corVarName]])
  pl <- pl +
    annotate("text", x=x_range[1], hjust=0, y=y_range[2], size=5, label=paste0("italic(R)^2==", formatC(adj.r.sq, format="fg", digits=2)), parse=T) +
    annotate("text", x=x_range[1], hjust=0, y=y_range[1]*0.1+y_range[2]*0.9, size=5, label=paste0("italic(P)==", formatC(p, format="fg", digits=2)), parse=T)
  print(pl)
  return(pl)
}

#' @export
#' @rdname SingleTCRAnalysis
#' @name univariateFeatureSelect
univariateFeatureSelect <- function(dt_univ, aaIndexIDSet, sig=0.05, bestOnly=T){
  d <- dt_univ %>%
    dplyr::filter(PValue<sig) %>%
    dplyr::arrange(PValue)
  if(bestOnly) d <- dplyr::distinct(d, AAIndexID, .keep_all=T)
  grep(paste0(paste0("sCPP_", aaIndexIDSet, "_"), collapse="|"), d$Feature, value=T)
}

#' @export
#' @rdname SingleTCRAnalysis
#' @name multivariateRegressionPlot
multivariateRegressionPlot <- function(dt_deltaG_sCPP, dt_deltaG_sCPP_outlier=NULL, corVarNames){
  full.model <- lm(
    as.formula(paste0("DeltaG~", paste0(corVarNames, collapse="+"))),
    data=dt_deltaG_sCPP
  )
  step.model <- MASS::stepAIC(full.model, direction="both", trace=F)
  summary.step.model <- summary(step.model)
  cat("The result of stepwise elimination.\n")
  print(summary.step.model)
  adj.r.sq <- summary.step.model$adj.r.squared
  f.stat <- summary.step.model$fstatistic
  p <- pf(f.stat[1], f.stat[2], f.stat[3], lower.tail=F)
  stat <- data.table::data.table("AdjRSq"=adj.r.sq, "PValue"=p)
  if(!is.null(dt_deltaG_sCPP_outlier)){
    dt_deltaG_sCPP <- dplyr::bind_rows(
      dplyr::mutate(dt_deltaG_sCPP, Type="MVR"),
      dplyr::mutate(dt_deltaG_sCPP_outlier, Type="Outlier")
    )
  }else{
    dt_deltaG_sCPP <- dplyr::mutate(dt_deltaG_sCPP, Type="MVR")
  }
  dt_deltaG_sCPP <- dplyr::mutate(dt_deltaG_sCPP, Prediction=predict(step.model, newdata=dt_deltaG_sCPP))
  pl <- ggpubr::ggscatter(
    dplyr::filter(dt_deltaG_sCPP, Type=="MVR"),
    x="DeltaG", y="Prediction", xlab="Delta G (kcal/mol)", ylab="Multivariate regression",
    color="black", shape=21, size=3,
    add="reg.line",
    add.params=list(color="blue", fill="lightgray"),
    conf.int=T
  )
  if(!is.null(dt_deltaG_sCPP_outlier)){
    pl <- pl + geom_point(data=dplyr::filter(dt_deltaG_sCPP, Type=="Outlier"),
                          aes(x=DeltaG, y=Prediction),
                          shape=21, fill="firebrick1", size=3)
  }
  pl <- pl + ggpubr::theme_pubr(base_size=16)
  x_range <- range(dt_deltaG_sCPP$DeltaG)
  y_range <- range(dt_deltaG_sCPP$Prediction)
  pl <- pl +
    annotate("text", x=x_range[1], hjust=0, y=y_range[2], size=5, label=paste0("italic(R)^2==", formatC(stat$AdjRSq, format="fg", digits=2)), parse=T) +
    annotate("text", x=x_range[1], hjust=0, y=y_range[1]*0.2+y_range[2]*0.8, size=5, label=paste0("italic(P)==", formatC(stat$PValue, format="fg", digits=2)), parse=T)
  dt_deltaG_sCPP_vif <- car::vif(step.model)
  dt_deltaG_sCPP_vif <- data.table::data.table("Var"=names(dt_deltaG_sCPP_vif), "VIF"=formatC(dt_deltaG_sCPP_vif, format="f", digits=2))
  tbl <- ggpubr::ggtexttable(
    dt_deltaG_sCPP_vif, rows=NULL, theme=ggpubr::ttheme("minimal", base_size=12)
  )
  pl <- ggpubr::ggarrange(pl, tbl, nrow=1, ncol=2, align="v")
  print(pl)
  return(pl)
}

#' @export
#' @rdname SingleTCRAnalysis
#' @name multivariateRegressionPrediction
multivariateRegressionPrediction <- function(dt_deltaG_sCPP, dt_deltaG_sCPP_new=NULL, corVarNames){
  full.model <- lm(
    as.formula(paste0("DeltaG~", paste0(corVarNames, collapse="+"))),
    data=dt_deltaG_sCPP
  )
  step.model <- MASS::stepAIC(full.model, direction="both", trace=F)
  summary.step.model <- summary(step.model)
  cat("The result of stepwise elimination.\n")
  print(summary.step.model)
  if(is.null(dt_deltaG_sCPP_new)){
    return(dplyr::mutate(dt_deltaG_sCPP, Prediction=predict(step.model, newdata=dt_deltaG_sCPP)))
  }else{
    return(dplyr::mutate(dt_deltaG_sCPP_new, Prediction=predict(step.model, newdata=dt_deltaG_sCPP_new)))
  }
}
