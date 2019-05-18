#' Utility functions for creating diagnostic plots for binary classifiers.
#'
#' @param trueClass A factor of true class labels.
#' @param predProb A numeric vector of predicted scores etc.
#' @param groups A vector indicating strata. Disable by setting NA.
#' @param colors A set of colors. The length should be the same with the number of groups.
#' @param plotTypes A character vector indicating which types of plots should be generated.
#' @export
#' @rdname Utility_Classifier
#' @name Utility_Classifier
compute_auc <- function(trueClass, predProb){
  roc_obj <- pROC::roc(response=trueClass, predictor=predProb)
  aucSet <- as.numeric(pROC::ci.auc(roc_obj, method="delong"))*100
  aucSet <- aucSet[c(2,1,3)]
  names(aucSet) <- c("AUC", "95%CI_Lo", "95%CI_Up")
  return(aucSet)
}

#' @export
#' @rdname Utility_Classifier
#' @name Utility_Classifier
rocPlot <- function(trueClass, predProb, groups=NA, colors=NA){
  ## Data
  dt <- data.table::data.table("Truth"=trueClass, "Probability"=predProb, "Group"=groups)
  lev <- levels(trueClass) ## c("Negative","Positive")
  
  ## ROC
  roc_dt <- data.table::rbindlist(lapply(split(dt, by="Group"), function(d){
    pr_obj <- precrec::evalmod(scores=d$"Probability", labels=d$"Truth")
    pr_obj <- autoplot(pr_obj, "ROC")
    data.table::as.data.table(pr_obj$data)[,x:=x*100][,y:=y*100][,"Group":=unique(d$"Group")]
  }))
  auc_dt <- data.table::rbindlist(lapply(split(dt, by="Group"), function(d){
    aucSet <- compute_auc(trueClass=d$"Truth", predProb=d$"Probability")
    data.table::data.table(Group=unique(d$"Group"), AUC=paste0("AUC: ", format(aucSet["AUC"]*100, digits=3), "% [", format(aucSet["95%CI_Lo"]*100, digits=3), "%-", format(aucSet["95%CI_Up"]*100, digits=3),"%]"))
  }))
  data.table::setorder(auc_dt, "Group")
  if(identical(groups, NA)){
    plt <- ggplot(roc_dt, aes(x=x, y=y)) +
      geom_line(size=1.5) +
      geom_abline(slope=1, intercept=0, linetype="dotted") +
      annotate("text", x=60, y=25, label=auc_dt$AUC, parse=F, size=4.5)
  }else{
    if(identical(colors, NA)) colors <- ggsci::pal_d3()(dplyr::n_distinct(groups))
    plt <- ggplot(roc_dt, aes(x=x, y=y)) +
      geom_line(aes(color=Group), size=1.5) +
      geom_abline(slope=1, intercept=0, linetype="dotted") +
      scale_color_manual(values=colors, guide=F)
    cols <- matrix(rep.int(colors, times=2), nrow=nrow(auc_dt), ncol=2)
    tt <- gridExtra::ttheme_minimal(base_size=12,
                                    core=list(fg_params=list(col=cols), bg_params=list(col=NA)),
                                    rowhead=list(bg_params=list(col=NA)),
                                    colhead=list(bg_params=list(col=NA)))
    plt <- plt +
      annotation_custom(gridExtra::tableGrob(auc_dt, rows=NULL, cols=NULL, theme=tt),
                        xmin=50, xmax=80, ymin=5, ymax=30)
  }
  plt <- plt +
    scale_x_continuous(name="FPR(%)", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    scale_y_continuous(name="TPR(%)", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    ggtitle("ROC") + ggpubr::theme_pubr(base_size=16)
  return(plt)
}

#' @export
#' @rdname Utility_Classifier
#' @name Utility_Classifier
prcPlot <- function(trueClass, predProb, groups=NA, colors=NA){
  ## Data
  dt <- data.table::data.table("Truth"=trueClass, "Probability"=predProb, "Group"=groups)
  lev <- levels(trueClass) ## c("Negative","Positive")
  
  ## Precision-recall curve
  prc_dt <- data.table::rbindlist(lapply(split(dt, by="Group"), function(d){
    pr_obj <- precrec::evalmod(scores=d$"Probability", labels=d$"Truth")
    pr_obj <- autoplot(pr_obj, "PRC")
    data.table::as.data.table(pr_obj$data)[,x:=x*100][,y:=y*100][,"Group":=unique(d$"Group")]
  }))
  if(identical(groups, NA)){
    plt <- ggplot(prc_dt, aes(x=x, y=y)) +
      geom_line(size=1.5) +
      geom_abline(slope=-1, intercept=100, linetype="dotted")
  }else{
    if(identical(colors, NA)) colors <- ggsci::pal_d3()(dplyr::n_distinct(groups))
    plt <- ggplot(prc_dt, aes(x=x, y=y)) +
      geom_line(aes(color=Group), size=1.5) +
      geom_abline(slope=-1, intercept=100, linetype="dotted") +
      scale_color_manual(values=colors)
  }
  plt <- plt +
    scale_x_continuous(name="Recall(%)", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    scale_y_continuous(name="Precision(%)", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    ggtitle("Precision-Recall") + ggpubr::theme_pubr(base_size=16) +
    theme(legend.justification=c(0, 0), legend.position=c(0.1, 0.1), legend.title=element_blank())
  return(plt)
}

#' @export
#' @rdname Utility_Classifier
#' @name Utility_Classifier
calibPlot <- function(trueClass, predProb, groups=NA, colors=NA){
  ## Data
  dt <- data.table::data.table("Truth"=trueClass, "Probability"=predProb, "Group"=groups)
  lev <- levels(trueClass) ## c("Negative","Positive")
  
  ## Calibrarion chart
  format_calib_dt <- function(trueClass, predProb){
    bucket_array <- seq(0.0, 1.0, by=0.1)
    positive_in_band <- function(bucket){
      in_bucket_indicator <- predProb >= bucket_array[bucket] & predProb < bucket_array[bucket+1]
      bucket_size <- sum(in_bucket_indicator)
      positive <- sum(trueClass[in_bucket_indicator]==lev[2])
      return(qbeta(c(llb=0.025, lb=0.25, y=0.5, ub=0.75, uub=0.965), 0.5+positive, 0.5+bucket_size-positive))
    }
    calib_dt <- data.table::data.table(bucket=1:10, percentage=5+bucket_array[1:10]*100, blb=bucket_array[1:10], bub=bucket_array[(1:10) + 1])
    calib_dt <- cbind(calib_dt, 100*t(sapply(calib_dt$bucket, positive_in_band)))
    return(calib_dt)
  }
  calib_dt <- data.table::rbindlist(lapply(split(dt, by="Group"), function(d){format_calib_dt(d$"Truth", d$"Probability")[,"Group":=unique(d$"Group")]}))
  if(identical(groups, NA)){
    plt <- ggplot(calib_dt, aes(x=percentage, y=y)) +
      geom_ribbon(aes(ymin=llb, ymax=uub), alpha=0.2) +
      geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.4) +
      geom_abline(slope=1, intercept=0, linetype="dotted")
  }else{
    if(identical(colors, NA)) colors <- ggsci::pal_d3()(dplyr::n_distinct(groups))
    plt <- ggplot(calib_dt, aes(x=percentage, y=y)) +
      geom_ribbon(aes(ymin=llb, ymax=uub, fill=Group), alpha=0.2) +
      geom_ribbon(aes(ymin=lb, ymax=ub, fill=Group), alpha=0.4) +
      geom_abline(slope=1, intercept=0, linetype="dotted") +
      scale_fill_manual(values=colors)
  }
  plt <- plt +
    scale_x_continuous(name="Predicted probability (%)", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    scale_y_continuous(name="Smoothed true probability (%)", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    ggtitle("Calibration") + ggpubr::theme_pubr(base_size=16) +
    theme(legend.justification=c(1, 0), legend.position=c(0.9, 0.1), legend.title=element_blank())
  return(plt)
}

#' @export
#' @rdname Utility_Classifier
#' @name Utility_Classifier
cumGainPlot <- function(trueClass, predProb, groups=NA, colors=NA){
  ## Data
  dt <- data.table::data.table("Truth"=trueClass, "Probability"=predProb, "Group"=groups)
  lev <- levels(trueClass) ## c("Negative","Positive")
  dt$Truth <- factor(dt$Truth, levels=rev(lev)) ## necessary for caret::lift
  
  ## Cumulative gain chart
  gain_dt <- data.table::rbindlist(lapply(split(dt, by="Group"), function(d){
    gain_obj <- caret::lift(Truth~Probability, data=d)
    data.table::as.data.table(gain_obj$data)[,"Group":=unique(d$"Group")]
  }))
  xyline_dt <- data.table::rbindlist(lapply(split(gain_dt, by="Group"), function(d){
    v <- 80
    window <- 5
    x <- d$CumTestedPct
    y <- d$CumEventPct
    res <- data.table::data.table(CumEventPct=v, CumTestedPct=NA, Group=unique(d$Group))
    for(i in seq(along=v)){
      nearest <- which.min((y - v[i])^2)
      index <- max(1, nearest - window):min(length(y), nearest + window)
      res$CumTestedPct[i] <-
        if(length(index) > 2){
          approx(y[index], x[index], xout=v[i])$y
        }else{
          NA
        }
    }
    return(res)
  }))
  if(identical(groups, NA)){
    plt <- ggplot(gain_dt, aes(x=CumTestedPct, y=CumEventPct)) +
      geom_line(size=1.5) +
      geom_abline(slope=1, intercept=0, linetype="dotted")
    plt <- plt +
      geom_segment(data=xyline_dt, aes(x=CumTestedPct, y=CumEventPct, xend=CumTestedPct, yend=0)) +
      geom_segment(data=xyline_dt, aes(x=CumTestedPct, y=CumEventPct, xend=0, yend=CumEventPct))
  }else{
    if(identical(colors, NA)) colors <- ggsci::pal_d3()(dplyr::n_distinct(groups))
    plt <- ggplot(gain_dt, aes(x=CumTestedPct, y=CumEventPct)) +
      geom_line(aes(color=Group), size=1.5) +
      geom_abline(slope=1, intercept=0, linetype="dotted") +
      scale_color_manual(values=colors)
    plt <- plt +
      geom_segment(data=xyline_dt, aes(x=CumTestedPct, y=CumEventPct, xend=CumTestedPct, yend=0, color=Group)) +
      geom_segment(data=xyline_dt, aes(x=CumTestedPct, y=CumEventPct, xend=0, yend=CumEventPct, color=Group))
  }
  plt <- plt +
    scale_x_continuous(name="% of data examined", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    scale_y_continuous(name="% of targets found", breaks=seq(0, 100, 20), limits=c(0, 102), expand=c(0, 0)) +
    ggtitle("Cumulative Gain") + ggpubr::theme_pubr(base_size=16) +
    theme(legend.justification=c(1, 0), legend.position=c(0.9, 0.1), legend.title=element_blank())
  return(plt)
}

#' @export
#' @rdname Utility_Classifier
#' @name Utility_Classifier
classifierDiagnosticPlots <- function(trueClass, predProb, groups=NA, colors=NA,
                                      plotTypes=c("ROC","PRC","Calibration","CumGain")){
  plotList <- list()
  if("ROC" %in% plotTypes) plotList$"ROC" <- rocPlot(trueClass, predProb, groups, colors)
  if("PRC" %in% plotTypes) plotList$"PRC" <- prcPlot(trueClass, predProb, groups, colors)
  if("Calibration" %in% plotTypes) plotList$"Calibration" <- calibPlot(trueClass, predProb, groups, colors)
  if("CumGain" %in% plotTypes) plotList$"CumGain" <- cumGainPlot(trueClass, predProb, groups, colors)
  plot_comb <- ggpubr::ggarrange(plotlist=plotList, ncol=length(plotList), nrow=1)
  print(plot_comb)
  return(plotList)
}
