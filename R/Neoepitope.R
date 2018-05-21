#' Neoepitope analysis.
#'
#' \code{Neoepitope_BurdenDT} returns a neoepitope summary datatable using the thresholds provided. \cr
#' \code{Neoepitope_PValueDT} calculates P-values by the log-rank test with varying thresholds for the neoepitope burden defined by the specified threshold of the neoepitope dissimilarity index. If no threshold is provided, neoepitope dissimilarity indices were aggregated by patients without thresholding. \cr
#' \code{Neoepitope_PValueDT_Batch} does the jobs for a set of thresholds of the neoepitope dissimilarity index. \cr
#' \code{Neoepitope_PValueLineChart} generates a line chart of threshold values and P-values.\cr
#' \code{Neoepitope_KMPlot} generates a Kaplan-Mayer curve with the thresholds provided.
#'
#' @param dt_neoepitope A datatable of neoepitopes that contains the following columns: "Dataset","Sample","Months","Status", "ImmunogenicityScore.WT", and, "ImmunogenicityScore.MT".
#' @param dt_neoepitope_burden A datatable returned by \code{Neoepitope_BurdenDT}.
#' @param dt_pvalue A datatable returned by \code{Neoepitope_PValueDT_Batch}.
#' @param thr A threshold of neoepitope dissimilarity index. Can be set as "none" to disabled. In that case, neoepitope dissimilarity indices were aggregated by patients without thresholding.
#' @param thrSet A set of thresholds of neoepitope dissimilarity index.
#' @param thr.ne A threshold of neoepitope burden or aggregated neoepitope dissimilarity index per patient.
#' @param coreN The number of cores to be used for parallelization. Set \code{NULL} to disable parallelization.
#' @export
#' @rdname Neoepitope
#' @name Neoepitope
Neoepitope_BurdenDT <- function(dt_neoepitope, thr="none", thr.ne=1){
  dt_neoepitope[,Neoepitope:=abs(ImmunogenicityScore.WT - ImmunogenicityScore.MT)/ImmunogenicityScore.WT]
  if(is.numeric(thr)){
    dt_neoepitope[,Neoepitope:=Neoepitope>=thr]
    if(sum(dt_neoepitope$Neoepitope)==0){
      cat("No neoepitopes were retained! Try different thresholds.", sep="")
      return(NULL)
    }
    neoepitopeCount <- sum(dt_neoepitope$Neoepitope)
    cat(formatC(neoepitopeCount/nrow(dt_neoepitope)*100, digits=1, format="f", drop0trailing=F), "% of neoepitope candidates were retained after filtering by neoepitope dissimilarity index.", sep="")
  }
  dt_neoepitope <- data.table::dcast.data.table(dt_neoepitope, Dataset+Sample+Months+Status~., value.var="Neoepitope", fun=sum)
  colnames(dt_neoepitope) <- c("Dataset","Sample","Months","Status","NeoepitopeBurden")
  dt_neoepitope[,NeoepitopeBurdenGroup:=factor(dplyr::if_else(dt_neoepitope$NeoepitopeBurden>=thr.ne, 1, 0), levels=c(0, 1), labels=c("Low","High"))]
  return(dt_neoepitope)
}

#' @export
#' @rdname Neoepitope
#' @name Neoepitope
Neoepitope_PValueDT <- function(dt_neoepitope, thr="none"){
  dt_neoepitope[,Neoepitope:=abs(ImmunogenicityScore.WT - ImmunogenicityScore.MT)/ImmunogenicityScore.WT]
  if(is.numeric(thr)){
    dt_neoepitope[,Neoepitope:=Neoepitope>=thr]
    if(sum(dt_neoepitope$Neoepitope)==0){
      dt_pval <- data.table::data.table("PValue"=NA, "NeoepitopeBurdenThreshold"=NA)
      dt_pval[,"Threshold":=thr]
      return(dt_pval)
    }
  }
  dt_neoepitope <- data.table::dcast.data.table(dt_neoepitope, Dataset+Sample+Months+Status~., value.var="Neoepitope", fun=sum)
  colnames(dt_neoepitope) <- c("Dataset","Sample","Months","Status","NeoepitopeBurden")
  if(is.numeric(thr)){
    neoepitope_thrSet <- unique(round(seq(from=min(dt_neoepitope$NeoepitopeBurden), to=max(dt_neoepitope$NeoepitopeBurden), by=1)))
  }else{
    neoepitope_thrSet <- seq(from=min(dt_neoepitope$NeoepitopeBurden), to=max(dt_neoepitope$NeoepitopeBurden), length.out=1000)
  }
  logrankPValue <- function(thr.ne){
    dt <- data.table::copy(dt_neoepitope)
    dt[,NeoepitopeBurdenGroup:=dplyr::if_else(NeoepitopeBurden>=thr.ne, 1, 0)]
    if(dplyr::n_distinct(dt$NeoepitopeBurdenGroup)==1){
      return(data.table::data.table("PValue"=NA, "NeoepitopeBurdenThreshold"=thr.ne))
    }else{
      diff <- survival::survdiff(survival::Surv(Months, Status)~NeoepitopeBurdenGroup, data=dt)
      p <- stats::pchisq(diff$chisq, length(diff$n)-1, lower.tail=F)
      return(data.table::data.table("PValue"=p, "NeoepitopeBurdenThreshold"=thr.ne))
    }
  }
  dt_pval <- data.table::rbindlist(lapply(neoepitope_thrSet, logrankPValue))
  dt_pval[,"Threshold":=thr]
  return(dt_pval)
}

#' @export
#' @rdname Neoepitope
#' @name Neoepitope
Neoepitope_PValueDT_Batch <- function(dt_neoepitope, thrSet=seq(0, 1, length.out=1000), coreN=parallel::detectCores()){
  if(is.null(coreN)){
    dt_pvalue <- foreach::foreach(thr=thrSet, .export="Neoepitope_PValueDT", .packages="data.table", .inorder=F) %do% {
      Neoepitope_PValueDT(dt_neoepitope, thr)
    } %>% data.table::rbindlist()
    return(dt_pvalue)
  }else{
    doParallel::registerDoParallel(coreN)
    dt_pvalue <- foreach::foreach(thr=thrSet, .export="Neoepitope_PValueDT", .packages="data.table", .inorder=F) %dopar% {
      Neoepitope_PValueDT(dt_neoepitope, thr)
    } %>% data.table::rbindlist()
    doParallel::stopImplicitCluster()
    return(dt_pvalue)
  }
}

#' @export
#' @rdname Neoepitope
#' @name Neoepitope
Neoepitope_PValueLineChart <- function(dt_pvalue){
  if(identical(unique(dt_pvalue$"Threshold"), "none")){
    g <- ggplot(dt_pvalue, aes_string(x="NeoepitopeBurdenThreshold", y="PValue")) + geom_line()
  }else{
    dt_pvalue_plot <- data.table::dcast.data.table(dt_pvalue[stats::complete.cases(dt_pvalue),], Threshold~., value.var="PValue", fun=min)
    colnames(dt_pvalue_plot) <- c("Threshold", "PValue")
    g <- ggplot(dt_pvalue_plot, aes_string(x="Threshold", y="PValue")) + geom_line()
  }
  g + xlab("Threshold") +
    scale_y_log10(name="P-value", breaks=c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1), labels=scales::comma) +
    ggpubr::theme_pubr(base_size=16)
}

#' @export
#' @rdname Neoepitope
#' @name Neoepitope
Neoepitope_KMPlot <- function(dt_neoepitope_burden){
  surv <- survival::survfit(survival::Surv(Months, Status)~NeoepitopeBurdenGroup, data=dt_neoepitope_burden)
  res <- survminer::ggsurvplot(
    fit=surv,
    data=dt_neoepitope_burden,
    palette=c("grey50", "firebrick1"), linetype="solid", conf.int=T, size=0.75,
    font.x=16, font.y=16, font.tickslab=12, font.legend=16,
    pval=T, pval.method=T, log.rank.weights="1", pval.size=6, pval.coord=c(0.1, 0.1), pval.method.coord=c(0.1, 0.2),
    legend.title="Neoepitope burden", legend.labs=c("Low","High"),
    risk.table=T, risk.table.y.text.col=T, tables.y.text=F
  ) + xlab("Time [month]")
  res$plot <- res$plot + ggpubr::rremove("x.title") + ggpubr::rremove("x.text")
  res$table <- res$table + ggpubr::rremove("y.title")
  ggpubr::ggarrange(res$plot, res$table, nrow=2, ncol=1, heights=c(1, 0.35), align="v")
}
