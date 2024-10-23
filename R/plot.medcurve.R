#' Plot the results of the function MediationCurve
#' 
#' This function plots one of the estimates of interest from the output of the 
#' \code{MediationCurve} function as either a smoothed curve or a series of 
#' straight lines connecting the point estimates. Point x-axis values are 
#' plotted at the mean of each treatment 'bin'. 
#' 
#' @param med.results This is an output from the function \code{MediationCurve} 
#' with the class \code{'medcurve'}.
#' @param plot.est This is the estimate of interest to be plotted. Options 
#' include \code{"ACME"} for the average causal mediation effect, \code{"ADE"} 
#' for the average direct effect, \code{"Total"} for the total effect, and 
#' \code{"Proportion"} for the proportion of the effect mediated (i.e., 
#' ACME/Total Effect). This defaults to \code{"ACME"}.
#' @param addci This is a logical value specifying whether a ribbon of 
#' 95% confidence intervals should be wrapped around the plotted curve or line. 
#' This defaults to \code{TRUE}.
#' @param smoothplot This is a logical value specifying whether the main line 
#' and 95% confidence interval ribbon (if \code{addci = TRUE}) connecting the 
#' point estimates should be a smoothed loess fit if \code{TRUE} or a series of 
#' straight lines connecting the point estimates if \code{FALSE}. Note that in 
#' the case of a few unevenly spaced out values with sharp curves, a smoothed 
#' fit may not go through all of the point estimates. This defaults to 
#' \code{TRUE}.
#' @param addpts This is a logical value specifying whether to add points for 
#' the point estimates of the estimate of interest. This can be valuable for 
#' assessing how much a smoothed fit may be deviating from the actual point 
#' estimates. This defaults to \code{FALSE}.
#' @param addptci This is a logical value specifying whether to add error bars 
#' for the 95% confidence intervals of the point estimates of the estimate of 
#' interest. This can be valuable for assessing how much a smoothed confidence 
#' interval ribbon fit may be deviating from the actual point estimates. This 
#' defaults to \code{FALSE}.
#' @param plotrug This is a logical value specifying whether to add a rug plot 
#' along the x-axis showing observed treatment values from the data used to 
#' fit the models. If \code{TRUE}, a vector of treatment values must be 
#' provided for the argument \code{rugvalues}. This defaults to \code{FALSE}.
#' @param rugvalues This is a numeric vector of the treatment values to be 
#' included in the rug plot if \code{plotrug = TRUE}. This defaults to 
#' \code{NULL}. If \code{plotrug = TRUE} and \code{rugvalues} is NULL, an error 
#' will be returned. 
#' 
#' @return \code{plot.medcurve} returns the specified plot of the class 
#' \code{'ggplot'}.
#' 
#' @import ggplot2
#' @export plot.medcurve
#' 
#' @examples
#' \dontrun{
#' #Smoothed ACME fit
#' plot(medres)
#' 
#' #Non-smoothed ADE fit with rug plot
#' plot(medres, "ADE", smoothplot = F, plotrug = T, rugvalues = meddata$treat)
#' 
#' #Smoothed Total Effect fit with point estimate points and CIs
#' plot(medres, "Total", addpts = T, addptci = T)
#' }
#' 
plot.medcurve <- function(med.results, plot.est = "ACME", addci = T, 
                          smoothplot = T, addpts = F, addptci = F, 
                          plotrug = F, rugvalues = NULL){
  if(!plot.est %in% c("ACME", "ADE", "Total", "Proportion")){
    stop(paste0("plot.est must be one of 'ACME', 'ADE', 'Total', ", 
                "or 'Proportion'"))
  }
  if(plotrug & is.null(rugvalues)){
    stop(paste0("If plotrug = TRUE, you must supply a vector of numeric values",
                " for rugvalues."))
  }
  plotdat <- med.results
  plotdat$gridmean <- rowMeans(plotdat[c("control.value", "treat.value")])
  if(plot.est == "ACME"){
    plotdat[, c("est", "lci", "uci")] <- plotdat[
      , c("d.avg", "d.avg.ci.2.5", "d.avg.ci.97.5")]
  } else if(plot.est == "ADE"){
    plotdat[, c("est", "lci", "uci")] <- plotdat[
      , c("z.avg", "z.avg.ci.2.5", "z.avg.ci.97.5")]
  } else if(plot.est == "Proportion"){
    plotdat[, c("est", "lci", "uci")] <- plotdat[
      , c("d.avg", "d.avg.ci.2.5", "d.avg.ci.97.5")]
    plot.est <- "Proportion Mediated"
  } else {
    plotdat[, c("est", "lci", "uci")] <- plotdat[
      , c("tau.coef", "tau.ci.2.5", "tau.ci.97.5")]
    plot.est <- "Total Effect"
  }
  gg1 <- ggplot() + theme_bw() +
    geom_hline(yintercept = 0)
  if(addci){
    if(smoothplot){
      gg1 <- gg1 + stat_smooth(data = plotdat, aes(
        x = gridmean, y = est, ymin = predict(loess(lci ~ gridmean)),
        ymax = predict(loess(uci ~ gridmean))), geom = "ribbon",
        alpha = 0.5, fill = "grey50", show.legend = F)
    } else {
      gg1 <- gg1 + geom_ribbon(data = plotdat, aes(x = gridmean, y = est,
                               ymin = lci, ymax = uci),
                               alpha = 0.5, fill = "grey50", show.legend = F)
    }
  }
  if(smoothplot){
    gg1 <- gg1 + geom_smooth(data = plotdat, aes(x = gridmean, y = est), 
                             se = F, formula = 'y ~ x', method = 'loess') 
  } else {
    gg1 <- gg1 + geom_line(data = plotdat, aes(x = gridmean, y = est), 
                           linewidth = 1, linetype = 1, color = "#3366FF") 
  }
  gg1 <- gg1 + xlab("Treatment") + ylab(plot.est)
  if(addpts) gg1 <- gg1 + geom_point(data = plotdat, 
                                     aes(x = gridmean, y = est), size = 2)
  if(addptci) gg1 <- gg1 + geom_errorbar(data = plotdat, aes(
    x = gridmean, y = est, ymin = lci, ymax = uci))
  if(plotrug){
    if(is.null(rugvalues)) stop(paste0(
      "If plotrug is TRUE, a vector of treatment values to be rug plotted ",
      "along the x-axis must be provided as the argument rugvalues."))
    rugdat <- data.frame(treat = rugvalues)
    gg1 <- gg1 + geom_rug(data = rugdat, aes(x = treat))
  }
  return(gg1)
}
