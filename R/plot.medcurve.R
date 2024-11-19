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
#' @param addci This is a logical value specifying whether 95% confidence 
#' intervals should be wrapped around the mean estimate curve or line. 
#' This defaults to \code{TRUE}.
#' @param smoothplot This is a logical value specifying whether the main line 
#' and 95% confidence interval ribbon (if \code{addci = TRUE}) connecting the 
#' point estimates should be a smoothed x-spline fit 
#' (\code{\link[ggalt]{stat_xspline}}) if \code{TRUE} or a series 
#' of straight lines connecting the point estimates if \code{FALSE}. The 
#' x-spline fit forces the curve to go through each observed point. This 
#' defaults to \code{TRUE}.
#' @param addpts This is a logical value specifying whether to add points for 
#' the point estimates of the estimate of interest. This defaults to 
#' \code{FALSE}.
#' @param addptci This is a logical value specifying whether to add error bars 
#' for the 95% confidence intervals of the point estimates of the estimate of 
#' interest. This defaults to \code{FALSE}.
#' @param plotrug This is a logical value specifying whether to add a rug plot 
#' along the x-axis showing observed treatment values from the data used to 
#' fit the models. If \code{TRUE}, a vector of treatment values must be 
#' provided for the argument \code{rugvalues}. This defaults to \code{FALSE}.
#' @param rugvalues This is a numeric vector of the treatment values to be 
#' included in the rug plot if \code{plotrug = TRUE}. This defaults to 
#' \code{NULL}. If \code{plotrug = TRUE} and \code{rugvalues} is NULL, an error 
#' will be returned. 
#' @param xsplineshape This is a value between -1 and 1 controlling the 
#' smoothing of the x-spline (see \code{\link[ggalt]{stat_xspline}}) smooth. 
#' The x-spline approximates a Catmull-Rom spline at -1, has sharp corners at 0, 
#' and approximates a cubic b-spline at 1. This defaults to -0.5, which has the 
#' x-spline smoothly curving directly through all points.
#' @param ribbon This is a logical value specifying whether the confidence 
#' interval curves (if \code{addci = TRUE}) will be plotted as curves (if 
#' \code{FALSE}) or as a filled-in area specified by \code{geom_ribbon} (if 
#' \code{TRUE}). This values defaults to \code{TRUE}. In some cases (e.g., 
#' sometimes with \code{plot.est = "Proportion"}), the confidence interval 
#' curves can be sufficiently different lengths the x-values of the CI curve 
#' boundaries will be different lengths and not coincide enough to merge to form 
#' the dataset used by geom_ribbon to fill in the ribbon area. If that's the 
#' case, a warning will print and the argument \code{'ribbon'} will be ignored.
#' @param addcurve This is a logical value specifying whether a curve linking 
#' the mediation estimate values (and the CI values if \code{addci = TRUE}) 
#' should be drawn. One may want to plot just the points and CIs without linking 
#' them with a curve line, in which case one should set this argument to 
#' \code{FALSE} and specify \code{addpts = TRUE} and \code{addptci = TRUE}. This 
#' defaults to \code{TRUE}.
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
#' plot.medcurve(medres)
#' 
#' #Non-smoothed ADE fit with rug plot
#' plot.medcurve(medres, "ADE", smoothplot = F, plotrug = T, rugvalues = meddata$treat)
#' 
#' #Smoothed Total Effect fit with point estimate points and CIs
#' plot.medcurve(medres, "Total", addpts = T, addptci = T)
#' }
#' 
plot.medcurve <- function(med.results, plot.est = "ACME", addci = T, 
                          smoothplot = T, addpts = F, addptci = F, 
                          plotrug = F, rugvalues = NULL, 
                          xsplineshape = -0.5, ribbon = T, addcurve = T){
  if(!plot.est %in% c("ACME", "ADE", "Total", "Proportion")){
    stop(paste0("plot.est must be one of 'ACME', 'ADE', 'Total', ", 
                "or 'Proportion'"))
  }
  if(plotrug & is.null(rugvalues)){
    stop(paste0("If plotrug = TRUE, you must supply a vector of numeric values",
                " for rugvalues."))
  }
  if(!addcurve & !addpts & !addptci){
    stop(paste0("At least one of addcurve, addpts, and addptci must be TRUE."))
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
      , c("n.avg", "n.avg.ci.2.5", "n.avg.ci.97.5")]
    plot.est <- "Proportion Mediated"
  } else {
    if("tau.coef" %in% names(plotdat)) names(plotdat) <- gsub("tau.coef", "tau", 
                                                              names(plotdat))
    plotdat[, c("est", "lci", "uci")] <- plotdat[
      , c("tau", "tau.ci.2.5", "tau.ci.97.5")]
    plot.est <- "Total Effect"
  }
  gg1 <- ggplot() + theme_bw() +
    geom_hline(yintercept = 0)
  if(addci){
    if(smoothplot){
      if(ribbon){
        tmp <- ggplot(data = plotdat, aes(x = gridmean)) +
          stat_xspline(aes(y = lci),
                       spline_shape = xsplineshape) +
          stat_xspline(aes(y = uci),
                       spline_shape = xsplineshape)
        tmp <- ggplot_build(tmp)
        if(nrow(tmp$data[[1]]) != nrow(tmp$data[[2]])){
          warning(paste0("The CI curve lengths differ too much to merge the ",
                         "curves into the same data frame to build the ",
                         "geom_ribbon, and so the argument 'ribbon' will be ",
                         "ignored and curve lines will be provided instead."))
          gg1 <- gg1 + stat_xspline(
            data = tmp$data[[1]], aes(x = x, y = y), spline_shape = -0.5, 
            color = "grey50", size = 1, show.legend = F) + 
            stat_xspline(
              data = tmp$data[[2]], aes(x = x, y = y), spline_shape = -0.5, 
              color = "grey50", size = 1, show.legend = F)
        } else {
          tmp <- data.frame(x = tmp$data[[1]]$x, ymin = tmp$data[[1]]$y,
                            ymax = tmp$data[[2]]$y)
          gg1 <- gg1 + geom_ribbon(data = tmp, aes(
            x = x, ymin = ymin, ymax = ymax), fill = "grey50", alpha = 0.5, 
            show.legend = F)
        }
      } else {
        gg1 <- gg1 + stat_xspline(
          data = plotdat, aes(x = gridmean, y = lci), spline_shape = -0.5, 
          color = "grey50", size = 1, show.legend = F) + 
          stat_xspline(
            data = plotdat, aes(x = gridmean, y = uci), spline_shape = -0.5, 
            color = "grey50", size = 1, show.legend = F)
      }
      
    } else {
      gg1 <- gg1 + geom_ribbon(data = plotdat, aes(x = gridmean, y = est,
                               ymin = lci, ymax = uci),
                               alpha = 0.5, fill = "grey50", show.legend = F)
    }
  }
  if(addcurve){
    if(smoothplot){
      gg1 <- gg1 + stat_xspline(
        data = plotdat, aes(x = gridmean, y = est), spline_shape = -0.5, 
        color = "#3366FF", size = 1)
    } else {
      gg1 <- gg1 + geom_line(data = plotdat, aes(x = gridmean, y = est), 
                             linewidth = 1, linetype = 1, color = "#3366FF") 
    }
  }
  gg1 <- gg1 + xlab("Treatment") + ylab(plot.est)
  if(addpts) gg1 <- gg1 + geom_point(data = plotdat, 
                                     aes(x = gridmean, y = est), size = 2)
  if(addptci) gg1 <- gg1 + geom_errorbar(data = plotdat, aes(
    x = gridmean, y = est, ymin = lci, ymax = uci))
  if(plotrug){
    if(is.null(rugvalues)) stop(paste0(
      "If plotrug is TRUE, a vector of treatment values to be rug plotted ",
      "along the x-axis must be provided as the argument 'rugvalues'."))
    rugdat <- data.frame(treat = rugvalues)
    gg1 <- gg1 + geom_rug(data = rugdat, aes(x = treat))
  }
  return(gg1)
}
