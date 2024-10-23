#' Examine causal mediation effects across a sequence of continuous treatment values
#' 
#' This function extends the functionality of the 
#' \code{\link[mediation]{mediate}} function for causal mediation analysis as 
#' defined by 
#' \href{https://imai.fas.harvard.edu/research/files/BaronKenny.pdf}{Imai et al. 2010} 
#' for modelling linear and nonlinear mediation effects in the case of a 
#' continuous treatment. This function is largely inspired by the 
#' \href{https://causal-curve.readthedocs.io/en/latest/intro.html}{causal-curve} 
#' Python package by Roni Kobrosly. 
#' 
#' A separate plotting function is provided as described in 
#' \code{\link[DrewDayRFunctions]{plot.medcurve}}.
#' 
#' @param medmodel This is a mediation model of any class acceptable by the 
#' \code{\link[mediation]{mediate}} function, namely 'lm', 'polr', 'bayespolr', 
#' 'glm', 'bayesglm', 'gam', 'rq', 'survreg', or 'merMod'. This model should 
#' include the treatment as an independent variable and the mediator as a 
#' dependent variable. 
#' @param outmodel This is a mediation model of any class acceptable by the 
#' \code{\link[mediation]{mediate}} function as explained above. This model 
#' should include the treatment and the mediator as independent variables and 
#' the outcome as a dependent variable.
#' @param treatname This is a character value for the column name of the 
#' treatment variable as it was called in the formula for \code{medmodel} and 
#' \code{outmodel}.
#' @param medname This is a character value for the column name of the 
#' mediator variable as it was called in the formula for \code{medmodel} and 
#' \code{outmodel}.
#' @param treat_values An optional numeric vector of treatment values to 
#' iterate across. This defaults to \code{NULL}, meaning that a vector of 
#' treatment values is created using the 'grid' arguments below. If a vector is 
#' provided for \code{treatvalues}, those 'grid' arguments are ignored.
#' @param gridlen A length for the vector of treatment values to be generated. 
#' This will be passed to the \code{length.out} argument in 
#' \code{\link[base]{seq}}. This defaults to \code{9}.
#' @param gridlowerq A lower quantile for the beginning of the sequence of 
#' treatment values. This defaults to \code{0.01}, meaning that the first value 
#' is the 1st percentile of the treatment values. 
#' @param gridupperq An upper quantile for the end of the sequence of 
#' treatment values. This defaults to \code{0.99}, meaning that the last value 
#' is the 99th percentile of the treatment values. 
#' @param qgrid This is a logical value determining whether the sequence of 
#' treatment values should be of quantiles if \code{TRUE} or of evenly spaced 
#' values on the original scale if \code{FALSE}. Here is an example if 
#' \code{gridlen = 5}, the treatment is normally distributed with a mean of 0 
#' and an SD of 1, \code{gridlowerq = 0.025}, and \code{gridupperq = 0.975}. 
#' If \code{qgrid = TRUE}, then the grid would be the 2.5th, 26.25th, 50th, 
#' 73.75th, and 97.5th percentiles of the treatment values, which in this case 
#' would be -1.96, -0.64, 0.00, 0.64, and 1.96. If \code{qgrid = FALSE}, these 
#' values would instead be evenly spaced between -1.96 and 1.96, so they would 
#' be -1.96, -0.98, 0, 0.98, 1.96. This defaults to \code{FALSE} because the 
#' smoothing at the edges of the curve tends to more severely miss the actual 
#' point estimates when not having an evenly spaced set of values. 
#' @param progress A logical value as to whether to include a progress bar 
#' showing the percentage of the mediation functions across the treatment values 
#' that has been completed so far. This function can take some time to run, so 
#' this can be pretty useful to see how far the function has progressed. This 
#' defaults to \code{TRUE}. 
#' @param ... These are additional arguments to be passed to the 
#' \code{\link[mediation]{mediate}} function. One important argument is 
#' \code{boot}, which is a logical value determining whether bootstrap CIs are 
#' calculated. This defaults to \code{FALSE}, but it is recommended to set it 
#' to \code{TRUE}. The number of bootstrap iterations is controlled with the 
#' argument \code{sims}, which defaults to 1000.
#' 
#' @returns \code{MediationCurve} returns a data frame containing the output of 
#' the \code{\link[mediation]{mediate}} function iterations. This data frame 
#' includes the following columns:
#' \item{control.value, treat.value}{These are the 'control' and 'treated' 
#' values for each pair of treatment values.}
#' \item{X, X.ci.2.5, X.ci.97.5, X.p}{Here 'X' is a placeholder. These are the 
#' estimate of interest, lower 95% CI, upper 95% CI, and p-value, 
#' respectively. These columns are included for all the estimates listed below.}
#' \item{d0, d1, d.avg}{Point estimates for the average causal mediation effect 
#' (ACME) under the control condition, treatment condition, and the average 
#' between the two. Note that if a GAM is used for the outcome model, these 
#' values will all be equal.}
#' \item{z0, z1, z.avg}{The same as above but for the average direct effect 
#' (ADE).}
#' \item{tau}{The total effect of the treatment on the outcome.}
#' \item{n0, n1, n.avg}{The same as for the ACME and ADE, but this is the 
#' proportion mediated (i.e., the ACME divided by the total effect).}
#' 
#' @import mediation
#' @export MediationCurve
#' 
#' @examples
#' # Linear treatment -> outcome association, nonlinear treatment -> mediator 
#' #  and mediator -> outcome association
#' 
#' \dontrun{
#' # Make a simple data frame
#' set.seed(90)
#' meddata <- data.frame(treat = rnorm(500))
#' meddata$medhat <- -2* meddata$treat^2
#' set.seed(101)
#' meddata$med <- meddata$medhat + rnorm(500)
#' meddata$outhat <- meddata$treat + -2* meddata$med^2
#' set.seed(11)
#' meddata$out <- meddata$outhat + rnorm(500)
#' 
#' #Run the mediation and outcome models
#' medmodel <- mgcv::gam(med ~ s(treat), data = meddata)
#' outmodel <- mgcv::gam(out ~ treat + s(med), data = meddata)
#' 
#' # Mediation with bootstrap CIs
#' medres <- MediationCurve(medmodel, outmodel, "treat", "med", boot = T)
#' 
#' plot.medres(medres)
#' }
#' 
#' @references 
#' Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal 
#' Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), 
#' pp. 309-334.
#' 
#' Kobrosly, R. W., (2020). causal-curve: A Python Causal Inference Package to 
#' Estimate Causal Dose-Response Curves. Journal of Open Source Software, 5(52), 
#' 2523, https://doi.org/10.21105/joss.02523.
#' 
MediationCurve <- function(medmodel, outmodel, treatname, medname, 
                           treat_values = NULL, gridlen = 9, gridlowerq = 0.01, 
                           gridupperq = 0.99, qgrid = F, progress = T, ...){
  rescols <- c("control.value", "treat.value", 
               c(sapply(c("d0", "d1", "z0", "z1", "n0", "n1"), 
                        function(x) paste0(x, c("", ".ci", ".p")))), 
               paste0("tau", c(".coef", ".ci", ".p")), 
               c(sapply(c("d", "z", "n"), 
                        function(x) paste0(x, c(".avg", ".avg.ci", ".avg.p")))))
  rescolslong <- c("control.value", "treat.value", 
                   c(sapply(c("d0", "d1", "z0", "z1", "n0", "n1"), 
                            function(x) paste0(x, c(
                              "", ".ci.2.5%", ".ci.97.5%", ".p")))), 
                   paste0("tau", c(".coef", ".ci.2.5%", ".ci.97.5%", ".p")), 
                   c(sapply(c("d", "z", "n"), 
                            function(x) paste0(x, c(
                              ".avg", ".avg.ci.2.5%", ".avg.ci.97.5%", 
                              ".avg.p")))))
  
  resdat <- as.data.frame(matrix(NA, gridlen - 1, length(rescolslong)))
  names(resdat) <- rescolslong
  if(!is.null(treat_values)){
    treatspread <- treat_values
    gridlen <- length(treat_values)
  } else {
    if(qgrid){
      treatspread <- quantile(outmodel$model[, treatname], 
                              probs = seq(gridlowerq, gridupperq, 
                                          length.out = gridlen))
    } else {
      treatspread <- seq(quantile(outmodel$model[, treatname], probs = gridlowerq), 
                         quantile(outmodel$model[, treatname], probs = gridupperq),
                         length.out = gridlen)
    }
  }
  
  if(progress){
    pb <- txtProgressBar(min = 1, max = gridlen - 1, initial = 1, style = 3)
  }
  for(ii in 1:(gridlen - 1)){
    tempmed <- suppressMessages(mediation::mediate(
      medmodel, outmodel, treat = treatname, 
      mediator = medname, control.value = treatspread[ii], 
      treat.value = treatspread[ii + 1], ...))
    if(progress) setTxtProgressBar(pb, ii)
    resdat[ii, ] <- unlist(tempmed[rescols])
    rm(tempmed)
  };rm(ii)
  if(progress) close(pb)
  names(resdat) <- gsub("%", "", names(resdat), fixed = T)
  names(resdat) <- gsub("tau.coef", "tau", names(resdat))
  resdat <- resdat[, c(1:10, 31:34, 11:18, 35:38, 27:30, 19:26, 
                       39:ncol(resdat))]
  class(resdat) <- c("medcurve", class(resdat))
  return(resdat)
}

