#' Visualize continuous by continuous interactions
#' 
#' \code{InteractionCoefPlot} produces a plot of how the coefficient between 
#' one variable in a bivariate continuous interaction and the outcome changes 
#' over levels of the other continuous variable. This is a useful visualization 
#' of bivariate interactions between two continuous variables. 
#' 
#' @param model A model of class \code{\link[stats]{lm}} or \code{\link[stats]{glm}}.
#' @param data The data frame used to generate \code{model}.
#' @param ixterm The column name of the interacting variable in the data frame 
#' \code{data}. 
#' @param multiplier A multiplier for the coefficients of interest (e.g., an IQR 
#' of the predictor of interest). Defaults to 1.
#' @param coeftransform An optional transformation function for the coefficients 
#' (e.g., \code{tenfoldperc} if the predictor of interest is on a log10 scale 
#' where \code{tenfoldperc <- function(x) ((10^x) - 1) * 100)}. Defaults to 
#' \code{NULL}.
#' @param predname An optional new name for the variable \code{pred}. Defaults to 
#' \code{NULL}.
#' @param ixname An optional new name for \code{ixterm}. Defaults to \code{NULL}.
#' @param outname An optional new name for the outcome variable. Defaults to 
#' \code{NULL}.
#' @param title An optional plot title. Defaults to \code{NULL}.
#' @param fillcolor The color for the shading of the 95% CI region. Defaults to 
#' \code{"grey50"}.
#' @param autotitle An option to automatically generate a plot title. Defaults to 
#' \code{FALSE}.
#' @param addpvallab An option to add a label to the plot that shows the 
#' interaction term p-value. Defaults to \code{FALSE}.
#' @param labsize The label text size if \code{addpvallab} is \code{TRUE}. Defaults
#' to 4.
#' @param lengthout The number of levels of \code{ixterm} generated with the 
#' \code{length.out} argument for the \code{\link[base]{seq}} function, defaults 
#' to 50.
#' @param shadebysig An option to shade the CI region areas that don't overlap 
#' the null a darker tint of \code{fillcolor}. Defaults to \code{FALSE}.
#' @param otherix An option to include a second plot in which \code{pred} and 
#' \code{ixterm} are flipped. Defaults to \code{FALSE}.
#' @param robust This is an option to use robust sandwich errors for the CIs
#' (see \code{\link[sandwich]{vcovHC}} for details on HC estimators). Defaults to 
#' \code{TRUE}. This should be set to \code{TRUE} if model is a robust Poisson 
#' model for RR estimation.
#' @param HC This is a character value defining the robust estimator used if 
#' \code{robust = TRUE} (see \code{\link[sandwich]{vcovHC}} for details on HC 
#' estimators). This defaults to \code{"HC0"}.
#' @param logistic An option for if model is a logistic \code{\link[stats]{glm}}. 
#' Defaults to \code{FALSE}. This option sets the null line to 1, exponentiates 
#' all the coefficients to bring them to an odds ratio scale, and changes the y 
#' label to "OR" instead of "Coefficient".
#' @param xlab10exp An option to change the x-axis labels to be of the form 
#' "10^x" (see \code{\link[scales]{label_math}}). This is useful if your 
#' predictor of interest (or also the \code{ixterm}
#' variable if \code{otherix} is \code{TRUE}) is on a logarithmic scale.
#' 
#' @details 
#' Plot axis labels, fonts, etc. can be changed with standard ggplot options. 
#' For example, if the model is a robust Poisson GLM, you can set 
#' \code{logistic = TRUE} when running this function and saving it into an 
#' object that in this example I will call \code{plotlist} and then change the 
#' "OR" label to "RR" by running \code{plotlist$MainGGplot + ylab("Main Variable 
#' RR")} or something similar.
#' 
#' @return \code{InteractionCoefPlot} returns a list of: 
#' \item{MainGGplot}{The interaction plot}
#' \item{OtherGGplot}{The other interaction plot if \code{otherix = TRUE}}
#' \item{MainMatrix}{The matrix of values used to generate the main plot}
#' \item{OtherMatrix}{The matrix of values used to generate the other plot if 
#' \code{otherix = TRUE}}
#' 
#' @import lmtest sandwich ggplot2 msm scales
#' @export InteractionCoefPlot
#' 
#' @examples
#' 
#' set.seed(1)
#' exdat <- as.data.frame(mvtnorm::rmvnorm(500, mean=rep(0,5)))
#' names(exdat) <- paste0("X", 1:ncol(exdat))
#' exdat$X1X2 <- exdat$X1 * exdat$X2
#' exdat$yhat <- c(-1 + as.matrix(exdat) %*% c(0.5, 1.5, 0, -2, 2, 1))
#' set.seed(2)
#' exdat$y <- exdat$yhat + rnorm(500,0,5)
#' exlm <- lm(y ~ X1 * X2 + X3 + X4 + X5, exdat)
#' explot <- InteractionCoefPlot(exlm, exdat, "X1", "X2", addpvallab=T, 
#' shadebysig=T, robust=F, otherix=T)
#' cowplot::plot_grid(explot$MainGGplot, explot$OtherGGplot, ncol=1)
#' 

InteractionCoefPlot <- function(model, data, pred, ixterm, multiplier = 1, 
                                coeftransform = NULL, predname = NULL, 
                                ixname = NULL, outname = NULL, title = NULL, 
                                fillcolor = NULL, autotitle = F, addpvallab = F, 
                                labsize = 4, lengthout = 50, shadebysig = F, 
                                otherix = F, robust = T, HC = "HC0", 
                                logistic = F, xlab10exp = F){
  if(is.null(outname)) outname <- names(model$model)[1]
  if(is.null(predname)) predname <- pred
  if(is.null(ixname)) ixname <- ixterm
  if(is.null(fillcolor)) fillcolor <- "grey50"
  predcoefdat <- data.frame(Outcome = outname, Predictor = predname, 
                            Interaction = ixname, IxLevel = seq(
                              min(data[, ixterm], na.rm = T), 
                              max(data[, ixterm], na.rm = T), 
                              length.out = lengthout))
  ixcoeflabel <- ifelse(paste0(pred, ":", ixterm) %in% names(coef(model)), 
                      paste0(pred, ":", ixterm), paste0(ixterm, ":", pred))
  predcoefdat$Beta <- coef(model)[pred] + predcoefdat$IxLevel * coef(model)[ixcoeflabel]
  modcoef <- coef(model)
  if(robust){
    vc <- vcovHC(model, type = HC)
  } else {
    vc <- vcov(model)
  }

  predcoefdat$SE <- sapply(predcoefdat$IxLevel, function(x){
    msm::deltamethod(formula(paste0("~x", which(names(modcoef) == pred), " * ", multiplier, " + x", 
                                    which(names(modcoef) == ixcoeflabel), " * ", x, " * ", multiplier)), 
                     modcoef, vc)
  })

  predcoefdat$LCI <- predcoefdat$Beta - (predcoefdat$SE * 1.96)
  predcoefdat$UCI <- predcoefdat$Beta + (predcoefdat$SE * 1.96)
  if(!is.null(coeftransform)){predcoefdat[, c("Beta", "LCI", "UCI")] <- 
    lapply(predcoefdat[, c("Beta", "LCI", "UCI")], coeftransform)}
  if(logistic){
    predcoefdat[, c("Beta", "LCI", "UCI")] <- lapply(predcoefdat[
      , c("Beta", "LCI", "UCI")], function(x) exp(x))
    ylabaux <- " OR"
    yint <- 1
  } else {
    ylabaux <- " Coefficient"
    yint <- 0
  }

  rugdat <- data.frame(RugX = data[, ixterm])
  if(shadebysig){
    lowshadedat <- predcoefdat[predcoefdat$UCI < 0 & predcoefdat$LCI < 0, ]
    highshadedat <- predcoefdat[predcoefdat$LCI > 0 & predcoefdat$UCI > 0, ]
    g1 <- ggplot(predcoefdat, aes(x = IxLevel, y = Beta)) + geom_line() + theme_bw() + 
      geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, fill = fillcolor) + 
      geom_hline(yintercept = yint, linetype = "dashed") + 
      geom_ribbon(data = lowshadedat, aes(x = IxLevel, ymin = LCI, ymax = UCI), 
                  inherit.aes = F, alpha = 0.3, fill = fillcolor) + 
      geom_ribbon(data = highshadedat, aes(x = IxLevel, ymin = LCI, ymax = UCI), 
                  inherit.aes = F, alpha = 0.3, fill = fillcolor) + 
      geom_rug(data = rugdat, aes(x = RugX), inherit.aes = F) + 
      ylab(paste0(predname, ylabaux)) + xlab(ixname)
  } else {
    g1 <- ggplot(predcoefdat, aes(x = IxLevel, y = Beta)) + geom_line() + theme_bw() + 
      geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, fill = fillcolor) + 
      geom_hline(yintercept = yint, linetype = "dashed") + 
      geom_rug(data = rugdat, aes(x = RugX), inherit.aes = F) + 
      ylab(paste0(predname, ylabaux)) + xlab(ixname)
  }

  if(is.null(title)&autotitle){
    g1 <- g1 + ggtitle(paste0("Change in ", outname, " Per ", multiplier, 
                          " Unit Increase in\n", predname, " Over Levels of ", ixname)) + 
      theme(plot.title = element_text(hjust = 0.5))
  } else if(!is.null(title)){
    g1 <- g1 + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  }
  if(addpvallab){
    if(robust){
      robustcoef <- robustse(model, HCtype = HC)
      templab <- paste0("Interaction p = ", round(robustcoef[
        paste0(pred, ":", ixterm), 4], 3))
      if(robustcoef[paste0(pred, ":", ixterm), 4] < 0.05) templab <- paste0(
        templab, " * ")
    } else {
      templab <- paste0("Interaction p = ", round(summary(model)$coef[
        paste0(pred, ":", ixterm), 4], 3))
      if(summary(model)$coef[paste0(pred, ":", ixterm), 4] < 0.05){
        templab <- paste0(templab, " * ")
      }
    }

    g1 <- g1 + annotate(geom = "text", x = mean(range(data[, ixterm])), 
                        y = max(predcoefdat$UCI), 
                        label = templab, size = labsize, hjust = 0.5, vjust = 1)
    if(xlab10exp){
      g1 <- g1 + scale_x_continuous(labels = scales::label_math(10^.x))
    }
  }
  retlist <- list(GGplot = g1, Matrix = predcoefdat)

  if(otherix){
    othercoefdat <- data.frame(Outcome = outname, Predictor = ixname, Interaction = predname, 
                            IxLevel = seq(min(data[, pred], na.rm = T), 
                                        max(data[, pred], na.rm = T), 
                                        length.out = lengthout))
    othercoefdat$Beta <- coef(model)[ixterm] + othercoefdat$IxLevel * coef(model)[
      paste0(pred, ":", ixterm)]
    modcoef <- coef(model)
    if(robust){
      vc <- vcovHC(model, type = HC)
    } else {
      vc <- vcov(model)
    }
    othercoefdat$SE <- sapply(othercoefdat$IxLevel, function(x) {
      msm::deltamethod(formula(paste0("~x", which(names(modcoef) == ixterm), 
                                      " * ", multiplier, " + x", 
                                      which(names(modcoef) == paste0(
                                        pred, ":", ixterm)), " * ", x, " * ", multiplier)), 
                       modcoef, vc)
    })
    othercoefdat$LCI <- othercoefdat$Beta - (othercoefdat$SE * 1.96)
    othercoefdat$UCI <- othercoefdat$Beta + (othercoefdat$SE * 1.96)

    if(!is.null(coeftransform)){othercoefdat[, c("Beta", "LCI", "UCI")] <- 
      lapply(othercoefdat[, c("Beta", "LCI", "UCI")], coeftransform)}
    if(logistic){
      othercoefdat[, c("Beta", "LCI", "UCI")] <- lapply(othercoefdat[
        , c("Beta", "LCI", "UCI")], function(x) exp(x))
    }
    if(logistic){othercoefdat[, c("Beta", "LCI", "UCI")] <- 
      lapply(othercoefdat[, c("Beta", "LCI", "UCI")], function(x) exp(x))}
    otherrugdat <- data.frame(RugX = data[, pred])
    if(shadebysig){
      lowshadedat <- othercoefdat[othercoefdat$UCI < 0 & othercoefdat$LCI < 0, ]
      highshadedat <- othercoefdat[othercoefdat$LCI > 0 & othercoefdat$UCI > 0, ]
      g2 <- ggplot(othercoefdat, aes(x = IxLevel, y = Beta)) + geom_line() + theme_bw() + 
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, fill = fillcolor) + 
        geom_hline(yintercept = yint, linetype = "dashed") + 
        geom_ribbon(data = lowshadedat, aes(x = IxLevel, ymin = LCI, ymax = UCI), 
                    inherit.aes = F, alpha = 0.3, fill = fillcolor) + 
        geom_ribbon(data = highshadedat, aes(x = IxLevel, ymin = LCI, ymax = UCI), 
                    inherit.aes = F, alpha = 0.3, fill = fillcolor) + 
        geom_rug(data = otherrugdat, aes(x = RugX), inherit.aes = F) + 
        ylab(paste0(ixname, ylabaux)) + xlab(predname)
    } else {
      g2 <- ggplot(othercoefdat, aes(x = IxLevel, y = Beta)) + geom_line() + theme_bw() + 
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, fill = fillcolor) + 
        geom_hline(yintercept = yint, linetype = "dashed") + 
        geom_rug(data = otherrugdat, aes(x = RugX), inherit.aes = F) + 
        ylab(paste0(ixname, ylabaux)) + xlab(predname)
    }

    if(is.null(title)&autotitle){
      g2 <- g2 + ggtitle(paste0("Change in ", outname, " Per ", multiplier, 
                            " Unit Increase in\n", predname, 
                            " Over Levels of ", ixname)) + 
        theme(plot.title = element_text(hjust = 0.5))
    } else if(!is.null(title)){
      g2 <- g2 + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    }
    if(addpvallab){
      if(robust){
        templab <- paste0("Interaction p = ", round(robustcoef[paste0(pred, ":", ixterm), 4], 3))
        if(robustcoef[paste0(pred, ":", ixterm), 4]<0.05) templab <- paste0(templab, " * ")
      } else {
        templab <- paste0("Interaction p = ", round(summary(model)$coef[
          paste0(pred, ":", ixterm), 4], 3))
        if(summary(model)$coef[paste0(pred, ":", ixterm), 4]<0.05) templab <- paste0(templab, " * ")
      }

      g2 <- g2 + annotate(geom = "text", x = mean(range(data[, pred])), y = max(othercoefdat$UCI), 
                      label = templab, size = labsize, hjust = 0.5, vjust = 1)
    }
    if(xlab10exp){
      g2 <- g2 + scale_x_continuous(labels = scales::label_math(10^.x))
    }
    retlist <- list(MainGGplot = g1, OtherGGplot = g2, MainMatrix = predcoefdat, 
                  OtherMatrix = othercoefdat)
  }
  return(retlist)
}
