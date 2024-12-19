#' Run multiple linear regressions and return results
#' 
#' \code{GLMResults} is a complicated helper function to handle all sorts of 
#' repeated linear regression analyses. It outputs a table of key regression 
#' output including optional diagnostic test values, a list of the model 
#' objects, and a coefficient forest plot or plots. Supported regression types 
#' include linear regression (\code{\link[stats]{lm}}) and GLMs 
#' (\code{\link[stats]{glm}} with any \code{\link[stats]{family}} or 
#' \code{\link[MASS]{glm.nb}}) or Firth's logistic 
#' regression using \code{\link[logistf]{logistf}}. Furthermore, multiple 
#' imputation through MICE is supported for all but the \code{logistf} models, 
#' and all binary contrast coefficients (e.g., \code{level 2 - level 1}, 
#' \code{level 3 - level 2}, etc.) are returned for categorical predictors and 
#' interaction terms (the latter using 
#' \code{\link[DrewDayRFunctions]{ContrastCoefficients}}). Additional options 
#' include estimating robust sandwich errors, incorporating continuous or 
#' categorical bivariate multiplicative interactions, an option to get k-fold 
#' cross validation RMSE or accuracy, and a sensitivity test where the model 
#' is run after omitting high leverage points as determined by Cook's distance.
#' 
#' @param prednames A character vector of the column names of predictors of 
#' interest you'd like to evaluate in separate regressions. These column names 
#' can refer to character, factor, or numeric columns in the data frame.
#' @param outnames A character vector of the column names of outcome variables 
#' you'd like to evaluate in separate regressions. These column names can refer 
#' to character, factor, or numeric columns in the data frame. For each 
#' \code{outnames}, \code{\link[stats]{lm}} linear regressions will be 
#' run if the referenced column is numeric, or a \code{\link[stats]{glm}} 
#' GLM with a \code{\link[stats]{binomial}} family will be run with the 
#' link function specified in the input \code{binomlink}. This does not yet 
#' support characters or factors with >2 unique values (i.e., it can only 
#' accommodate continuous or binary outcome variables).
#' @param covnames Either a character vector of column names for all the 
#' covariates you'd like to include in each regression or a list of character 
#' vectors equal in length to the length of the unique combination of prednames 
#' and \code{outnames}. This defaults to \code{NULL}, meaning no covariates are 
#' included. This can be a list in order to control for separate covariates for 
#' each unique combination of \code{prednames} and \code{outnames}, in which 
#' case the list must be of length \code{length(prednames) * length(outnames)}. 
#' Note that in ordering this list the function loops through the 
#' \code{prednames} and then the \code{outnames} (e.g., Outcome 1 - Predictor 1, 
#' Outcome 1 - Predictor 2, Outcome 1 - Predictor 3, 
#' Outcome 2 - Predictor 1, etc.), and so you should order the list of 
#' covariates accordingly. These can include spline terms, namely
#' \code{\link[splines]{ns}} natural splines or \code{\link[splines]{bs}} 
#' b-splines as defined by the 'splines' package (e.g., 
#' \code{"ns(Year, df = 3)"}). These can also include interaction terms denoted 
#' by the ":" separator , though the user should make sure to include the main 
#' effects too (e.g., 
#' \code{covnames = c("Income", "HouseholdSize", "Income:HouseholdSize")}).
#' @param Data A data frame object of class "data.frame". All variables 
#' defined in \code{prednames}, \code{outnames}, \code{covnames}, and 
#' \code{ixterm} should be present in this data frame.
#' @param logout If \code{TRUE}, this will log-transform each outcome variable 
#' (the log base will be \code{logbaseout}) prior to running the regression. 
#' This can be either a single value or a logical vector of length equal to the 
#' length of \code{outnames}. This defaults to \code{FALSE}.
#' @param logpred If \code{TRUE}, this will log-transform each predictor 
#' variable (the log base will be \code{logbasepred}) prior to running the 
#' regression. This can be either a single value or a logical vector of length 
#' equal to the length of \code{prednames}. This defaults to \code{FALSE}.
#' @param logbasepred The base for the log-transformation of one or more of the 
#' predictor variables. This defaults to 10.
#' @param logbaseout The base for the log transformation of the outcome 
#' variables. Defaults to \code{exp(1)} (i.e., natural log transformation).
#' @param Outtitle Defines the title used for the column of outcome variables 
#' in the table and on the plot. This defaults to \code{"Outcome"}.
#' @param Predtitle Defines the title used for the column of predictor variables 
#' in the table and on the plot. This defaults to \code{"Exposure"}.
#' @param ixterm A column name for an optional interaction term for the 
#' \code{prednames} predictors of interest (I abbreviate 'interaction' as 'ix'). 
#' This can refer to either a character, factor, or numeric-class column in the 
#' \code{Data} data frame. All relevant interaction and main effect coefficients 
#' will be output. This defaults to \code{NULL}.
#' @param Firth If \code{TRUE} and if a given \code{outnames} variable is 
#' binary, this will perform a Firth's logistic regression 
#' (\code{\link[logistf]{logistf}}) instead of a GLM 
#' (\code{\link[stats]{glm}}). This defaults to \code{FALSE}.
#' @param altprednames This is a vector of character values for alternate names
#' for the predictor variable column names specified in \code{prednames} to 
#' replace those column names in the output tables and plots. This defaults to 
#' \code{NULL}.
#' @param altoutnames This is a vector of character values for alternate names
#' for the outcome variable column names specified in \code{outnames} to 
#' replace those column names in the output tables and plots. This defaults to 
#' \code{NULL}.
#' @param altixname This is a character value for an alternate name for the 
#' interaction variable column name specified in \code{ixterm} to replace that 
#' column name in the output tables and plots. This defaults to \code{NULL}.
#' @param binomfam The family function (see \code{\link[stats]{family}}) to 
#' use for GLMs of binary or count \code{outnames}. Options include 
#' \code{\link[stats]{binomial}}, \code{\link[stats]{poisson}}, 
#' \code{\link[stats]{quasibinomial}}, or \code{\link[stats]{quasipoisson}}. 
#' Refer to the help pages for those families to see options for link functions. 
#' If the Poisson family is used here, one should make sure to set 
#' \code{robust = TRUE} to perform a robust Poisson regression for estimating 
#' risk ratios. This defaults to \code{binomial(link = "logit")}. 
#' @param integerascount If \code{TRUE}, all outcome variables specified in 
#' \code{outnames} are treated as count variables if they're of class 
#' \code{"integer"}, thereby running a GLM for those outcomes with the family 
#' defined by \code{countfam}. If \code{FALSE}, all outcome variables of class 
#' \code{"integer"} are treated as continuous, meaning that a \code{lm} is run 
#' for those outcomes.
#' @param countfam The family function (see \code{\link[stats]{family}}) to 
#' use for GLMs of count \code{outnames} if \code{integerascount} is 
#' \code{TRUE}. Options include \code{\link[stats]{poisson}} or
#' \code{\link[stats]{quasipoisson}}. Refer to the help pages for those 
#' families to see options for link functions. An additional option is to 
#' provide the character value \code{"negbin"}, which indicates that a 
#' negative binomial regression (\code{\link[MASS]{glm.nb}}) will be used with 
#' the log link function. Note that zero-inflated models require some 
#' additional model specification and are currently not supported in this 
#' function. This defaults to \code{poisson("log")}. 
#' @param robust If \code{TRUE}, this calculates and returns robust CIs and 
#' p-values based on the \code{\link[sandwich]{vcovHC}} function instead of the 
#' default CIs and p-values. This defaults to \code{TRUE}.
#' @param robustdf If \code{FALSE}, this sets \code{\link[lmtest]{coefci}} to 
#' calculate z-score-based robust confidence intervals if \code{robust} is 
#' \code{TRUE}. If \code{TRUE}, t test-based robust confidence intervals are 
#' instead calculated based on the residual degrees of freedom in each 
#' regression if \code{robust} is \code{TRUE}. This defaults to \code{FALSE}, 
#' which is the default for the \code{df} parameter in 
#' \code{\link[lmtest]{coefci}}.
#' @param HCtype This is the formula used for the \code{\link[sandwich]{vcovHC}} 
#' function calculation of robust CIs and p-values. This defaults to 
#' \code{"HC0"}.
#' @param predspline If \code{TRUE}, this transforms each predictor 
#' variable in \code{prednames} that is not categorical into a spline term using 
#' the 'splines' package. A coefficient for each degree of freedom is returned. 
#' This defaults to \code{FALSE}. 
#' @param splinetype This specifies the type of spline function to apply to the 
#' continuous predictor variables in \code{prednames} if 
#' \code{predspline = TRUE}. Options include \code{"ns"} for natural splines 
#' (\code{\link[splines]{ns}}) or \code{"bs"} for b-splines 
#' (\code{\link[splines]{bs}}). This defaults to \code{"ns"}.
#' @param predsplinedf This specifies the number of degrees of freedom for the 
#' predictor variable splines if \code{predspline = TRUE}. This defaults to 3.
#' @param coeftrans_cont This is an optional function to transform the 
#' coefficient estimates for models with continuous outcomes (\code{lm} models) 
#' into another scale to be included in the results table and to be plotted. 
#' For example, if the outcome is on a log base 2 scale, one could set 
#' \code{coeftrans_cont = function(x) 2^x} to include coefficients on the 
#' original outcome scale in the results table and in the plots. If 
#' \code{coeftrans_cont = NULL} (the default), then only coefficients on the 
#' original scale will be included in the results table and plotted for all 
#' continuous outcomes.
#' @param coeftrans_binom This is an optional function to transform the 
#' coefficient estimates for models with binary outcomes into another scale to 
#' be included in the results table and to be plotted. An example would be if 
#' the outcome is binary and a logistic regression is run, one could set 
#' \code{coeftrans_binom = exp} to get coefficients on the odds ratio scale. 
#' If \code{coeftrans_binom = NULL}, then only coefficients on the original 
#' scale will be included in the results table and plotted for all 
#' binary outcomes. This defaults to \code{exp}.
#' @param coeftrans_count This is an optional function to transform the 
#' coefficient estimates for models with count outcomes into another scale to 
#' be included in the results table and to be plotted. An example would be if 
#' the outcome is count data and a Poisson regression is run, one could set 
#' \code{coeftrans_count = exp} to get coefficients on the incidence ratio 
#' scale. If \code{coeftrans_count = NULL}, then only coefficients on the 
#' original scale will be included in the results table and plotted for all 
#' count outcomes. This defaults to \code{exp}.
#' @param KFCV If \code{TRUE}, this performs k-fold cross validation 
#' using the 'caret' package. If a given model is a linear regression 
#' (\code{\link[stats]{lm}}), the RMSE will be reported in the tables. If a 
#' model is a GLM (\code{\link[stats]{glm}}), the accuracy will be reported 
#' in the tables. The train function input method is 
#' \code{trainControl(method  =  "repeatedcv", number  =  KF, repeats  =  50)}. 
#' This will not be run if the model is of class \code{logistf} or if MICE is 
#' performed. This defaults to \code{FALSE}. 
#' @param KF This is an integer value of the number of k-folds to perform if 
#' \code{KFCV} is \code{TRUE}. This defaults to 10 for 10-fold cross-validation.
#' @param facetcol This defines the number of facet_wrap columns for the 
#' coefficient plots (facetted by outcome variable). If \code{NULL}, this will 
#' equal the number of unique \code{outnames} values. This defaults to 
#' \code{NULL}.
#' @param covix This character vector of covariate column names allows for 
#' covariates to also have interaction terms with the \code{ixterm} interaction 
#' variable. If one wishes to have all covariates have interaction terms with 
#' \code{ixterm}, the character vector for \code{covix} should be the same as 
#' that for \code{covnames}. This defaults to \code{NULL}, meaning that there 
#' are no interactions with any of the covariates.
#' @param ixpred If \code{TRUE}, an interaction term is included between each of 
#' the \code{prednames} predictors of interest and \code{ixterm}. This can be 
#' set to \code{FALSE} to allow for interactions only with the covariates 
#' defined in \code{covix}. If \code{covix} is not \code{NULL} and 
#' \code{ixpred} is \code{TRUE}, there will be interaction terms with both the 
#' predictor of interest and the covariates. This defaults to \code{TRUE}.
#' @param extradiag If \code{TRUE}, this will include extra diagnostic 
#' information in the tables. For linear regressions (\code{\link[stats]{lm}}), 
#' this will include the number of Cook's distances > 0.5 and the 
#' maximum Cook's distance, both of which define high leverage points, as well 
#' as a p-value for the heteroskedasticity of model residuals 
#' (\code{\link[lmtest]{bptest}}). For GLMs (\code{\link[stats]{glm}}), this 
#' will include output from 
#' \href{https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html}{DHARMa package tests}, 
#' including the number of Cook's distances > 0.5, the maximum Cook's distance, 
#' a p-value testing the heteroskedasticity of residuals from 
#' \code{\link[DHARMa]{testQuantiles}}, the minimum p-value from 
#' \code{\link[DHARMa]{testQuantiles}}, a p-value for a test of uniformity 
#' (\code{\link[DHARMa]{testUniformity}}), a p-value for outliers 
#' (\code{\link[DHARMa]{testOutliers}}), and a p-value for dispersion 
#' (\code{\link[DHARMa]{testDispersion}}). If count outcome variables are 
#' included, an additional p-value for zero inflation 
#' (\code{\link[DHARMa]{testZeroInflation}}) will be included in the 
#' results. Note that both the \code{lm} and \code{glm} forms of the p-value 
#' for heteroskedasticity will be called \code{"Heterosk.p"} in the model 
#' output table. If the model is of class \code{logistf} or if MICE is used, no 
#' extra diagnositics will be performed. This defaults to \code{TRUE}.
#' @param leverage.test If \code{TRUE}, this runs all regressions while 
#' excluding "high leverage" observations as defined by the user. This defaults 
#' to \code{FALSE}.
#' @param leverage.cutoff This is the Cook's distance cutoff above which 
#' observations will be excluded if \code{leverage.test} is \code{TRUE}. 
#' This defaults to 0.2.
#' @param leverage.meancutoff If \code{leverage.cutoff} is set to \code{NULL}, 
#' this is used as a multiplier by the mean Cook's distance to get a new cutoff 
#' above which observations will be excluded if \code{leverage.test} is 
#' \code{TRUE}. For example, if set to 2, the cutoff will be set to the mean 
#' Cook's distance times 2. If both \code{leverage.cutoff} and 
#' \code{leverage.meancutoff} are set to \code{NULL}, the cutoff will be set 
#' to the mean Cook's distance times 4. Defaults to \code{NULL}.
#' @param mice If TRUE, this performs multiple imputation (MICE using 
#' \code{\link[mice]{mice}}) on the user-specified variables and pools model 
#' estimates across the imputed datasets. This defaults to \code{FALSE}.
#' @param micevars A character vector of the names of variables to be imputed 
#' using MICE. This defaults to \code{NULL}, in which case it will include all 
#' predictor, outcome, and covariate variables that have missing values.
#' @param miceiter The number of multiple imputations to perform if \code{mice} 
#' is set to \code{TRUE}. This value will be applied to both the \code{m} and 
#' \code{maxit} arguments of \code{\link[mice]{mice}}. This defaults to 10.
#' @param returnplot If \code{TRUE}, this function will also return a 
#' \code{ggplot} plot object summarizing the results, but it will not if this 
#' is set to \code{FALSE}. This is a forest plot of the coefficients of 
#' interest. If there is an interaction term, a list of \code{ggplot} objects 
#' labelled \code{"Marginal"} and \code{"Interaction"} will be returned, which 
#' showcase the marginal main effect coefficients and the interaction term 
#' coefficients, respectively. If multiple outcome types are included in 
#' \code{outnames}, a list of separate plots for each of these types will be 
#' returned. This defaults to \code{TRUE}.
#' @param colorbypred If \code{TRUE}, the forest plot will be colored by 
#' each predictor variable if there is no interaction or if that interaction is 
#' either continuous or categorical with only two levels. This is especially 
#' useful when multiple discrete predictors, each with several factor levels, 
#' have been included in the function. This defaults to \code{FALSE}.
#' @param colorpal A palette function for coloring the forest plot by predictor 
#' variables if \code{colorbypred = TRUE}. If \code{NULL}, the default ggplot2 
#' palette (\code{\link[scales]{hue_pal}}) will be used. This defaults to 
#' \code{NULL}. 
#' @param log10yaxis If \code{TRUE}, the y-axis of each plot is transformed to 
#' a log base 10 scale with labels \code{10^x}. This is useful if coefficients 
#' are on a multiplicative scale and have a wide range. This defaults to 
#' \code{FALSE}.
#' @param yl The y-axis title to apply to the forest plot(s). This defaults to 
#' \code{"Coefficient"}.
#' 
#' @details
#' \code{GLMResults} was built over several years and so has incorporated all 
#' sorts of functionality related to linear regressions that I've found useful. 
#' This is my main workhorse function for most analyses calling for regressions.
#' 
#' @return \code{GLMResults} returns a list of:
#' \item{Table}{A table of all the key variables. Coefficient, LCI, and UCI 
#' refer to the coefficient for the predictor of interest and the lower and 
#' upper 95% confidence intervals, respectively. CoefTrans, etc. are those 
#' estimates transformed in some way as specified by the user if 
#' \code{coeftrans} is not \code{NULL}, for instance by being exponentiated 
#' or transformed into percent changes.}
#' \item{ModelList}{A list of all the \code{lm}, \code{glm}, \code{glm.nb}, or 
#' \code{logistf} model objects run for each unique combination of predictor and 
#' outcome variable names (i.e., \code{prednames} and \code{outnames}).}
#' \item{GGplot}{A forest ggplot of all predictor of interest coefficients and 
#' confidence intervals or a list of those ggplots.} 
#' 
#' @import caret lmtest sandwich logistf splines ggplot2 DHARMa pbapply cowplot
#' @importFrom purrr map
#' @importFrom MASS glm.nb
#' 
#' @export GLMResults
#' 
#' @examples
#' # Note that all examples simply showcase the use of different types of 
#' # variables with built-in example datasets, not necessarily the best methods 
#' # to analyze those particular datasets.
#' 
#' # OLS linear regression with a simple set of outcomes, predictors, and 
#' #  covariates without using robust errors
#' 
#' data("mtcars")
#' carpreds <- c("disp", "cyl", "drat", "wt")
#' carouts <- c("mpg", "hp")
#' carcovars <- c("am", "vs", "gear")
#' 
#' lmres <- GLMResults(prednames = carpreds, outnames = carouts, 
#'                    covnames = carcovars, Data = mtcars, robust = F)
#' View(lmres$Table)
#' lmres$GGplot
#' 
#' # Probit regression with an interaction term on the predictor of interest and
#' #  on the covariate "gear" while also using robust errors (the default).
#' 
#' data("infert")
#' infert$anyspont <- ifelse(infert$spontaneous > 0, 1, 0)
#' infert$induced <- as.factor(infert$induced)
#' glmres <- GLMResults(c("induced", "age"), c("case", "anyspont"), 
#' c("education"), infert, ixterm = "parity", covix = "education", 
#' binomfam = binomial("probit"), altprednames = c("Induced", "Age"), 
#' altoutnames = c("Infertile", "Spontaneous"), altixname = "Parity", 
#' extradiag = T)
#'                    
#' # Visualize continuous * continuous interaction
#' interact.plot <- DrewDayRFunctions::InteractionCoefPlot(
#' glmres$LMlist$`anyspont~age`, infert, "age", "parity", 
#' addpvallab = T, shadebysig = T)
#' interact.plot$GGplot
#' 
#' # Robust Poisson regression with a binary outcome
#' glmres <- GLMResults(carpreds[-1], c("vs", "am"), c("gear", "carb"), mtcars, 
#'                      binomfam  =  poisson("log"),  
#'                      altprednames  =  c("Displacement", "Cylinders", 
#'                      "Axle Ratio", "Weight"), 
#'                      altoutnames  =  c("Engine Type", "Transmission Type"))
#' glmres$GGplot
#' 
#' # Mix of continuous outcomes, binary outcomes with logistic GLMs, 
#' # and count outcomes with negative binomial GLMs, all with continuous and 
#' # categorical predictors and a categorical interaction term.
#' infert$parity <- as.integer(infert$parity)
#' set.seed(6)
#' infert$contout <- infert$age * 0.1 + rnorm(nrow(infert))
#' glmres <- GLMResults(c("age", "induced"), c("parity", "anyspont", "contout"), 
#' NULL, infert, ixterm = "education", countfam = "negbin", integerascount = T, 
#' extradiag = F, coeftrans_binom = NULL, altprednames  =  c("Age", "Induced"), 
#' altoutnames  =  c("Parity", "Spontaneous", "Continuous Outcome"), 
#' altixname  =  "Education")
#' if(require(cowplot)){
#' plotlist <- lapply(glmres$GGplot, 
#' function(x) cowplot::plot_grid(plotlist = x, nrow = 1))
#' cowplot::plot_grid(plotlist = plotlist, ncol = 1, labels = names(plotlist), 
#' hjust = 0)
#' }
#' 
#' # MICE example
#' # note that robust errors are not attainable when using MICE, so I use 
#' # robust  =  F
#' mtcars_miss <- mtcars
#' set.seed(101)
#' mtcars_miss[sample(1:nrow(mtcars_miss), 3), "cyl"] <- NA
#' set.seed(104)
#' mtcars_miss[sample(1:nrow(mtcars_miss), 4), "wt"] <- NA        
#' 
#' lmres <- GLMResults(prednames  =  carpreds, outnames  =  carouts, 
#'                    covnames  =  carcovars, Data  =  mtcars, mice  =  T, 
#'                    micevars  =  c("cyl", "wt"), robust  =  F, extradiag = F)
#' lmres$GGplot  
#'             
GLMResults <- function(prednames, outnames, covnames = NULL, Data, logout = F, 
                       logpred = F, logbasepred = 10, logbaseout = exp(1), 
                       Outtitle = "Outcome", Predtitle = "Exposure", 
                       ixterm = NULL, Firth = F, altprednames = NULL, 
                       altoutnames = NULL, altixname = NULL, 
                       binomfam = binomial("logit"), integerascount = F,
                       countfam = poisson("log"), robust = T, robustdf = F, 
                       HCtype = "HC0", predspline = F, splinetype = "ns", 
                       predsplinedf = 3, coeftrans_cont = NULL, 
                       coeftrans_binom = exp, coeftrans_count = exp, 
                       KFCV = F, KF = 10, facetcol = NULL, covix = NULL, 
                       ixpred = T, extradiag = T, 
                       leverage.test = F, leverage.cutoff = 0.2, 
                       leverage.meancutoff = NULL, mice = F, micevars = NULL, 
                       miceiter = 10, returnplot = T, colorbypred = F, 
                       colorpal = NULL, log10yaxis = F, yl = "Coefficient"){
  #checks and warnings
  if(!is.null(ixterm)){
    if(predspline) stop(paste0("This function is not yet set up to accommodate 
                               spline interactions. Please set 'ixterm' to NULL 
                               if 'predspline' is TRUE."))
  }
  if(is.null(ixterm)){
    if(!is.null(altixname)) altixname <- NULL
    if(!is.null(covix)) covix <- NULL
  }

  if(robust & Firth){
    warning(paste0("This function is not yet set up to extract robust sandwich",
                   " errors from logistf objects. The 'robust' argument will",
                   " be ignored."))
    robust <- F
  }
  if(Firth & extradiag){
    warning(paste0("This function is not yet set up to perform diagnostics on",
                   " logistf objects. The 'extradiag' argument will be",
                   " ignored."))
    extradiag <- F
  }
  if(mice & extradiag){
    warning(paste0("This function is not yet set up to perform diagnostics on",
                   " lists of models run on MICE-generated datasets. The ",
                   "'extradiag' argument will be ignored."))
    extradiag <- F
  }
  if(Firth & mice){
    warning(paste0("The mice package is not set up to pool estimates across",
                " logistf model objects, and so 'mice' is being ignored here,",
                " and complete case Firth's logistic regressions will be ",
                "performed instead."))
    mice <- F
  }
  if(robust & mice){
    warning(paste0("It is unclear how to obtain robust errors from the list of",
                   " models produced by MICE, so the 'robust' argument", 
                   " is being ignored."))
    robust <- F
  }
  if(is.null(ixterm) & predspline){
    if(!all(class(Data[, prednames]) %in% c("numeric", "integer", "double"))) 
      stop(paste0("All prednames columns must be of class 'numeric', 'integer',",
                  " or 'double' if predspline is TRUE."))
  }
  if(predspline && !splinetype %in% c("ns", "bs")){
    stop(paste0("The argument 'splinetype' must be one of 'ns' or 'bs' if ",
                "the argument 'predspline' is TRUE."))
  }
  
  #if families are entered as functions, set them to families
  if(is.function(countfam)) countfam <- countfam()
  if(is.function(binomfam)) binomfam <- binomfam()
  
  #family check
  nbbool <- is.character(countfam) && countfam == "negbin"
  if(!nbbool && class(countfam) != "family"){
    warning(paste0("An unknown family was provided for the argument",
                   " 'countfam'. This argument is instead being set to ", 
                   "the default value of poisson('log')."))
    countfam <- poisson("log")
  }
  if(class(binomfam) != "family"){
    warning(paste0("An unknown family was provided for the argument",
                   " 'binomfam'. This argument is instead being set to ", 
                   "the default value of binomial('logit')."))
    countfam <- binomial("logit")
  }
  
  #set character predictor columns to factor
  if(length(prednames) > 1) predclasses <- sapply(Data[, prednames], class) else 
    predclasses <- class(Data[, prednames])
  if(any(predclasses == "character")){
    for(cl in which(predclasses == "character")){
      Data[, prednames[cl]] <- as.factor(Data[, prednames[cl]])
    }; rm(cl)
  }
  if(class(Data[, ixterm]) == "character") Data[, ixterm]<- 
      as.factor(Data[, ixterm])
  
  lmlist <- list() #empty list for model objects

  if(length(logout) == 1){
    logout2 <- rep(logout, length(outnames))
  } else {
    logout2 <- logout
  }

  if(length(logpred) == 1){
    logpred2 <- rep(logpred, length(prednames))
  } else {
    logpred2 <- logpred
  }
  
  #get outcome variable types
  outtype <- sapply(outnames, function(x) getvartype(x, Data, integerascount))

  for(i in 1:length(outnames)){
    for(j in 1:length(prednames)){
      lmnum <- (1 + (j - 1)) + ((i - 1) * length(prednames))
      # Alter outcome and predictor expressions
      if(logout2[i] & outtype[i] == "continuous"){
        outexpr <- paste0("log(", outnames[i], ", base = ", logbaseout, ")")
      } else {
        outexpr <- outnames[i]
      }

      if(logpred2[j] & predspline && predclasses[j] != "factor"){
        predexpr <- paste0(splinetype, "(log(", prednames[j], ", base = ", 
                           logbasepred, "), ", predsplinedf, ")")
      } else if(logpred2[j]  && predclasses[j] != "factor"){
        predexpr <- paste0("log(", prednames[j], ", base = ", logbasepred, ")")
      } else if(predspline && predclasses[j] != "factor"){
        predexpr <- paste0(splinetype, "(", prednames[j], ", ", 
                           predsplinedf, ")")
      } else {
        predexpr <- prednames[j]
      }
      
      #Build the formula
      if(is.null(covnames)){
        if(is.null(ixterm)){
          if(mice){
            ccvars <- c(outnames[i], prednames[j])
            if(is.null(micevars)){micevars <- ccvars[which(sapply(Data[
              , ccvars], function(x) any(is.na(x))))]}
            ccDat <- Data[complete.cases(Data[, ccvars[ - which(ccvars%in%micevars)]]), ]
            micedata <- mice(ccDat[, ccvars], m = miceiter, maxit = miceiter, printFlag = F)
          } else {
            ccDat <- Data[complete.cases(Data[, c(outnames[i], prednames[j])]), ]
          }
          nobs <- dim(ccDat)[1]
          form1 <- formula(paste0(outexpr, "~", predexpr))
        } else {
          if(mice){
            ccvars <- c(outnames[i], prednames[j], ixterm)
            if(is.null(micevars)){micevars <- ccvars[which(sapply(Data[
              , ccvars], function(x) any(is.na(x))))]}
            ccDat <- Data[complete.cases(Data[, ccvars[
              -which(ccvars %in% micevars)]]), ]
            micedata <- mice(ccDat[, ccvars], m = miceiter, 
                             maxit = miceiter, printFlag = F)
          } else {
            ccDat <- Data[complete.cases(Data[, c(
              outnames[i], prednames[j], ixterm)]), ]
          }
          nobs <- dim(ccDat)[1]
          form1 <- formula(paste0(outexpr, "~", predexpr, " * ", ixterm))
        }
      } else {
        if(is.list(covnames)){
          if(length(covnames) != (length(prednames) * length(outnames))){
            stop(paste0("Error: covnames list must have a length of ", 
                        length(prednames) * length(outnames)))
          }
          mycovnames <- covnames[[lmnum]]
        } else {
          mycovnames <- covnames
        }
        if(any(grepl("(", mycovnames, fixed = T)) | any(
          grepl(":", mycovnames, fixed = T))){
          mycovnames_df <- gsub(".*\\(", "", mycovnames)
          mycovnames_df <- gsub(",.*", "", mycovnames_df)
          mycovnames_df <- gsub("\\)", "", mycovnames_df)
          mycovnames_df <- mycovnames_df[-which(grepl(
            ":", mycovnames_df, fixed = T))]
        } else if(any(grepl(":", mycovnames, fixed = T))){
          mycovnames_df <- mycovnames[-which(grepl(":", mycovnames, fixed = T))]
        } else if(any(grepl("(", mycovnames, fixed = T))){
          mycovnames_df <- gsub(".*\\(", "", mycovnames)
          mycovnames_df <- gsub(",.*", "", mycovnames_df)
          mycovnames_df <- gsub("\\)", "", mycovnames_df)
        } else {
          mycovnames_df <- mycovnames
        }
        if(is.null(ixterm) & is.null(covix)){
          if(mice){
            ccvars <- c(outnames[i], prednames[j], mycovnames_df)
            if(is.null(micevars)){micevars <- ccvars[which(sapply(Data[
              , ccvars], function(x) any(is.na(x))))]}
            ccDat <- Data[complete.cases(Data[, ccvars[
              -which(ccvars%in%micevars)]]), ]
            micedata <- mice(ccDat[, ccvars], m = miceiter, 
                             maxit = miceiter, printFlag = F)
          } else {
            ccDat <- Data[complete.cases(Data[, c(outnames[i], prednames[j], 
                                                  mycovnames_df)]), ]
          }
          nobs <- dim(ccDat)[1]
          form1 <- formula(paste0(outexpr, " ~ ", predexpr, " + ", paste(
            mycovnames, collapse = " + ")))
        } else if(is.null(covix)){
          if(mice){
            ccvars <- unique(c(outnames[i], prednames[j], mycovnames_df, ixterm))
            if(is.null(micevars)){micevars <- ccvars[which(sapply(Data[
              , ccvars], function(x) any(is.na(x))))]}
            ccDat <- Data[complete.cases(Data[, ccvars[
              -which(ccvars%in%micevars)]]), ]
            micedata <- mice(ccDat[, ccvars], m = miceiter, 
                             maxit = miceiter, printFlag = F)
          } else {
            ccDat <- Data[complete.cases(Data[, c(outnames[i], prednames[j], 
                                                  mycovnames_df, ixterm)]), ]
          }
          nobs <- dim(ccDat)[1]
          form1 <- formula(paste0(outexpr, " ~ ", predexpr, " * ", ixterm, " + ", 
                                paste(mycovnames, collapse = " + ")))
        } else {
          if(mice){
            ccvars <- unique(c(outnames[i], prednames[j], mycovnames_df, 
                               ixterm))
            if(is.null(micevars)) micevars <- ccvars[which(sapply(Data[
              , ccvars], function(x) any(is.na(x))))]
            ccDat <- Data[complete.cases(Data[, ccvars[
              -which(ccvars %in% micevars)]]), ]
            micedata <- mice(ccDat[, ccvars], m = miceiter, maxit = miceiter, 
                             printFlag = F)
          } else {
            ccDat <- Data[complete.cases(Data[, c(outnames[i], prednames[j], 
                                                  mycovnames_df, ixterm)]), ]
          }
          nobs <- dim(ccDat)[1]
          covixaux <- rep("", length(mycovnames))
          covixaux[which(mycovnames %in% covix)] <- paste0(" * ", ixterm)
          for(cv in 1:length(mycovnames)){
            mycovnames[cv] <- paste0(mycovnames[cv], covixaux[cv])
          }
          if(ixpred == T){
            form1 <- formula(paste0(outexpr, " ~ ", predexpr, " * ", ixterm, " + ", 
                                  paste(mycovnames, collapse = " + ")))
          } else {
            form1 <- formula(paste0(outexpr, " ~ ", predexpr, " + ", 
                                  paste(mycovnames, collapse = " + ")))
          }
        }
      } #end formula building
      
      #build model
      if(outtype[i] == "binary"){
        if(binomfam$family %in% c("poisson", "quasipoisson")){
          if(is.character(Data[, outnames[i]])) Data[, outnames[i]] <- 
              as.factor(Data[, outnames[i]])
          if(is.factor(Data[, outnames[i]])){
            levels(Data[, outnames[i]]) <- c("0", "1")
            Data[, outnames[i]] <- as.numeric(as.character(Data[, outnames[i]]))
          }
        }
        if(Firth){ #logistf
          lm1 <- logistf(form1, data = Data, na.action = na.exclude)
          if(binomfam$family != "binomial"|binomfam$link != "logit"){
            warning(paste0("Note that only binomfam  =  binomial('logit')",
                           "  is accepted for Firth's logistic", 
                           " regressions, so the alternative binomfam is being",
                           " ignored while the argument 'Firth' is TRUE."))
          }
          if(leverage.test){
            warning(paste0("It's unclear how to get Cook's distance for a",
                           " 'logistf' object. Ignoring leverage.test while ",
                           "the argument 'Firth' is TRUE."))
          }
        } else { #binary outcome glm
          if(mice){
            lm1 <- list()
            tempdatalist <- purrr::map(1:micedata$m, 
                                       function(x) mice::complete(micedata, x))
            lm1 <- lapply(tempdatalist, function(x) glm(
              form1, x, family = binomfam, na.action = na.exclude))
          } else {
            lm1 <- tryCatch(glm(form1, data = Data, family = binomfam, 
                                na.action = na.exclude), 
                            error = function(e) {
                              warning(e)
                              warning(
                                paste0("Maybe you should try a different ",
                                       "binomfam family and/or link function."))
                              stop(paste0("Halting run at predname ", 
                                          prednames[j], " and outname ", 
                                          outnames[i]))
                            })
          }
          if(leverage.test){
            omithighcooks(lm1, mice, leverage.cutoff, 
                          leverage.meancutoff)
          } #end if leverage.test
        } #end if binary glm else
      } else if(outtype[i] == "count"){ #count
        if(is.character(countfam) && countfam == "negbin"){ #glm.nb
          if(mice){
            lm1 <- list()
            tempdatalist <- purrr::map(1:micedata$m, 
                                       function(x) mice::complete(micedata, x))
            lm1 <- lapply(tempdatalist, function(x) glm.nb(
              form1, data = x, na.action = na.exclude))
          } else {
            lm1 <- glm.nb(form1, data = Data, na.action = na.exclude)
          }
          if(leverage.test){
            omithighcooks(lm1, mice, leverage.cutoff, 
                          leverage.meancutoff)
          }
        } else { #count outcome glm
          if(mice){
            lm1 <- list()
            tempdatalist <- purrr::map(1:micedata$m, 
                                       function(x) mice::complete(micedata, x))
            lm1 <- lapply(tempdatalist, function(x) glm(
              form1, x, family = countfam, na.action = na.exclude))
          } else {
            lm1 <- tryCatch(glm(form1, data = Data, family = countfam, 
                       na.action = na.exclude), 
                       error = function(e) {
                         warning(e)
                         warning(paste0("Maybe you should try a different ",
                                        "countfam family and/or link function."))
                         stop(paste0("Halting run at predname ", prednames[j], 
                                     " and outname ", outnames[i]))
                       })
          }
          if(leverage.test){
            omithighcooks(lm1, mice, leverage.cutoff, 
                          leverage.meancutoff)
          }
        } #end if binary glm else
      } else { #lm
        if(mice){
          tempdatalist <- purrr::map(1:micedata$m, 
                                     function(x) mice::complete(micedata, x))
          lm1 <- lapply(tempdatalist, function(x) lm(
            form1, x, na.action = na.exclude))
        } else {
          lm1 <- lm(form1, data = Data, na.action = na.exclude)
        }

        if(leverage.test){
          omithighcooks(lm1, mice, leverage.cutoff, 
                        leverage.meancutoff)
        }
      } #end model building
      
      #Begin results collation
      lmlist[[lmnum]] <- lm1
      names(lmlist)[[lmnum]] <- paste0(outnames[i], "~", prednames[j])
      if(mice) resData <- tempdatalist else resData <- Data
      modcoef <- ContrastCoefficients( #DrewDayRFunctions::
        lm1, resData, prednames[j], ixterm = ixterm, robust = robust, 
        HCtype = HCtype, usedf = robustdf)
      mc_ncol <- ncol(modcoef)
      modcoef <- as.data.frame(modcoef)
      modcoef$Outcome <- outnames[i]
      modcoef$Predictor <- prednames[j]
      modcoef$Variable <- rownames(modcoef)
      rownames(modcoef) <- NULL
      modcoef$N <- nobs
      modcoef <- modcoef[, c("Outcome", "Predictor", "Variable", "N", 
                             names(modcoef)[1:mc_ncol])]
      
      #Additional transformations of coefficients
      if(outtype[i] == "continuous" && !is.null(coeftrans_cont)){
        if(!is.function(coeftrans_cont)){
          warning(paste0("'coeftrans_cont' must be a function, but the ", 
                         "provided value is not a function, and therefore ",
                         "'coeftrans_cont' will be ignored."))
          coeftrans_cont <- NULL
        } else {
          modcoef$TransCoef <- coeftrans_cont(modcoef$Coefficient)
          modcoef$TransLCI <- coeftrans_cont(modcoef$LCI)
          modcoef$TransUCI <- coeftrans_cont(modcoef$UCI)
        }
      } else if(outtype[i] == "binary" && !is.null(coeftrans_binom)){
        if(!is.function(coeftrans_binom)){
          warning(paste0("'coeftrans_binom' must be a function, but the ", 
                         "provided value is not a function, and therefore ",
                         "'coeftrans_binom' will be ignored."))
          coeftrans_binom <- NULL
        } else {
          modcoef$TransCoef <- coeftrans_binom(modcoef$Coefficient)
          modcoef$TransLCI <- coeftrans_binom(modcoef$LCI)
          modcoef$TransUCI <- coeftrans_binom(modcoef$UCI)
        }
      } else if(outtype[i] == "count" && !is.null(coeftrans_count)){
        if(!is.function(coeftrans_count)){
          warning(paste0("'coeftrans_count' must be a function, but the ", 
                         "provided value is not a function, and therefore ",
                         "'coeftrans_count' will be ignored."))
          coeftrans_count <- NULL
        } else {
          modcoef$TransCoef <- coeftrans_count(modcoef$Coefficient)
          modcoef$TransLCI <- coeftrans_count(modcoef$LCI)
          modcoef$TransUCI <- coeftrans_count(modcoef$UCI)
        }
      }
      
      #Add diagnostic values
      if(mice){
        modcoef$AIC <- NA
      } else if(any(class(lm1) == "logistf")){
        extractAIC.logistf <- function(fit,  scale,  k = 2,  ...){
          dev <-   - 2  *  (fit$loglik['null'] - fit$loglik['full'])
          AIC <-  dev + k * fit$df
          return(AIC)
        }
        modcoef$AIC <- extractAIC.logistf(lm1)
      } else {
        modcoef$AIC <- AIC(lm1)
      }
      if(extradiag){
        if(!mice & !any(class(lm1) == "glm")){
          cooks <- cooks.distance(lm1)
          modcoef[, "N.Cooks.D_>_0.5"] <- length(which(cooks>0.5))
          modcoef[, "MaxCooks.D"] <- max(cooks, na.rm = T)
          modcoef[, "Heterosk.p"] <- 
            lmtest::bptest(lm1)$p.value
        } else if(!mice & any(class(lm1) == "glm")) {
          cooks <- cooks.distance(lm1)
          altmod <- glm(form1, data = lm1$model, 
                        family = lm1$family)
          simout <- DHARMa::simulateResiduals(altmod)
          quanttest <- DHARMa::testQuantiles(simout, plot = F)
          outliertest <- DHARMa::testOutliers(simout, plot = F, 
                                              type = 'bootstrap')
          uniftest <- DHARMa::testUniformity(simout, plot = F)
          disptest <- DHARMa::testDispersion(simout, plot = F)
          modcoef[, "N.Cooks.D_>_0.5"] <- length(which(cooks>0.5))
          modcoef[, "MaxCooks.D"] <- max(cooks, na.rm = T)
          modcoef[, "Heterosk.p"] <- quanttest$p.value
          modcoef[, "MinQuantileSplinePval"] <- min(quanttest$pvals)
          modcoef[, "UnifTestPval"] <- uniftest$p.value
          modcoef[, "OutlierTestPval"] <- outliertest$p.value
          modcoef[, "DispTestPval"] <- disptest$p.value
          if(outtype[i] == "count"){
            zeroinfltest <- DHARMa::testZeroInflation(simout, plot = F)
            modcoef[, "ZeroInflTestPval"] <- zeroinfltest$p.value
          }
        }
      }
      
      #Optional KFCV for lm, glm, or glm.nb
      if(KFCV){
        if(class(lm1) == "logistf" | mice){
          warning(paste0("This function is not set up to perform k-fold ",
                         "cross validation for list of linear models from",
                         " MICE-generated datasets or for 'logistf' ",
                         "model objects. Therefore, 'KFCV' is being ignored."))
        } else if(any(class(lm1) == "negbin")){
          t1 <- train(form1, ccDat, trControl  =  trainControl(
            method = "repeatedcv", number = KF, repeats = 50), method = "glm.nb")
          modcoef$KFCV <- t1$results$Accuracy
        } else if(any(class(lm1) == "glm")){
          if(outtype[i] == "count"){
            t1 <- train(form1, ccDat, 
                        trControl  =  trainControl(method = "repeatedcv", 
                                                   number = KF, repeats = 50), 
                        method = "glm", family = countfam)
          } else {
            t1 <- train(formula(paste0("factor(", outexpr, ")~", 
                                       as.character(form1)[3])), ccDat, 
                        trControl  =  trainControl(method = "repeatedcv", 
                                                   number = KF, repeats = 50), 
                        method = "glm", family = binomfam)
          }
          modcoef$KFCV <- t1$results$Accuracy
        } else if(any(class(lm1) == "lm")){
          t1 <- train(form1, ccDat, trControl  =  trainControl(
            method = "repeatedcv", number = KF, repeats = 50), method = "lm")
          modcoef$KFCV <- t1$results$RMSE
        } else {
          warning(paste0("An unrecognized model type is present. KFCV is not ",
                         "being performed."))
        }
      }
      #End results collation
      
      #Add temp data frame to main results
      if(i == 1 & j == 1) Resultsmat <- modcoef else {
        if(length(setdiff(names(modcoef), names(Resultsmat))) > 0){
          Resultsmat[, setdiff(names(modcoef), names(Resultsmat))] <- NA
        } 
        if(length(setdiff(names(Resultsmat), names(modcoef))) > 0){
          modcoef[, setdiff(names(Resultsmat), names(modcoef))] <- NA
        }
        Resultsmat <- rbind(Resultsmat, modcoef[, match(names(Resultsmat), 
                                                        names(modcoef))])
      }
    } #end j loop
  } #end i loop
  
  
  Resultsmat$Significant <- ifelse(Resultsmat$'p-value' < 0.05, "Yes", "")
  Resultsmat[, 1] <- factor(Resultsmat[, 1], levels = outnames)
  if(!is.null(altoutnames)) levels(Resultsmat[, 1]) <- altoutnames
  Resultsmat[, 2] <- factor(Resultsmat[, 2], levels = prednames)
  if(!is.null(altprednames)){
    levels(Resultsmat[, 2]) <- altprednames
    for(ll in 1:length(prednames)){
      Resultsmat$Variable <- gsub(prednames[ll], altprednames[ll], 
                                  Resultsmat$Variable)
    }; rm(ll)
  }
  if(!is.null(altixname)) Resultsmat$Variable <- 
    gsub(ixterm, altixname, Resultsmat$Variable)
  
  Resmatout <- Resultsmat
  names(Resmatout)[1:2] <- c(Outtitle, Predtitle)
  
  #Plotting
  if(is.null(coeftrans_cont)) coeftrans_cont <- function(x) x
  if(is.null(coeftrans_binom)) coeftrans_binom <- function(x) x
  if(is.null(coeftrans_count)) coeftrans_count <- function(x) x
  if(returnplot){
    if(length(unique(outtype)) == 1){
      if(all(outtype == "continuous")){
        horint <- coeftrans_cont(0)
      } else if(all(outtype == "binary")){
        horint <- coeftrans_binom(0)
      } else if(all(outtype == "count")){
        horint <- coeftrans_count(0)
      } else horint <- 0
      gg1 <- GLMResults_plot(Resultsmat, horint, facetcol, yl, colorbypred, 
                             colorpal, log10yaxis, Predtitle)
    } else {
      retplotlist <- list()
      for(mm in unique(outtype)){
        tmpplotdat <- Resultsmat[which(Resultsmat$Outcome %in% outnames[
          which(outtype == mm)]), ]
        if("TransCoef" %in% names(tmpplotdat) && all(
          is.na(tmpplotdat$TransCoef))){
          tmpplotdat$TransCoef <- tmpplotdat$Coefficient
          tmpplotdat$TransLCI <- tmpplotdat$LCI
          tmpplotdat$TransUCI <- tmpplotdat$UCI
        }
        if(mm == "continuous"){
          horint <- coeftrans_cont(0)
        } else if(mm == "binary"){
          horint <- coeftrans_binom(0)
        } else if(mm == "count"){
          horint <- coeftrans_count(0)
        } else horint <- 0
        retplotlist[[mm]] <- GLMResults_plot(
          tmpplotdat, horint, facetcol, yl, colorbypred, colorpal, log10yaxis, 
          Predtitle)
      }; rm(mm)
      listnamepos <- match(c("continuous", "binary", "count"), 
                           names(retplotlist))
      listnamepos <- listnamepos[which(!is.na(listnamepos))]
      retplotlist <- retplotlist[listnamepos]
      names(retplotlist) <- capfirst(names(retplotlist))
      gg1 <- retplotlist
    }
    
    retlist <- list(Resmatout, lmlist, gg1)
    names(retlist) <- c("Table", "ModelList", "GGplot")
  } else {
    retlist <- list(Resmatout, lmlist)
    names(retlist) <- c("Table", "ModelList")
  }

  return(retlist)
}
