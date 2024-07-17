#' Run multiple linear regressions and return results
#' 
#' \code{GLMResults} is a complicated helper function I've piecemeal added to 
#' throughout the years to handle all sorts of repeated linear regression analyses.
#' It outputs a table of key regression output including optional diagnostic test 
#' values, a list of the model objects, and a coefficient forest plot. Supported 
#' regression types include linear regression (\code{\link[stats]{lm}}) and 
#' binomial GLMs (\code{\link[stats]{glm}} with \code{\link[stats]{binomial}} 
#' family and link \code{"logit"} for logistic regressions or Firth's logistic 
#' regression using \code{\link[logistf]{logistf}}, \code{"log"} for binomial 
#' regressions, or \code{"probit"} for probit regressions, in addition to 
#' robust Poisson regressions for estimating risk ratios for binary 
#' outcomes where the family is \code{\link[stats]{poisson}}, the link function
#' is \code{"log"}, and robust errors are used).
#' 
#' @param prednames A character vector of the column names of predictors of 
#' interest you'd like to evaluate in separate regressions. These column names 
#' can refer to character, factor, or numeric columns in the data frame.
#' @param outnames A character vector of the column names of outcome variables 
#' you'd like to evaluate in separate regressions. These column names can refer 
#' to character, factor, or numeric columns in the data frame. For each 
#' \code{outnames}, \code{\link[stats]{lm}\code{()}} linear regressions will be 
#' run if the referenced column is numeric, or a \code{\link[stats]{glm}}\code{()} 
#' GLM with a \code{\link[stats]{binomial}}\code{()} family will be run with the 
#' link function specified in the input \code{binomlink}. This does not yet support 
#' characters or factors with >2 unique values (i.e., it can only accommodate 
#' continuous or binary outcome variables).
#' @param covnames Either a character vector of column names for all the 
#' covariates you'd like to include in each regression or a list of character 
#' vectors equal in length to the length of the unique combination of prednames 
#' and \code{outnames}. Defaults to \code{NULL}, meaning no covariates are 
#' included. This can be a list in order to control for separate covariates for 
#' each unique combination of \code{prednames} and \code{outnames}, in which 
#' case the list must be of length \code{length(prednames) * length(outnames)}. 
#' Note that in ordering this list the function loops through the \code{prednames} 
#' and then the \code{outnames} (e.g., Outcome 1 - Predictor 1, Outcome 1 - Predictor 2, 
#' Outcome 1 - Predictor 3, Outcome 2 - Predictor 1, etc.), and so you should 
#' order the list of covariates accordingly. These can include spline terms, such as
#' \code{\link[splines]{ns}} natural splines or \code{\link[splines]{bs}} 
#' b-splines as defined by the 'splines' package (e.g., \code{"ns(Year, df=3)"}).
#' These can also include interaction terms denoted by the ":" separator , though 
#' you should make sure to include the main effects too (e.g., \code{c("Income", 
#' "HouseholdSize", "Income:HouseholdSize")}).
#' @param Data A data frame object of class "data.frame".
#' @param logout If \code{TRUE}, this will log-transform each outcome variable 
#' (the log base will be \code{logbaseout}) prior to running the regression and then 
#' return both the raw coefficients as well as converted coefficients to percent 
#' changes using the formula percent change = (exp(raw coefficient)-1) x 100. 
#' Can be either a single value or a logical vector of length equal to the 
#' length of \code{outnames}. Defaults to \code{FALSE}.
#' @param logpred If \code{TRUE}, this will log-transform each predictor variable 
#' (the log base will be \code{logbasepred}) prior to running the regression and 
#' then return the raw coefficients. Can be either a single value or a logical 
#' vector of length equal to the length of \code{prednames}. Defaults to \code{FALSE}.
#' @param logbasepred The base for the log-transformation of one or more of the 
#' predictor variables. Defaults to 10.
#' @param logbaseout The base for the log transformation of the outcome variables. 
#' Defaults to \code{exp(1)} (i.e., natural log transformation).
#' @param Outtitle Defines the title used for the column of outcome variables 
#' in the table and on the plot. Defaults to \code{"Outcome"}.
#' @param Predtitle Defines the title used for the column of predictor variables 
#' in the table and on the plot. Defaults to \code{"Exposure"}.
#' @param MyMult The table will output raw and transformed coefficients and CIs. 
#' This multiplier is an optional scalar to multiply the raw coefficients for an 
#' estimate of the coefficient for something other than a 1-unit increase. 
#' Defaults to 1, and can be either a number or the character \code{"IQR"}, in 
#' which case the raw coefficient will be multiplied by the interquartile range 
#' of the predictor variable.
#' @param ixterm A column name for an optional interaction term for the 
#' \code{prednames} predictors of interest (I abbreviate interaction as 'ix'). 
#' Can refer to either a character, factor, or numeric-class column in the 
#' \code{Data} data frame. All relevant interaction and main coefficients will 
#' be output. Defaults to \code{NULL}.
#' @param Firth If \code{TRUE} and if a given \code{outnames} variable is binary, 
#' this will perform a Firth's logistic regression (\code{\link[logistf]{logistf}}) 
#' instead of a logistic regression (\code{\link[stats]{glm}}).Defaults to 
#' \code{FALSE}.
#' @param altprednames This is a vector of character values for alternate names
#' for the predictor variable column names specified in \code{prednames} to 
#' replace those column names in the output tables and plots. Defaults to \code{NULL}.
#' @param altoutnames This is a vector of character values for alternate names
#' for the outcome variable column names specified in \code{outnames} to 
#' replace those column names in the output tables and plots. Defaults to \code{NULL}.
#' @param altixname This is a character value for an alternate name for the 
#' interaction variable column name specified in \code{ixterm} to replace that 
#' column name in the output tables and plots. Defaults to \code{NULL}.
#' @param binomlink The link function to use for the family specified in \code{binfunc} 
#' for GLMs of binary \code{outnames}. Can be \code{"logit"}, \code{"log"}, or 
#' \code{"probit"}. Defaults to \code{"logit"}.
#' @param binfunc The family to use when \code{outnames} is binary. Can be 
#' \code{\link[stats]{binomial}}, \code{\link[stats]{poisson}}, or 
#' \code{\link[stats]{quasibinomial}}. Defaults to \code{\link[stats]{binomial}}. 
#' If the Poisson family is used, one should make sure to set \code{binomlink} to
#' \code{"log"} and \code{robust} to \code{TRUE}.
#' @param robust If \code{TRUE}, this calculates and returns robust CIs and 
#' p-values based on the \code{\link[sandwich]{vcovHC}} function instead of the 
#' default CIs and p-values. Defaults to \code{TRUE}.
#' @param robustdf If \code{FALSE}, this sets \code{\link[lmtest]{coefci}} to 
#' calculate z-score-based robust confidence intervals if \code{robust} is 
#' \code{TRUE}. If \code{TRUE}, t test-based robust confidence intervals are instead 
#' calculated based on the residual degrees of freedom in each regression if 
#' \code{robust} is \code{TRUE}. Defaults to \code{FALSE}, which is the default 
#' for the \code{df} parameter in \code{\link[lmtest]{coefci}}.
#' @param HCtype This is the formula used for the \code{\link[sandwich]{vcovHC}} 
#' function calculation of robust CIs and p-values. Defaults to \code{"HC0"}.
#' @param predspline: If \code{TRUE}, this transforms each predictor of interest 
#' into a natural spline term using \code{\link[splines]{ns}} from the '
#' splines' package. A coefficient for each degree of freedom is returned. 
#' Defaults to \code{FALSE}.
#' @param predsplinedf This specifies the number of degrees of freedom for the 
#' predictor variable natural splines if \code{predspline} is \code{TRUE}. 
#' Defaults to 3.
#' @param plotExpBeta If \code{TRUE}, this will plot exp(coefficient) values) 
#' in the coefficient forest plot. Defaults to \code{FALSE}.
#' @param plotPercChange If \code{TRUE}, this will plot percent changes (i.e., 
#' it will plot (exp(coefficient)-1) x 100 values) in the coefficient forest plot.
#' Defaults to \code{FALSE}. If logistic regressions were run and both 
#' \code{plotExpBeta} and \code{plotPercChange} are \code{FALSE}, raw coefficients 
#' will be plotted.
#' @param LOOCV If \code{TRUE}, this performs leave one out cross validation 
#' using the 'caret' package. If a given model is a linear regression 
#' (\code{\link[stats]{lm}}), the RMSE will be reported in the tables. If a model 
#' is a GLM (\code{\link[stats]{glm}}), the accuracy will be reported in the 
#' tables. The train function input method is \code{trainControl(method = 
#' "repeatedcv", number = 10, repeats = 50)}. Defaults to \code{FALSE}. 
#' @param facetcol This defines the number of facet_wrap columns for the 
#' coefficient plots (facetted by outcome variable). If \code{NULL}, this will 
#' equal the number of unique \code{outnames} values. Defaults to \code{NULL}.
#' @param covix This character vector of covariate column names allows for 
#' covariates to also have interaction terms with the \code{ixterm} interaction 
#' variable. If one wishes to have all covariates have interaction terms with 
#' \code{ixterm}, the character vector for \code{covix} should be the same as 
#' that for \code{covnames}. Defaults to \code{NULL}, meaning that there are no 
#' interactions with any of the covariates.
#' @param ixpred If \code{TRUE}, an interaction term is included between each of 
#' the \code{prednames} predictors of interest and \code{ixterm}. This can be 
#' set to \code{FALSE} to allow for interactions only with the covariates defined 
#' in \code{covix}. If \code{covix} is not \code{NULL} and \code{ixpred} is 
#' \code{TRUE}, there will be interaction terms with both the predictor of 
#' interest and the covariates. Defaults to \code{TRUE}.
#' @param extradiag If \code{TRUE}, this will include extra diagnostic 
#' information in the tables. For linear regressions (\code{\link[stats]{lm}}), 
#' this will include a p-value for the heteroskedasticity of model residuals 
#' (\code{\link[lmtest]{bptest}}), the number of Cook's distances > 0.5, and the 
#' maximum Cook's distance, the latter two of which define high leverage points.
#' For GLMs (\code{\link[stats]{glm}}), this will include output from DHARMa 
#' package tests, including a heteroskedasticity of residuals p-value from 
#' \code{\link[DHARMa]{testQuantiles}}, the minimum p-value from 
#' \code{\link[DHARMa]{testQuantiles}}, the number of Cook's distances > 0.5, 
#' the maximum Cook's distance, a p-value for a test of uniformity 
#' (\code{\link[DHARMa]{testUniformity}}), a p-value for outliers 
#' (\code{\link[DHARMa]{testOutliers}}), and a p-value for dispersion 
#' (\code{\link[DHARMa]{testDispersion}}). Defaults to \code{TRUE}.
#' @param leverage.test If \code{TRUE}, this runs all regressions while 
#' excluding "high leverage" observations as defined by the user. Defaults to 
#' \code{FALSE}.
#' @param leverage.cutoff This is the Cook's distance cutoff above which 
#' observations will be excluded if \code{leverage.test} is \code{TRUE}. 
#' Defaults to 0.2.
#' @param leverage.meancutoff If \code{leverage.cutoff} is set to \code{NULL}, 
#' this is used as a multiplier by the mean Cook's distance to get a new cutoff 
#' above which observations will be excluded if \code{leverage.test} is \code{TRUE}. 
#' For example, if set to 2, the cutoff will be set to the mean Cook's distance 
#' times 2. If both \code{leverage.cutoff} and \code{leverage.meancutoff} are set 
#' to \code{NULL}, the cutoff will be set to the mean Cook's distance times 4. 
#' Defaults to \code{NULL}.
#' @param post.power If \code{TRUE}, this repeatedly simulates data with similar 
#' structure to the actual data using \code{\link[SimMultiCorrData]{rcorrvar}} 
#' and then uses a user-specified predictor of interest coefficient to calculate 
#' a simulated outcome variable (all covariate coefficients are drawn from the 
#' original model estimates), runs regressions on these simulated data frames, 
#' collects the iterated predictor of interest coefficient p-values, and 
#' calculates how many of those p-values are <0.05. This is a form of posterior 
#' power test for different effect sizes of the predictor of interest association 
#' with the outcome variable. Defaults to \code{FALSE}.
#' @param effect.size This is the user-specified effect size for which the 
#' post-hoc power is tested if \code{post.power} is set to \code{TRUE}. Defaults 
#' to 0.5.
#' @param nsim This is the number of simulations for the post-hoc power test if 
#' \code{post.power} is set to \code{TRUE}. Defaults to 1E3.
#' @param mice If TRUE, this performs multiple imputation (MICE using 
#' \code{\link[mice]{mice}}) on the user-specified variables and pools model 
#' estimates across the imputed datasets. Defaults to \code{FALSE}.
#' @param micevars A character vector of the names of variables to be imputed 
#' using MICE. Defaults to \code{NULL}, in which case it will include all 
#' predictor, outcome, and covariate variables that have missing values.
#' @param miceiter The number of multiple imputations to perform if \code{mice} 
#' is set to \code{TRUE}. Defaults to 10.
#' 
#' @details
#' \code{GLMResults} was built over several years and so has incorporated all 
#' sorts of functionality related to linear regressions that I've found useful. 
#' This is my main workhorse function for most analyses calling for regressions.
#' 
#' @return \code{GLMResults} returns a list of:
#' \item{Matrix}{A table of all the key variables. Beta, LCI, and UCI refer to 
#' the coefficient for the predictor of interest and the lower and upper 95% 
#' confidence intervals, respectively. Beta.Transform, etc. are those estimates 
#' transformed in some way as specified by the user, for instance by being 
#' multiplied by the IQR, exponentiated or transformed into percent changes (done
#' by default when the outcome is binary), or some other multiplication specified 
#' by the user.}
#' \item{GGplot}{A forest ggplot of all predictor of interest coefficients and CIs.} 
#' \item{LMlist}{A list of all the lm, glm, or logistf objects run for each unique 
#' combination of predictor and outcome variable names (i.e., \code{prednames} 
#' and \code{outnames}).}
#' 
#' @import caret lmtest sandwich logistf splines ggplot2 emmeans SimMultiCorrData DHARMa pbapply cowplot
#' @export GLMResults
#' 
#' @examples
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
#' 
#' # Probit regression with an interaction term on the predictor of interest and
#' #  on the covariate "gear" while also using robust errors (the default)
#' 
#' glmres <- GLMResults(carpreds, c("vs", "am"), c("gear", "carb"), mtcars, 
#'                      ixterm = "qsec", covix = "gear", binomlink = "probit", 
#'                      altprednames = c("Displacement", "Cylinders", "Axle Ratio", 
#'                      "Weight"), altoutnames = c("Engine Type", "Transmission Type"), 
#'                      altixname = "Quarter Mile Time")
#' 
#' # Robust Poisson regression
#' glmres <- GLMResults(carpreds[-1], c("vs", "am"), c("gear", "carb"), mtcars, 
#'                      binomlink = "log", binfunc = poisson,  
#'                      altprednames = c("Displacement", "Cylinders", "Axle Ratio", 
#'                      "Weight"), altoutnames = c("Engine Type", "Transmission Type"))
#' 
#' # MICE example
#' # note that robust errors are not attainable when using MICE, so I use robust = F
#' 
#' mtcars_miss <- mtcars
#' set.seed(101)
#' mtcars_miss[sample(1:nrow(mtcars_miss), 3), "cyl"] <- NA
#' set.seed(104)
#' mtcars_miss[sample(1:nrow(mtcars_miss), 4), "wt"] <- NA        
#' 
#' lmres <- GLMResults(prednames = carpreds, outnames = carouts, 
#'                    covnames = carcovars, Data = mtcars, mice = T, 
#'                    micevars = c("cyl", "wt"), robust = F)
#'             


GLMResults <- function(prednames, outnames, covnames = NULL, Data, logout = F, 
                       logpred = F, logbasepred = 10, logbaseout = exp(1), 
                       Outtitle = "Outcome", Predtitle = "Exposure", MyMult = 1, 
                       ixterm = NULL, Firth = F, altprednames = NULL, 
                       altoutnames = NULL, altixname = NULL, binomlink = "logit", 
                       binfunc = binomial, robust = T, robustdf = F, 
                       HCtype = "HC0", predspline = F, predsplinedf = 3, 
                       plotExpBeta = F, plotPercChange = F, LOOCV = F, 
                       facetcol = NULL, covix = NULL, ixpred = T, extradiag = T, 
                       leverage.test = F, leverage.cutoff = 0.2, 
                       leverage.meancutoff = NULL, post.power = F, 
                       effect.size = 0.5, nsim = 1E3, mice = F, micevars = NULL, 
                       miceiter = 10){
  robustse <- function(x, HCtype="HC0",usedf=robustdf) {
    mod1 <- coeftest(x, vcov = function(x) vcovHC(x,type=HCtype))
    if(!usedf){
      cis<-coefci(x, vcov=function(x) vcovHC(x,type=HCtype))
    } else {
      cis<-coefci(x, vcov = function(x) vcovHC(x,type=HCtype),
                  df=x$df.residual)
    }
    mod1<-cbind(mod1,cis)
    return(mod1)
  }
  confint.qt<-function(beta,se,DF,IQR=1,level=0.95){
    ci.lower<-(beta*IQR)-((se*IQR)*qt(((level/2)+0.5), DF))
    ci.upper<-(beta*IQR)+((se*IQR)*qt(((level/2)+0.5), DF))
    CIs<-data.frame(CIL=ci.lower,CIU=ci.upper)
    return(CIs)
  }
  if(!is.null(ixterm)&predspline){
    stop(paste0("Code not yet set up to accommodate spline interactions. ",
                "Please set 'ixterm' to NULL if 'predspline' is TRUE."))
  }

  if(robust&Firth){
    stop(paste0("Code not yet set up to extract robust sandwich SEs from logistf",
                " objects. Please set 'robust' to FALSE if 'Firth' is TRUE."))
  }
  if(Firth&extradiag){
    stop(paste0("Code not yet set up to perform GLM diagnostics on logistf",
                " objects. Please set 'extradiag' to FALSE if 'Firth' is TRUE."))
  }

  predclasses<-sapply(Data[,prednames],class)
  if(any(predclasses=="character")){
    for(cl in which(predclasses=="character")){
      Data[,prednames[cl]]<-as.factor(Data[,prednames[cl]])
    }; rm(cl)
  }

  if(any(predclasses%in%c("factor"))){
    predlevellengths<-sapply(Data[,prednames],function(x) length(levels(x)))
    predlevelrows<-ifelse(predlevellengths%in%c(0,1),1,predlevellengths-1)
    predrowlength<-sum(predlevelrows)
  } else {
    predlevelrows<-rep(1,length(prednames))
    predrowlength<-length(prednames)
  }

  if(is.null(ixterm)&predspline!=F){
    Resultsmat<-as.data.frame(matrix(NA,predsplinedf*predrowlength*length(outnames),12))
  } else if(!is.null(ixterm)){
    nixlvl<-length(levels(Data[,ixterm]))
    if(class(Data[,ixterm])%in%c("factor","character")){
      if(class(Data[,ixterm])=="character"){Data[,ixterm]<-as.factor(Data[,ixterm])}
      if(is.null(ncol(combn(levels(Data[,ixterm]),2)))) {
        ixtermlength<-nixlvl+(nixlvl-1)+1
      } else {
        ixtermlength<-nixlvl+(nixlvl-1)+
          ncol(combn(levels(Data[,ixterm]),2))
      }
    } else {
      ixtermlength<-3
    }

    Resultsmat<-as.data.frame(matrix(NA,ixtermlength*predrowlength*length(outnames),
                                     12))
  } else {
    Resultsmat<-as.data.frame(matrix(NA,predrowlength*length(outnames),11))
  }

  lmlist<-list()

  if(extradiag){Extramat<-data.frame(Heterosk.p=rep(NA,nrow(Resultsmat)))}

  if(length(logout)==1){
    logout2<-rep(logout,length(outnames))
  } else {
    logout2<-logout
  }

  if(length(logpred)==1){
    logpred2<-rep(logpred,length(prednames))
  } else {
    logpred2<-logpred
  }

  for(i in 1:length(outnames)){
    isbinary <- length(unique(Data[which(!is.na(Data[,outnames[i]])),outnames[i]]))==2
    
    for(j in 1:length(prednames)){
      if(is.null(ixterm)&predspline==T){
        rowstart<-(1+(j-1))+((i-1)*predsplinedf*predrowlength)+((predsplinedf-1)*(j-1))
        rowend<-(predsplinedf+(j-1))+((i-1)*predsplinedf*predrowlength)+ 
          ((predsplinedf-1)*(j-1))
        rownum<-rowstart:rowend
        lmnum<-(1+(j-1))+((i-1)*length(prednames))
      } else if(!is.null(ixterm)){
        rowstart<-(1+(j-1))+((i-1)*ixtermlength*predrowlength)+((ixtermlength-1)*(j-1))
        rowend<-(ixtermlength+(j-1))+((i-1)*ixtermlength*predrowlength)+ 
          ((ixtermlength-1)*(j-1))
        rownum<-rowstart:rowend
        lmnum<-(1+(j-1))+((i-1)*predrowlength) 
      } else if(any(sapply(Data[,prednames],FUN=function(x) class(x)=="factor"))) {
        rowstarts<-cumsum(c(1,predlevelrows[-length(predlevelrows)]))+((i-1)*predrowlength)
        rowends<-cumsum(predlevelrows)+((i-1)*predrowlength)
        rownum<-rowstarts[j]:rowends[j]
        lmnum<-(1+(j-1))+((i-1)*length(prednames))
      } else {
        rownum<-(1+(j-1))+((i-1)*length(prednames))
        lmnum<-rownum
      }

      if(logout2[i]==T){
        outexpr<-paste0("log(",outnames[i],",base=",logbaseout,")")
      } else {
        outexpr<-outnames[i]
      }

      if(logpred2[j]==T&predspline==T){
        predexpr<-paste0("ns(log(",prednames[j],",base=",logbasepred,"),",predsplinedf,")")
      } else if(logpred2[j]==T){
        predexpr<-paste0("log(",prednames[j],",base=",logbasepred,")")
      } else if(predspline==T){
        predexpr<-paste0("ns(",prednames[j],",",predsplinedf,")")
      } else {
        predexpr<-prednames[j]
      }

      if(is.null(covnames)){
        if(is.null(ixterm)){
          if(mice){
            ccvars<-c(outnames[i],prednames[j])
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j])]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr))
        } else {
          if(mice){
            ccvars<-c(outnames[i],prednames[j],ixterm)
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],ixterm)]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm))
        }
      } else {
        if(is.list(covnames)){
          if(length(covnames)!=(length(prednames)*length(outnames))){
            stop(paste0("Error: covnames list must have a length of ",
                        length(prednames)*length(outnames)))
          }
          mycovnames<-covnames[[lmnum]]
        } else {
          mycovnames<-covnames
        }
        if(any(grepl("(",mycovnames,fixed=T))&any(grepl(":",mycovnames,fixed=T))){
          mycovnames_df<-gsub(".*\\(","",mycovnames)
          mycovnames_df<-gsub(",.*","",mycovnames_df)
          mycovnames_df<-mycovnames_df[-which(grepl(":",mycovnames_df,fixed=T))]
        } else if(any(grepl(":",mycovnames,fixed=T))){
          mycovnames_df<-mycovnames[-which(grepl(":",mycovnames,fixed=T))]
        } else if(any(grepl("(",mycovnames,fixed=T))){
          mycovnames_df<-gsub(".*\\(","",mycovnames)
          mycovnames_df<-gsub(",.*","",mycovnames_df)
        } else {
          mycovnames_df<-mycovnames
        }
        if(is.null(ixterm)&is.null(covix)){
          if(mice){
            ccvars<-c(outnames[i],prednames[j],mycovnames_df)
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovnames_df)]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr,"+",paste(mycovnames,collapse="+")))
        } else if(is.null(covix)){
          if(mice){
            ccvars<-unique(c(outnames[i],prednames[j],mycovnames_df,ixterm))
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovnames_df,ixterm)]),]
          }
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm,"+",
                                paste(mycovnames,collapse="+")))
        } else {
          if(mice){
            ccvars<-unique(c(outnames[i],prednames[j],mycovnames_df,ixterm))
            if(is.null(micevars)){micevars<-ccvars[which(sapply(Data[,ccvars],function(x) any(is.na(x))))]}
            ccDat<-Data[complete.cases(Data[,ccvars[-which(ccvars%in%micevars)]]),]
            micedata<-mice(ccDat[,ccvars],m=miceiter,maxit=miceiter,printFlag=F)
          } else {
            ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovnames_df,ixterm)]),]
          }
          nobs<-dim(ccDat)[1]
          covixaux<-rep("",length(mycovnames))
          covixaux[which(mycovnames%in%covix)]<-paste0("*",ixterm)
          for(cv in 1:length(mycovnames)){
            mycovnames[cv]<-paste0(mycovnames[cv],covixaux[cv])
          }
          if(ixpred==T){
            form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm,"+",
                                  paste(mycovnames,collapse="+")))
          } else {
            form1<-formula(paste0(outexpr,"~",predexpr,"+",
                                  paste(mycovnames,collapse="+")))
          }
        }
      }
      
      if(isbinary){
        if(binfunc(binomlink)$family=="poisson"){
          if(binomlink!="log"){
            warning(paste0("The only acceptable link function for the Poisson ",
                           "family is 'log', so 'binomlink' is being changed to",
                           " 'log'."))
            binomlink <- "log"
          }
          if(is.character(Data[,outnames[i]])) Data[,outnames[i]]<-as.factor(Data[,outnames[i]])
          if(is.factor(Data[,outnames[i]])){
            levels(Data[,outnames[i]])<-c("0","1")
            Data[,outnames[i]]<-as.numeric(as.character(Data[,outnames[i]]))
          }
        }
        if(Firth==T){
          lm1<-logistf(form1,data=Data,na.action=na.exclude)
          if(binfunc(binomlink)$family!="binomial"|binomlink!="logit"){
            warning(paste0("Note that only binfunc = binomial and binomlink = ",
                           "'logit' are accepted for Firth's logistic regressions,",
                           " so those specified parameters are being ignored",
                           " while Firth is TRUE."))
          }
          if(leverage.test==T){
            warning(paste0("It's unclear how to get Cook's distance for a 'logistf'",
                           " object. Ignoring leverage.test while Firth is TRUE."))
          }
        } else {
          if(mice){
            lm1<-list()
            tempdatalist<-list()
            for(mim in 1:micedata$m){
              tempdatalist[[mim]]<-ccDat
              for(mip in 1:length(micevars)){
                tempdatalist[[mim]][match(rownames(micedata$imp[[micevars[mip]]]),rownames(ccDat)),
                                    micevars[mip]]<-
                  micedata$imp[[micevars[mip]]][,mim]
              }; rm(mip)
              lm1[[mim]]<-glm(form1,data=tempdatalist[[mim]],family=binfunc(link=binomlink),na.action=na.exclude)
            }; rm(mim)
          } else {
            lm1<-glm(form1,data=Data,family=binfunc(link=binomlink),na.action=na.exclude)
          }
          if(leverage.test==T){
            if(mice){
              for(mim in 1:length(lm1)){
                cooks<-cooks.distance(lm1[[mim]])
                if(!is.null(leverage.cutoff)){
                  threshold<-leverage.cutoff
                } else if(!is.null(leverage.meancutoff)){
                  threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
                } else {threshold<-mean(cooks,na.rm=T)*4}
                lm1<-glm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                         family=binfunc(link=binomlink),na.action=na.exclude)
              }
            } else {
              cooks<-cooks.distance(lm1)
              if(!is.null(leverage.cutoff)){
                threshold<-leverage.cutoff
              } else if(!is.null(leverage.meancutoff)){
                threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
              } else {threshold<-mean(cooks,na.rm=T)*4}
              lm1<-glm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                       family=binfunc(link=binomlink),na.action=na.exclude)
            }
          }
          if(robust==T){
            if(mice){
              warning("Unclear how to get pooled robust SEs from imputed GLMs; ignoring robust command")
            } else {
              robustcoef<-robustse(lm1,HCtype=HCtype)
            }
          }
        }

      } else {
        if(mice){
          lm1<-list()
          tempdatalist<-list()
          for(mim in 1:micedata$m){
            tempdatalist[[mim]]<-ccDat
            for(mip in 1:length(micevars)){
              tempdatalist[[mim]][match(rownames(micedata$imp[[micevars[mip]]]),rownames(ccDat)),
                                  micevars[mip]]<-micedata$imp[[micevars[mip]]][,mim]
            }; rm(mip)
            lm1[[mim]]<-lm(form1,data=tempdatalist[[mim]],na.action=na.exclude)
          }; rm(mim)
        } else {
          lm1<-lm(form1,data=Data,na.action=na.exclude)
        }

        if(leverage.test==T){
          if(mice){
            for(mim in 1:length(lm1)){
              cooks<-cooks.distance(lm1[[mim]])
              if(!is.null(leverage.cutoff)){
                threshold<-leverage.cutoff
              } else if(!is.null(leverage.meancutoff)){
                threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
              } else {threshold<-mean(cooks,na.rm=T)*4}
              lm1<-lm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                       na.action=na.exclude)

            }
          } else {
            cooks<-cooks.distance(lm1)
            if(!is.null(leverage.cutoff)){
              threshold<-leverage.cutoff
            } else if(!is.null(leverage.meancutoff)){
              threshold<-mean(cooks,na.rm=T)*leverage.meancutoff
            } else {
              threshold<-mean(cooks,na.rm=T)*4
            }
            lm1<-lm(form1,data=Data[-c(as.numeric(names(cooks[which(cooks>threshold)]))),],
                    na.action=na.exclude)
          }
        }
        if(robust==T){
          if(mice){
            warning("Unclear how to do robust SEs with MICE, so ignoring 'robust' command")
          } else {
            robustcoef<-robustse(lm1,HCtype=HCtype)
          }
        }
      }
      lmlist[[lmnum]]<-lm1
      names(lmlist)[[lmnum]]<-paste0(outnames[i],"~",prednames[j])
      Resultsmat[rownum,1]<-outnames[i]
      Resultsmat[rownum,3]<-nobs
      if(mice){
        Resultsmat[rownum,2]<-as.character(summary(pool(lm1))$term[2:(1+length(rownum))])
      } else if(Firth){
        Resultsmat[rownum,2]<-lm1$terms[2:(1+length(rownum))]
      } else {
        Resultsmat[rownum,2]<-dimnames(summary(lm1)$coef)[[1]][2:(1+length(rownum))]
      }
      if(!is.null(ixterm)){
        Resultsmat[rownum,2]<-prednames[j]
        if(mice) {
          micesumm<-summary(pool(lm1),conf.int=T)
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            Resultsmat[rownum,12]<-c(
              paste(micesumm$term[2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(micesumm$term[2],"*",comb1[2,]," - ",comb1[1,])
            )
            miceixsumm<-as.data.frame(matrix(NA,ixtermlength,6))
            names(miceixsumm)<-c("beta","SE","DF","pval","CIL","CIU")
            miceixsumm[1,]<-micesumm[2,c(2:3,5:8)]
            miceixsumm[(1+nixlvl:((2*nixlvl)-2)),]<-
              micesumm[c(3:(2+(nixlvl-1))),c(2:3,5:8)]
            miceixsumm[(2*nixlvl):
                         ((2*nixlvl)+((nixlvl)-2)),]<-
              micesumm[grep(paste0(prednames[j],":"),micesumm$term),c(2:3,5:8)]
            newData<-ccDat
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              newlm1<-list()
              tempdatalist<-list()
              for(mim in 1:micedata$m){
                tempdatalist[[mim]]<-newData
                for(mip in 1:length(micevars)){
                  tempdatalist[[mim]][match(rownames(micedata$imp[[micevars[mip]]]),rownames(newData)),
                                      micevars[mip]]<-micedata$imp[[micevars[mip]]][,mim]
                }; rm(mip)
                if(any(class(lm1[[1]])=="glm")){
                  newlm1[[mim]]<-glm(form1,data=tempdatalist[[mim]],family=binfunc(link=binomlink),
                                     na.action=na.exclude)
                } else {
                  newlm1[[mim]]<-lm(form1,data=tempdatalist[[mim]],na.action=na.exclude)
                }
              }; rm(mim)

              newmicesumm<-summary(pool(newlm1),conf.int=T)
              miceixsumm[1+rp,]<-newmicesumm[2,c(2:3,5:8)]

              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )
                if(is.null(covix)){
                  newtakerows<-c((nrow(newmicesumm)-nixlvl+2):
                                   (nrow(newmicesumm)-(nixlvl-1)+nix2take))
                } else {
                  newimpixrows<-which(grepl(paste0(prednames[j],":"),newmicesumm$term))
                  newtakerows<-newimpixrows[1:nix2take]
                }

                miceixsumm[newfillrows,]<- newmicesumm[newtakerows,c(2:3,5:8)]
              }
            }; rm(rp)
            betalist<-miceixsumm$beta
            selist<-miceixsumm$SE
            dflist<-miceixsumm$DF
            pvallist<-miceixsumm$pval
            cilist<-miceixsumm[,c("CIL","CIU")]
          } else {
            betalist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),2])
            pvallist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),6])
            selist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),3])
            dflist<-c(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),5])
            cilist<-cbind(micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),7],
                          micesumm[c(2:3,grep(paste0(prednames[j],":"),micesumm$term)),8])
            Resultsmat[rownum,12]<-
              as.character(micesumm$term[c(2:3,grep(paste0(prednames[j],":"),
                                                          micesumm$term))])
          }
        } else if(any(class(lm1)=="logistf")){
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            if(is.null(ncol(comb1))) comb1<-as.matrix(comb1,nrow=length(comb1))
            Resultsmat[rownum,12]<-c(
              paste(names(lm1$coefficients)[2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(lm1$terms[2],"*",comb1[1,]," - ",comb1[2,])
            )
            betalist<-rep(NA,ixtermlength)
            betalist[1]<-lm1$coefficients[2]
            betalist[(1+nixlvl:((2*nixlvl)-2))]<-
              lm1$coefficients[c(3:(2+(nixlvl-1)))]
            betalist[(2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))]<-
              -1*lm1$coefficients[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )]

            selist<-rep(NA,ixtermlength)
            selist[1]<-sqrt(diag(lm1$var))[2]
            selist[(1+nixlvl:((2*nixlvl)-2))]<-
              sqrt(diag(lm1$var))[c(3:(2+(nixlvl-1)))]
            selist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              sqrt(diag(lm1$var))[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )]

            cilist<-matrix(NA,ixtermlength,2)
            cilist[1,]<-c(lm1$ci.lower[2],lm1$ci.upper[2])
            cilist[c(1+nixlvl:((2*nixlvl)-2)),]<-
              cbind(lm1$ci.lower[c(3:(2+(nixlvl-1)))],
                    lm1$ci.upper[c(3:(2+(nixlvl-1)))])

            cilist[c((2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))),]<-
              cbind(-1*lm1$ci.upper[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )],
              -1*lm1$ci.lower[c(
                (length(lm1$coefficients)-(nixlvl-2)):
                  length(lm1$coefficients)
              )])

            pvallist<-rep(NA,ixtermlength)
            pvallist[1]<-lm1$prob[2]
            pvallist[(1+nixlvl):((2*nixlvl)-1)]<-
              lm1$prob[c(3:(2+(nixlvl-1)))]
            pvallist[(2*nixlvl):
                       ((2*nixlvl)+(nixlvl-2))]<-
              lm1$prob[c(
                (length(lm1$prob)-(nixlvl-2)):
                  length(lm1$prob)
              )]

            newData<-Data
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              newlm1<-logistf(form1,data=newData,family=binomial(link=binomlink),na.action=na.exclude)
              betalist[1+rp]<-newlm1$coefficients[2]
              selist[1+rp]<-sqrt(diag(newlm1$var))[2]
              pvallist[1+rp]<-newlm1$prob[2]
              cilist[1+rp,]<-cbind(newlm1$ci.lower[2],
                                   newlm1$ci.upper[2])

              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newtakerows<-c((length(lm1$coefficients)-nixlvl+2):
                                 (length(lm1$coefficients)-(nixlvl-1)+nix2take))
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )

                betalist[newfillrows]<- -1*newlm1$coefficients[newtakerows]
                selist[newfillrows]<-sqrt(diag(newlm1$var))[newtakerows]
                pvallist[newfillrows]<-newlm1$prob[newtakerows]
                cilist[newfillrows,]<-cbind(-1*newlm1$ci.upper[newtakerows],
                                            -1*newlm1$ci.lower[newtakerows])
              }
            }
          } else {
            betalist<-c(lm1$coefficients[c(2:3,length(lm1$coefficients))])
            pvallist<-c(lm1$prob[c(2:3,length(lm1$prob))])
            selist<-c(sqrt(diag(lm1$var))[c(2:3,nrow(lm1$var))])
            cilist<-cbind(lm1$ci.lower[c(2:3,length(lm1$prob))],
                          lm1$ci.upper[c(2:3,length(lm1$prob))])
          }
        } else if (robust==T) {
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            Resultsmat[rownum,12]<-c(
              paste(dimnames(robustcoef)[[1]][2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(dimnames(robustcoef)[[1]][2],"*",comb1[2,]," - ",comb1[1,])
            )
            betalist<-rep(NA,ixtermlength)
            betalist[1]<-robustcoef[2,1]
            betalist[(1+nixlvl:((2*nixlvl)-2))]<-
              robustcoef[c(3:(2+(nixlvl-1))),1]
            betalist[(2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))]<-
              robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),1]
            selist<-rep(NA,ixtermlength)
            selist[1]<-robustcoef[2,2]
            selist[(1+nixlvl:((2*nixlvl)-2))]<-
              robustcoef[c(3:(2+(nixlvl-1))),2]
            selist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),2]
            cilist<-matrix(NA,ixtermlength,2)
            cilist[1,]<-c(robustcoef[2,5],robustcoef[2,6])
            cilist[1+nixlvl:((2*nixlvl)-2),]<-
              cbind(robustcoef[c(3:(2+(nixlvl-1))),5],
                    robustcoef[c(3:(2+(nixlvl-1))),6])

            cilist[c((2*nixlvl):((2*nixlvl)+((nixlvl)-2))),]<-
              cbind(robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),5],
                    robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),6])

            pvallist<-rep(NA,ixtermlength)
            pvallist[1]<-robustcoef[2,4]
            pvallist[(1+nixlvl:((2*nixlvl)-2))]<-
              robustcoef[c(3:(2+(nixlvl-1))),4]
            pvallist[(2*nixlvl):((2*nixlvl)+(nixlvl-2))]<-
              robustcoef[grep(paste0(prednames[j],":"),dimnames(robustcoef)[[1]]),4]
            

            newData<-Data
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              if(any(class(lm1)=="glm")){
                newlm1<-glm(form1,data=newData,
                            family=binfunc(link=binomlink),na.action=na.exclude)
              } else {
                newlm1<-lm(form1,data=newData,na.action=na.exclude)
              }

              newrobustcoef<-robustse(newlm1,HCtype=HCtype)
              betalist[1+rp]<-newrobustcoef[2,1]
              selist[1+rp]<-newrobustcoef[2,2]
              pvallist[1+rp]<-newrobustcoef[2,4]
              cilist[1+rp,]<-cbind(newrobustcoef[2,5],
                                   newrobustcoef[2,6])

              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )
                if(is.null(covix)){
                  newtakerows<-c((nrow(newrobustcoef)-nixlvl+2):
                                   (nrow(newrobustcoef)-(nixlvl-1)+nix2take))
                } else {
                  newimpixrows<-which(grepl(paste0(prednames[j],":"),rownames(newrobustcoef)))
                  newtakerows<-newimpixrows[1:nix2take]
                }

                betalist[newfillrows]<- newrobustcoef[newtakerows,1]
                selist[newfillrows]<-newrobustcoef[newtakerows,2]
                pvallist[newfillrows]<-newrobustcoef[newtakerows,4]
                cilist[newfillrows,]<-cbind(newrobustcoef[newtakerows,5],
                                            newrobustcoef[newtakerows,6])
              }
            }
          } else {
            betalist<-c(robustcoef[c(2:3,grep(paste0(prednames[j],":"),
                                              dimnames(robustcoef)[[1]])),1])
            pvallist<-c(robustcoef[c(2:3,grep(paste0(prednames[j],":"),
                                              dimnames(robustcoef)[[1]])),4])
            selist<-c(robustcoef[c(2:3,grep(paste0(prednames[j],":"),
                                            dimnames(robustcoef)[[1]])),2])
            cilist<-cbind(robustcoef[c(2:3,grep(paste0(prednames[j],":"),
                                                dimnames(robustcoef)[[1]])),5],
                          robustcoef[c(2:3,grep(paste0(prednames[j],":"),
                                                dimnames(robustcoef)[[1]])),6])
            Resultsmat[rownum,12]<-
              dimnames(summary(lm1)$coef)[[1]][c(2:3,grep(paste0(prednames[j],":"),
                                                          dimnames(robustcoef)[[1]]))]
          }
        } else { #if robust==F
          if(class(Data[,ixterm])=="factor"){
            comb1<-combn(levels(Data[,ixterm]),2)
            Resultsmat[rownum,12]<-c(
              paste(rownames(summary(lm1)$coef)[2],levels(Data[,ixterm]),sep="|"),
              paste0(comb1[2,1:(nixlvl-1)]," - ",comb1[1,1:(nixlvl-1)]),
              paste0(rownames(summary(lm1)$coef)[2],"*",comb1[2,]," - ",comb1[1,])
            )
            betalist<-rep(NA,ixtermlength)
            betalist[1]<-lm1$coefficients[2]
            betalist[(1+nixlvl:((2*nixlvl)-2))]<-
              lm1$coefficients[c(3:(2+(nixlvl-1)))]
            betalist[(2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))]<-
              lm1$coefficients[grep(paste0(prednames[j],":"),names(lm1$coef))]

            selist<-rep(NA,ixtermlength)
            selist[1]<-summary(lm1)$coef[2,2]
            selist[(1+nixlvl:((2*nixlvl)-2))]<-
              summary(lm1)$coef[c(3:(2+(nixlvl-1))),2]
            selist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              summary(lm1)$coef[grep(paste0(prednames[j],":"),names(lm1$coef)),2]

            cilist<-matrix(NA,ixtermlength,2)
            lm1.ci<-suppressMessages(confint(lm1,names(lm1$coef)[c(
              2,c(3:(2+(nixlvl-1))),grep(paste0(prednames[j],":"),names(lm1$coef)))]))
            cilist[1,]<-lm1.ci[1,]
            cilist[c(1+nixlvl:((2*nixlvl)-2)),]<-
              lm1.ci[c(2:(1+(nixlvl-1))),]
            cilist[c((2*nixlvl):
                       ((2*nixlvl)+((nixlvl)-2))),]<-
              lm1.ci[grep(paste0(prednames[j],":"),rownames(lm1.ci)),]

            pvallist<-rep(NA,ixtermlength)
            pvallist[1]<-summary(lm1)$coef[2,4]
            pvallist[(1+nixlvl:((2*nixlvl)-2))]<-
              summary(lm1)$coef[c(3:(2+(nixlvl-1))),4]
            pvallist[(2*nixlvl):
                     ((2*nixlvl)+((nixlvl)-2))]<-
              summary(lm1)$coef[grep(paste0(prednames[j],":"),names(lm1$coef)),4]

            newData<-Data
            for(rp in 1:(nixlvl-1)){
              newData[,ixterm]<-factor(newData[,ixterm],
                                       levels=levels(newData[,ixterm])[
                                         c(2:(length(levels(newData[,ixterm]))),1)])
              if(any(class(lm1)=="glm")){
                newlm1<-glm(form1,data=newData,family=binfunc(link=binomlink),na.action=na.exclude)
              } else {
                newlm1<-lm(form1,data=newData,na.action=na.exclude)
              }

              betalist[1+rp]<-newlm1$coefficients[2]
              selist[1+rp]<-summary(newlm1)$coef[2,2]
              pvallist[1+rp]<-summary(newlm1)$coef[2,4]
              cilist[1+rp,]<-suppressMessages(confint(newlm1,2:3))[1,]

              nix2take<-(nixlvl-1)-rp
              if(nix2take>0){
                newfillrows<-c(
                  (((3*nixlvl)-1)+((nix2take+1)*(rp-1))):
                    ((((3*nixlvl)-1)+((nix2take+1)*(rp-1)))+(nix2take-1))
                )

                betalist[newfillrows]<- newlm1$coefficients[grep(paste0(prednames[j],":"),names(newlm1$coef))]
                selist[newfillrows]<-summary(newlm1)$coef[grep(paste0(prednames[j],":"),names(newlm1$coef)),2]
                pvallist[newfillrows]<-summary(newlm1)$coef[grep(paste0(prednames[j],":"),names(newlm1$coef)),4]
                cilist[newfillrows,]<-suppressMessages(confint(newlm1),grep(paste0(prednames[j],":"),names(newlm1$coef)))
              }
            }
          } else {
            Resultsmat[rownum,12]<-
              dimnames(summary(lm1)$coef)[[1]][c(2:3,grep(paste0(prednames[j],":"),
                                                          dimnames(summary(lm1)$coef)[[1]]))]
          }
        }

      } else if(predspline==T){
        Resultsmat[rownum,12]<-dimnames(summary(lm1)$coef)[[1]][2:(predsplinedf+1)]
      }

      if(MyMult=="IQR"){
        if(logpred2[j]==TRUE){
          Mult<-as.numeric((IQR(Data[,prednames[j]],na.rm=T)/
                              summary(Data[,prednames[j]])[2])*100)
        } else {
          Mult<-IQR(Data[,prednames[j]],na.rm=T)
        }
      } else {
        Mult<-MyMult
      }

      if(!is.null(ixterm)){
        if(mice|any(class(lm1)=="logistf")|robust|class(Data[,ixterm])=="factor"){
          Resultsmat[rownum,4]<-betalist
          Resultsmat[rownum,5:6]<-cilist
          Resultsmat[rownum,7]<-pvallist
        } else {
          mysumrows<-c(2:3,nrow(summary(lm1)$coef))
          betalist<-summary(lm1)$coef[mysumrows,1]
          selist<-summary(lm1)$coef[mysumrows,2]
          cilist<-suppressMessages(confint(lm1))[mysumrows,1:2]
          pvallist<-summary(lm1)$coef[mysumrows,4]
          Resultsmat[rownum,4:7]<-c(betalist,cilist,pvallist)
        }
      } else { #below if ixterm is NULL
        if(predspline==T){
          mysumrows<-c(2:(predsplinedf + 1))
        } else if(class(Data[,prednames[j]])=="factor"){
          mysumrows<-c(2:(1+length(rownum)))
        } else {
          mysumrows<-2
        }

        if(mice){
          micesumm<-summary(pool(lm1),conf.int=T)
          betalist<-micesumm[mysumrows,2]
          selist<-micesumm[mysumrows,3]
          cilist<-micesumm[mysumrows,7:8]
          dflist<-micesumm[mysumrows,"df"]
          pvallist<-micesumm[mysumrows,6]
        } else if(any(class(lm1)=="logistf")){
          betalist<-lm1$coefficients[2]
          selist<-sqrt(diag(lm1$var))[2]
          cilist<-cbind(lm1$ci.lower[2],lm1$ci.upper[2])
          pvallist<-lm1$prob[2]
        } else if(robust==T){
          betalist<-robustcoef[mysumrows,1] 
          selist<-robustcoef[mysumrows,2]
          cilist<-robustcoef[mysumrows,5:6]
          pvallist<-robustcoef[mysumrows,4]
        } else {
          betalist<-summary(lm1)$coef[mysumrows,1]
          selist<-summary(lm1)$coef[mysumrows,2]
          cilist<-confint.lm(lm1)[mysumrows,1:2]
          pvallist<-summary(lm1)$coef[mysumrows,4]
        }

        Resultsmat[rownum,4:7]<-c(betalist,cilist,pvallist)
      }

      if(mice){
        if(logpred2[j]==T){
          IQRBeta<-betalist*log(1+(Mult/100),base=logbasepred)
          IQRCIs<-confint.qt(betalist,selist,dflist,log(1+(Mult/100),
                                                        base=logbasepred))
        } else {
          IQRBeta<-betalist*Mult
          IQRCIs<-confint.qt(betalist,selist,dflist,Mult)
        }
      } else if(any(class(lm1)=="logistf")|any(class(lm1)=="glm")&robust==F){
        if(Mult!=1){
          warning("The multiplying factor for beta was ignored because it's 
          unclear how to get profile penalized log likelihood CIs for logistf or 
          profile log likelihood CIs for glm for a change in predictor not equal to 1.")
        }
        IQRBeta<-betalist
        IQRCIs<-cilist
      } else {
        if(logpred2[j]==T){
          IQRBeta<-betalist*log(1+(Mult/100),base=logbasepred)
          IQRCIs<-confint.qt(betalist,selist,summary(lm1)$df[2],
                             log(1+(Mult/100),base=logbasepred))
        } else {
          IQRBeta<-betalist*Mult
          IQRCIs<-confint.qt(betalist,selist,summary(lm1)$df[2],Mult)
        }
      }
      
      if(plotPercChange){
        Resultsmat[rownum,8]<-(exp(IQRBeta)-1)*100
        Resultsmat[rownum,9:10]<-(exp(IQRCIs)-1)*100
      } else if(plotExpBeta){
        Resultsmat[rownum,8]<-exp(IQRBeta)
        Resultsmat[rownum,9:10]<-exp(IQRCIs)
      } else if(logout2[i]){
        Resultsmat[rownum,8] <- logbaseout^IQRBeta
        Resultsmat[rownum,9:10] <- logbaseout^IQRCIs
      } else {
        Resultsmat[rownum,8]<-IQRBeta
        Resultsmat[rownum,9:10]<-IQRCIs
      }
      
      
      if(LOOCV==T){
        if(isbinary){
          t1<-train(form1,ccDat,trControl = trainControl(method="repeatedcv",number=10,
                                                         repeats=50),method="lm")
          Resultsmat[rownum,11]<-t1$results$RMSE
        } else if(any(class(lm1)=="glm")){
          t1<-train(formula(paste0("factor(",outexpr,")~",as.character(form1)[3])),ccDat,
                    trControl = trainControl(method="repeatedcv",number=10,repeats=50),
                    method="glm",family=binfunc(link=binomlink))
          Resultsmat[rownum,11]<-t1$results$Accuracy
        } else {
          Resultsmat[rownum,11]<-NA
        }
      } else {
        Resultsmat[rownum,11]<-NA
      }
      if(mice){
        Resultsmat[rownum,"AIC"]<-NA
      } else if(any(class(lmlist[[lmnum]])=="logistf")){
        extractAIC.logistf<-function(fit, scale, k=2, ...){
          dev<- -2 * (fit$loglik['null']-fit$loglik['full'])
          AIC<- dev+k*fit$df
          edf<- fit$df
          return(AIC)
        }
        Resultsmat[rownum,"AIC"]<-extractAIC.logistf(lmlist[[lmnum]])
      } else {
        Resultsmat[rownum,"AIC"]<-AIC(lmlist[[lmnum]])
      }
      if(extradiag){
        if(mice){
          warning("Extra diagnostics not attainable for imputed models; ignoring extradiag.")
        } else if(!any(class(lmlist[[lmnum]])=="glm")){
          cooks<-cooks.distance(lmlist[[lmnum]])
          Extramat[rownum,"Heterosk.p"]<-lmtest::bptest(lmlist[[lmnum]])$p.value
          Extramat[rownum,"N.Cooks.D_>_0.5"]<-length(which(cooks>0.5))
          Extramat[rownum,"MaxCooks.D"]<-max(cooks,na.rm=T)
        } else {
          cooks<-cooks.distance(lmlist[[lmnum]])
          altmod<-glm(form1,data=lmlist[[lmnum]]$model,family=lmlist[[lmnum]]$family)
          simout<-DHARMa::simulateResiduals(altmod)
          quanttest<-DHARMa::testQuantiles(simout,plot=F)
          outliertest<-DHARMa::testOutliers(simout,plot=F,type='bootstrap')
          uniftest<-DHARMa::testUniformity(simout,plot=F)
          disptest<-DHARMa::testDispersion(simout,plot=F)
          Extramat[rownum,"Heterosk.p"]<-quanttest$p.value
          Extramat[rownum,"MinQuantileSplinePval"]<-min(quanttest$pvals)
          Extramat[rownum,"N.Cooks.D_>_0.5"]<-length(which(cooks>0.5))
          Extramat[rownum,"MaxCooks.D"]<-max(cooks,na.rm=T)
          Extramat[rownum,"UnifTestPval"]<-uniftest$p.value
          Extramat[rownum,"OutlierTestPval"]<-outliertest$p.value
          Extramat[rownum,"DispTestPval"]<-disptest$p.value
        }
      }
      if(post.power==T){
        Xmat<-model.matrix(form1,ccDat)[,-1]
        Xmatbinary<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==T)
        Xmatcont<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==F)
        Xmat<-Xmat[,c(Xmatcont,Xmatbinary)]
        Xmatbinary<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==T)
        Xmatcont<-which(apply(Xmat,2,function(x) setequal(unique(x),c(0,1)))==F)
        M<-apply(Xmat,2,calc_moments)
        catmarginals<-as.list(1-M[1,Xmatbinary])
        support <- list() # default support will be generated inside simulation
        seeds<-sample(1:(nsim*10),nsim)
        powersimfunc<-function(seed=seeds){
          if(is.null(ixterm)){
            myvc<-cov2cor(cov(Xmat))
            invisible(
              Simmat<-rcorrvar(nrow(Xmat),length(Xmatcont),length(Xmatbinary),means = M[1,Xmatcont],
                               vars =  (M[2,Xmatcont])^2, skews = M[3,Xmatcont], skurts = M[4,Xmatcont],
                               fifths = M[5,Xmatcont], sixths = M[6,Xmatcont], marginal = catmarginals,
                               rho = myvc, seed = seed)
            )
            newX<-cbind(Simmat$continuous_variables,Simmat$ordinal_variables-1)
            names(newX)<-dimnames(Xmat)[[2]]
          } else {
            myvc<-cov2cor(cov(Xmat[,-length(Xmatcont)]))
            invisible(
              Simmat<-rcorrvar(nrow(Xmat),length(Xmatcont)-1,length(Xmatbinary),
                               means = M[1,Xmatcont[-length(Xmatcont)]],
                               vars =  (M[2,Xmatcont[-length(Xmatcont)]])^2,
                               skews = M[3,Xmatcont[-length(Xmatcont)]],
                               skurts = M[4,Xmatcont[-length(Xmatcont)]],
                               fifths = M[5,Xmatcont[-length(Xmatcont)]],
                               sixths = M[6,Xmatcont[-length(Xmatcont)]], marginal = catmarginals,
                               rho = myvc, seed = seed)
            )
            newX<-cbind(Simmat$continuous_variables,Simmat$ordinal_variables-1)
            names(newX)<-dimnames(Xmat)[[2]][-length(Xmatcont)]
            newX[,dimnames(Xmat)[[2]][length(Xmatcont)]]<-
              newX[,dimnames(Xmat)[[2]][1]]*newX[,dimnames(Xmat)[[2]][length(Xmatcont)+1]]
            newX<-newX[,c(1:(length(Xmatcont)-1),ncol(newX),length(Xmatcont):(ncol(newX)-1))]
          }
          betas<-lm1$coef[names(newX)]+sapply(summary(lm1)$coef[names(newX),2],
                                              function(x) rnorm(1,0,x))
          #keep general beta structure but introduce noise with SD = beta SE
          if(is.null(ixterm)){
            betas[1]<-effect.size
          } else {
            betas[length(betas)]<-effect.size
          }
          my.y<-as.matrix(newX)%*%betas+lm1$coef[1]+rnorm(nrow(newX),sd(lm1$residuals)) #keep former residual SD
          newDat<-cbind(my.y,newX)
          names(newDat)[1]<-outnames[i]
          names(newDat)<-gsub(" ","_",names(newDat))
          names(newDat)<-gsub("$","",names(newDat),fixed=T)
          names(newDat)<-gsub(":","_",names(newDat),fixed=T)
          names(newDat)<-gsub("+","_",names(newDat),fixed=T)
          names(newDat)<-gsub("-","_",names(newDat),fixed=T)
          names(newDat)<-gsub("<","_",names(newDat),fixed=T)
          names(newDat)<-gsub(">","_",names(newDat),fixed=T)
          newlm<-lm(formula(paste0(outnames[i],"~",paste(names(newDat)[-1],collapse="+"))),newDat)
          if(robust==T){
            newcoefs<-robustse(newlm)
            rets<-newcoefs[2,c(4)]
          } else {
            rets<-summary(newlm)$coef[2,c(4)]
          }
          return(rets)
        }
        pboptions(type="timer")
        pvalvec<-pbsapply(seeds,powersimfunc)
        pow<-length(which(pvalvec<0.05))/length(pvalvec)
        Resultsmat[rownum,"Posthoc_Power"]<-pow
      }
    } #end j loop
  } #end i loop
  if(post.power==T){
    if(!is.null(ixterm)|predspline==T){
      names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                           "Beta.Transform","LCI.Transform","UCI.Transform",
                           "K-fold RMSE","Variable","AIC","Posthoc_Power")
    } else {
      names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                           "Beta.Transform","LCI.Transform","UCI.Transform",
                           "K-fold RMSE","AIC","Posthoc_Power")
    }
  } else {
    if(!is.null(ixterm)|predspline==T){
      names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                           "Beta.Transform","LCI.Transform","UCI.Transform",
                           "K-fold RMSE","Variable","AIC")
    } else {
      names(Resultsmat)<-c(Outtitle,Predtitle,"N","Beta","LCI","UCI","p-value",
                           "Beta.Transform","LCI.Transform","UCI.Transform",
                           "K-fold RMSE","AIC")
    }
  }

  if(!is.null(ixterm)|predspline==T){
    if(post.power==T){
      Resultsmat<-Resultsmat[,c(1:2,12,3:11,13:14)]
    } else {
      Resultsmat<-Resultsmat[,c(1:2,12,3:11,13)]
    }
  }
  if(!LOOCV){
    Resultsmat$'K-fold RMSE'<-NULL
    Resultsmat$'K-fold Accuracy'<-NULL
  }
  Resultsmat$Significant<-ifelse(Resultsmat$'p-value'<0.05,"Yes","")
  Resultsmat<-Resultsmat[which(apply(Resultsmat,1,function(x) !all(is.na(x)))),]
  if(extradiag){
    if(isbinary){
      names(Extramat)[which(names(Extramat)=="Heterosk.p")]<-"QuantileSplinePval"
    }
    Resultsmat<-cbind(Resultsmat,Extramat)
  }
  Resmatout<-Resultsmat

  names(Resultsmat)[2]<-"Predictor"

  Resultsmat[,1]<-factor(Resultsmat[,1],levels=outnames)
  if(!is.null(altoutnames)) levels(Resultsmat[,1])<-altoutnames
  Resultsmat[,2]<-factor(Resultsmat[,2],levels=prednames)
  if(!is.null(altprednames)) levels(Resultsmat[,2])<-altprednames
  if(!is.null(altixname)) Resultsmat$Variable<-gsub(ixterm,altixname,Resultsmat$Variable)

  if(MyMult=="IQR"){
    Multname<-"IQR"
  } else {
    Multname<-Mult
  }

  if(is.null(facetcol)){
    facetcol<-length(outnames)
  }

  if(!is.null(ixterm)&class(Data[,ixterm])=="factor"){
    plotData<-Resultsmat[grep("|",Resultsmat$Variable,fixed=T),]
    plotData$Interaction_Level<-gsub(".*[|]","",plotData$Variable)
    plotData$star<-ifelse(plotData$Significant=="Yes","*","")
    plotData.ix<-Resultsmat[grep("*",Resultsmat$Variable,fixed=T),]
    plotData.ix$Interaction_Level<-gsub(".*[*]","",plotData.ix$Variable)
    plotData.ix$star<-ifelse(plotData.ix$Significant=="Yes","*","")
  } else if(!is.null(ixterm)) {
    plotData<-Resultsmat[-grep(":",Resultsmat$Variable,fixed=T),]
    plotData$Interaction_Level<-plotData$Variable
    plotData$star<-ifelse(plotData$Significant=="Yes","*","")
    plotData.ix<-Resultsmat[grep(":",Resultsmat$Variable,fixed=T),]
    plotData.ix$Interaction_Level<-plotData.ix$Variable
    plotData.ix$star<-ifelse(plotData.ix$Significant=="Yes","*","")
  } else {
    plotData<-Resultsmat
    plotData$star<-ifelse(plotData$Significant=="Yes","*","")
  }


  if(predspline==F){
    if(all(logout)|isbinary) horint<-1 else horint<-0
    if(all(logout)|isbinary|plotExpBeta){
      yl <- "Exp(Beta)"
    } else if(plotPercChange) yl <- "% Change" else yl <- "Beta"
    
    if(is.null(ixterm)){
      gg1<-ggplot(data=plotData,aes(x=Predictor,y=Beta.Transform,
                                    color=Significant))+
        geom_errorbar(aes(ymin=LCI.Transform,ymax=UCI.Transform),
                      width=0.4,linewidth=1)+
        geom_point()+geom_hline(aes(yintercept=horint))+
        geom_text(aes(y=UCI.Transform,vjust=-0.5,
                      label=star),show.legend=F,size=6)+
        facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
        ylab(yl)+xlab(paste0(Predtitle))+
        scale_color_manual(values=c("black","red"),guide="none")+
        theme(axis.text.x=element_text(angle=45,hjust=1))
    } else {
      if(class(Data[,ixterm])=="factor"){
        gg1.x<-ggplot(data=plotData,aes(x=Predictor,y=Beta.Transform,
                                        color=Interaction_Level))+
          geom_errorbar(aes(ymin=LCI.Transform,ymax=UCI.Transform),
                        width=0.4,linewidth=1,position=position_dodge(0.75))+
          geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=horint))+
          geom_text(aes(y=UCI.Transform,label=star),vjust=-0.5,
                    position=position_dodge(0.75),show.legend=F,size=6)+
          facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
          ylab("Exp(Beta)")+xlab(paste0(Predtitle))+ggtitle("Marginal Coefficients")+
          scale_color_brewer(palette="Set1",name="Interaction Level")+
          theme(axis.text.x=element_text(angle=45,hjust=1),
                legend.position = "bottom",
                plot.title=element_text(hjust=0.5))
        gg2<-ggplot(data=plotData.ix,
                    aes(x=Predictor,y=Beta.Transform,color=Interaction_Level))+
          geom_errorbar(aes(ymin=LCI.Transform,ymax=UCI.Transform),
                        width=0.4,linewidth=1,position=position_dodge(0.75))+
          geom_point(position=position_dodge(0.75))+geom_hline(aes(yintercept=horint))+
          geom_text(aes(y=UCI.Transform,label=star),vjust=-0.5,
                    position=position_dodge(0.75),show.legend=F,size=6)+
          facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
          ylab(yl)+xlab(paste0(Predtitle))+ggtitle("Interaction Coefficients")+
          scale_color_brewer(palette="Dark2",name="Interaction")+
          theme(axis.text.x=element_text(angle=45,hjust=1),
                legend.position = "bottom",
                plot.title=element_text(hjust=0.5))
        gg1<-cowplot::plot_grid(gg1.x,gg2,nrow=2)
      } else {
        gg1<-ggplot(data=plotData.ix,
                    aes(x=Predictor,y=Beta.Transform,color=Significant))+
          geom_errorbar(aes(ymin=LCI.Transform,ymax=UCI.Transform),
                        width=0.4,linewidth=1)+
          geom_point()+geom_hline(aes(yintercept=horint))+
          geom_text(aes(y=UCI.Transform,label=star),vjust=-0.5,
                    show.legend=F,size=6)+
          facet_wrap(formula(paste0(".~",Outtitle)),scales="free",ncol=facetcol)+
          ylab(yl)+xlab(paste0(Predtitle))+ggtitle("Interaction Coefficients")+
          scale_color_manual(values=c("black","red"),guide="none")+
          theme(axis.text.x=element_text(angle=45,hjust=1),
                legend.position = "bottom",
                plot.title=element_text(hjust=0.5))
      }
      
    } 
  } else {
    warning(paste0("No plots are being output, as it was unclear how best to",
                   " plot results when the predictor is a spline."))
    gg1<-list()
  }

  retlist<-list(Resmatout,gg1,lmlist)
  names(retlist)<-c("Matrix","GGplot","LMlist")

  return(retlist)
}
