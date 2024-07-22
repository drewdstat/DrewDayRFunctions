#' Run multiple GAMs and return results
#' 
#' \code{GAMResults} runs GAM models with smooths on a series of predictor 
#' variables of interest for all combinations of predictor and outcome variables, 
#' allowing for categorical and continuous interaction terms. A blog post with 
#' some nice visualizations and explanations of what GAMs are can be
#' \href{https://ecogambler.netlify.app/blog/interpreting-gams/}{found here}.
#'
#' @param prednames A character vector of the predictor variables.
#' @param outnames A character vector of the outcome variables. The function 
#' currently allows for numeric or binary ("factor" or "character" class) outcome
#' variables.
#' @param covnames A character vector of all covariates or a list of character 
#' vectors of covariates of the same length as \code{length(prednames) * 
#' length(outnames)} so that each unique combination of prednames and outnames 
#' has a defined set of covariates. Note that in ordering this list the function 
#' loops through the \code{prednames} and then the \code{outnames} (e.g., 
#' Outcome 1 - Predictor 1, Outcome 1 - Predictor 2, Outcome 1 - Predictor 3, 
#' Outcome 2 - Predictor 1, etc.), and so you should order the list of covariates 
#' accordingly. These covariates will be included as linear (i.e., not smooth) 
#' terms in the model. These can include spline terms, such as
#' \code{\link[splines]{ns}} natural splines or \code{\link[splines]{bs}} 
#' b-splines as defined by the 'splines' package (e.g., \code{"ns(Year, df=3)"}).
#' These can also include interaction terms denoted by the ":" separator , though 
#' you should make sure to include the main effects too (e.g., \code{c("Income", 
#' "HouseholdSize", "Income:HouseholdSize")}).
#' @param Data A data frame containing all prednames, outnames, and covnames as 
#' columns.
#' @param logout A TRUE/FALSE value for whether to natural log-transform each 
#' outcome variable. Defaults to FALSE.
#' @param logpred A TRUE/FALSE value for whether to natural log-transform each
#' predictor variable. Defaults to FALSE.
#' @param Outtitle A single character value for the name of the column containing
#' outcome variables in the data frame of results to be returned. Defaults to 
#' "Outcome".
#' @param Predtitle A single character value for the name of the column containing
#' predictor variables in the data frame of results to be returned. Defaults to
#' "Exposure".
#' @param ixterm A single character value for the name of an interacting variable.
#' This variable can be either categorical, in which case the column should be 
#' of a character or factor class, or it can be continuous.
#' @param smooth.select A TRUE/FALSE value for whether to add an extra penalty to
#' each smooth term so that it can be penalized to zero. Defaults to FALSE, and 
#' this corresponds to the 'select' parameter in the \code{\link[mgcv]{gam}} function.
#' @param bs_s A character value wrapped in single quotes reflecting the basis 
#' spline to be used in the \code{\link[mgcv]{s}} smooths, which will be used as 
#' the smooth terms in all models except if a continuous * continous interaction 
#' is included (i.e., if \code{ixterm} is a continuous variable, in which case 
#' \code{\link[mgcv]{ti}} smooths are used. Defaults to \code{"'tp'"}, indicating 
#' thin-plate regression splines (TPRS) will be used for the basis spline (see 
#' \code{?tprs}). This is also the default basis spline for \code{\link[mgcv]{s}}.
#' @param K_s An integer reflecting the \code{k} parameter for \code{\link[mgcv]{s}}
#' smooth terms. The meaning of this term is that the maximum effective degrees 
#' of freedom (EDF) of each smooth is this value minus one (see \code{?choose.k}). 
#' This defaults to -1, which is the default value for \code{\link[mgcv]{s}}. 
#' Setting K_s to -1 is the equivalent of not specifying a value, in which case 
#' the 'k' value default value for \code{\link[mgcv]{s}} is 10 (i.e., maximum 
#' EDF=9) when there is a single term in the smooth (e.g., \code{s(X1)}) and 
#' \code{bs="tp"}. For your information, the default \code{k} when using 
#' two-term smooths (e.g., \code{s(X1, X2)}) for \code{\link[mgcv]{s}} when 
#' \code{bs="tp"} is 30 and for >=3-term smooths this is 110 (see ?tprs). 
#' Note that increasing this value increases the subspace of functions (can be 
#' thought of as a ceiling on possible EDF values), and therefore higher 
#' \code{K} values can lead to higher EDF (see \code{?choose.k}).
#' @param bs_ti A character value wrapped in single quotes reflecting the basis 
#' spline to be used in the \code{\link[mgcv]{ti}} tensor product interaction 
#' smooths, which will be used only if a continuous * continous interaction 
#' is included (i.e., if \code{ixterm} is a continuous variable). Defaults to 
#' \code{"'cr'"}, indicating cubic regression splines will be used for the basis 
#' spline (see \code{?smooth.construct.cr.smooth.spec}).
#' @param K_ti An integer reflecting the \code{k} parameter for 
#' \code{\link[mgcv]{ti}}. This defaults to NA, which is also the default \code{k}
#' for \code{\link[mgcv]{ti}}, and which means default \code{k} values will be used.
#' The documentation shown with \code{?smooth.construct.cr.smooth.spec} suggests 
#' that the default K is 10, but in practice it appears that the default \code{k}
#' is actually 5 based on the \code{bsdim} values output in my own experiments. 
#' This means 4 knots in the univariate \code{\link{ti}} terms and 4 each for the 
#' bivariate \code{ti(var1, var2)} terms. 
#' @param na.action This is the \code{\link[stats]{na.action}} that will be 
#' passed to \code{\link[mgcv]{gam}}. Defaults to \code{\link[stats]{na.exclude}}. 
#' @param plotresid A TRUE/FALSE value indicating whether the GAM plots showing
#' the smooth fit and CIs should also include points showing the model residuals 
#' on the plot (see \code{\link[mgcViz]{l_points}}). Defaults to FALSE.
#' @param quantiletrimcheck A TRUE/FALSE value indicating whether a certain 
#' quantile of one or both of the predictor and/or outcome variables should be
#' trimmed as a sensitivity analysis. Defaults to FALSE. If TRUE, the data will 
#' be trimmed at the top \code{trimperc} percentiles. The purpose of this check
#' is to make sure that the observed association curves don't change radically
#' when removing edge points considering that splines can be at times sensitive
#' to data at the edges of the distribution.
#' @param trimdim A character value or vector of \code{"x"}, \code{"y"}, or 
#' \code{c("x", "y")} indicating whether to trim the predictor (x) or outcome (y)
#' values when \code{quantiletrimcheck} is TRUE. Defaults to \code{c("x", "y")}.
#' @param trimperc An integer or numeric value indicating total percentage of the 
#' data to be trimmed from the high and low ends of the data distribution if 
#' \code{quantiletrimcheck} is TRUE. Defaults to 5, meaning that everything below 
#' the 2.5th percentile and above the 97.5th percentile will be removed before 
#' running the GAM models.
#' 
#' @details
#' \code{GAMResults} bases its models on the assumption that there will be only 
#' one univariate smooth term in each model for the predictor of interest and that
#' utilizes the smooth function \code{\link[mgcv](s)} with a user-defined basis spline
#' and "k" parameter. This is akin to \code{gam(y ~ s(x, bs = bs_s, k = K_s) + z, ...)} 
#' where y is an outcome variable, x is a predictor variable of interest, and z 
#' is a covariate. If an outcome variable is binary, a logistic GAM is run. 
#' 
#' In the case of a continuous interacting variable (i.e., \code{ixterm}), this 
#' function utilizes \code{\link[mgcv](ti)} for univariate tensor product 
#' interaction terms for both the predictor of interest and interaction term 
#' (similar to main effects in a linear model with a multiplicative interaction 
#' term), as well as a bivariate tensor product interaction for those two variables. 
#' This is the same as \code{gam(y ~ ti(x, bs = bs_ti, k = K_ti) + ti(iota, bs = bs_ti, 
#' k = K_ti) + ti(x, iota, bs = bs_ti, k = K_ti) + z, ...)} where iota is the interacting 
#' variable. This is implemented as it is in the 
#' \href{https://pennlinc.github.io/ModelArray/reference/generator_gamFormula_continuousInteraction.html}{ModelArray R package} 
#' and as succinctly described on this linguistics professor's 
#' \href{https://janhove.github.io/posts/2017-06-26-continuous-interactions/}{blog post}.
#' 
#' In the case of a categorical interacting variable (i.e., \code{ixterm}), two 
#' separate models are run to obtain all the relevant parameter estimates. The 
#' first model is structured as \code{gam(y ~ s(x, bs = bs_s, k = K_s, by = iota) 
#' + iota + z)}, while the second is \code{gam(y ~ s(x, bs = bs_s, k = K_s, by = 
#' ordered(iota)) + s(x, bs = bs_s, k = K_s) + iota + z)}. From the first model
#' we are able to obtain the factor level-specific curves, and from the second
#' model we can obtain curves of the differences between those level-specific 
#' curves. Note that theoretically the curve for the referent level produced by 
#' the second term "s(x)" in the second model should theoretically be equivalent
#' to the referent level curve produced by the first term in the first model, but 
#' in practice they tend to have minuscule difference in EDF, exact fitted curve 
#' values, etc., but it's close enough I find it convenient to grab terms from
#' these two models as if they were equivalent. Further explanation of both models can 
#' be found at 
#' \href{https://stats.stackexchange.com/questions/416991/gam-factor-smooth-interaction-include-main-effect-smooth}{this CrossValidated blog post}. 
#' Examples of the second model are shown in \href{https://fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/}{this blog post}
#' and in example papers like \href{https://www.nature.com/articles/s44161-022-00131-8}{Zhernakova et al. 2022}.
#'  
#' @return \code{GAMResults} returns a list of:
#' \item{Matrix}{A data frame of important results from the GAM models.}
#' \item{Plots}{A list of ggplot plots of class \code{c("plotSmooth", "gg")} showing 
#' the predictor of interest curves plotted against predicted outcome values. 
#' These plots are based on the mgcViz package plotting functions \code{l_ciPoly}, 
#' \code{l_fitLine}, and \code{l_rug}.}
#' \item{GAMlist}{A list of all GAM models stored as objects of class "gam".}
#' \item{IxGAMlist}{If ixterm is not NULL and is a categorical variable, the 
#' additional GAM models that were run to calculate the interaction term 
#' coefficients are stored here as a list of "gam" class objects.}
#' 
#' @import mgcv mgcViz dplyr
#' @export GAMResults
#' 
#' @examples
#' # The first example involves analyzing several GAMs with no interaction terms.
#' 
#' set.seed(100)
#' gamdat1 <- data.frame(X1 = rnorm(500), X2 = rnorm(500), Z = rnorm(500,0,2))
#' 
#' gamdat1$yhat1 <- exp(sin(gamdat1$X1)^0.5) + gamdat1$Z * 0.75 - 2
#' gamdat1$yhat2 <- exp(cos(gamdat1$X2)^2) + gamdat1$Z * 2 + 3
#' 
#' set.seed(102)
#' gamdat1$y1 <- gamdat1$yhat1 + rnorm(500)
#' set.seed(103)
#' gamdat1$y2 <- gamdat1$yhat2 + rnorm(500)
#' 
#' gamres1 <- GAMResults(prednames = c("X1", "X2"), outnames = c("y1", "y2"), 
#'                       covnames = "Z", Data = gamdat1)
#' 
#' gridPrint(grobs = gamres1$Plots, ncol = 2)
#' 
#' # This next example involves a factor-smooth interaction
#' 
#' set.seed(100)
#' gamdat2 <- data.frame(X = rnorm(500), Z = rnorm(500,0,2), 
#'                       S = as.factor(sample(c("M","F"), 500, replace = T)))
#' gamdat2$yhat <- ifelse(gamdat2$S == "M", sin(gamdat2$X) * 3, 
#'                        sin(gamdat2$X) * 0.5)
#' gamdat2$yhat <- gamdat2$yhat + gamdat2$Z * 2
#' set.seed(101)
#' gamdat2$y <- gamdat2$yhat + rnorm(500)
#' 
#' gamres2 <- GAMResults("X", "y", "Z", gamdat2, ixterm = "S")
#' 
#' gridPrint(grobs = gamres2$Plots$`y~X`)
#' 


GAMResults <- function(prednames, outnames, covnames = NULL, Data, logout = F, 
                     logpred=F, Outtitle = "Outcome", Predtitle = "Exposure",
                     ixterm = NULL, smooth.select = F, bs_s = "'tp'", K_s = -1,
                     bs_ti = "'cr'", K_ti = NA, na.action = na.exclude, 
                     plotresid = F, quantiletrimcheck = F, trimdim = c("x", "y"), 
                     trimperc = 5){
  if(!is.null(ixterm)){
    if(class(Data[,ixterm])%in%c("factor","character")){
      nixlvl<-length(levels(Data[,ixterm]))
      if(class(Data[,ixterm])=="character"){Data[,ixterm]<-as.factor(Data[,ixterm])}
      ixtermlength<-nixlvl+ncol(combn(levels(Data[,ixterm]),2))
      K <- K_s
      bs <- bs_s
    } else {
      ixtermlength <- 3
      K <- K_ti
      bs <- bs_ti
    }
    Resultsmat<-as.data.frame(matrix(NA,ixtermlength*length(prednames)*length(outnames),
                                     12))
  } else {
    Resultsmat<-as.data.frame(matrix(NA,length(prednames)*length(outnames),11))
    K <- K_s
    bs <- bs_s
  }
  
  lmlist<-list()
  if(!is.null(ixterm)&&class(Data[,ixterm])=="factor"){
    lmlist2<-list()
  }
  
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
    for(j in 1:length(prednames)){
      if(!is.null(ixterm)){
        rowstart<-(1+(j-1))+((i-1)*ixtermlength*length(prednames))+
          ((ixtermlength-1)*(j-1))
        rowend<-(ixtermlength+(j-1))+((i-1)*ixtermlength*length(prednames))+
          ((ixtermlength-1)*(j-1))
        rownum<-rowstart:rowend
        lmnum<-(1+(j-1))+((i-1)*length(prednames))
      } else {
        rownum<-(1+(j-1))+((i-1)*length(prednames))
        lmnum<-rownum
      }
      
      if(logout2[i]==T){
        outexpr<-paste0("log(",outnames[i],")")
      } else {
        outexpr<-outnames[i]
      }
      
      if(logpred2[j]==T){
        predexpr<-paste0("log(",prednames[j],")")
      } else {
        predexpr<-prednames[j]
      }
      
      if(is.null(covnames)){
        if(is.null(ixterm)){
          ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j])]),]
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr))
          gamform1<-formula(paste0(outexpr,"~s(",predexpr,",bs=",bs,",k=",K,")"))
        } else {
          ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],ixterm)]),]
          nobs<-dim(ccDat)[1]
          form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm))
          if(class(Data[,ixterm])=="factor"){
            gamform1<-
              formula(paste0(outexpr,"~s(",predexpr,",by=",ixterm,",bs=",bs,
                             ",k=",K,")+", ixterm))
            gamform2<-
              formula(paste0(outexpr,"~s(",predexpr,",by=ordered(",ixterm,
                             "),bs=",bs,",k=",K,")+s(",predexpr,",bs=",
                             bs,",k=",K,")+", ixterm))
          } else {
            gamform1<-formula(paste0(outexpr,"~ti(",predexpr,",bs=",bs,")+ti(",
                                     ixterm, ",bs=",bs,")+ti(",predexpr,",",
                                     ixterm,",bs=",bs,")"))
          }
          
        }
      } else {
        if(is.list(covnames)){
          if(length(covnames)!=(length(prednames)*length(outnames))){
            stop(paste0("Error: covnames must have a length of ", 
                        length(prednames)*length(outnames)))
          }
          mycovnames<-covnames[[lmnum]]
        } else {
          mycovnames<-covnames
        }
        
        if(any(grepl("(",mycovnames,fixed=T))&any(grepl(":",mycovnames,fixed=T))){
          mycovcols<-gsub(".*\\(","",mycovnames)
          mycovcols<-gsub(",.*","",mycovcols)
          mycovcols<-mycovcols[-which(grepl(":",mycovcols,fixed=T))]
        } else if(any(grepl(":",mycovnames,fixed=T))){
          mycovcols<-mycovnames[-which(grepl(":",mycovnames,fixed=T))]
        } else if(any(grepl("(",mycovnames,fixed=T))){
          mycovcols<-gsub(".*\\(","",mycovnames)
          mycovcols<-gsub(",.*","",mycovcols)
        } else {
          mycovcols<-mycovnames
        }
        
        if(is.null(ixterm)){
          ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovcols)]),]
          form1<-formula(paste0(outexpr,"~",predexpr,"+",paste(mycovnames,collapse="+")))
          gamform1<-formula(paste0(outexpr,"~s(",predexpr,",bs=",bs,", k=",K,")+",
                                   paste(mycovnames,collapse="+")))
        } else {
          ccDat<-Data[complete.cases(Data[,c(outnames[i],prednames[j],mycovcols,ixterm)]),]
          form1<-formula(paste0(outexpr,"~",predexpr,"*",ixterm,"+",
                                paste(mycovnames,collapse="+")))
          if(class(Data[,ixterm])=="factor"){
            gamform1<-
              formula(paste0(outexpr,"~s(",predexpr,",by=",ixterm,",bs=",bs,
                             ",k=",K,")+", ixterm,"+",paste(mycovnames,collapse="+")))
            gamform2<-
              formula(paste0(outexpr,"~s(",predexpr,",by=ordered(",ixterm,
                             "),bs=",bs,",k=",K,")+s(",predexpr,",bs=",bs,
                             ",k=",K,")+",ixterm,"+",paste(mycovnames,collapse="+")))
          } else {
            gamform1<-formula(paste0(outexpr,"~ti(",predexpr,",bs=",bs,",k=",
            K,")+ti(",ixterm,",bs=",bs,",k=",K,")+ti(",
                                     predexpr,",",ixterm,",bs=",bs,",k=",K,")+",
                                     paste(mycovnames,collapse="+")))
          } 
        }
      }
      
      if(quantiletrimcheck){
        trimperc<-trimperc/2
        if("x" %in% trimdim){
          Data<-Data[-which(Data[,prednames[j]]<quantile(Data[,prednames[j]],trimperc/100,na.rm=T)|
                              Data[,prednames[j]]>quantile(Data[,prednames[j]],
                                                           1-(trimperc/100),na.rm=T)),]
          ccDat<-ccDat[-which(ccDat[,prednames[j]]<quantile(ccDat[,prednames[j]],trimperc/100)|
                              ccDat[,prednames[j]]>quantile(ccDat[,prednames[j]],
                                                           1-(trimperc/100))),]
          if(!is.null(ixterm)&&class(Data[,ixterm])!="factor"){
            Data<-Data[-which(Data[,ixterm]<quantile(Data[,ixterm],trimperc/100,na.rm=T)|
                                Data[,ixterm]>quantile(Data[,ixterm],
                                                             1-(trimperc/100),na.rm=T)),]
            ccDat<-ccDat[-which(ccDat[,ixterm]<quantile(ccDat[,ixterm],trimperc/100)|
                                  ccDat[,ixterm]>quantile(ccDat[,ixterm],
                                                                1-(trimperc/100))),]
          }
        }
        if("y" %in% trimdim){
          if(length(unique(Data[which(!is.na(Data[, outnames[i]])), outnames[i]])) == 2){
            warning(paste0("'y' trim dimension does not make sense for a binary ", 
                           "outcome. Ignoring 'y' trim dimension."))
          } else {
            Data<-Data[-which(Data[,outnames[i]]<quantile(Data[,outnames[i]],trimperc/100)|
                                Data[,outnames[i]]>quantile(Data[,outnames[i]],
                                                            1-(trimperc/100))),]
            ccDat<-ccDat[-which(ccDat[,outnames[i]]<quantile(ccDat[,outnames[i]],trimperc/100)|
                                  ccDat[,outnames[i]]>quantile(ccDat[,outnames[i]],
                                                               1-(trimperc/100))),]
          }
        }
      }
      
      if(length(unique(Data[which(!is.na(Data[, outnames[i]])), outnames[i]])) == 2){
        lm1<-glm(form1,data=Data,na.action=na.action,family="binomial")
        g1<-gam(gamform1,data=Data,select=smooth.select,na.action=na.action,family=binomial)
        g1.nona<-gam(gamform1,data=ccDat,select=smooth.select,na.action=na.action,
                     family=binomial)
        if(!is.null(ixterm)&&class(Data[,ixterm])=="factor"){
          g2<-gam(gamform2,data=Data,select=smooth.select,na.action=na.action,
                  family=binomial)
          g2.nona<-gam(gamform2,data=ccDat,select=smooth.select,na.action=na.action,
                       family=binomial)
          lmlist2[[lmnum]]<-g2.nona
        }
      } else {
        lm1<-lm(form1,data=Data,na.action=na.action)
        g1<-gam(gamform1,data=Data,select=smooth.select,na.action=na.action)
        g1.nona<-gam(gamform1,data=ccDat,select=smooth.select,na.action=na.action)
        if(!is.null(ixterm)&&class(Data[,ixterm])=="factor"){
          g2<-gam(gamform2,data=Data,select=smooth.select,na.action=na.action)
          g2.nona<-gam(gamform2,data=ccDat,select=smooth.select,na.action=na.action)
          lmlist2[[lmnum]]<-g2.nona
        }
      }
      lmlist[[lmnum]]<-g1.nona
      names(lmlist)[[lmnum]]<-paste0(outnames[i],"~",prednames[j])
      Resultsmat[rownum,1]<-outnames[i]
      Resultsmat[rownum,2]<-prednames[j]
      Resultsmat[rownum,3]<-nrow(ccDat)
      if(!is.null(ixterm)){
        if(class(Data[,ixterm])=="factor"){
          myvalmat<-cbind(c(summary(g1)$s.table[,1],
                            summary(g2)$s.table[1,1]),
                          c(summary(g1)$s.table[,4],
                            summary(g2)$s.table[1,4]),
                          c(k.check(g1.nona)[,1],
                            k.check(g2.nona)[1,1]),
                          c(k.check(g1.nona)[,3],
                            k.check(g2.nona)[1,3]),
                          c(k.check(g1.nona)[,4],
                            k.check(g2.nona)[1,4]),
                          summary(g1)$dev.expl*100
          )
          Resultsmat[rownum,12]<-
            c(dimnames(summary(g1)$s.table)[[1]],dimnames(summary(g2)$s.table)[[1]][1])
        } else {
          Resultsmat[rownum,12]<-
            dimnames(summary(g1)$s.table)[[1]][1:3]
        }
      }
      
      if(!is.null(ixterm)&&class(Data[,ixterm])!="factor"){
        mysumrows<-c(1:3)
      } else {
        mysumrows<-1
      }
      if(!is.null(ixterm)&&class(Data[,ixterm])=="factor"){
        Resultsmat[rownum,4:9]<-myvalmat
      } else if(!is.null(ixterm)&&class(Data[,ixterm])!="factor"){
        Resultsmat[rownum,4:9]<-cbind(summary(g1)$s.table[mysumrows,c(1,4)],
                                      k.check(g1.nona)[mysumrows,c(1,3:4)],
                                      (summary(g1)$dev.expl*100))
      } else {
        Resultsmat[rownum,4:9]<-c(summary(g1)$s.table[mysumrows,c(1,4)],
                                  k.check(g1.nona)[mysumrows,c(1,3:4)],
                                  (summary(g1)$dev.expl*100))
      }
      Resultsmat[rownum,10]<-AIC(lm1,g1)[1,2]-AIC(lm1,g1)[2,2]
      if(any(class(lm1)=="glm")){
        Resultsmat[rownum,11]<-anova(lm1,g1,test="Chisq")[2,5]
      } else {
        Resultsmat[rownum,11]<-anova(lm1,g1)[2,6]
      }
      
      Resultsmat[rownum,"AIC"]<-g1$aic
      
    } #j
  } #i
  
  
  if(!is.null(ixterm)){
    names(Resultsmat)<-c(Outtitle,Predtitle,"N","EDF","Smooth p-value","K'","K-index",
                         "K check p-value","%Dev Expl","deltaAIC","ANOVA p-value",
                         "Variable","AIC")
    Resultsmat<-Resultsmat[,c(1:2,12,3:11,13)]
  } else {
    names(Resultsmat)<-c(Outtitle,Predtitle,"N","EDF","Smooth p-value","K'","K-index",
                         "K check p-value","%Dev Expl","deltaAIC","ANOVA p-value","AIC")
  }
  Resultsmat$Significant<-ifelse(Resultsmat$'Smooth p-value'<0.05,"Yes","")
  Resmatout<-Resultsmat
  
  names(Resultsmat)[2]<-"Predictor"
  
  if(!is.null(ixterm)){
    if(class(Data[,ixterm])=="factor"){
      gg1<-list()
      for(i in 1:length(outnames)){
        for(j in 1:length(prednames)){
          ggnum<-(1+(j-1))+((i-1)*length(prednames))
          gg1[[ggnum]]<-list()
          names(gg1)[ggnum]<-paste0(outnames[i],"~",prednames[j])
          summtable<-summary(lmlist[[ggnum]])$s.table
          summtable2<-summary(lmlist2[[ggnum]])$s.table
          thistitle<-list(gg1=list(),gg2=list())
          for(k in 1:length(levels(Data[,ixterm]))){
            thistitle$gg1[[k]]<-
              paste0("EDF=",round(summtable[k,1],2),", p=",signif(summtable[k,4],3))
          }; rm(k)
          for(k in 1:length(grep("ordered",rownames(summtable2)))){
            thistitle$gg2[[k]]<-
              paste0("EDF=",round(summtable2[k,1],2),", p=",signif(summtable2[k,4],3))
          }; rm(k)
          thisylab<-lmlist[[ggnum]]$formula[[2]]
          allviz<-getViz(lmlist[[ggnum]])
          allviz2<-getViz(lmlist2[[ggnum]])
          totalplotlength<-length(levels(Data[,ixterm]))+length(grep("ordered",rownames(summtable2)))
          multplotpalette<-RColorBrewer::brewer.pal(totalplotlength,"Set1")
          for(k in 1:totalplotlength){
            if(k <= length(levels(Data[,ixterm]))){
              gg1[[ggnum]][[k]]<-plot(sm(allviz,k))+
                geom_hline(yintercept=0)+
                l_ciPoly(fill=multplotpalette[k],alpha=0.3)+
                l_fitLine(colour=multplotpalette[k],size=0.75)+
                l_rug(mapping = aes(x=x),alpha=0.8)+theme_bw()+
                ylab(thisylab)+ggtitle(thistitle$gg1[[k]])+
                xlab(rownames(summtable)[k])+
                theme(plot.title = element_text(hjust=0.5,size=15),
                      axis.title=element_text(size=12),
                      axis.text=element_text(size=11))
              names(gg1[[ggnum]])[k]<-rownames(summtable)[k]
            } else {
              gg1[[ggnum]][[k]]<-plot(sm(allviz2,k-length(levels(Data[,ixterm]))))+
                geom_hline(yintercept=0)+
                l_ciPoly(fill=multplotpalette[k],alpha=0.3)+
                l_fitLine(colour=multplotpalette[k],size=0.75)+
                l_rug(mapping = aes(x=x),alpha=0.8)+theme_bw()+
                ylab(thisylab)+ggtitle(thistitle$gg2[[k-length(levels(Data[,ixterm]))]])+
                xlab(rownames(summtable2)[k-length(levels(Data[,ixterm]))])+
                theme(plot.title = element_text(hjust=0.5,size=15),
                      axis.title=element_text(size=12),
                      axis.text=element_text(size=11))
              names(gg1[[ggnum]])[k]<-rownames(summtable2)[
                grep("ordered",rownames(summtable2))][k-length(levels(Data[,ixterm]))]
            }
          }; rm(k)
        }
      }
    } else {
      gg1<-list()
      for(i in 1:length(outnames)){
        for(j in 1:length(prednames)){
          ggnum<-(1+(j-1))+((i-1)*length(prednames))
          summtable<-summary(lmlist[[ggnum]])$s.table
          thistitle1<-paste0("EDF=",round(summtable[1,1],2),", p=",signif(summtable[1,4],3))
          thistitle2<-paste0("EDF=",round(summtable[2,1],2),", p=",signif(summtable[2,4],3))
          thistitle3<-paste0("EDF=",round(summtable[3,1],2),", p=",signif(summtable[3,4],3))
          thisylab<-lmlist[[ggnum]]$formula[[2]]
          allviz<-getViz(lmlist[[ggnum]])
          mygg1<-plot(sm(allviz,1))+
            geom_hline(yintercept=0)+
            l_ciPoly(fill="blue",alpha=0.3)+
            l_fitLine(colour="blue",size=0.75)+
            l_rug(mapping = aes(x=x),alpha=0.8)+theme_bw()+
            ylab(thisylab)+ggtitle(thistitle1)+
            theme(plot.title = element_text(hjust=0.5,size=15),
                  axis.title=element_text(size=12),
                  axis.text=element_text(size=11))
          thistitle<-paste0("EDF=",round(summtable[2,1],2),", p=",signif(summtable[2,4],3))
          mygg2<-plot(sm(allviz,2))+
            geom_hline(yintercept=0)+
            l_ciPoly(fill="red",alpha=0.3)+
            l_fitLine(colour="red",size=0.75)+
            l_rug(mapping = aes(x=x),alpha=0.8)+theme_bw()+
            ylab(thisylab)+ggtitle(thistitle2)+
            theme(plot.title = element_text(hjust=0.5,size=15),
                  axis.title=element_text(size=12),
                  axis.text=element_text(size=11))
          thistitle<-paste0("EDF=",round(summtable[3,1],2),", p=",signif(summtable[3,4],3))
          mygg3<-plot(sm(allviz,3))+
            l_fitRaster() + l_fitContour() + l_points() + l_rug() + theme_bw() +
            ggtitle(thistitle3)+
            theme(plot.title = element_text(hjust=0.5,size=15),
                  legend.position="bottom",
                  axis.title=element_text(size=12),
                  axis.text=element_text(size=11))
          gg1[[ggnum]]<-list(mygg1,mygg2,mygg3)
          names(gg1)[ggnum]<-paste0(outnames[i],"~",prednames[j])
        }
      }
    }
    
  } else {
    gg1<-list()
    for(i in 1:length(outnames)){
      for(j in 1:length(prednames)){
        ggnum<-(1+(j-1))+((i-1)*length(prednames))
        summtable<-summary(lmlist[[ggnum]])$s.table
        thistitle<-paste0("EDF=",round(summtable[1,1],2),", p=",signif(summtable[1,4],3))
        thisylab<-lmlist[[ggnum]]$formula[[2]]
        
        allviz<-getViz(lmlist[[ggnum]])
        if(plotresid){
          gg1[[ggnum]]<-plot(sm(allviz,1))+
            geom_hline(yintercept=0)+
            l_ciPoly(fill="blue",alpha=0.3)+
            l_fitLine(colour="blue",size=0.75)+
            l_rug(mapping = aes(x=x, y=y),alpha=0.8)+
            l_points(shape=19,size=1,alpha=0.1)+theme_classic()+
            ylab(thisylab)+ggtitle(thistitle)+
            theme(plot.title = element_text(hjust=0.5,size=15),
                  axis.title=element_text(size=12),
                  axis.text=element_text(size=11))
        } else {
          gg1[[ggnum]]<-plot(sm(allviz,1))+
            geom_hline(yintercept=0)+
            l_ciPoly(fill="blue",alpha=0.3)+
            l_fitLine(colour="blue",size=0.75)+
            l_rug(mapping = aes(x=x),alpha=0.8)+theme_bw()+
            ylab(thisylab)+ggtitle(thistitle)+
            theme(plot.title = element_text(hjust=0.5,size=15),
                  axis.title=element_text(size=12),
                  axis.text=element_text(size=11))
        }
      } #end of plot j loop
    }#end of plot i loop
  }
  
  if(!is.null(ixterm)&&class(Data[,ixterm])=="factor"){
    names(lmlist2)<-names(lmlist)
    retlist<-list(Resmatout,gg1,lmlist,lmlist2)
    names(retlist)<-c("Matrix","Plots","GAMlist","IxGAMlist")
  } else {
    retlist<-list(Resmatout,gg1,lmlist)
    names(retlist)<-c("Matrix","Plots","GAMlist")
  }
  return(retlist)
}
