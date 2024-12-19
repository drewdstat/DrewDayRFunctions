#' Extract all binary contrast coefficients of interest from a regression model 
#' 
#' \code{ContrastCoefficients} is a function designed to extract all the binary 
#' combinations of contrasts for a categorical (factor) predictor of interest 
#' from a linear model. The typical output of linear model only returns the 
#' coefficients for categorical variables in the context of a single referent 
#' level. For example, if a predictor variable of interest is a 4-level 
#' categorical treatment variable with factor levels 
#' \code{c("placebo", "drug1", "drug2", "drug3")}, a linear model will only 
#' return the "drug1 - placebo", "drug2 - placebo", and "drug3 - placebo" binary 
#' contrast coefficients. \code{ContrastCoefficients} will return all binary 
#' contrast coefficients for all possible referent levels, meaning the 
#' aforementioned coefficients will be returned, as well the coefficients for 
#' "drug2 - drug1", "drug 3 - drug 1", and "drug 3 - drug 2" contrasts. 
#' Furthermore, \code{ContrastCoefficients} can also provide all contrasts for 
#' an interaction term as well. This works for when both the predictor of 
#' interest and the interacting variable are categorical or when either the 
#' predictor of interest or the interacting variable are categorical. If only a 
#' continuous predictor of interest or both a continuous predictor and 
#' interacting variable are specified, the function will simply return the 
#' coefficients of interest.
#' 
#' The way this function obtains these additional coefficients is by iteratively 
#' updating the referent level for the predictor of interest and, if specified, 
#' the interacting variable and rerunning the linear models. Other methods may 
#' instead try to calculate these coefficients from the initial model 
#' coefficients, but the iterative method has the benefit of uniformly applying 
#' the same error estimation procedures for each coefficient. Some helpful 
#' internal functions used in this functions include 
#' \code{DrewDayRFunctions:::getcombos} to find all n combinations within one or 
#' between two categorical variables and 
#' \code{DrewDayRFunctions:::getmodelcoefs}, which extracts all model 
#' coefficients within a given model (i.e., without providing all contrast 
#' combinations) for the specified predictor and/or interacting variable of 
#' interest for the specified model. 
#' 
#' This function is currently set up to support \code{\link[stats]{lm}}, 
#' \code{\link[stats]{glm}}, \code{\link[MASS]{glm.nb}}, 
#' and \code{\link[logistf]{logistf}} (i.e., Firth's 
#' logistic regression) model objects, as well as a \code{list} 
#' or \code{\link[mice]{mira}} class list of models produced by the 
#' multiple imputation function \code{\link[mice]{mice}}, wherein each model is 
#' either a \code{lm}, \code{glm}, or \code{glm.nb} model.
#' 
#' @param mod This is a model object of class \code{\link[stats]{lm}}, 
#' \code{\link[stats]{glm}}, \code{negbin} (see \code{\link[MASS]{glm.nb}}), 
#' \code{\link[logistf]{logistf}}, or a \code{list} 
#' or \code{\link[mice]{mira}} class list of models, where in each model is a 
#' \code{lm}, \code{glm}, or \code{negbin} object. If \code{mod} is a 
#' \code{list} or \code{mira} object, estimates will be pooled using the 
#' \code{\link[mice]{pool}} function.
#' @param Data This is a data frame used to run the model \code{mod}. This 
#' should be of class \code{data.frame} unless \code{mod} is of class 
#' \code{list} or \code{mira}, in which case \code{Data} should be a list of 
#' datasets of class \code{list} or \code{\link[mice]{mids}}.
#' @param pred This is a character value defining the column name of the 
#' predictor of interest in the model \code{mod} and the data frame \code{Data}. 
#' @param ixterm This is an optional character value defining the column name 
#' of the interacting variable included in the model \code{mod} and present in 
#' the data frame \code{Data}. This defaults to \code{NULL}, meaning that no 
#' interacting variable contrasts will be output. It is necessary for there to 
#' be an interaction with the variable \code{ixterm} in the model \code{mod} for 
#' the function to work. Also, the \code{pred} variable should be on the left 
#' side of the interaction term while \code{ixterm} is on the right side (i.e., 
#' \code{pred:ixterm} and not \code{ixterm:pred}).
#' @param robust If \code{TRUE}, this calculates and returns robust CIs and 
#' p-values based on the \code{\link[sandwich]{vcovHC}} function instead of the 
#' default CIs and p-values. This defaults to \code{TRUE}, and is only applied 
#' if \code{mod} is of class \code{'lm'} or \code{'glm'}. 
#' @param HCtype This is the formula used for the \code{\link[sandwich]{vcovHC}} 
#' function calculation of robust CIs and p-values. Defaults to \code{"HC0"}.
#' @param usedf If \code{FALSE}, this sets \code{\link[lmtest]{coefci}} to 
#' calculate z-score-based robust confidence intervals if \code{robust} is 
#' \code{TRUE}. If \code{TRUE}, t test-based robust confidence intervals are 
#' instead calculated based on the residual degrees of freedom in each 
#' regression if \code{robust} is \code{TRUE}. This defaults to \code{FALSE}, 
#' which is the default for the \code{df} argument in 
#' \code{\link[lmtest]{coefci}}.
#' @param savecontrastmodels If \code{TRUE}, this will output the additional 
#' iterated models (beyond the original model \code{mod}) as a list. This 
#' defaults to \code{FALSE}.
#' 
#' @return If \code{savecontrastmodels} is \code{FALSE}, 
#' \code{ContrastCoefficients} will simply return a matrix of the coefficients, 
#' their error and confidence interval estimates, and their p-values. Otherwise, 
#' \code{ContrastCoefficients} will return a list of the following items:
#' \item{Coefficients}{This is the matrix of coefficients and related values 
#' returned by this function regardless.}
#' \item{ContrastModels}{This is the list of iterated model objects. If 
#' \code{mod} was a \code{list} or \code{mira} object, this will be a list of 
#' lists of model objects.}
#' 
#' @import logistf MASS sandwich
#' @export ContrastCoefficients
#' 
#' @details
#' \code{ContrastCoefficients} returns model coefficients for all unique binary 
#' combinations of categorical predictor variables when there is no interaction 
#' term specified in the model and the predictor variable of interest is 
#' categorical. When there is an interaction term, \code{ContrastCoefficients} 
#' returns all marginal coefficients for all unique binary contrast 
#' coefficients at each referent level of the interaction term are reported as 
#' "marginal coefficients". 
#' 
#' As an example, we will define Model 1 with a 
#' continuous predictor "age" and a two-level interaction term "sex" and Model 2 
#' with a 3-level categorical predictor "treatment" and the same two-level 
#' categorical interaction term "sex". If \code{allmarginalcontrasts = TRUE}, 
#' then the Model 1 coefficients will include age given sex = male (
#' "age|sex = male"), "age|sex = female", the interacting variable coefficient 
#' "sex (female - male)", and the interaction term coefficient 
#' "age:sex (female - male)". For Model 2, the coefficients would include 
#' 1. all binary treatment contrasts at each sex 
#' referent level (i.e., "treatment (2 - 1)|sex = male", 
#' "treatment (2 - 1)|sex = female", "treatment (3 - 1)|sex = male", 
#' "treatment (3 - 1)|sex = female", etc.), 2. the sex 
#' coefficient "sex (female - male)", and 3. all unique binary contrast 
#' interaction term coefficients (i.e., "treatment (2 - 1):sex (female - male)", 
#' "treatment (3 - 1):sex (female - male)", & 
#' "treatment (3 - 2):sex (female - male)"). 
#' 
#' Note that this function does not support greater than two-way interactions
#' (e.g., three-way interactions or higher), and it also does not support 
#' interactions that do not include the main effects of both the predictor of 
#' interest and the interacting variable (i.e., formulae with 
#' \code{pred * ixterm} or \code{pred + ixterm + pred:ixterm} are okay, but 
#' not those with \code{pred + pred:ixterm} or \code{ixterm + pred:ixterm}).
#' 
#' @examples
#' # Note that these examples with the infert dataset do not reflect any real 
#' # hypotheses and are instead just to demonstrate the function with example 
#' # continuous and categorical variables.
#' 
#' data("infert")
#' set.seed(5)
#' infert$fac4 <- as.factor(sample(1:4, nrow(infert), replace = T))
#' 
#' #A categorical predictor without any interaction terms
#' mod.1 <- glm(case ~ fac4 + education + age, data = infert, family = binomial)
#' ContrastCoefficients(mod.1, infert, "fac4")
#' 
#' #A categorical predictor with a categorical interaction term
#' mod.1 <- glm(case ~ fac4*education + age, data = infert, family = binomial)
#' ContrastCoefficients(mod.1, infert, "fac4", "education")
#' 
#' #A categorical predictor with a continuous interaction term
#' mod.1 <- glm(case ~ fac4*age + education, data = infert, family = binomial)
#' ContrastCoefficients(mod.1, infert, "fac4", "age")
#' 
#' #A continuous predictor with a categorical interaction term
#' mod.1 <- glm(case ~ age*fac4 + education, data = infert, family = binomial)
#' ContrastCoefficients(mod.1, infert, "age", "fac4")
#' 
#' #Missing data imputation with mice & a categorical by categorical interaction
#' infert$education_missing <- infert$education
#' set.seed(3)
#' infert[sample(1:nrow(infert), 50), "education_missing"] <- NA
#' infert.mice <- mice::mice(infert[, c(
#' "education_missing", "case", "fac4", "age")], maxit = 3, m = 3, seed = 20, 
#' printFlag = F)
#' # maxit and m should be higher, but this is just a fast-running example
#' mod.1 <- with(data = infert.mice, 
#' exp = glm(case ~ fac4*education_missing + age, family = binomial))
#' ContrastCoefficients(mod.1, infert.mice, "fac4", "education_missing")
#' 
ContrastCoefficients <- function(mod, Data, pred, ixterm = NULL, robust = T, 
                              HCtype = "HC0", usedf = F, savecontrastmodels = F){
  #check variable classes and test for pred/ixterm in mod
  if(!any(c("glm", "lm", "logistf", "mira", "list") %in% class(mod))) {
    stop(paste0("'mod' must be of class 'glm', 'lm', 'glm.nb', 'logistf', ",
                "'list', or 'mira'."))
  }
  if(any(c("list", "mira") %in% class(mod))){
    miceFlag <- T
    if(class(Data) == "mids"){
      DataList <- purrr::map(1:Data$m, function(x) mice::complete(Data, x))
    } else DataList <- Data
    Data <- DataList[[1]] #take the first data.frame to fit checks below
    if(any("mira" %in% class(mod))){
      origform <- formula(mod$analyses[[1]])
      origfam <- mod$analyses[[1]]$family
      modclasses <- sapply(mod$analyses, 
                           function(x) any(c("glm", "lm", "negbin") %in% class(x)))
    } else {
      origform <- mod[[1]]$formula
      origfam <- mod[[1]]$family
      modclasses <- sapply(mod, function(x) any(c(
        "glm", "lm") %in% class(x)))
    }
    if(!all(modclasses)){
      stop(paste0("If mod is a list of models or a mira object, each model ",
                  "within those lists must be of class 'glm', 'lm', or 'negbin'."))
    }
  } else miceFlag <- F
  
  tmpcoef <- tryCatch(getmodelcoefs(
    mod, pred = pred, ixterm = ixterm, robust = robust, HCtype = HCtype, 
    usedf = usedf), error = function(e) {
      warning(e)
      NULL})
  if(is.null(tmpcoef)){
    stop(paste0("It appears that pred and/or ixterm was misspecified. ", 
                "Please make sure that either pred (if ixterm is NULL) or both",
                " pred and ixterm are both in the model mod and that there is",
                " an interaction term with both pred and ixterm included."))
    tstnm <- NULL
  } else tstnm <- rownames(tmpcoef)[grepl(pred, rownames(tmpcoef))][1]
  if(is.null(tstnm) || length(tstnm) == 0){
    stop(paste0("You specified a pred that does not appear in the model mod.",
                " Please either provide a model with that pred variable or ",
                "specify a different pred that is in mod."))
  }
  if(class(Data[, pred]) == "character") Data[, pred] <- 
      as.factor(Data[, pred])

  if(!is.null(ixterm)){
    if(!any(grepl(":", rownames(tmpcoef)))){
      stop(paste0("An ixterm was specified, but there does not appear to be ", 
                  "any interactions with ixterm in the model mod. Please ",
                  "either provide a model with interaction terms with ixterm ",
                  "or set ixterm to NULL (its default)."))
    }
    tstnm <- rownames(tmpcoef)[grepl(pred, rownames(tmpcoef)) & 
                                 grepl(ixterm, rownames(tmpcoef))][1]
    if(length(tstnm) == 0 | is.na(tstnm)){
      stop(paste0("An ixterm was specified that does not ", 
                  "appear to be either present at all in the model mod or is ",
                  "not part of an interaction term. Please either provide a ", 
                  "model with the variable ixterm being part of one of its ",
                  "interaction terms or specify a different ixterm."))
    }
    tst <- strsplit(tstnm, split = ":")[[1]]
    if(grepl(pred, tst[2]) | grepl(ixterm, tst[1])){
      stop(paste0("If specifying an interaction term, one must make sure that ", 
                  "pred is on the left side of the interaction term and ixterm", 
                  " is on the right side (i.e., pred:ixterm and not ",
                  "ixterm:pred). Please make sure the order of pred and ixterm", 
                  " are correct."))
    }
    if(class(Data[, ixterm]) == "character") Data[, ixterm] <- 
      as.factor(Data[, ixterm])
  }
  
  if(!is.null(ixterm)){
    varclasses <- sapply(Data[, c(pred, ixterm)], class)
  } else varclasses <- class(Data[, pred])
  
  #Get combination labels
  changecombolabel <- function(x){
    tmp <- strsplit(x, split = "|", fixed = T)
    sapply(tmp, function(x) paste(x[2], x[1], sep = " - "))
  }
  if(is.factor(Data[, pred])){
    predcombos <- changecombolabel(getcombos(Data[, pred])[[1]])
    predlvls <- levels(Data[, pred])
  } else predcombos <- pred
  if(!is.null(ixterm) && class(Data[, ixterm]) == "factor"){
    ixcombos <- changecombolabel(getcombos(Data[, ixterm])[[1]])
    ixlvls <- levels(Data[, ixterm])
    predixcombos <- getcombos(Data[, pred], Data[, ixterm])$xycomb
    if(is.factor(Data[, pred])){
      predixcombos[, ] <- lapply(predixcombos, 
                                 function(x) changecombolabel(as.character(x)))
      predixcombos[, 1] <- paste0(paste(pred, predixcombos[, 1], sep = " ("), ")")
      predixcombos[, 2] <- paste0(paste(ixterm, predixcombos[, 2], sep = " ("), ")")
      predixcombos <- paste(predixcombos[, 1], predixcombos[, 2], sep = ":")
    } else {
      predixcombos <- paste0(paste(ixterm, changecombolabel(predixcombos), 
                                   sep = " ("), ")")
      predixcombos <- paste(pred, predixcombos, sep = ":")
    }
  }
  
  #set initial model list and coefficient table
  if(savecontrastmodels) contrastmodels <- list()
  modcoef <- getmodelcoefs(mod, pred = pred, ixterm = ixterm, robust = robust, 
                           HCtype = HCtype, usedf = usedf)
  
  #begin conditional loops
  if(length(varclasses) == 2 && all(varclasses == "factor")){ 
    #get rownames for modcoef
    rownames(modcoef)[which(grepl(pred, rownames(modcoef)) & 
                              grepl(ixterm, rownames(modcoef)))] <- 
      predixcombos[which(grepl(paste(escapespecialchars(
        predcombos[grep(paste0(levels(Data[, pred])[1], "$"), predcombos)]), 
        collapse = "|"), predixcombos) & 
          grepl(paste(escapespecialchars(ixcombos[grep(paste0(levels(Data[
            , ixterm])[1], "$"), ixcombos)]), collapse = "|"), predixcombos))]
    rownames(modcoef)[which(grepl(pred, rownames(modcoef)) & 
                              !grepl(paste0(ixterm, "|\\:"), rownames(modcoef)))] <- 
      paste0(paste0(paste(pred, predcombos[grep(paste0(escapespecialchars(levels(Data[
        , pred])[1]), "$"), predcombos)], sep = " ("), ")"), "|",
        ixterm, " = ", levels(Data[, ixterm])[1])
    rownames(modcoef)[which(!grepl(pred, rownames(modcoef)) & 
                              grepl(ixterm, rownames(modcoef)))] <- 
      paste0(paste(ixterm, ixcombos[grep(paste0(escapespecialchars(
        levels(Data[, ixterm])[1]), "$"), ixcombos)], 
        sep = " ("), ")")
    
    #set up new data and model with new referent levels for pred and ixterm
    for(il in 1:length(ixlvls)){
      #predterm releveling to original levels
      if(il > 2){
        #this allows cycling of ixterm referents to be saved but to reset the 
        # cycling of the pred referents
        if(miceFlag){
          tmpdat <- newData
          for(dl in 1:length(tmpdat)){
            tmpdat[[dl]][, pred] <- factor(tmpdat[[dl]][, pred], 
                                           levels = levels(Data[, pred]))
          }; rm(dl)
          newData <- tmpdat
        } else {
          tmpdat <- newData
          tmpdat[, pred] <- factor(tmpdat[, pred], levels = levels(Data[, pred]))
          newData <- tmpdat
        }
      } else {
        if(miceFlag) newData <- DataList else newData <- Data
      }
      #ixterm releveling
      if(il > 1){
        if(miceFlag){
          for(dl in 1:length(newData)){
            newData[[dl]][, ixterm] <- factor(
              newData[[dl]][, ixterm], levels = levels(newData[[dl]][, ixterm])[
                c(2:(length(levels(newData[[dl]][, ixterm]))), 1)])
          }; rm(dl)
        } else {
          newData[, ixterm] <- factor(
            newData[, ixterm], levels = levels(newData[, ixterm])[
              c(2:(length(levels(newData[, ixterm]))), 1)])
        }
        getmaineff <- T
      } else getmaineff <- F
      
      if(length(predlvls) > 2){
        if(il == 1) plcycle <- 1:(length(predlvls) - 2) else 
          plcycle <- 1:(length(predlvls) - 1)
        
        for(pl in plcycle){
          if(!(il > 1 & pl == 1)){
            #pred releveling
            if(miceFlag){
              for(dl in 1:length(newData)){
                newData[[dl]][, pred] <- factor(
                  newData[[dl]][, pred], levels = levels(newData[[dl]][, pred])[
                    c(2:(length(levels(newData[[dl]][, pred]))), 1)])
              }; rm(dl)
            } else {
              newData[, pred] <- factor(
                newData[, pred], levels = levels(newData[, pred])[
                  c(2:(length(levels(newData[, pred]))), 1)])
            }
          }
          if(miceFlag){
            newpredlvls <- levels(newData[[1]][, pred])
            newixlvls <- levels(newData[[1]][, ixterm])
          } else {
            newpredlvls <- levels(newData[, pred])
            newixlvls <- levels(newData[, ixterm])
          }
          
          #run new model
          newmod <- getnewmod(mod, newData, miceFlag)
          
          if(savecontrastmodels){ #save temporary model if requested
            contrastmodels[[paste0(pred, " ref=", newpredlvls[
              1], ", ", ixterm, " ref=", newixlvls[1])]] <- newmod
          }
          ncoef2take <- (length(predlvls) - 1) - pl 
          #To do: Change above to simpler length(predlvls) - pl so always get all marginals
          nixcoef2take <- (length(ixlvls)) - il
          if(il > 1) ncoef2take <- ncoef2take + 1
          if(ncoef2take > 0){
            tmp <- getmodelcoefs(newmod, pred = pred, ixterm = ixterm, 
                                 robust = robust, HCtype = HCtype, usedf = usedf)
            
            tmp1 <- tmp[which(grepl(pred, rownames(tmp)) &
                                !grepl(ixterm, rownames(tmp))), ]
            if(!is.matrix(tmp1)){ #set back to matrix if it's just 1 row
              tmp1 <- matrix(tmp1, nrow = 1)
              colnames(tmp1) <- colnames(tmp)
            }
            tmp1 <- tmp1[1:ncoef2take, ]
            if(!is.matrix(tmp1)){
              tmp1 <- matrix(tmp1, nrow = 1)
              colnames(tmp1) <- colnames(tmp)
            }
            rownames(tmp1) <- paste0(
              paste0(paste(pred, predcombos[
                grep(paste0(escapespecialchars(newpredlvls[1]), "$"),
                     predcombos)], sep = " ("), ")"), "|", ixterm, " = ",
              newixlvls[1])
            modcoef <- rbind(modcoef, tmp1)
            if(nixcoef2take > 0){
              tmp2 <- tmp[which(grepl(ixterm, rownames(tmp)) & 
                                  !grepl(pred, rownames(tmp))), ]
              if(!is.matrix(tmp2)){
                tmp2 <- matrix(tmp2, nrow = 1)
                colnames(tmp2) <- colnames(tmp)
              }
              tmp2 <- tmp2[1:nixcoef2take, ]
              if(!is.matrix(tmp2)){
                tmp2 <- matrix(tmp2, nrow = 1)
                colnames(tmp2) <- colnames(tmp)
              }
              rownames(tmp2) <- paste0(paste(ixterm, ixcombos[grep(paste0(
                escapespecialchars(newixlvls[1]), "$"), ixcombos)], 
                sep = " ("), ")")
              
              ixlvl2keep <- newixlvls[2:(length(ixlvls) - (il - 1))]
              tmp3 <- tmp[which(grepl(paste(escapespecialchars(rownames(tmp)[
                1:ncoef2take]), collapse = "|"), rownames(tmp)) & 
                  grepl(paste(paste0(ixterm, escapespecialchars(ixlvl2keep)), 
                              collapse = "|"), rownames(tmp))), ]
              if(!is.matrix(tmp3)){
                tmp3 <- matrix(tmp3, nrow = 1)
                colnames(tmp3) <- colnames(tmp)
              }
              rownames(tmp3) <- predixcombos[which(grepl(paste(escapespecialchars(
                predcombos[grep(paste0(newpredlvls[1], "$"), predcombos)]), 
                collapse = "|"), predixcombos) & 
                  grepl(paste(escapespecialchars(ixcombos[grep(paste0(newixlvls[
                    1], "$"), ixcombos)]), collapse = "|"), predixcombos))]
              if(getmaineff & pl == 1){
                modcoef <- rbind(modcoef, tmp2)
              }
              modcoef <- rbind(modcoef, tmp3)
            }
          }
        }; rm(pl)
      } else {
        #if predictor has 2 levels, there's no need to cycle pred levels
        if(il == 1) next else {
          if(miceFlag){
            newixlvls <- levels(newData[[1]][, ixterm])
          } else {
            newixlvls <- levels(newData[, ixterm])
          }
          #run new model
          newmod <- getnewmod(mod, newData, miceFlag)
          
          if(savecontrastmodels){ #save temporary model if requested
            contrastmodels[[paste0(pred, " ref=", levels(Data[, pred])[
              1], ", ", ixterm, " ref=", newixlvls[1])]] <- newmod
          }
          
          nixcoef2take <- (length(ixlvls)) - il
          
          tmp <- getmodelcoefs(newmod, pred = pred, ixterm = ixterm, 
                               robust = robust, HCtype = HCtype, usedf = usedf)
          
          tmp1 <- tmp[which(grepl(pred, rownames(tmp)) &
                              !grepl(ixterm, rownames(tmp))), ]
          if(!is.matrix(tmp1)){ #set back to matrix if it's just 1 row
            tmp1 <- matrix(tmp1, nrow = 1)
            colnames(tmp1) <- colnames(tmp)
          }
          rownames(tmp1) <- paste0(
            paste0(paste(pred, predcombos[
              grep(paste0(escapespecialchars(levels(Data[, pred])[1]), "$"),
                   predcombos)], sep = " ("), ")"), "|", ixterm, " = ",
            newixlvls[1])
          modcoef <- rbind(modcoef, tmp1)
          if(nixcoef2take > 0){
            tmp2 <- tmp[which(grepl(ixterm, rownames(tmp)) & 
                                !grepl(pred, rownames(tmp))), ]
            if(!is.matrix(tmp2)){
              tmp2 <- matrix(tmp2, nrow = 1)
              colnames(tmp2) <- colnames(tmp)
            }
            tmp2 <- tmp2[1:nixcoef2take, ]
            if(!is.matrix(tmp2)){
              tmp2 <- matrix(tmp2, nrow = 1)
              colnames(tmp2) <- colnames(tmp)
            }
            rownames(tmp2) <- paste0(paste(ixterm, ixcombos[grep(paste0(
              escapespecialchars(newixlvls[1]), "$"), ixcombos)], 
              sep = " ("), ")")
            
            ixlvl2keep <- newixlvls[2:(length(ixlvls) - (il - 1))]
            tmp3 <- tmp[which(grepl(paste(escapespecialchars(rownames(tmp)[
              1:ncoef2take]), collapse = "|"), rownames(tmp)) & 
                grepl(paste(paste0(ixterm, escapespecialchars(ixlvl2keep)), 
                            collapse = "|"), rownames(tmp))), ]
            if(!is.matrix(tmp3)){
              tmp3 <- matrix(tmp3, nrow = 1)
              colnames(tmp3) <- colnames(tmp)
            }
            rownames(tmp3) <- predixcombos[which(grepl(paste(escapespecialchars(
              predcombos[grep(paste0(newpredlvls[1], "$"), predcombos)]), 
              collapse = "|"), predixcombos) & 
                grepl(paste(escapespecialchars(ixcombos[grep(paste0(newixlvls[
                  1], "$"), ixcombos)]), collapse = "|"), predixcombos))]
            if(getmaineff & pl == 1){
              modcoef <- rbind(modcoef, tmp2)
            }
            modcoef <- rbind(modcoef, tmp3)
          }
        }
      } #end else if length(levels(Data[, pred]) == 2)
    }; suppressWarnings(rm(il, tmp, tmp1, tmp2, tmp3, newmod))
    
    #reorder
    predlabs <- paste0(paste(pred, predcombos, sep = " ("), ")")
    margixlabs <- paste0(paste(ixterm, "=", ixlvls))
    marglabs <- c()
    for(ml in margixlabs){
      tl <- paste0(predlabs, "|", ml)
      marglabs <- c(marglabs, tl)
    }; rm(ml, tl)
    modcoef <- modcoef[c(marglabs, 
                         paste0(paste(ixterm, ixcombos, sep = " ("), ")"), 
                         predixcombos), ]
  } else if(length(varclasses) == 2 && varclasses[2] == "factor" & 
            varclasses[1] != "factor"){
    #make rownames for modcoef
    rownames(modcoef)[which(grepl(pred, rownames(modcoef)) & 
                              grepl(ixterm, rownames(modcoef)))] <- 
      predixcombos[which(grepl(paste(escapespecialchars(
        predcombos[grep(paste0(levels(Data[, pred])[1], "$"), predcombos)]), 
        collapse = "|"), predixcombos) & 
          grepl(paste(escapespecialchars(ixcombos[grep(paste0(levels(Data[
            , ixterm])[1], "$"), ixcombos)]), collapse = "|"), predixcombos))]
    rownames(modcoef)[which(grepl(pred, rownames(modcoef)) & 
                              !grepl(ixterm, rownames(modcoef)))] <- 
      paste0(pred, "|", ixterm, " = ", levels(Data[, ixterm])[1])
    rownames(modcoef)[which(!grepl(pred, rownames(modcoef)) & 
                              grepl(ixterm, rownames(modcoef)))] <- 
      paste0(paste(ixterm, ixcombos[grep(paste0(escapespecialchars(
        levels(Data[, ixterm])[1]), "$"), ixcombos)], 
        sep = " ("), ")")
    
    #set up new data and model with new referent levels for ixterm
    if(miceFlag) newData <- DataList else newData <- Data
    for(il in 1:length(ixlvls)){ 
      if(miceFlag){
        for(dl in 1:length(newData)){
          newData[[dl]][, ixterm] <- factor(
            newData[[dl]][, ixterm], levels = levels(newData[[dl]][, ixterm])[
              c(2:(length(levels(newData[[dl]][, ixterm]))), 1)])
        }; rm(dl)
        newixlvls <- levels(newData[[1]][, ixterm])
      } else {
        newData[, ixterm] <- factor(
          newData[, ixterm], levels = levels(newData[, ixterm])[
            c(2:(length(levels(newData[, ixterm]))), 1)])
        newixlvls <- levels(newData[, ixterm])
      }
      
      #run new model
      newmod <- getnewmod(mod, newData, miceFlag)
      
      if(savecontrastmodels){ #save temporary model if requested
        contrastmodels[[paste0(ixterm, " ref=", 
                               newixlvls[1])]] <- newmod
      }
      nixcoef2take <- (length(ixlvls) - 1) - il
      if((nixcoef2take + 1) > 0){
        tmp <- getmodelcoefs(newmod, pred = pred, ixterm = ixterm, 
                             robust = robust, HCtype = HCtype, usedf = usedf)
        tmp1 <- tmp[which(grepl(pred, rownames(tmp)) & 
                            !grepl(ixterm, rownames(tmp))), ]
        if(!is.matrix(tmp1)){
          tmp1 <- matrix(tmp1, nrow = 1)
          colnames(tmp1) <- colnames(tmp)
        }
        rownames(tmp1) <- paste0(pred, "|", ixterm, " = ",
                                 levels(newData[, ixterm])[1])
        if(nixcoef2take > 0){
          tmp2 <- tmp[which(grepl(ixterm, rownames(tmp)) & 
                              !grepl(pred, rownames(tmp))), ]
          if(!is.matrix(tmp2)){
            tmp2 <- matrix(tmp2, nrow = 1)
            colnames(tmp2) <- colnames(tmp)
          }
          tmp2 <- tmp2[1:nixcoef2take, ]
          if(!is.matrix(tmp2)){
            tmp2 <- matrix(tmp2, nrow = 1)
            colnames(tmp2) <- colnames(tmp)
          }
          rownames(tmp2) <- paste0(paste(ixterm, ixcombos[grep(paste0(
            escapespecialchars(newixlvls[1]), "$"), ixcombos)], 
            sep = " ("), ")")
          
          ixlvl2keep <- newixlvls[2:(length(ixlvls) - il)]
          tmp3 <- tmp[which(grepl(escapespecialchars(pred), rownames(tmp)) & 
                              grepl(paste(paste0(ixterm, escapespecialchars(ixlvl2keep)), 
                                          collapse = "|"), rownames(tmp))), ]
          if(!is.matrix(tmp3)){
            tmp3 <- matrix(tmp3, nrow = 1)
            colnames(tmp3) <- colnames(tmp)
          }
          rownames(tmp3) <- predixcombos[grep(paste(escapespecialchars(
            ixcombos[grep(paste0(newixlvls[1], "$"), ixcombos)]), 
            collapse = "|"), predixcombos)]
          modcoef <- rbind(modcoef, tmp1, tmp2, tmp3)
        } else modcoef <- rbind(modcoef, tmp1)
      }
    }; suppressWarnings(rm(il, tmp, tmp2, tmp3))
    #reorder
    modcoef <- modcoef[c(paste0(predcombos, "|", ixterm, " = ", levels(Data[
      , ixterm])), paste0(paste(ixterm, ixcombos, sep = " ("), ")"), 
      predixcombos), ]
  } else if(length(varclasses) > 1 && varclasses[2] != "factor" & 
            varclasses[1] == "factor"){
    rownames(modcoef)[which(grepl(pred, rownames(modcoef)) & 
                              !grepl(ixterm, rownames(modcoef)))] <- 
      paste0(paste(pred, predcombos[grep(paste0(escapespecialchars(levels(Data[
        , pred])[1]), "$"), predcombos)], sep = " ("), ")")
    rownames(modcoef)[which(grepl(pred, rownames(modcoef)) & 
                              grepl(ixterm, rownames(modcoef)))] <- 
      paste(paste0(paste(pred, predcombos[grep(paste0(escapespecialchars(levels(Data[
        , pred])[1]), "$"), predcombos)], sep = " ("), ")"), ixterm, sep = ":")
    
    if(length(levels(Data[, pred])) > 2){
      if(miceFlag) newData <- DataList else newData <- Data
      
      for(pl in 1:(length(predlvls) - 2)){
        if(miceFlag){
          for(dl in 1:length(newData)){
            newData[[dl]][, pred] <- factor(
              newData[[dl]][, pred], levels = levels(newData[[dl]][, pred])[
                c(2:(length(levels(newData[[dl]][, pred]))), 1)])
          }; rm(dl)
          newpredlvls <- levels(newData[[1]][, pred])
        } else {
          newData[, pred] <- factor(
            newData[, pred], levels = levels(newData[, pred])[
              c(2:(length(levels(newData[, pred]))), 1)])
          newpredlvls <- levels(newData[, pred])
        }
        
        newmod <- getnewmod(mod, newData, miceFlag)
        
        if(savecontrastmodels){
          contrastmodels[[paste0(pred, " ref=", 
                                 newpredlvls[1])]] <- newmod
        }
        ncoef2take <- (length(predlvls) - 1) - pl
        if(ncoef2take > 0){
          tmp <- getmodelcoefs(newmod, pred = pred, ixterm = ixterm, 
                               robust = robust, HCtype = HCtype, usedf = usedf) 
          tmp1 <- tmp[which(grepl(pred, rownames(tmp)) & 
                              !grepl(ixterm, rownames(tmp))), ]
          if(!is.matrix(tmp1)){
            tmp1 <- matrix(tmp1, nrow = 1)
            colnames(tmp1) <- colnames(tmp)
          }
          tmp1 <- tmp1[1:ncoef2take, ]
          if(!is.matrix(tmp1)){
            tmp1 <- matrix(tmp1, nrow = 1)
            colnames(tmp1) <- colnames(tmp)
          }
          rownames(tmp1) <- paste0(paste(pred, predcombos[
            grep(paste0(escapespecialchars(newpredlvls[1]), "$"), 
                 predcombos)], sep = " ("), ")")
          tmp3 <- tmp[which(grepl(paste(escapespecialchars(rownames(tmp)[
            1:ncoef2take]), collapse = "|"), rownames(tmp)) & 
              grepl(paste0(":", escapespecialchars(ixterm)), 
                    rownames(tmp))), ]
          if(!is.matrix(tmp3)){
            tmp3 <- matrix(tmp3, nrow = 1)
            colnames(tmp3) <- colnames(tmp)
          }
          rownames(tmp3) <- paste(rownames(tmp1), ixterm, sep = ":")
          modcoef <- rbind(modcoef, tmp1, tmp3)
        }
      }; rm(pl)
      modcoef <- modcoef[c(paste0(paste(pred, predcombos, sep = " ("), ")"), 
                           ixterm, 
                           paste(paste0(paste(pred, predcombos, sep = " ("), ")"), 
                                 ixterm, sep = ":")), ]
    } 
  } else if(class(Data[, pred]) == "factor" & is.null(ixterm)){
    rownames(modcoef) <- 
      paste0(paste(pred, predcombos[grep(paste0(escapespecialchars(levels(Data[
        , pred])[1]), "$"), predcombos)], sep = " ("), ")")
    if(miceFlag) newData <- DataList else newData <- Data
    for(pl in 1:(length(predlvls) - 2)){
      if(miceFlag){
        for(dl in 1:length(newData)){
          newData[[dl]][, pred] <- factor(
            newData[[dl]][, pred], levels = levels(newData[[dl]][, pred])[
              c(2:(length(levels(newData[[dl]][, pred]))), 1)])
        }; rm(dl)
        newpredlvls <- levels(newData[[1]][, pred])
      } else {
        newData[, pred] <- factor(
          newData[, pred], levels = levels(newData[, pred])[
            c(2:(length(levels(newData[, pred]))), 1)])
        newpredlvls <- levels(newData[, pred])
      }
      
      #run new model
      newmod <- getnewmod(mod, newData, miceFlag)
      
      if(savecontrastmodels){
        contrastmodels[[paste0(pred, " ref=", 
                               newpredlvls[1])]] <- newmod
      }
      ncoef2take <- (length(predlvls) - 1) - pl
      if(ncoef2take > 0){
        tmp <- getmodelcoefs(newmod, pred = pred, ixterm = ixterm, 
                             robust = robust, HCtype = HCtype, usedf = usedf) 
        tmp1 <- tmp[grep(pred, rownames(tmp)), ]
        if(!is.matrix(tmp1)){
          tmp1 <- matrix(tmp1, nrow = 1)
          colnames(tmp1) <- colnames(tmp)
        }
        tmp1 <- tmp[1:ncoef2take, ]
        if(!is.matrix(tmp1)){
          tmp1 <- matrix(tmp1, nrow = 1)
          colnames(tmp1) <- colnames(tmp)
        }
        rownames(tmp1) <- paste0(paste(pred, predcombos[
          grep(paste0(escapespecialchars(newpredlvls[1]), "$"), 
               predcombos)], sep = " ("), ")")
        modcoef <- rbind(modcoef, tmp1)
      }
    }; rm(pl)
  } else {
    if(savecontrastmodels){
      if(is.null(ixterm)) modtitle <- pred else 
        modtitle <- paste0(pred, ", ", ixterm)
      contrastmodels[[modtitle]] <- newmod
    }
  } 
  #end conditional loops
  
  #output
  if(savecontrastmodels){
    return(list(Coefficients = modcoef, ContrastModels = contrastmodels))
  } else return(modcoef)
}
