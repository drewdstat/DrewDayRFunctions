#Get a character string of all counts and percentages (e.g., "10 (5%)")
countperc <- function(x, roundplace = 2, na.rm = T){
  if(na.rm) summ <- summary(x[which(!is.na(x))]) else summ <- summary(x)
  perc <- prop.table(summ)
  summ <- paste0(summ, " (", round(perc*100, roundplace), "%)")
  names(summ) <- names(perc)
  return(summ)
}

#Continuous variable summary statistics with the option of comparing by a 
# categorical variable with 2 or more levels with t-tests or ANOVAs
myvarcompsummary_cat_tab <- 
  function(catnames, Data, compvar = NULL, catvarnames = NULL, kable = T, 
           roundplace = 2, fontsize = 12, totalcol = F, nmissing = F, sigtest = T, 
           allcomplete = F,  scroll = T,  scrollwidth = "100%", 
           scrollheight = "200px"){
  library(rcompanion)
  
  if(!is.null(compvar)){
    if(!is.factor(Data[, compvar])) Data[, compvar] <- as.factor(Data[, compvar])
    ncomplevel <- length(levels(Data[, compvar]))
    if(allcomplete){
      responses <- lapply(Data[complete.cases(Data[, catnames]), catnames], 
                        function(x) dimnames(table(x, Data[complete.cases(Data[
                          , catnames]), compvar]))[[1]])
    } else {
      responses <- lapply(Data[, catnames], 
                        function(x) dimnames(table(x, Data[, compvar]))[[1]])
    }
    if(any(sapply(responses, length)!= 2)){
      newmat <- fastDummies::dummy_cols(Data[, c(catnames, compvar)], 
                                      select_columns = catnames[which(sapply(
                                        responses, length)!= 2)], 
                                      ignore_na = T)[, -which(sapply(
                                        responses, length)!= 2)] 
    } else {
      newmat <- Data[, c(catnames, compvar)]
    }
    newmat2 <- newmat[, -which(names(newmat) == compvar)]
    if(is.list(responses)){responses <- unlist(responses)}
    pwdf <- as.data.frame(t(sapply(newmat[, -which(names(newmat) == compvar)], 
                                 function(x) pairwiseNominalIndependence(
                                   table(x, newmat[, compvar]), fisher = F, 
                                   gtest = F, chisq = T, method = "none"))))
    summncol <- 3 + ncomplevel
    if(totalcol){summncol <- summncol + 1}
    summtab <- as.data.frame(matrix(NA, 2 * nrow(pwdf), summncol))
    summtab[, 1] <- rep(dimnames(pwdf)[[1]], each = 2)
    for(i in 1:ncol(newmat2)){
      tab1 <- table(newmat2[, i], newmat[, compvar])
      perctab1 <- t(prop.table(t(tab1), 1) * 100)
      tab1 <- matrix(paste0(tab1, " (", round(perctab1, roundplace), "%)"), 
                     nrow(tab1), ncol(tab1))
      dimnames(tab1) <- dimnames(perctab1)
      summtab[(1 + ((i-1) * 2)):(i * 2), 2] <- dimnames(tab1)[[1]]
      summtab[(1 + ((i-1) * 2)):(i * 2), 3:(2 + ncol(tab1))] <- tab1
      if(totalcol){
        tab2 <- table(newmat2[, i])
        perctab2 <- prop.table(t(tab2), 1)*100
        tab2 <- paste0(tab2, " (", round(perctab2, roundplace), "%)")
        summtab[(1 + ((i - 1) * 2)):(i * 2), summncol] <- tab2
      }
    }
    if(totalcol){
      summtab[, (summncol - 1)] <- rep(unlist(pwdf$p.Chisq), each = 2)
      stars <- ifelse(summtab[, (summncol-1)] < 0.05, "*", "")
      summtab[, (summncol - 1)] <- paste0(summtab[, (summncol-1)], stars)
      names(summtab) <- c("Variable", "Group", 
                        paste0(compvar, " = ", levels(Data[, compvar])), 
                        "Chisq p-value", "Total")
    } else {
      summtab[, summncol] <- rep(unlist(pwdf$p.Chisq), each = 2)
      stars <- ifelse(summtab[, summncol]<0.05, "*", "")
      summtab[, summncol] <- paste0(summtab[, summncol], stars)
      names(summtab) <- c("Variable", "Group", 
                        paste0(compvar, " = ", levels(Data[, compvar])), 
                        "Chisq p-value")
    }
    if(!sigtest) summtab <- summtab[, -grep("p-value", names(summtab))]
    if(kable){
      K1 <- kable(summtab, "html") %>%
        kable_styling(font_size = fontsize, full_width  =  F) %>%
        collapse_rows(columns = c(1, ncol(summtab)), valign = "top") %>%
        column_spec(1, bold = T)
      if(scroll) K1 <- K1 %>% scroll_box(width = scrollwidth, height = scrollheight)
      print(K1)
    } else {
      return(summtab)
    }
  } else {
    if(length(catnames) == 1){
      responses <- list()
      if(nmissing) responses[[1]] <- summary(as.factor(Data[
        which(!is.na(Data[, catnames])), catnames])) else 
          responses[[1]] <- summary(as.factor(Data[, catnames]))
    } else if(allcomplete){
      if(nmissing){
        responses <- c(lapply(Data[complete.cases(Data[, catnames]), catnames], 
                                     function(x) summary(as.factor(x))))
      } else {
        responses <- c(lapply(Data[complete.cases(Data[, catnames]), catnames], 
                                   function(x) summary(as.factor(x[which(!is.na(x))]))))
      }
    } else {
      if(nmissing){
        responses <- c(lapply(Data[, catnames], function(x) summary(as.factor(x))))
      } else {
        responses <- c(lapply(Data[, catnames], function(x) summary(as.factor(x[which(!is.na(x))]))))
      }
    }
    
    summtabcol <- 3
    summtab <- as.data.frame(matrix(NA, 0, summtabcol))
    if(!is.null(catvarnames)) catnames <- catvarnames
    for(ll in 1:length(responses)){
      temptab <- as.data.frame(matrix(NA, sapply(responses, length)[ll], 3))
      temptab[, 1] <- catnames[ll]
      temptab[, 2] <- names(responses[[ll]])
      temptab[, 3] <- responses[[ll]]
      temptab[, 4] <- 100 * (responses[[ll]]/sum(responses[[ll]]))
      temptab[, 3] <- paste0(temptab[, 3], " (", round(temptab[, 4], roundplace), "%)")
      temptab[, 4] <- NULL
      summtab <- rbind(summtab, temptab)
    };rm(ll, temptab)
    
    names(summtab) <- c("Variable", "Category", "N (%)")
    if(kable){
      K1 <- kable(summtab[, -1], "html") %>%
        kable_styling(font_size = fontsize, full_width  =  F) %>%
        pack_rows(index = table(forcats::fct_inorder(summtab[, 1]))) %>%
        column_spec(1, bold = T)
      
      if(scroll) K1 <- K1 %>% scroll_box(width = scrollwidth, height = scrollheight)
      print(K1)
    } else {
      return(summtab)
    }
  }
}

# Continuous variable summary statistics with a t-test or ANOVA if a 
# categorical comparison variable is provided
summtable_tt <- function(vars, Data, varnames = NULL, LODvars = NULL, 
                         compvar = NULL, newcomplvls = NULL, roundplace = 2, 
                         padj = F, padjmethod = "BH", nonpar = F, perm = F, 
                         signif = F){
  if(signif){
    roundfunc <- function(x, roundplace){
      if(!any(c("POSIXt", "Date") %in% class(x))) signif(x, roundplace) else x
    }
  } else {
    roundfunc <- function(x, roundplace){
      if(!any(c("POSIXt", "Date") %in% class(x))) round(x, roundplace) else x
    }
  }
  if(is.null(varnames)){
    varnames <- vars
  }
  if(!is.null(compvar)){
    if(!is.null(newcomplvls)){
      Data[, compvar] <- factor(as.character(Data[, compvar]), levels = newcomplvls)
    } else {
      if(class(Data[, compvar]) == "character"){
        Data[, compvar] <- as.factor(Data[, compvar])
      }
    }
    nlevels <- length(levels(Data[, compvar]))
    if(!is.null(LODvars)){
      summtable1 <- data.frame(matrix(NA, length(vars), 7))
      summtable2 <- data.frame(matrix(NA, length(vars)*nlevels, 7))
    } else {
      summtable1 <- data.frame(matrix(NA, length(vars), 6))
      summtable2 <- data.frame(matrix(NA, length(vars)*nlevels, 6))
    }
    for(i in 1:length(vars)){
      summtable1[i, 1] <- varnames[i]
      summtable1[i, 2] <- "Total"
      summtable1[i, 4] <- paste0(roundfunc(mean(Data[, vars[i]], na.rm = T), roundplace), " (", 
                              roundfunc(sd(Data[, vars[i]], na.rm = T), roundplace), ")")
      summtable1[i, 5] <- paste0(roundfunc(median(Data[, vars[i]], na.rm = T), roundplace), " (", 
                              roundfunc(min(Data[, vars[i]], na.rm = T), roundplace), " - ", 
                              roundfunc(max(Data[, vars[i]], na.rm = T), roundplace), ")")
      nobs <- length(which(is.na(Data[, vars[i]]) == F)) #nmiss
      summtable1[i, 3] <- nobs
      if(nonpar){
        if(nlevels == 2){
          tt <- tryCatch(wilcox.test(formula(paste0(vars[i], "~", compvar)), Data), 
                       error = function(e) NULL)
          if(is.null(tt)) ttp <- NA else ttp <- tt$p.value
        } else {
          tt <- tryCatch(kruskal.test(formula(paste0(vars[i], "~", compvar)), Data), 
                       error = function(e) NULL)
          if(is.null(tt)) ttp <- NA else ttp <- tt$p.value
        }
      } else if(perm){
        ttp <- perm.ttest(vars[i], NULL, compvar, Data)
      } else {
        if(nlevels == 2){
          tt <- tryCatch(t.test(formula(paste0(vars[i], "~", compvar)), Data), 
                       error = function(e) NULL)
          if(is.null(tt)) ttp <- NA else ttp <- tt$p.value
        } else {
          tt <- tryCatch(aov(formula(paste0(vars[i], "~", compvar)), Data), 
                       error = function(e) NULL)
          if(is.null(tt)) ttp <- NA else ttp <- summary(tt)[[1]][["Pr(>F)"]][1]
        }
      }
      
      summtable1[i, 6] <- signif(ttp, roundplace)
      if(!is.null(LODvars)){
        nLOD <- length(which(Data[, LODvars[i]] == 1))
        summtable1[i, 7] <- paste0(nLOD, " (", round((nLOD/nobs)*100, roundplace), "%)")
      }
      for(j in 1:nlevels){
        rownum <- (1+(j-1))+((i-1)*nlevels)
        thislevel <- levels(Data[, compvar])[j]
        summtable2[rownum, 1] <- varnames[i]
        summtable2[rownum, 2] <- thislevel
        summtable2[rownum, 4] <- paste0(roundfunc(mean(Data[Data[, compvar] == thislevel, vars[i]], na.rm = T), roundplace), " (", 
                                     roundfunc(sd(Data[Data[, compvar] == thislevel, vars[i]], na.rm = T), roundplace), ")")
        summtable2[rownum, 5] <- paste0(roundfunc(median(Data[Data[, compvar] == thislevel, vars[i]], na.rm = T), roundplace), " (", 
                                     roundfunc(min(Data[Data[, compvar] == thislevel, vars[i]], na.rm = T), roundplace), " - ", 
                                     roundfunc(max(Data[Data[, compvar] == thislevel, vars[i]], na.rm = T), roundplace), ")")
        nobs <- length(which(is.na(Data[Data[, compvar] == thislevel, vars[i]]) == F))#nmiss
        summtable2[rownum, 3] <- nobs
        if(nlevels > 2){
          if(nonpar){
            tt <- tryCatch(wilcox.test(Data[which(Data[, compvar] == levels(Data[, compvar])[j]), vars[i]], 
                                Data[which(Data[, compvar]!= levels(Data[, compvar])[j]), vars[i]]), 
                         error = function(e) NULL)
            if(is.null(tt)) ttp <- NA else ttp <- tt$p.value
          } else if(perm){
            ttp <- perm.ttest(Data[which(Data[, compvar] == levels(Data[, compvar])[j]), vars[i]], 
                           Data[which(Data[, compvar]!= levels(Data[, compvar])[j]), vars[i]])
          } else {
            tt <- tryCatch(t.test(Data[which(Data[, compvar] == levels(Data[, compvar])[j]), vars[i]], 
                                Data[which(Data[, compvar]!= levels(Data[, compvar])[j]), vars[i]]), 
                         error = function(e) NULL)
            if(is.null(tt)) ttp <- NA else ttp <- tt$p.value
          }
        }
        summtable2[rownum, 6] <- signif(ttp, roundplace)
        if(!is.null(LODvars)){
          nLOD2 <- length(which(Data[Data[, compvar] == thislevel, LODvars[i]] == 1))
          LODperc <- round((nLOD2/nrow(Data[which(Data[, compvar] == thislevel&
                                            !is.na(Data[, vars[i]])), ]))*100, roundplace)
          summtable2[rownum, 7] <- paste0(nLOD2, " (", LODperc, "%)")
        }
      }
    }
    summtable <- rbind(summtable1, summtable2)
    if(nlevels>2){testname <- "ANOVA/t-test"} else {testname <- "t-test"}
    if(!is.null(LODvars)){
      names(summtable) <- c("Variable", "Group", "N", "Mean (SD)", "Median (Range)", 
                          paste0(testname, " p-value"), "N <LOD (%)")
    } else {
      names(summtable) <- c("Variable", "Group", "N", "Mean (SD)", "Median (Range)", 
                          paste0(testname, " p-value"))
    }
    summtable$Variable <- factor(summtable$Variable, levels = varnames)
    summtable <- summtable[with(summtable, order(Variable)), ]
    summtable$Group <- factor(summtable$Group, levels = unique(summtable$Group))
    if(padj){
      if(nlevels == 2){
        pagg <- summtable[which(summtable$Group == "Total"), c("Variable", "t-test p-value")]
        pagg$padj <- stats::p.adjust(pagg[, ncol(pagg)], method = padjmethod)
        summtable <- merge(summtable, pagg[, c(1, ncol(pagg))], by = "Variable", all.x = T)
      } else {
        pagg <- summtable[, c(1:2, ncol(summtable))]
        pagg[which(pagg$Group!= "Total"), "padj"] <- 
          stats::p.adjust(pagg[which(pagg$Group!= "Total"), grep("p-value", names(pagg))], method = padjmethod)
        pagg[which(pagg$Group == "Total"), "padj"] <- 
          stats::p.adjust(pagg[which(pagg$Group == "Total"), grep("p-value", names(pagg))], method = padjmethod)
        summtable <- merge(summtable, pagg[, c(1:2, ncol(pagg))], by = c("Variable", "Group"), all.x = T)
      }
      summtable <- summtable[with(summtable, order(Variable, Group)), ]
      names(summtable)[ncol(summtable)] <- "Adjusted p-value"
      stars <- ifelse(summtable[, "Adjusted p-value"]<0.05, "**", 
                    ifelse(summtable[, paste0(testname, " p-value")]<0.05, "*", ""))
      summtable[, paste0(testname, " p-value")] <- 
        paste0(signif(summtable[, paste0(testname, " p-value")], roundplace), stars)
      summtable[, "Adjusted p-value"] <- paste0(signif(summtable[, "Adjusted p-value"], roundplace), stars)
      if(!is.null(LODvars)){summtable <- summtable[, c(1:3, 7, 4:6, 8)]}
    } else {
      stars <- ifelse(summtable[, 6]<0.05, "*", "")
      summtable[, 6] <- paste0(summtable[, 6], stars)
      if(!is.null(LODvars)){summtable <- summtable[, c(1:3, 7, 4:6)]}
    }
  } else {
    if(!is.null(LODvars)){
      summtable <- data.frame(matrix(NA, length(vars), 5))
    } else {
      summtable <- data.frame(matrix(NA, length(vars), 4))
    }
    for(i in 1:length(vars)){
      summtable[i, 1] <- varnames[i]
      summtable[i, 3] <- paste0(roundfunc(mean(Data[, vars[i]], na.rm = T), roundplace), " (", 
                             roundfunc(sd(Data[, vars[i]], na.rm = T), roundplace), ")")
      summtable[i, 4] <- paste0(roundfunc(median(Data[, vars[i]], na.rm = T), roundplace), " (", 
                             roundfunc(min(Data[, vars[i]], na.rm = T), roundplace), " - ", 
                             roundfunc(max(Data[, vars[i]], na.rm = T), roundplace), ")")
      summtable[i, 2] <- length(which(is.na(Data[, vars[i]]) == F))
      if(!is.null(LODvars)){
        nLOD <- length(which(Data[, LODvars[i]] == 1))
        summtable[i, 5] <- paste0(nLOD, " (", round((nLOD/length(which(is.na(Data[, vars[i]]) == F)))*100, roundplace), "%)")
      }
    }
    if(!is.null(LODvars)){
      names(summtable) <- c("Variable", "N", "Mean (SD)", "Median (Range)", 
                          "N <LOD (%)")
    } else {
      names(summtable) <- c("Variable", "N", "Mean (SD)", "Median (Range)")
    }
  }
  return(summtable)
}

#Plot of correlation coefficients
mypearson <- function(contnames, catnames = NULL, Data = tides_ndphxcc, plot = T, exclude = F, cutoff = 0.4, 
                    corindex = NULL, usemethod = "pearson", plottitle = "", fontsize = 0.6){
  if(!is.null(catnames)){
    mydums <- predict(dummyVars(formula(paste0("~", paste(catnames, collapse = "+"))), 
                              Data), Data)
    corprepmat <- cbind(Data[, contnames], mydums)
  } else {
    corprepmat <- Data[, contnames]
  }
  cor1 <- cor(corprepmat, use = "pairwise.complete.obs", method = usemethod)
  if(exclude == T){
    cor_index <- which(unlist(lapply(as.data.frame(cor1), 
                                   function(x) max(abs(x[x%nin%c(-1, 1)]))))>cutoff)
  }
  if(plot == T){
    col1 <- colorRampPalette(c("red", "grey", "blue"))
    if(exclude == T){
      corrplot::corrplot.mixed(cor1[cor_index, cor_index], upper = "ellipse", lower = "number", 
                               lower.col = col1(20), tl.pos = "lt", tl.col = "black", tl.cex = fontsize, 
                               number.cex = fontsize)
    } else if(!is.null(corindex)){
      corrplot::corrplot.mixed(cor1[corindex[[1]], corindex[[2]]], upper = "ellipse", lower = "number", 
                               lower.col = col1(20), tl.pos = "lt", tl.col = "black", tl.cex = fontsize, 
                               number.cex = fontsize)
    } else {
      corrplot::corrplot.mixed(cor1, upper = "ellipse", lower = "number", 
                               lower.col = col1(20), tl.pos = "lt", tl.col = "black", tl.cex = fontsize, 
                               number.cex = fontsize)
    }
  } else {
    if(exclude == T){
      return(cor1[cor_index, cor_index])
    } else if(!is.null(corindex)){
      return(cor1[corindex[[1]], corindex[[2]]])
    } else {
      return(cor1)
    }
  }
}

#Continuous variable summary statistics
myvarcompsummary <- function(colnames, Data, compvar = NULL, outname = "Score", 
                           compname = "Version", footnote = NULL, kable = T){
  myskew <- function(x){
    if(length(x[!is.na(x)])<3){
      skew <- NA 
    } else {
      skew <- e1071::skewness(x, na.rm = T, type = 2)
    }
    return(skew)
  }
  mykurt <- function(x){
    if(length(x[!is.na(x)])<3){
      kurt <- NA 
    } else {
      kurt <- e1071::kurtosis(x, na.rm = T, type = 2)
    }
    return(kurt)
  }
  nacountperc <- function(x) {paste0(length(which(is.na(x) == T)), " (",     
                                   round((length(which(is.na(x) == T))/length(x))*100, 2), "%)")}
  
  if(is.null(compvar)){
    Outtable <- as.data.frame(matrix(NA, length(colnames), 13))
    Outtable[, 1] <- colnames
    Outtable[, 2] <- "Combined"
    aggform1 <- formula(paste0("cbind(", paste(colnames, collapse = ", "), ")~1"))
    Outtable[, 3] <- unlist(aggregate(aggform1, Data, FUN = function(x) length(x[!is.na(x)]), na.action = NULL))
    for(i in 1:length(colnames)){
      aggform2 <- formula(paste0(colnames[i], "~1"))
      totagg1 <- as.data.frame(aggregate(aggform2, Data[!is.na(Data[, colnames[i]]), ], 
                                       FUN = summary, na.action = NULL))
      Outtable[i, 4:9] <- unlist(totagg1)[1:6]
    }
    Outtable[, 10] <- unlist(aggregate(aggform1, Data, FUN = function(x) sd(x, na.rm = T), 
                                    na.action = NULL))
    Outtable[, 11] <- unlist(aggregate(aggform1, Data, FUN = myskew, na.action = NULL))
    Outtable[, 12] <- unlist(aggregate(aggform1, Data, FUN = mykurt, na.action = NULL))
    Outtable[, 13] <- unlist(aggregate(aggform1, Data, FUN = nacountperc, na.action = NULL))
  } else {
    ncomplevels <- length(levels(as.factor(Data[, compvar])))
    Outtable <- as.data.frame(matrix(NA, length(colnames)*ncomplevels, 11))
    Outtable[, 1] <- rep(colnames, each = ncomplevels)
    Outtable[, 2] <- rep(levels(as.factor(Data[, compvar])), times = length(colnames))
    aggform1 <- formula(paste0("cbind(", paste(colnames, collapse = ", "), ")~", compvar))
    totagg1 <- aggregate(aggform1, Data, FUN = function(x) summary(x)[1:6], na.action = NULL)
    Outtable[, 3] <- unlist(aggregate(aggform1, Data, FUN = function(x) length(x[!is.na(x)]), 
                                   na.action = NULL))[(1+ncomplevels):
                                                      (ncomplevels*length(colnames)+ncomplevels)]
    for(i in 1:length(colnames)){
      for(j in 1:ncomplevels){
        rownum <- (1+(j-1))+((i-1)*ncomplevels)
        Outtable[rownum, 4:9] <- totagg1[[i+1]][j, ]
      }
    }
    Outtable[, 10] <- as.numeric(as.character(
      unlist(aggregate(aggform1, Data, FUN = function(x) sd(x, na.rm = T), 
                                    na.action = NULL))[(1+ncomplevels):
                                                       (ncomplevels*length(colnames)+ncomplevels)]))
    Outtable[, 11] <- as.numeric(as.character(unlist(aggregate(aggform1, Data, FUN = myskew, 
                                    na.action = NULL))[(1+ncomplevels):
                                                       (ncomplevels*length(colnames)+ncomplevels)]))
    Outtable[, 12] <- as.numeric(as.character(unlist(aggregate(aggform1, Data, FUN = mykurt, 
                                    na.action = NULL))[(1+ncomplevels):
                                                       (ncomplevels*length(colnames)+ncomplevels)]))
    Outtable[, 13] <- unlist(aggregate(aggform1, Data, FUN = nacountperc, 
                                    na.action = NULL))[(1+ncomplevels):
                                                       (ncomplevels*length(colnames)+ncomplevels)]
  }
  
  
  Outtable[, 4:12] <- round(Outtable[, 4:12], 2)
  names(Outtable) <- c(outname, compname, "N", "Min", "1st Q", "Median", "Mean", "3rd Q", "Max", 
                     "SD", "Skew", "Kurtosis", "# Missing (%)")
  Outtable <- Outtable[!is.nan(Outtable$Mean), ]
  rownames(Outtable) <- NULL
  if(is.null(compvar)){
    Outtable <- Outtable[, -2]
  }
  
  if(kable == T){
    if(is.null(footnote)){
      K <- kable(Outtable, "html") %>%
        kable_styling(font_size = 16) %>%
        column_spec(1, bold = T) %>%
        collapse_rows(columns = 1, valign = "top")
    } else {
      K <- kable(Outtable, "html") %>%
        kable_styling(font_size = 16) %>%
        column_spec(1, bold = T) %>%
        collapse_rows(columns = 1, valign = "top") %>%
        footnote(general = footnote)
    }
    return(list(Table = Outtable, Kable = K))
  } else {
    return(Outtable)
  }
}

#Boxplots of variables compared by levels of a categorical variable
subset_bplot <- function(grepname = NULL, Data, compvar = "Trimester", compname = "Trimester", 
                       comp = T, colnames = NULL, meancolor = "purple", palette = NULL, matchaxes = F, 
                       starvec = NULL){
  if(is.null(colnames)){
    if(is.null(grepname)){
      meltDat <- Data
      if(comp){
        meltDat <- meltDat[, -which(names(meltDat) == compvar)]
        colnames <- setdiff(names(Data), compvar)
      } else {
        colnames <- names(Data)
      }
    } else {
      meltDat <- Data[, grep(grepname, names(Data))]
    }
  } else {
    if(is.null(grepname)){
      meltDat <- Data[, colnames]
    } else {
      meltDat <- Data[, c(which(names(Data)%in%colnames), grep(grepname, names(Data)))]
    }
  }
  
  names(meltDat) <- gsub("log10|_SGadj|_cradj", "", names(meltDat))
  subData_long <- reshape2::melt(meltDat)
  subData_long[, 3] <- "Combined"
  names(subData_long)[3] <- "Comparison"
  
  if(comp == T){
    greppattern <- paste(c(grepname, colnames, compvar), collapse = "|")
    if(any(grepl("(", names(Data), fixed = T))){
      greppattern <- gsub("(", "\\(", greppattern, fixed = T)
      greppattern <- gsub(")", "\\)", greppattern, fixed = T)
    }
    subData_longer <- reshape2::melt(Data[, c(grep(greppattern, 
                                                names(Data)))], id.vars = compvar)
    subData_longer$variable <- gsub("log10|_SGadj|_cradj", "", subData_longer$variable)
    names(subData_longer)[1] <- "Comparison"
    
    for(i in levels(subData_long$variable)){
      myarg <- length(which(is.na(subData_long[subData_long$variable == i, "value"]) == T)) >= 
        min(unlist(summary(as.factor(Data[!is.na(Data[, compvar]), compvar]))))
      if(myarg == TRUE){
        subData_long[subData_long$variable == i, "value"] <- NA
      } 
    }
    subData_longest <- rbind(subData_long[, c(3, 1:2)], subData_longer)
    subData_longest$Comparison <- factor(subData_longest$Comparison, 
                                       levels = c(levels(subData_longer$Comparison), "Combined"))
    if(is.null(palette)){
      palette <- RColorBrewer::brewer.pal(length(unique(subData_longest$Comparison)), "Set1")
    }
    if(!is.null(starvec)){
      extraDat <- data.frame(variable = names(meltDat), stars = starvec, yval = NA)
      for(i in names(meltDat)){
        extraDat[extraDat$variable == i, "yval"] <- 
          max(subData_longest[subData_longest$variable == i, "value"])+
          diff(range(subData_longest[subData_longest$variable == i, "value"]))/10
      }
      subData_longest$stars <- 
        extraDat[match(subData_longest$variable, extraDat$variable), "stars"]
      subData_longest$staryval <- 
        extraDat[match(subData_longest$variable, extraDat$variable), "yval"]
      if(matchaxes){
        gg1 <- ggplot(data = subData_longest, aes(x = variable, y = value, fill = Comparison))+
          geom_boxplot()+theme_bw()+
          stat_summary(geom = "point", color = meancolor, size = 3, 
                       shape = 18, fun = mean, position = position_dodge(0.75))+
          geom_text(aes(label = stars, y = staryval), size = 6, color = "black")+
          scale_fill_manual(values = palette, name = compname)+
          facet_grid(~variable, scales = "free_x", space = "free_x")+
          theme(legend.position = "bottom", 
                axis.text.x = element_text(vjust = 0, angle = 90, size = 12), 
                axis.title = element_blank(), 
                strip.background  =  element_blank(), 
                strip.text = element_blank(), 
                legend.text = element_text(size = 18), 
                legend.title = element_text(size = 18), 
                axis.text.y = element_text(angle = 90, hjust = 0.5))
      } else {
        gg1 <- ggplot(data = subData_longest, aes(x = variable, y = value, fill = Comparison))+
          geom_boxplot()+theme_bw()+
          stat_summary(geom = "point", color = meancolor, size = 3, 
                       shape = 18, fun = mean, position = position_dodge(0.75))+
          geom_text(aes(label = stars, y = staryval), size = 8, color = "black")+
          scale_fill_manual(values = palette, name = compname)+
          facet_wrap(~variable, scales = "free", ncol = length(unique(subData_longest$variable)))+
          theme(legend.position = "bottom", 
                axis.text.x = element_text(vjust = 0, angle = 90, size = 12), 
                axis.title = element_blank(), 
                strip.background  =  element_blank(), 
                strip.text = element_blank(), 
                legend.text = element_text(size = 18), 
                legend.title = element_text(size = 18), 
                axis.text.y = element_text(angle = 90, hjust = 0.5))
      }
    } else {
      if(matchaxes){
        gg1 <- ggplot(data = subData_longest, aes(x = variable, y = value, fill = Comparison))+
          geom_boxplot()+theme_bw()+
          stat_summary(geom = "point", color = meancolor, size = 3, 
                       shape = 18, fun = mean, position = position_dodge(0.75))+
          scale_fill_manual(values = palette, name = compname)+
          facet_grid(~variable, scales = "free_x", space = "free_x")+
          theme(legend.position = "bottom", 
                axis.text.x = element_text(vjust = 0, angle = 90, size = 12), 
                axis.title = element_blank(), 
                strip.background  =  element_blank(), 
                strip.text = element_blank(), 
                legend.text = element_text(size = 18), 
                legend.title = element_text(size = 18), 
                axis.text.y = element_text(angle = 90, hjust = 0.5))
      } else {
        gg1 <- ggplot(data = subData_longest, aes(x = variable, y = value, fill = Comparison))+
          geom_boxplot()+theme_bw()+
          stat_summary(geom = "point", color = meancolor, size = 3, 
                       shape = 18, fun = mean, position = position_dodge(0.75))+
          scale_fill_manual(values = palette, name = compname)+
          facet_wrap(~variable, scales = "free", ncol = length(unique(subData_longest$variable)))+
          theme(legend.position = "bottom", 
                axis.text.x = element_text(vjust = 0, angle = 90, size = 12), 
                axis.title = element_blank(), 
                strip.background  =  element_blank(), 
                strip.text = element_blank(), 
                legend.text = element_text(size = 18), 
                legend.title = element_text(size = 18), 
                axis.text.y = element_text(angle = 90, hjust = 0.5))
      }
    }
    
    
  } else {
    subData_longest <- subData_long
    subData_longest <- subData_longest[!is.na(subData_longest$value), ]
    if(is.null(palette)){
      palette <- RColorBrewer::brewer.pal(length(unique(subData_longest$Comparison)), "Set1")
    }
    if(matchaxes){
      gg1 <- ggplot(data = subData_longest, aes(x = variable, y = value, fill = Comparison))+
        geom_boxplot()+theme_bw()+
        stat_summary(geom = "point", color = meancolor, size = 3, 
                     shape = 18, fun = mean, position = position_dodge(0.75))+
        scale_fill_manual(values = palette)+
        facet_grid(~variable, scales = "free_x", space = "free_x")+
        theme(legend.position = "none", 
              axis.text.x = element_text(vjust = 0, angle = 90, size = 12), 
              axis.title = element_blank(), 
              strip.background  =  element_blank(), 
              strip.text = element_blank(), 
              axis.text.y = element_text(angle = 90, hjust = 0.5))
    } else {
      gg1 <- ggplot(data = subData_longest, aes(x = variable, y = value, fill = Comparison))+
        geom_boxplot()+theme_bw()+
        stat_summary(geom = "point", color = meancolor, size = 3, 
                     shape = 18, fun = mean, position = position_dodge(0.75))+
        scale_fill_manual(values = palette)+
        facet_wrap(~variable, scales = "free", ncol = ncol(meltDat))+
        theme(legend.position = "none", 
              axis.text.x = element_text(vjust = 0, angle = 90, size = 12), 
              axis.title = element_blank(), 
              strip.background  =  element_blank(), 
              strip.text = element_blank(), 
              axis.text.y = element_text(angle = 90, hjust = 0.5))
    }
    
  }
  
  return(gg1)
}

#Simple summary output for categorical variables
myvarcompsummary_cat <- function(catnames, Data, catvarnames = NULL, fontsize = 12, 
                                 omitmissing = F){
  if(is.null(catvarnames)) catvarnames <- catnames
  for(i in 1:length(catnames)){
    if(omitmissing) tempDat <- Data[complete.cases(Data[, catnames[i]]), 
                                    ] else tempDat <- Data
    if(!is.factor(tempDat[, catnames[i]])){
      tempDat[, catnames[i]] <- as.factor(tempDat[, catnames[i]])}
    summ <- summary(tempDat[, catnames[i]])
    kablemat <- matrix(NA, 2, length(summ))
    kablemat[1, ] <- names(summ)
    for(j in 1:length(names(summ))){
      kablemat[2, j] <- paste0(unlist(summ)[j], " (",  round((unlist(summ)[j]/sum(
        unlist(summ)))*100, 2), "%)")
    }
    K1 <- kable(kablemat, "html", caption = paste0(catvarnames[i])) %>%
      kable_styling(font_size = fontsize, full_width  =  F) %>%
      row_spec(1, bold = T)
    print(K1)
  }
}

#Quantile transform a numeric vector into q quantiles (e.g., quartiles if q=4)
myquantile <- function(x, q, labelnumbers = T){
  mybreaks <- quantile(x, probs = seq(0, 1, by = 1/q))
  if(labelnumbers){
    myquants <- cut(x, breaks = mybreaks, labels = 1:q)
  } else {
    myquants <- cut(x, breaks = mybreaks)
  }
  return(myquants)
}

#Perform a permutation test form of a Student's t-test
perm.ttest <- function(x, y = NULL, compvar = NULL, Data, niter = 10000){
  if(!is.null(y)){
    if(is.character(x)){x <- Data[, x]}
    if(is.character(y)){y <- Data[, y]}
  } else {
    if(!is.factor(Data[, compvar])) Data[, compvar] <- as.factor(Data[, compvar])
    y <- Data[, x][which(Data[, compvar] == levels(Data[, compvar])[2])]
    x <- Data[, x][which(Data[, compvar] == levels(Data[, compvar])[1])]
  }
  origdiff <- mean(y, na.rm = T)-mean(x, na.rm = T)
  comb <- c(x, y)
  permcomb <- matrix(NA, length(comb), niter)
  permcomb <- apply(permcomb, 2, function(u) sample(comb))
  meandiff <- function(x, y) mean(y, na.rm = T) - mean(x, na.rm = T)
  assign_contrast <- function(combined){
    a1 <- combined[1:length(x)]
    a2 <- combined[(length(x)+1):length(combined)]
    meandiff(a1, a2)
  }
  permdiffs <- apply(permcomb, 2, assign_contrast)
  permp <- length(which(abs(permdiffs) > abs(origdiff)))/niter
  return(permp)
}

#Output coefficients and robust standard errors and p-values
robustse <- function(x, HCtype = "HC0", usedf = F){
  mod1 <- lmtest::coeftest(x, vcov = function(x) vcovHC(x, type = HCtype))
  if(!usedf){
    cis<-lmtest::coefci(x, vcov=function(x) sandwich::vcovHC(x, type = HCtype))
  } else {
    cis<-lmtest::coefci(x, vcov = function(x) sandwich::vcovHC(x, type = HCtype),
                df = x$df.residual)
  }
  mod1<-cbind(mod1, cis)
  return(mod1)
}

#calculate a t-distribution-based set of confidence intervals
confint_qt <- function(beta, se, DF, IQR = 1, level = 0.95){
  ci.lower<-(beta * IQR) - ((se * IQR) * qt(((level/2) + 0.5), DF))
  ci.upper<-(beta * IQR) + ((se * IQR) * qt(((level/2) + 0.5), DF))
  CIs<-data.frame(CIL = ci.lower, CIU = ci.upper)
  return(CIs)
}

#Capitalize the 1st letter of either just the 1st word or each word in a string
capfirst <- function(x, eachword = F) {
  if(eachword){
    s1 <- strsplit(x, " ")[[1]]
    paste0(toupper(substring(s1, 1, 1)), substring(s1, 2), collapse = " ")
  } else {
    paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  }
}

#The combination formula
nCr <- function(n, r) factorial(n)/(factorial(r) * factorial(n - r))

#Get all combos
getcombos <- function(x, y = NULL, r = 2){
  if(is.character(x)) x <- as.factor(x)
  if(is.null(levels(x))) xcomb <- "No combos" else {
    tmp <- tryCatch(combn(levels(x), r, simplify = F), error = function(e) NULL)
    if(is.null(tmp)) xcomb <- "No combos" else xcomb <- sapply(
      tmp, function(ww) paste(ww, collapse = "|"))
  }
  if(!is.null(y)){
    if(is.character(y)) y <- as.factor(y)
    if(is.null(levels(y))){
      ycomb <- "No combos"
      xycomb <- xcomb
    } else {
      tmp2 <- tryCatch(combn(levels(y), r, simplify = F), 
                       error = function(e) NULL)
      if(is.null(tmp2)){
        ycomb <- "No combos"
        xycomb <- xcomb} else {
        ycomb <- sapply(tmp2, function(ww) paste(ww, collapse = "|"))
        if(length(xcomb) ==1 && xcomb == "No combos"){
          xycomb <- ycomb
        } else {
          xycomb <- expand.grid(X = xcomb, Y = ycomb)
        }
      }
    }
  }
  if(is.null(y)) retlist <- list(xcomb = xcomb) else 
    retlist <- list(xcomb = xcomb, ycomb = ycomb, xycomb = xycomb)
  return(retlist)
}

#Find the N closest values in a vector
findClosest <- function(x, n){
  x <- sort(x)
  x[seq.int(which.min(diff(x, lag = n - 1L)), length.out = n)]
}

#Remove all non-numeric characters from vector and make it as.numeric
numerize <- function(x) as.numeric(gsub("[^0-9.-]", "", x))

#Remove all columns from a data frame that are all NA, NaN, or user-specified
# missing characters (by default this includes the empty string "")
dropnacols <- function(dat, misschars = c("")){
  dropcols <- which(sapply(dat, function(x) all(is.na(x)|is.nan(
    as.numeric(x))|as.character(x) %in% misschars)))
  return(dat[, -dropcols])
}

#Get path to example data files
findexampledata <- function(path = NULL){
  if(is.null(path)){ 
    dir(system.file("exdata", package = "DrewDayRFunctions"))
  } else {
    system.file("exdata", path, package = "DrewDayRFunctions", mustWork = T)
  }
}

#The following rewrites ggalt::stat_xspline so that you don't have to download 
# that package's annoying proj4 dependency, which can cause download errors. 
stat_xspline <- function (mapping = NULL, data = NULL, geom = "line", 
                          position = "identity", na.rm = TRUE, show.legend = NA, 
                          inherit.aes = TRUE, spline_shape = -0.25, 
                          open = TRUE, rep_ends = TRUE, ...){
  layer(stat = StatXspline, data = data, mapping = mapping, 
        geom = geom, position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, params = list(
          spline_shape = spline_shape, open = open, na.rm = na.rm, 
          rep_ends = rep_ends, ...))
}
StatXspline <- ggproto("StatXspline", Stat, required_aes = c("x", "y"),
                       compute_group = function(self, data, scales, params,
                                                spline_shape=-0.25, open=TRUE, 
                                                rep_ends=TRUE) {
                         tf <- tempfile(fileext=".png")
                         png(tf)
                         plot.new()
                         tmp <- graphics::xspline(data$x, data$y, spline_shape, open, 
                                        rep_ends, draw=FALSE, NA, NA)
                         invisible(dev.off())
                         unlink(tf)
                         data.frame(x=tmp$x, y=tmp$y)
                       })
geom_xspline <- function(mapping = NULL, data = NULL, stat = "xspline",
                         position = "identity", na.rm = TRUE, show.legend = NA,
                         inherit.aes = TRUE,
                         spline_shape=-0.25, open=TRUE, rep_ends=TRUE, ...) {
  layer(geom = GeomXspline,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      spline_shape = spline_shape,
      open = open,
      na.rm = na.rm,
      rep_ends = rep_ends,
      ...)
  )
}
GeomXspline <- ggproto("GeomXspline", GeomLine, required_aes = c("x", "y"),
                       default_aes = aes(colour = "black", linewidth = 0.5, 
                                         linetype = 1, alpha = NA))