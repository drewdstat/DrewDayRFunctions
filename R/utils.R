countperc <- function(x, roundplace = 2, na.rm = T){
  if(na.rm) summ <- summary(x[which(!is.na(x))]) else summ <- summary(x)
  perc <- prop.table(summ)
  summ <- paste0(summ, " (", round(perc*100, roundplace), "%)")
  names(summ) <- names(perc)
  return(summ)
}
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
                        function(x) dimnames(table(x, Data[complete.cases(Data[, catnames]), compvar]))[[1]])
    } else {
      responses <- lapply(Data[, catnames], 
                        function(x) dimnames(table(x, Data[, compvar]))[[1]])
    }
    if(any(sapply(responses, length)!= 2)){
      newmat <- fastDummies::dummy_cols(Data[, c(catnames, compvar)], 
                                      select_columns = catnames[which(sapply(responses, length)!= 2)], 
                                      ignore_na = T)[, -which(sapply(responses, length)!= 2)] 
    } else {
      newmat <- Data[, c(catnames, compvar)]
    }
    newmat2 <- newmat[, -which(names(newmat) =  = compvar)]
    if(is.list(responses)){responses <- unlist(responses)}
    pwdf <- as.data.frame(t(sapply(newmat[, -which(names(newmat) =  = compvar)], 
                                 function(x) pairwiseNominalIndependence(table(x, newmat[, compvar]), fisher = F, 
                                                                         gtest = F, chisq = T, method = "none"))))
    summncol <- 3+ncomplevel
    if(totalcol){summncol <- summncol+1}
    summtab <- as.data.frame(matrix(NA, 2*nrow(pwdf), summncol))
    summtab[, 1] <- rep(dimnames(pwdf)[[1]], each = 2)
    for(i in 1:ncol(newmat2)){
      tab1 <- table(newmat2[, i], newmat[, compvar])
      perctab1 <- t(prop.table(t(tab1), 1)*100)
      tab1 <- matrix(paste0(tab1, " (", round(perctab1, roundplace), "%)"), nrow(tab1), ncol(tab1))
      dimnames(tab1) <- dimnames(perctab1)
      summtab[(1+((i-1)*2)):(i*2), 2] <- dimnames(tab1)[[1]]
      summtab[(1+((i-1)*2)):(i*2), 3:(2+ncol(tab1))] <- tab1
      if(totalcol){
        tab2 <- table(newmat2[, i])
        perctab2 <- prop.table(t(tab2), 1)*100
        tab2 <- paste0(tab2, " (", round(perctab2, roundplace), "%)")
        summtab[(1+((i-1)*2)):(i*2), summncol] <- tab2
      }
    }
    if(totalcol){
      summtab[, (summncol-1)] <- rep(unlist(pwdf$p.Chisq), each = 2)
      stars <- ifelse(summtab[, (summncol-1)]<0.05, "*", "")
      summtab[, (summncol-1)] <- paste0(summtab[, (summncol-1)], stars)
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
    #responses <- c(unique(unlist(sapply(as.character(Data[
    #  complete.cases(Data[, catnames]), catnames]), unique))))
    #responses <- unique(c(sapply(Data[complete.cases(Data[, catnames]), catnames], unique)))
    if(allcomplete){
      if(nmissing){
        responses <- unique(c(lapply(Data[complete.cases(Data[, catnames]), catnames], function(x) summary(as.factor(x)))))
      } else {
        responses <- unique(c(lapply(Data[complete.cases(Data[, catnames]), catnames], 
                                   function(x) summary(as.factor(x[which(!is.na(x))])))))
      }
    } else {
      if(nmissing){
        responses <- unique(c(lapply(Data[, catnames], function(x) summary(as.factor(x)))))
      } else {
        responses <- unique(c(lapply(Data[, catnames], function(x) summary(as.factor(x[which(!is.na(x))])))))
      }
    }
    
    # if(nmissing) summtabcol <- 4 else 
    summtabcol <- 3
    summtab <- as.data.frame(matrix(NA, 0, summtabcol))
    if(!is.null(catvarnames)) catnames <- catvarnames
    for(ll in 1:length(responses)){
      temptab <- as.data.frame(matrix(NA, sapply(responses, length)[ll], 3))
      temptab[, 1] <- catnames[ll]
      temptab[, 2] <- names(responses[[ll]])
      temptab[, 3] <- responses[[ll]]
      temptab[, 4] <- 100*(responses[[ll]]/sum(responses[[ll]]))
      temptab[, 3] <- paste0(temptab[, 3], " (", round(temptab[, 4], roundplace), "%)")
      temptab[, 4] <- NULL
      # if(nmissing){
      #   temptab[, 4] <- paste0(length(which(is.na(Data[, catnames[ll]]))), " (", 
      #                       round(100*(length(which(is.na(Data[, catnames[ll]])))/
      #                                    length(Data[, catnames[ll]])), roundplace), "%)")
      # }
      summtab <- rbind(summtab, temptab)
    };rm(ll, temptab)
    # if(nmissing) names(summtab) <- c("Variable", "Category", "N (%)", "N Missing (%)") else{}
    names(summtab) <- c("Variable", "Category", "N (%)")
    # summtab <- as.data.frame(matrix(NA, length(catnames), (length(responses)+1)))
    # summtab[, 1] <- catnames
    # for(i in 1:length(responses)){
    #   lengths <- apply(Data[, catnames], 2, function(x) length(which(as.character(x) =  = responses[i])))
    #   summtab[, (i+1)] <- paste0(lengths, " (", round((lengths/nrow(Data))*100, roundplace), "%)")
    # }
    # names(summtab) <- c("Variable", as.character(responses))
    # if(nmissing){
    #   NMissing <- apply(Data[, catnames], 2, function(x) length(which(is.na(x))))
    #   summtab$Missing <- paste0(NMissing, " (", round((NMissing/nrow(Data))*100, roundplace), "%)")
    # }
    if(kable){
      K1 <- kable(summtab[, -1], "html") %>%
        kable_styling(font_size = fontsize, full_width  =  F) %>%
        pack_rows(index = table(forcats::fct_inorder(summtab[, 1]))) %>%
        column_spec(1, bold = T)
      print(K1)
    } else {
      return(summtab)
    }
  }
}
summtable_tt <- function(vars, Data, varnames = NULL, LODvars = NULL, compvar = NULL, newcomplvls = NULL, 
                       roundplace = 2, padj = F, padjmethod = "BH", nonpar = F, perm = F, signif = F){
  if(signif){
    roundfunc <- function(x, roundplace){signif(x, roundplace)}
  } else {
    roundfunc <- function(x, roundplace){round(x, roundplace)}
  }
  if(is.null(varnames)){
    varnames <- vars
  }
  if(!is.null(compvar)){
    if(!is.null(newcomplvls)){
      Data[, compvar] <- factor(as.character(Data[, compvar]), levels = newcomplvls)
    } else {
      if(class(Data[, compvar]) =  = "character"){
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
      nobs <- length(which(is.na(Data[, vars[i]]) =  = F)) #nmiss
      summtable1[i, 3] <- nobs
      #summtable1[i, 5] <- paste0(nmiss, " (", round((nmiss/nrow(Data))*100, roundplace), "%)")
      if(nonpar){
        if(nlevels =  = 2){
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
        if(nlevels =  = 2){
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
        nLOD <- length(which(Data[, LODvars[i]] =  = 1))
        summtable1[i, 7] <- paste0(nLOD, " (", round((nLOD/nobs)*100, roundplace), "%)")
      }
      for(j in 1:nlevels){
        rownum <- (1+(j-1))+((i-1)*nlevels)
        thislevel <- levels(Data[, compvar])[j]
        summtable2[rownum, 1] <- varnames[i]
        summtable2[rownum, 2] <- thislevel
        summtable2[rownum, 4] <- paste0(roundfunc(mean(Data[Data[, compvar] =  = thislevel, vars[i]], na.rm = T), roundplace), " (", 
                                     roundfunc(sd(Data[Data[, compvar] =  = thislevel, vars[i]], na.rm = T), roundplace), ")")
        summtable2[rownum, 5] <- paste0(roundfunc(median(Data[Data[, compvar] =  = thislevel, vars[i]], na.rm = T), roundplace), " (", 
                                     roundfunc(min(Data[Data[, compvar] =  = thislevel, vars[i]], na.rm = T), roundplace), " - ", 
                                     roundfunc(max(Data[Data[, compvar] =  = thislevel, vars[i]], na.rm = T), roundplace), ")")
        nobs <- length(which(is.na(Data[Data[, compvar] =  = thislevel, vars[i]]) =  = F))#nmiss
        summtable2[rownum, 3] <- nobs
        #summtable2[rownum, 5] <- paste0(nmiss, " (", 
        #                             round((nmiss/nrow(Data[Data[, compvar] =  = thislevel, ]))*100, roundplace), "%)")
        if(nlevels>2){
          if(nonpar){
            tt <- tryCatch(wilcox.test(Data[which(Data[, compvar] =  = levels(Data[, compvar])[j]), vars[i]], 
                                Data[which(Data[, compvar]!= levels(Data[, compvar])[j]), vars[i]]), 
                         error = function(e) NULL)
            if(is.null(tt)) ttp <- NA else ttp <- tt$p.value
          } else if(perm){
            ttp <- perm.ttest(Data[which(Data[, compvar] =  = levels(Data[, compvar])[j]), vars[i]], 
                           Data[which(Data[, compvar]!= levels(Data[, compvar])[j]), vars[i]])
          } else {
            tt <- tryCatch(t.test(Data[which(Data[, compvar] =  = levels(Data[, compvar])[j]), vars[i]], 
                                Data[which(Data[, compvar]!= levels(Data[, compvar])[j]), vars[i]]), 
                         error = function(e) NULL)
            if(is.null(tt)) ttp <- NA else ttp <- tt$p.value
          }
        }
        summtable2[rownum, 6] <- signif(ttp, roundplace)
        if(!is.null(LODvars)){
          nLOD2 <- length(which(Data[Data[, compvar] =  = thislevel, LODvars[i]] =  = 1))
          LODperc <- round((nLOD2/nrow(Data[which(Data[, compvar] =  = thislevel&
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
    #"N Missing (%)", 
    summtable$Variable <- factor(summtable$Variable, levels = varnames)
    summtable <- summtable[with(summtable, order(Variable)), ]
    summtable$Group <- factor(summtable$Group, levels = unique(summtable$Group))
    if(padj){
      if(nlevels =  = 2){
        pagg <- summtable[which(summtable$Group =  = "Total"), c("Variable", "t-test p-value")]
        pagg$padj <- stats::p.adjust(pagg[, ncol(pagg)], method = padjmethod)
        summtable <- merge(summtable, pagg[, c(1, ncol(pagg))], by = "Variable", all.x = T)
      } else {
        pagg <- summtable[, c(1:2, ncol(summtable))]#which(summtable$Group!= "Total")
        pagg[which(pagg$Group!= "Total"), "padj"] <- 
          stats::p.adjust(pagg[which(pagg$Group!= "Total"), grep("p-value", names(pagg))], method = padjmethod)
        pagg[which(pagg$Group =  = "Total"), "padj"] <- 
          stats::p.adjust(pagg[which(pagg$Group =  = "Total"), grep("p-value", names(pagg))], method = padjmethod)
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
      #nmiss <- length(which(is.na(Data[, vars[i]]) =  = T))
      summtable[i, 2] <- length(which(is.na(Data[, vars[i]]) =  = F))#paste0(nmiss, " (", nmiss/nrow(Data), "%)")
      if(!is.null(LODvars)){
        nLOD <- length(which(Data[, LODvars[i]] =  = 1))
        summtable[i, 5] <- paste0(nLOD, " (", round((nLOD/length(which(is.na(Data[, vars[i]]) =  = F)))*100, roundplace), "%)")
      }
    }
    #if(!is.null(LODvars)){summtable <- summtable[, c(1:2, 5, 3:4)]}
    if(!is.null(LODvars)){
      names(summtable) <- c("Variable", "N", "Mean (SD)", "Median (Range)", 
                          "N <LOD (%)")
    } else {
      names(summtable) <- c("Variable", "N", "Mean (SD)", "Median (Range)")
    }
  }
  return(summtable)
}
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
  if(exclude =  = T){
    cor_index <- which(unlist(lapply(as.data.frame(cor1), 
                                   function(x) max(abs(x[x%nin%c(-1, 1)]))))>cutoff)
  }
  if(plot =  = T){
    col1 <- colorRampPalette(c("red", "grey", "blue"))
    if(exclude =  = T){
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
    if(exclude =  = T){
      return(cor1[cor_index, cor_index])
    } else if(!is.null(corindex)){
      return(cor1[corindex[[1]], corindex[[2]]])
    } else {
      return(cor1)
    }
  }
}
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
  nacountperc <- function(x) {paste0(length(which(is.na(x) =  = T)), " (",     
                                   round((length(which(is.na(x) =  = T))/length(x))*100, 2), "%)")}
  
  if(is.null(compvar)){
    Outtable <- as.data.frame(matrix(NA, length(colnames), 13))
    Outtable[, 1] <- colnames
    Outtable[, 2] <- "Combined"
    aggform1 <- formula(paste0("cbind(", paste(colnames, collapse = ", "), ")~1"))
    #totagg1 <- as.data.frame(aggregate(aggform1, Data, FUN = summary, na.action = NULL))
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
        Outtable[rownum, 4:9] <- totagg1[[i+1]][j, ]#unlist(totagg1[j, 1+i][1:6]) 
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
  
  
  #Outtable[, 4:12] <- apply(Outtable[, 4:12], 2, function(x) ifelse(x>10&x<100, round(x, 2), signif(x, 3)))
  Outtable[, 4:12] <- round(Outtable[, 4:12], 2)
  names(Outtable) <- c(outname, compname, "N", "Min", "1st Q", "Median", "Mean", "3rd Q", "Max", 
                     "SD", "Skew", "Kurtosis", "# Missing (%)")
  Outtable <- Outtable[!is.nan(Outtable$Mean), ]
  rownames(Outtable) <- NULL
  if(is.null(compvar)){
    Outtable <- Outtable[, -2]
  }
  
  if(kable =  = T){
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
subset_bplot <- function(grepname = NULL, Data, compvar = "Trimester", compname = "Trimester", 
                       comp = T, colnames = NULL, meancolor = "purple", palette = NULL, matchaxes = F, 
                       starvec = NULL){
  if(is.null(colnames)){
    if(is.null(grepname)){
      meltDat <- Data
      if(comp){
        meltDat <- meltDat[, -which(names(meltDat) =  = compvar)]
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
  
  if(comp =  = T){
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
      myarg <- length(which(is.na(subData_long[subData_long$variable =  = i, "value"]) =  = T))> = 
        min(unlist(summary(as.factor(Data[!is.na(Data[, compvar]), compvar]))))
      if(myarg =  = TRUE){
        subData_long[subData_long$variable =  = i, "value"] <- NA
      } 
    }
    subData_longest <- rbind(subData_long[, c(3, 1:2)], subData_longer)
    #subData_longest <- subData_longest[!is.na(subData_longest$value), ]
    subData_longest$Comparison <- factor(subData_longest$Comparison, 
                                       levels = c(levels(subData_longer$Comparison), "Combined"))
    if(is.null(palette)){
      palette <- RColorBrewer::brewer.pal(length(unique(subData_longest$Comparison)), "Set1")
    }
    if(!is.null(starvec)){
      extraDat <- data.frame(variable = names(meltDat), stars = starvec, yval = NA)
      for(i in names(meltDat)){
        extraDat[extraDat$variable =  = i, "yval"] <- 
          max(subData_longest[subData_longest$variable =  = i, "value"])+
          diff(range(subData_longest[subData_longest$variable =  = i, "value"]))/10
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
myvarcompsummary_cat <- function(catnames, Data, catvarnames = NULL, fontsize = 12, omitmissing = F){
  if(is.null(catvarnames)) catvarnames <- catnames
  for(i in 1:length(catnames)){
    if(omitmissing) tempDat <- Data[complete.cases(Data[, catnames[i]]), ] else tempDat <- Data
    if(!is.factor(tempDat[, catnames[i]])){tempDat[, catnames[i]] <- as.factor(tempDat[, catnames[i]])}
    summ <- summary(tempDat[, catnames[i]])
    kablemat <- matrix(NA, 2, length(summ))
    kablemat[1, ] <- names(summ)
    for(j in 1:length(names(summ))){
      kablemat[2, j] <- paste0(unlist(summ)[j], " (",  round((unlist(summ)[j]/sum(unlist(summ)))*100, 2), "%)")
    }
    K1 <- kable(kablemat, "html", caption = paste0(catvarnames[i])) %>%
      kable_styling(font_size = fontsize, full_width  =  F) %>%
      #group_rows(catnames[i], 1, 2) %>%
      row_spec(1, bold = T)
    print(K1)
  }
}
myquantile <- function(x, q, labelnumbers = T){
  mybreaks <- quantile(x, probs = seq(0, 1, by = 1/q))
  if(labelnumbers){
    myquants <- cut(x, breaks = mybreaks, labels = 1:q)
  } else {
    myquants <- cut(x, breaks = mybreaks)
  }
  return(myquants)
}
perm.ttest <- function(x, y = NULL, compvar = NULL, Data, niter = 10000){
  if(!is.null(y)){
    if(is.character(x)){x <- Data[, x]}
    if(is.character(y)){y <- Data[, y]}
  } else {
    if(!is.factor(Data[, compvar])) Data[, compvar] <- as.factor(Data[, compvar])
    y <- Data[, x][which(Data[, compvar] =  = levels(Data[, compvar])[2])]
    x <- Data[, x][which(Data[, compvar] =  = levels(Data[, compvar])[1])]
  }
  origdiff <- mean(y, na.rm = T)-mean(x, na.rm = T)
  comb <- c(x, y)
  permcomb <- matrix(NA, length(comb), niter)
  permcomb <- apply(permcomb, 2, function(u) sample(comb))
  meandiff <- function(x, y) mean(y, na.rm = T)-mean(x, na.rm = T)
  assign_contrast <- function(combined){
    a1 <- combined[1:length(x)]
    a2 <- combined[(length(x)+1):length(combined)]
    meandiff(a1, a2)
  }
  permdiffs <- apply(permcomb, 2, assign_contrast)
  permp <- length(which(abs(permdiffs)>abs(origdiff)))/niter
  return(permp)
}