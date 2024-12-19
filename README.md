**DrewDayRFunctions** is an R package I created to share convenience functions I've created throughout my career to aid in epidemiologic and other statistical analyses. These include summary statistics functions, functions for rapidly performing linear and functional regressions for combinations of predictors and outcomes and returning relevant parameters, visualization functions, and functions for nonlinear mediation analysis. Currently this package contains the following functions:

  1. **GLMResults**: A function for taking combinations of predictor of interest, outcome, and covariate variables, performing linear or binary (logistic, log binomial, probit, robust Poisson, etc.) generalized linear regressions, and returning 1. a table of coefficient and diagnostic values, 2. a summary plot, and 3. a list of all model objects for any follow-up analyses. Options include adding categorical or continuous interaction terms, using robust errors, performing multiple imputation using MICE, performing post-hoc power tests, and rerunning models when excluding high-leverage points.
  2. **ContrastCoefficients**: A function to extract all the binary combinations of contrasts for a categorical (factor) predictor of interest and/or interaction term from a linear model. The typical output of linear model only returns the coefficients for categorical variables in the context of a single referent level. For example, if a predictor variable of interest is a 4-level categorical treatment variable with factor levels {"placebo", "drug1", "drug2", "drug3"}, a linear model will only return the "drug1 - placebo", "drug2 - placebo", and "drug3 - placebo" binary contrast coefficients. This function will return all binary contrast coefficients for all possible referent levels, meaning the aforementioned coefficients will be returned, as well the coefficients for "drug2 - drug1", "drug 3 - drug 1", and "drug 3 - drug 2" contrasts. Furthermore, ContrastCoefficients can also provide all contrasts for a categorical interaction term as well, including all marginal coefficients for the main effects (e.g., predictor|interacting variable level 1, predictor|interacting variable level 2, etc).
  3. **GAMResults**: A function for taking combinations of predictor of interest, outcome, and covariate variables, performing linear or binary generalized additive models, and returning 1. a table of parameter and diagnostic values, 2. summary plot(s), and 3. a list of all model objects. Options include adding continuous and categorical variable by smooth interactions, performing smooth selection, and rerunning models when excluding upper/lower quantiles of data to examine edge effects on smooth fits.
  4. **InteractionCoefPlot**: A function for visualizing continuous by continuous interactions from linear regression models. These plots show how the coefficient between one variable in a bivariate continuous interaction and the outcome changes over levels of the interacting variable.
  5. **GenericDataSummary**: This function takes a data frame and outputs an html Markdown document that contains an interactive table with a data dictionary if one is provided, an interactive table of summary statistics for the continuous variables as well as boxplots and histograms for those variables, and a searchable table showing frequencies for all categorical variables that also uses a bar plot to visualize frequencies. This is an easy way to summarize a newly obtained dataset or to send a general summary to someone who will be receiving a dataset.
  6. **GetExcelColors**: This function takes an Excel spreadsheet from a .xlsx file and creates a data frame where the values in the dataset are replaced with all the background colors in that Excel sheet. It also provides a table and heatmap plot showing how many cells have each color detected in the Excel spreadsheet. This function is useful for when important information is encoded in the background colors of an Excel file.
  7. **MediationCurve**: This function extends the functionality of the 'mediate' function from the 'mediation' R package for causal mediation analysis as defined by Imai et al. 2010 for modelling linear and nonlinear mediation effects in the case of a continuous treatment. This function is largely inspired by the causal-curve Python package by Roni Kobrosly. Mediation analysis is applied over a range of continuous treatment values to obtain a curve of the average causal mediation effect (ACME), average direct effect (ADE), total effect, and the proportion mediated. This is useful for characterizing nonlinear mediation effects.
