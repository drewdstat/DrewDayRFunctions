**DrewDayRFunctions** is an R package I created to share convenience functions I've created throughout my career to aid in epidemiologic and other statistical analyses. These include summary statistics functions, functions for rapidly performing linear and functional regressions for combinations of predictors and outcomes and returning relevant parameters, visualization functions, and functions for calculating statistical power. Currently this package contains the following functions:

  1. **GLMResults**: A function for taking combinations of predictor of interest, outcome, and covariate variables, performing linear or binary (logistic, log binomial, probit, robust Poisson, etc.) generalized linear regressions, and returning 1. a table of coefficient and diagnostic values, 2. a summary plot, and 3. a list of all model objects for any follow-up analyses. Options include adding categorical or continuous interaction terms, using robust errors, performing multiple imputation using MICE, performing post-hoc power tests, and rerunning models when excluding high-leverage points. 
  2. **GAMResults**: A function for taking combinations of predictor of interest, outcome, and covariate variables, performing linear or binary generalized additive models, and returning 1. a table of parameter and diagnostic values, 2. summary plot(s), and 3. a list of all model objects. Options include adding continuous and categorical variable by smooth interactions, performing smooth selection, and rerunning models when excluding upper/lower quantiles of data to examine edge effects on smooth fits.
  3. **InteractionCoefPlot**: A function for visualizing continuous by continuous interactions from linear regression models. These plots show how the coefficient between one variable in a bivariate continuous interaction and the outcome changes over levels of the interacting variable.
