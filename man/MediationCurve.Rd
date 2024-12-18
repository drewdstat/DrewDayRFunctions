% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MediationCurve.R
\name{MediationCurve}
\alias{MediationCurve}
\title{Examine causal mediation effects across a sequence of continuous treatment values}
\usage{
MediationCurve(
  medmodel,
  outmodel,
  treatname,
  medname,
  treat_values = NULL,
  vary_t0 = F,
  fixed_t0 = NULL,
  gridlen = 9,
  gridlowerq = 0.01,
  gridupperq = 0.99,
  qgrid = F,
  progress = T,
  boot = T,
  ...
)
}
\arguments{
\item{medmodel}{This is a mediation model of any class acceptable by the
\code{\link[mediation]{mediate}} function, namely 'lm', 'polr', 'bayespolr',
'glm', 'bayesglm', 'gam', 'rq', 'survreg', or 'merMod'. This model should
include the treatment as an independent variable and the mediator as a
dependent variable.}

\item{outmodel}{This is a mediation model of any class acceptable by the
\code{\link[mediation]{mediate}} function as explained above. This model
should include the treatment and the mediator as independent variables and
the outcome as a dependent variable.}

\item{treatname}{This is a character value for the column name of the
treatment variable as it was called in the formula for \code{medmodel} and
\code{outmodel}.}

\item{medname}{This is a character value for the column name of the
mediator variable as it was called in the formula for \code{medmodel} and
\code{outmodel}.}

\item{treat_values}{An optional numeric vector of treatment values to
iterate across. This defaults to \code{NULL}, meaning that a vector of
treatment values is created using the 'grid' arguments below. If a vector is
provided for \code{treatvalues}, those 'grid' arguments are ignored.}

\item{vary_t0}{This is a logical value determining whether the control
treatment t0 varies along with treated treatment value t1. There are two
theoretical approaches to exploring the mediation results of a continuous
treatment. One is suggested by
\href{https://imai.fas.harvard.edu/research/files/BaronKenny.pdf}{Imai et al. 2010}
and involves holding a "control" treatment value t0 at the same value while
iteratively increasing the "treated" treatment value t1 and calculating
mediation estimates each time. For an example sequence of treatment values
T1, T2, T3, ..., this approach would set t0 = T1 for each iteration, while
t1 = T2 for the first iteration, t1 = T3 for the second, and so on. The
second approach is used by Roni Kobrosly in his
\href{https://causal-curve.readthedocs.io/en/latest/intro.html}{causal-curve}
Python package, and it involves iteratively increasing both t0 and t1 along
the grid of treatment values. In the prior example, t0 = T1 and t1 = T2 for
the first iteration, then t0 = T2 and t1 = T3 for the second, and so on. When
\code{vary_t0 = FALSE}, the former Imai et al. approach is used, and if
\code{vary_t0 = TRUE}, the latter approach is used. This defaults to
\code{FALSE}, meaning that t0 remains fixed at a set value.}

\item{fixed_t0}{If \code{vary_t0 = FALSE}, this argument defines the fixed
control treatment level t0. This can be anywhere between the sequence of
values between \code{gridlowerq} and \code{gridupperq}, which if
\code{vary_t0 = FALSE}, will define t1 treatment values. If \code{NULL}
(the default), this will default to \code{gridlowerq}.}

\item{gridlen}{A length for the vector of treatment values to be generated.
This will be passed to the \code{length.out} argument in
\code{\link[base]{seq}}. This defaults to \code{9}.}

\item{gridlowerq}{A lower quantile for the beginning of the sequence of
treatment values. This defaults to \code{0.01}, meaning that the first value
is the 1st percentile of the treatment values.}

\item{gridupperq}{An upper quantile for the end of the sequence of
treatment values. This defaults to \code{0.99}, meaning that the last value
is the 99th percentile of the treatment values.}

\item{qgrid}{This is a logical value determining whether the sequence of
treatment values should be of quantiles if \code{TRUE} or of evenly spaced
values on the original scale if \code{FALSE}. Here is an example if
\code{gridlen = 5}, the treatment is normally distributed with a mean of 0
and an SD of 1, \code{gridlowerq = 0.025}, and \code{gridupperq = 0.975}.
If \code{qgrid = TRUE}, then the grid would be the 2.5th, 26.25th, 50th,
73.75th, and 97.5th percentiles of the treatment values, which in this case
would be -1.96, -0.64, 0.00, 0.64, and 1.96. If \code{qgrid = FALSE}, these
values would instead be evenly spaced between -1.96 and 1.96, so they would
be -1.96, -0.98, 0, 0.98, 1.96. This defaults to \code{FALSE} because the
smoothing at the edges of the curve tends to more severely miss the actual
point estimates when not having an evenly spaced set of values.}

\item{progress}{A logical value as to whether to include a progress bar
showing the percentage of the mediation functions across the treatment values
that has been completed so far. This function can take some time to run, so
this can be pretty useful to see how far the function has progressed. This
defaults to \code{TRUE}.}

\item{boot}{This is a logical value to be passed to the
\code{\link[mediation]{mediate}} function that determines whether bootstrap
CIs are calculated. This is recommended for GAM models, which are likely to
be used for \code{medmodel} and/or \code{outmodel} in this function, and so
this is set by default to \code{TRUE}.}

\item{...}{These are additional arguments to be passed to the
\code{\link[mediation]{mediate}} function. One important argument is the
number of bootstrap iterations for calculating the bootstrap CIs when
\code{boot = TRUE}. This is controlled with the argument \code{sims}, which
defaults to 1000.}
}
\value{
\code{MediationCurve} returns a data frame containing the output of
the \code{\link[mediation]{mediate}} function iterations. This data frame
includes the following columns:
\item{control.value, treat.value}{These are the 'control' and 'treated'
values for each pair of treatment values.}
\item{X, X.ci.2.5, X.ci.97.5, X.p}{Here 'X' is a placeholder. These are the
estimate of interest, lower 95\% CI, upper 95\% CI, and p-value,
respectively. These columns are included for all the estimates listed below.}
\item{d0, d1, d.avg}{Point estimates for the average causal mediation effect
(ACME) under the control condition, treatment condition, and the average
between the two. Note that if a GAM is used for the outcome model, these
values will all be equal.}
\item{z0, z1, z.avg}{The same as above but for the average direct effect
(ADE).}
\item{tau}{The total effect of the treatment on the outcome.}
\item{n0, n1, n.avg}{The same as for the ACME and ADE, but this is the
proportion mediated (i.e., the ACME divided by the total effect).}
}
\description{
This function extends the functionality of the
\code{\link[mediation]{mediate}} function for causal mediation analysis as
defined by
\href{https://imai.fas.harvard.edu/research/files/BaronKenny.pdf}{Imai et al. 2010}
for modelling linear and nonlinear mediation effects in the case of a
continuous treatment. This function is largely inspired by the
\href{https://causal-curve.readthedocs.io/en/latest/intro.html}{causal-curve}
Python package by Roni Kobrosly.
}
\details{
A separate plotting function is provided as described in
\code{\link[DrewDayRFunctions]{plot.medcurve}}.
}
\examples{
# Linear treatment -> outcome association, nonlinear treatment -> mediator 
#  and mediator -> outcome association

\dontrun{
# Make a simple data frame
set.seed(90)
meddata <- data.frame(treat = rnorm(500))
meddata$medhat <- -0.2* meddata$treat^2
set.seed(101)
meddata$med <- meddata$medhat + rnorm(500)
meddata$outhat <- meddata$treat + -0.2* meddata$med^2
set.seed(11)
meddata$out <- meddata$outhat + rnorm(500)

#Run the mediation and outcome models
medmodel <- mgcv::gam(med ~ s(treat), data = meddata)
outmodel <- mgcv::gam(out ~ treat + s(med), data = meddata)

#Mediation with t0 fixed at the approximate median
medres <- MediationCurve(medmodel, outmodel, "treat", "med", 
treat_values = c(seq(-2, -0.5, 0.5), seq(0.5, 2, 0.5)), 
fixed_t0 = 0)

plot.medcurve(medres)

# Mediation with varying t0 and t1
medres2 <- MediationCurve(medmodel, outmodel, "treat", "med", 
treat_values = seq(-2, 2, 0.5), vary_t0 = T)

plot.medcurve(medres2)
}

}
\references{
Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal
Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December),
pp. 309-334.

Kobrosly, R. W., (2020). causal-curve: A Python Causal Inference Package to
Estimate Causal Dose-Response Curves. Journal of Open Source Software, 5(52),
2523, https://doi.org/10.21105/joss.02523.
}
