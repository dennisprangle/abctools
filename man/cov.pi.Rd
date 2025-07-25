\name{cov.pi}
\alias{cov.pi}
\alias{cov.mc}
\alias{covstats.pi}
\alias{covstats.mc}
\title{
Coverage property diagnostics
}
\description{
  These functions produce diagnostic statistics for an ABC analysis to
  judge when the tolerance level is small enough to produce roughly
  no approximation error.  This is done by running analyses for many
  pseudo-observed test data sets and assessing whether the results satisfy the "coverage
  property" (roughly speaking: credible intervals have the claimed
  coverage levels).
}
\usage{
cov.pi(param, sumstat, testsets, tol, eps, diagnostics = c(), 
multicore = FALSE, cores, method = "rejection", nacc.min=20, ...)

cov.mc(index, sumstat, testsets, tol, eps, diagnostics = c(), 
multicore = FALSE, cores, method = "rejection", nacc.min=20, ...)

covstats.pi(raw, diagnostics = c("KS", "CGR"), nacc.min = 20)

covstats.mc(raw, index, diagnostics = c("freq", "loglik.binary", 
"loglik.multi"), nacc.min = 20) 
}
\arguments{
  \item{param}{
    A data frame of parameter values.  It must have the same number of
    rows as \code{sumstat} and contain numeric values only.
}
  \item{index}{
    A vector of model indices. Any value which can be converted to
    factor is ok (e.g. character or numeric entries).  It must have the
    same length as \code{nrow(sumstat)}.
}
  \item{sumstat}{
    A data frame of summary statistic values whose the ith row has been
    simulated using \code{param[i,]} or \code{index[i]}.
}
  \item{testsets}{
    A numeric vector giving the rows of \code{sumstat} to be used as
    pseudo-observed data to test the coverage property.
}
  \item{tol}{
    A vector of proportions of ABC acceptances which will be
    investigated.
  }
  \item{eps}{
    A vector of ABC thresholds which will be investigated.  These are
    used when \code{tol} is missing.  One of \code{eps} and \code{tol}
    must be supplied.
}
  \item{diagnostics}{
    A character vector containing diagnostics to be calculated.
    Allowable values for parameter inference are "KL" (Kullback-Leibler
    based test) or "CGR" (Cook, Gelman and Rubin test). Allowable values
    for model choice are "freq" (a separate frequency-based test for
    each model), "loglik.binary" (a separate \cr log-likelihood based test
    for each model) or "loglik.multi" (single log-likelihood based test).  If
    \code{diagnostics} is empty only raw results will be returned.
}
  \item{multicore}{
    Whether to use the \code{\link{parallel}} package to perform
    analyses of test datasets in parallel.
}
  \item{cores}{
    Number of cores to use when \code{multicore==TRUE}.
}
  \item{method}{
    Method used for ABC analysis.  The default is "rejection".  For
    alternatives see \code{\link[abc]{abc}} (parameter inference) or
    \code{\link[abc]{postpr}} (model choice).
}
  \item{nacc.min}{
    Minimum number of ABC acceptances required to compute diagnostics.
    See Values for details of how this is used.
}
  \item{\dots}{
    Extra arguments to be supplied to the function performing abc
    analysis i.e. \code{\link[abc]{abc}} (parameter inference) or
    \code{\link[abc]{postpr}} (model choice).
}
    \item{raw}{
    Raw output component from \code{cov.pi} or \code{cov.mc} for which
    diagnostics are to be calculated.
}
}

\details{
  These functions are intended to be applied as follows (i) random
  models/parameters are generated, data sets simulated for each and
  summary statistics calculated (ii) these are input to \code{cov.pi}
  (parameter inference) or \code{cov.mc} (model choice) which outputs raw
  results and diagnostics (see below) (iii) the output can be passed to
  \code{covstats.pi} or \code{covstats.mc} if further diagnostics are
  required later (or to find diagnostics for a subset of the pseudo-observed
  data).

  The \code{cov.pi} and \code{cov.mc} functions operate by performing
  many ABC analyses.  The user specifies which datasets amongst those
  simulated will used as pseudo-observed "test" data to be analysed.
  The results of each analysis are compared to the known
  model/parameters which produced the data to see whether they are
  consistent in a particular sense (i.e. if the coverage property is
  satisfied). Various diagnostics are provided to judge this easily, and
  determine what happens as the ABC threshold is varied.  Raw results
  are also returned which can be investigated in more detail.
  
  All ABC analyses use a rejection sampling algorithm implemented by the
  \code{abc} package.  The user may specify regression
  post-processing as part of this analysis.
}
\value{
  Output of \code{cov.pi} or \code{cov.mc} is a list of two data frames,
  \code{raw} and \code{diag}.  The \code{covstats.pi} and
  \code{covstats.mc} functions just return the latter data frame.

  For parameter inference, \code{raw} contains estimated cdfs (referred
  to as p0 estimates in Prangle et al 2013) of the true parameter values
  for each input configuration (i.e. for every tol/eps value at every
  test dataset).  \code{diag} is a data frame of tol/eps value,
  parameter name, diagnostic name and p-value.  Here the p-value relates
  to the test statistic used as a diagnostic.  It is NA if any analyses
  had fewer than \code{nacc.min} acceptances (Diagnostics based on a
  small number of acceptances can be misleading.)

  For model choice, \code{raw} contains estimated model weights for each input
  configuration, and \code{diag} is a data frame of tol/eps value, model,
  diagnostic name and p-value (NA under the same conditions as before.)

  In both cases, \code{raw} also reports the number of acceptances.  Note that
  \code{raw} contains p0 estimates/weights of NA if regression correction is
  requested but there are too few acceptances to compute it.
}
\references{
  Nunes, M. A. and Prangle, D. (2016) abctools: an R package for tuning
  approximate Bayesian computation analyses. \emph{The R Journal}
  \bold{7}, Issue 2, 189--205.
  
  Prangle D., Blum M. G. B., Popovic G., Sisson S. A. (2014) Diagnostic
  tools of approximate Bayesian computation using the coverage
  property. \emph{Australian and New Zealand Journal of Statistics}
  \bold{56}, Issue 4, 309--329.
}
\author{
Dennis Prangle
}

\seealso{
  \code{\link{mc.ci}} for a diagnostic plot of raw model choice results
  
  \code{\link[abc]{abc}} and \code{\link[abc]{postpr}} to perform ABC for a given dataset
}

\examples{
  ##The examples below are chosen to run relatively quickly (<5 mins)
  ##and do not represent recommended tuning choices.
  \dontrun{
  data(musigma2)
  library(ggplot2)
  ##Parameter inference example
  parameters <- data.frame(par.sim)
  sumstats <- data.frame(stat.sim)
  covdiag <- cov.pi(param=parameters, sumstat=sumstats, testsets=1:100, 
  tol=seq(0.1,1,by=0.1), diagnostics=c("KS"))

  #Plot of diagnostic results
  qplot(x=tol, y=pvalue, facets=.~parameter, data=covdiag$diag) 
  #Plot of raw results for tol=0.5
  qplot(x=mu, data=subset(covdiag$raw, tol==0.5)) 
  #Plot of raw results for tol=0.5
  qplot(x=sigma2, data=subset(covdiag$raw, tol==0.5)) 

  #Compute CGR statistic and plot
  cgrout <- covstats.pi(covdiag$raw, diagnostics="CGR") 
  qplot(x=tol, y=pvalue, facets=.~parameter, data=cgrout) 

  ##Model choice example, based on simple simulated data
  index <- sample(1:2, 1E4, replace=TRUE)
  sumstat <- ifelse(index==1, rnorm(1E4,0,1), rnorm(1E4,0,rexp(1E4,1)))
  sumstat <- data.frame(ss=sumstat)
  covdiag <- cov.mc(index=index, sumstat=sumstat, testsets=1:100, 
  tol=seq(0.1,1,by=0.1), diagnostics=c("freq"))
  qplot(x=tol, y=pvalue, data=covdiag$diag)
  llout <- covstats.mc(covdiag$raw, index=index, 
  diagnostics="loglik.binary")
  qplot(x=tol, y=pvalue, data=llout)
  mc.ci(covdiag$raw, tol=0.5, modname=1, modtrue=index[1:200])
  }
}

\keyword{ htest }
