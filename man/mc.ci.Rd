\name{mc.ci}
\alias{mc.ci}
\title{
  Diagnostic plots for model choice coverage output
}
\description{
  Plots credible interval estimates for raw model choice output from
  \code{\link{cov.mc}}.  This is used to investigate whether the coverage
  property holds and validate whether diagnostic statistics are acting
  as intended.
}
\usage{
  mc.ci(raw, tol, eps, modname, modtrue, nbins=5, 
bintype=c("interval", "quantile"), bw=FALSE, ...)
}

\arguments{
  \item{raw}{
    The \code{raw} item from the list output by \code{cov.mc}.
  }
  \item{tol}{
    The value of \code{tol} to test.
  }
  \item{eps}{
    The value of \code{eps} to test.  This is used when \code{tol} is
    missing.  One of \code{eps} and \code{tol} must be supplied.
  }
  \item{modname}{
    The name of the model to test.
  }
  \item{modtrue}{
    Vector containing the true models generating the pseudo-observed
    test data.  i.e. \code{modtrue[i]} is the model generating dataset i.    
  }
  \item{nbins}{
    Number of bins to display.
  }
  \item{bintype}{
    How to choose the bins (see Details).
  }
  \item{bw}{
    Whether to produce a black and white image.  Default is FALSE.
    Colour is used to make different bins stand out.
  }
  \item{...}{
    Additional plotting arguments - anything that can be used by \code{plot}.
  }
}

\details{
  This function provides a plot which can be used as an informal test of
  the model choice coverage hypothesis for a particular value of
  \code{eps} or \code{tol} and choice of model.  The plot is more
  flexible than the diagnostics, but not suitable as the basis of a
  formal test.

  For each pseudo-observed data set, the ABC probability that the model
  is \code{modname} is taken from \code{raw}, and the true model is
  taken from \code{modtrue}.  The probabilities are binned into
  \code{nbins} intervals, either of equal length or based on \code{nbins+1}
  equally spaced empirical quantiles.  The function estimates the
  observed probability of \code{modname} within each bin using Bayesian
  inference for a binomial proportion under a uniform prior.  The plot
  shows the mean and 95\% credible interval plotted against predicted
  probabilities.  Informally, the coverage property should be rejected
  if predicted values are too unlikely given the observed values.
}

\references{
  Prangle D., Blum M. G. B., Popovic G., Sisson S. A. (2013) Diagnostic
  tools of approximate Bayesian computation using the coverage
  property. \emph{(preprint)} \url{arxiv.org/abs/1301.3166}
}

\author{
Dennis Prangle
}

\seealso{
  \code{\link{cov.mc}} to produce the input for this function
}

\examples{
  ##The examples below are chosen to run relatively quickly (<5 mins)
  ##and do not represent recommended tuning choices.
  \dontrun{
  index <- sample(1:2, 1E4, replace=TRUE)
  sumstat <- ifelse(index==1, rnorm(1E4,0,1), rnorm(1E4,0,rexp(1E4,1)))
  sumstat <- data.frame(ss=sumstat)
  covdiag <- cov.mc(index=index, sumstat=sumstat, testsets=1:100, 
  tol=seq(0.1,1,by=0.1), diagnostics=c("freq"))
  mc.ci(covdiag$raw, tol=0.5, modname=1, modtrue=index[1:100])
  }
}

\keyword{ htest }
\keyword{ hplot }
