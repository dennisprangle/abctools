\name{saABC}
\alias{saABC}
\title{
  Summary statistic construction by semi-automatic ABC
}
\description{
  \code{saABC} fits parameter estimators based on simulated data to be
  used as summary statistics within ABC.  Fitting is by linear
  regression.  Some simple diagnostics are provided for assistance.
}
\usage{
saABC(theta, X, plot = TRUE)
}
\arguments{
  \item{theta}{
    A n x d matrix or data frame of simulated parameter values.
    \code{theta[i,j]} is the ith simulated value of parameter j.
}
  \item{X}{
    A n x p matrix or data frame of simulated data and/or associated
    transformations.  \code{X[i,]} is a vector of the data for parameter
    values \code{theta[i,]}.  A constant term should not be included.
    
}
  \item{plot}{
    When \code{plot==TRUE}, a plot of parameter values against fitted
    values is produced for each parameter as a side-effect.
}
}
\details{
  The semi-automatic ABC method of Fearnhead and Prangle (2012) is as
  follows:

  1) Simulate parameter vectors \eqn{\theta_i}{theta_i} and corresponding data sets \eqn{x_i}{x_i} for i=1,2,...,N.
  
  2) Use the simulations to fit an estimator of each parameter as a
  linear combination of f(x), where f(x) is a vector of
  transformations of x (including a constant term).
  
  3) Run ABC using these simulations.

  The \code{saABC} function automates step 2 of this process.  The user
  must supply simulated parameter values \code{theta} and corresponding
  f(x) values \code{x} (n.b. excluding the constant term).  The function
  returns weights for the linear combinations which can easily be used
  for step 3.  In particular, fitted weights are returned as a matrix
  of weights for the columns of \code{x} and a vector of constants.  The
  vector can usually be discarded, as it is not needed to find
  differences between summary statistics.

  The function also returns BIC values for each parameter so that the
  user can judge the quality of the fits, and compare different choices
  of f(x).  Diagnostic plots of supplied parameter values against fitted
  values are also optionally provided.  These are useful for exploratory
  purposes when there are a small number of parameters, but provide less
  protection from overfitting than BIC values.
}
\value{
  \item{B0 }{Vector of constant terms from fitted regressions.}
  \item{B }{Matrix of weights from fitted regressions.}
  \item{BICs }{Vector of BIC values for each fitted regression.}
}
\references{
Blum, M. G. B, Nunes, M. A., Prangle, D. and Sisson, S. A. (2013) A
comparative review of dimension reduction methods in approximate
Bayesian computation. \emph{Stat. Sci.} \bold{28}, Issue 2, 189--208.\cr\cr
Fearnhead, P. and Prangle, D. (2012) Constructing summary statistics for approximate Bayesian computation: 
semi-automatic approximate Bayesian
computation. \emph{J. R. Stat. Soc. B} \bold{74}, Part 3, 1--28.\cr\cr
Nunes, M. A. and Prangle, D. (2016) abctools: an R package for tuning
approximate Bayesian computation analyses. \emph{The R Journal}
\bold{7}, Issue 2, 189--205.
}
\author{
Dennis Prangle
}

\examples{
set.seed(1)
theta <- matrix(runif(2E3),ncol=2)
colnames(theta) <- c("Mean", "Variance")

X <- replicate(5, rnorm(1E3, theta[,1], theta[,2]))

saABC(theta, X)$BICs
saABC(theta, cbind(X, X^2))$BICs ##Variance parameter estimated better 
}
