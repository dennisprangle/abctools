\name{stage2}
\alias{stage2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{stage2}
\description{
Summary statistics selection for ABC inference using estimated posterior error. 
}
\usage{
stage2(obs, param, sumstats, obspar = NULL, init.best, dsets = 100, 
sumsubs = 1:ncol(sumstats), limit = length(sumsubs), do.only=NULL, 
do.err = FALSE, final.dens = FALSE,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{observed summary statistics.}
  \item{param}{matrix of simulated model parameter values.}
  \item{sumstats}{matrix of simulated summary statistics.}
  \item{obspar}{optional observed parameters (for use to assess simulation performance).}
  \item{init.best}{an initial estimate of the best summary statistics subset.  Can be either an index into the summaries combination table (see \code{\link{combmat}}) or a vector of indices into \code{1:nstats}.  See details.}
  \item{dsets}{the number of simulated datasets to treat as observed when estimating the posterior error.  See details.}
  \item{sumsubs}{an optional index into the summary statistics to limit summary selection to a specific subset of summaries.}
  \item{limit}{an optional integer indicating whether to limit summary selection to subsets of a maximum size.}
  \item{do.only}{an optional index into the summary statistics combination table.  Can be used to limit entropy calculations to certain summary statistics subsets only.}
  \item{do.err}{a boolean value indicating whether the simulation error should be returned.  Note: if \code{do.err=TRUE}, \code{obspar} must be supplied.}
  \item{final.dens}{a boolean value indicating whether the posterior sample should be returned.}
%  \item{verbose}{a boolean value indicating whether informative statements should be printed to screen.}
  \item{\dots}{any other optional arguments to the ABC inference procedure (e.g. arguments to the \code{abc} function).}
}
\details{The function uses the \code{init.best} set of summaries to determine the \code{dsets} simulated datasets which are closest (in Euclidean norm) to the observed dataset.  Since the model parameters generating the summary statistics are known for these simulated datasets, for each candidate subset of summary statistics, we can compute the error under ABC inference for each of these datasets.  The best subset of summary statistics is that which minimizes this (average) error over all \code{dsets} datasets.}
\value{
  A list with the following components:
  \item{best}{the best subset of statistics.}
  \item{closest}{the indices of the \code{dsets} simulated datasets closest to the oberved dataset as measured by the \code{init.best} subset of summaries.}
  \item{err}{simulation error (if \code{obspar} is supplied and \code{do.err=TRUE}).}
  \item{order}{the subsets considered during the algorithm (same as the input \code{do.only}.}
  \item{post.sample}{an array of dimension \code{nacc x npar x ndatasets} giving the posterior sample for each observed dataset.  Not returned if \code{final.dens=FALSE}.}
  \item{sumsubs}{an index into the subsets considered during the algorithm.}
}
\references{
Blum, M. G. B, Nunes, M. A., Prangle, D. and Sisson, S. A. (2013) A
comparative review of dimension reduction methods in approximate
Bayesian computation. \emph{Stat. Sci.} \bold{28}, Issue 2, 189--208.\cr\cr
Nunes, M. A. and Balding, D. J. (2010) On Optimal Selection of Summary
Statistics for Approximate Bayesian Computation.
\emph{Stat. Appl. Gen. Mol. Biol.} \bold{9}, Iss. 1, Art. 34.\cr\cr
Nunes, M. A. and Prangle, D. (2016) abctools: an R package for tuning
approximate Bayesian computation analyses. \emph{The R Journal}
\bold{7}, Issue 2, 189--205.
}
\author{Matt Nunes}
%\note{
%%  ~~further notes~~
%}

\section{Warning}{This function is computationally intensive due to its cyclic ABC inference procedure.}

%\seealso{
%% ~~objects to See Also as \code{\link{}}, ~~~
%}
\examples{

# load example data:

data(coal)
data(coalobs)

param<-coal[,2]
simstats<-coal[,5:8]

# use matrix below just in case to preserve dimensions.

obsstats<-matrix(coalobs[1,5:8],nrow=1)
obsparam<-matrix(coalobs[1,1])

\dontrun{
tmp<-stage2(obsstats, param, simstats, init.bes=c(1,3), dsets = 10) 
tmp$best
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
