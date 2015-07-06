\name{AS.select}
\alias{AS.select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary statistics selection using approximate sufficiency.
}
\description{
This function uses approximate sufficiency to assess subsets of summary statistics for ABC inference.
}
\usage{
AS.select(obs, param, sumstats, obspar=NULL, abcmethod=abc, 
grid=10, inturn=TRUE, limit=ncol(sumstats), allow.none=FALSE, 
do.err=FALSE, final.dens=FALSE, errfn=rsse, trace=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{(matrix of) observed summary statistics.}
  \item{param}{matrix of simulated model parameter values.}
  \item{sumstats}{matrix of simulated summary statistics.}
  \item{obspar}{optional observed parameters (for use to assess simulation performance).}
  \item{abcmethod}{a function to perform ABC inference, e.g. the \code{abc} function from package \emph{abc}.}
  \item{grid}{the number of bins into which to divide the posterior sample for the approximate sufficiency calculation.}
  \item{inturn}{a boolean value indicating whether "bad" statistics should be dropped and tested sequentially (\code{inturn=TRUE}) or all at the end (\code{inturn=FALSE}).}
  \item{limit}{an optional integer indicating whether to limit summary selection to subsets of a maximum size.}
  \item{allow.none}{a boolean values indicating whether an empty subset of statistics is considered in the selection procedure.}
  \item{do.err}{a boolean value indicating whether the simulation error should be returned.  Note: if \code{do.err=TRUE}, \code{obspar} must be supplied.}
%  \item{sumsubs}{an optional index into the summary statistics to limit summary selection to a specific subset of summaries.}
  \item{final.dens}{a boolean value indicating whether the posterior sample should be returned.}
  \item{errfn}{an error function to assess ABC inference performance.}
  \item{trace}{whether to show progress messages.}
  \item{\dots}{any other optional arguments to the ABC inference procedure (e.g. arguments to the \code{abc} function).}
}
\details{The summary selection procedure works by sequentially testing randomly chosen statistics for inclusion, using the ratio of ABC posterior samples to determine whether a statistic is added.  Since adding a statistic may result in a suboptimal subset of summaries, the included statistics are then individually dropped and retested, to determine whether a smaller subset of statistics is equally / more informative than the accepted set of statistics.}
\value{A list with the following components:
  \item{best}{the final subset of included statistics.}
  \item{err}{simulation error (if \code{obspar} is supplied and \code{do.err=TRUE}).}
  \item{post.sample}{an array of dimension \code{nacc x npar x ndatasets} giving the posterior sample for each observed dataset.  Not returned if \code{final.dens=FALSE}.}
}

\references{
Blum, M. G. B, Nunes, M. A., Prangle, D. and Sisson, S. A. (2013) A
comparative review of dimension reduction methods in approximate
Bayesian computation. \emph{Stat. Sci.} \bold{28}, Issue 2, 189--208.\cr\cr
Joyce, P. and P. Marjoram (2008) Approximately sufficient statistics and Bayesian computation. \emph{Stat. Appl. Gen. Mol. Biol.} \bold{7} Article 26.
}
\author{Matt Nunes}

\note{The approximate sufficiency techniques described here are only suitable for single parameters only.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AS.test}} 
}
\examples{

# load example data:

data(coal)
data(coalobs)

param<-coal[,2]
simstats<-coal[,4:6]

# use matrix below just in case to preserve dimensions.

obsstats<-matrix(coalobs[1,4:6],nrow=1)

# example of AS.select:

\dontrun{
tmp <-AS.select(obsstats, param, simstats, tol=.1, method="neuralnet",
nument=5, allow.none=FALSE, inturn=TRUE)

tmp$best
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
