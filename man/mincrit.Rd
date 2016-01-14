\name{mincrit}
\alias{mincrit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary statistics selection by minimizing a posterior sample measure. 
}
\description{
The function cycles through all possible subsets of summary statistics and computes a criterion from the posterior sample.  The subset which achieves the minimum is chosen as the most informative subset.
}
\usage{
mincrit(obs, param, sumstats, obspar = NULL, abcmethod = abc,
crit = nn.ent, sumsubs = 1:ncol(sumstats), limit=length(sumsubs), 
do.only = NULL, verbose = TRUE, do.crit = TRUE, do.err = FALSE, 
final.dens = FALSE, errfn = rsse, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{(matrix of) observed summary statistics.}
  \item{param}{matrix of simulated model parameter values.}
  \item{sumstats}{matrix of simulated summary statistics.}
  \item{obspar}{optional observed parameters (for use to assess simulation performance).}
  \item{abcmethod}{a function to perform ABC inference, e.g. the \code{abc} function from package \emph{abc}.}
  \item{crit}{a function to minimize to measure information from a posterior sample, e.g. \code{nn.ent}.}
  \item{sumsubs}{an optional index into the summary statistics to limit summary selection to a specific subset of summaries.}
  \item{limit}{an optional integer indicating whether to limit summary selection to subsets of a maximum size.}
  \item{do.only}{an optional index into the summary statistics combination table.  Can be used to limit entropy calculations to certain summary statistics subsets only.}
  \item{verbose}{a boolean value indicating whether informative statements should be printed to screen.}
%  \item{k}{the order of the nearest neighbour entropy if \code{do.ent=TRUE}.}
  \item{do.crit}{a boolean value indicating whether the measure on the posterior sample should be returned.}
  \item{do.err}{a boolean value indicating whether the simulation error should be returned.  Note: if \code{do.err=TRUE}, \code{obspar} must be supplied.}
  \item{final.dens}{a boolean value indicating whether the posterior sample should be returned.}
  \item{errfn}{an error function to assess ABC inference performance.}
  \item{\dots}{any other optional arguments to the ABC inference procedure (e.g. arguments to the \code{abc} function).}
}
\details{The function uses a criterion (e.g.sample entropy) as a proxy for information in a posterior sample.  The criterion for each possible subset of statistics is computed, and the best subset is judged as the one which minimises this vector of values.}
\value{
  A list with the following components:
  \item{best}{the best subset(s) of statistics.}
  \item{critvals}{the calculated criterion values (if \code{do.crit=TRUE}).}
  \item{err}{simulation error (if \code{obspar} is supplied and \code{do.err=TRUE}).}
  \item{order}{the subsets considered during the algorithm (same as the input \code{do.only}.}
  \item{post.sample}{an array of dimension \code{nacc x npar x ndatasets} giving the posterior sample for each observed dataset.  Not returned if \code{final.dens=FALSE}.}
  \item{sumsubs}{an index into the subsets considered during the algorithm.}
}
\references{
Nunes, M. A. and Balding, D. J. (2010) On Optimal Selection of Summary
Statistics for Approximate Bayesian Computation.
\emph{Stat. Appl. Gen. Mol. Biol.} \bold{9}, Iss. 1, Art. 34.

Nunes, M. A. and Prangle, D. (2016) abctools: an R package for tuning
approximate Bayesian computation analyses. \emph{The R Journal}
\bold{7}, Issue 2, 189--205.
}
\author{Matt Nunes}
%\note{
%%  ~~further notes~~
%}

\section{Warning}{These functions are computationally intensive due to the cyclic ABC inference procedure.}

\seealso{
\code{\link{nn.ent}}
}
\examples{

# load example data:

data(coal)
data(coalobs)

param<-coal[,2]
simstats<-coal[,4:6]

# use matrix below just in case to preserve dimensions.

obsstats<-matrix(coalobs[1,4:6],nrow=1)
obsparam<-matrix(coalobs[1,1])

# example of entropy minimization algorithm:

tmp <-mincrit(obsstats, param, simstats, tol=.01, method="rejection", 
do.crit=TRUE)

tmp$critvals

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
