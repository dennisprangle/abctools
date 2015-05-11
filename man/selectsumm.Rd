\name{selectsumm}
\alias{selectsumm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function for selecting summary statistics in ABC inference.
}
\description{The function implements functions which implement summary statistics selection methods.}
\usage{
selectsumm(obs, param, sumstats, obspar=NULL, ssmethod = mincrit, 
verbose = TRUE, final.dens = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{(matrix of) observed summary statistics.}
  \item{param}{matrix of simulated model parameter values.}
  \item{sumstats}{matrix of simulated summary statistics.}
  \item{obspar}{optional observed parameters (for use to assess simulation performance).}
  \item{ssmethod}{a function to perform summary statistics selection.  Current methods are \cr
	\code{AS.select} and \code{mincrit}.}
  \item{verbose}{a boolean value indicating whether informative statements should be printed to screen.}
  \item{final.dens}{a boolean value indicating whether the posterior sample should be returned.}
  \item{\dots}{any other optional arguments to the summary selection procedure.}
}

\details{
The function is essentially a wrapper for more specific summary selection methods, and is designed to be flexible for future additions and minimization criteria.  See the help files for each summary selection method for more details.
}
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
Blum, M. G. B, Nunes, M. A., Prangle, D. and Sisson, S. A. (2013) A
comparative review of dimension reduction methods in approximate
Bayesian computation. \emph{Stat. Sci.} \bold{28}, Issue 2, 189--208.\cr\cr
Fearnhead, P. and Prangle, D. (2012) Constructing summary statistics for approximate Bayesian computation: semi-automatic approximate Bayesian computation. \emph{J. R. Stat. Soc. B} \bold{74}, Part 3, 1--28.\cr\cr
Joyce, P. and P. Marjoram (2008) Approximately sufficient statistics and Bayesian computation. \emph{Stat. Appl. Gen. Mol. Biol.} \bold{7} Article 26.\cr\cr
Nunes, M. A. and Balding, D. J. (2010) On Optimal Selection of Summary Statistics for Approximate Bayesian Computation.  \emph{Stat. Appl. Gen. Mol. Biol.} \bold{9}, Iss. 1, Art. 34. \cr\cr
}
\author{Matt Nunes}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mincrit}}, \code{\link{AS.select}} 
}

\examples{

# load example data:

data(coal)
data(coalobs)

param<-coal[,2]
simstats<-coal[,4:6]

# use matrix below just in case to preserve dimensions.

obsstats<-matrix(coalobs[1,4:6],nrow=1)

tmp<-selectsumm(obsstats, param, simstats, ssmethod =AS.select, tol =.1,
method = "rejection", allow.none = FALSE, inturn = TRUE, hcorr = TRUE)

tmp$best

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
