\name{semiauto.abc}
\alias{semiauto.abc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs semi-automatic ABC based on summary statistics regression.}
\description{
Performs semi-automatic ABC based on summary statistics regression.
}
\usage{
semiauto.abc(obs, param, sumstats, obspar=NULL, abcmethod = abc, 
saprop = 0.5, abcprop = 0.5, overlap = FALSE, satr = list(), 
plot = FALSE, verbose = TRUE, do.err = FALSE, final.dens = FALSE, 
errfn = rsse, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{(matrix of) observed summary statistics.}
  \item{param}{matrix of simulated model parameter values.}
  \item{sumstats}{matrix of simulated summary statistics.}
  \item{obspar}{optional observed parameters (for use to assess simulation performance).}
  \item{abcmethod}{a function to perform ABC inference, e.g. the \code{abc} function from package \emph{abc}.}
  \item{saprop}{a proportion, denoting the proportion of simulated datasets with which to perform semi-automatic ABC regression.}
  \item{abcprop}{a proportion, denoting the proportion of simulated datasets with which to perform ABC using \code{abcmethod}.}
  \item{overlap}{a boolean value indicating whether the simulated datasets specified by \code{saprop} and \code{abcprop} are disjoint (\code{overlap=FALSE}) or not.}
  \item{satr}{a list of functions indicating transformations of the summary statistics \code{sumstats}.  These must be *suitable* functions, and must each return a vector, matrix or array with the number of elements being a multiple of the rows of \code{sumstats}.  See details and examples sections for more information}
  \item{plot}{When plot==TRUE, a plot of parameter values against fitted values is produced for each parameter as a side-effect.  This is most useful when the number of parameters is reasonably small.}
  \item{verbose}{a boolean value indicating whether informative statements should be printed to screen.}
  \item{do.err}{a boolean value indicating whether the simulation error should be returned.  Note: if \code{do.err=TRUE}, \code{obspar} must be supplied.}
  \item{final.dens}{a boolean value indicating whether the posterior sample should be returned.}
  \item{errfn}{an error function to assess ABC inference performance.}
  \item{\dots}{any other optional arguments to the ABC inference procedure (e.g. arguments to the \code{abc} function).}
}
\details{This function is essentially a wrapper for \code{saABC}.  See the details section of \code{saABC} for more details on the implementation.  The argument \code{satr} can be almost anything sensible in \code{function} form, see Examples section for example specifications.}
\value{ A list with the following components:
\item{err}{simulation error (if \code{obspar} is supplied and \code{do.err=TRUE}).}
\item{post.sample}{an array of dimension \code{nacc x npar x ndatasets} giving the posterior sample for each observed dataset.  Not returned if \code{final.dens=FALSE}.}
\item{sainfo}{A list with the following information about the semi-automatic ABC run:\cr \code{saprop}, \code{abcprop} ,\code{overlap}, \code{satr}.  See arguments for more details.}
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
\author{Matt Nunes and Dennis Prangle}

\seealso{
\code{\link{saABC}}, \code{\link{selectsumm}}
}

\section{Warning}{The argument \code{satr} must be supplied with valid functions.  Whilst there are checks, these are minimal, since doing sophisticated checks is quite difficult.

}

\examples{

data(coal)
data(coalobs)

param<-coal[,2]
simstats<-coal[,4:6]

# use matrix below just in case to preserve dimensions.

obsstats<-matrix(coalobs[1,4:6],nrow=1)
obsparam<-matrix(coalobs[1,1])

# perform semi-automatic ABC with summary statistics defined by 
# X, X^2,X^3,X^4:
# other alternative specifications for this could be:
# list(function(x){ cbind(x,x^2,x^3,x^4) })
# list(as.function(alist(x=,cbind(x,x^2,x^3)))) etc

tmp<-semiauto.abc(obsstats, param, simstats,tol=.01,method="rejection",
satr=list(function(x){outer(x,Y=1:4,"^")}))

tmp$sa.info

# both these functions may be problematic: 

\dontrun{
 tmp<-semiauto.abc(obsstats, param, simstats,tol=.01,method="rejection",
satr=list(unique,sum))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
