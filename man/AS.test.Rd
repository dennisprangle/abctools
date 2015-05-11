\name{AS.test}
\alias{AS.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for relative approximate sufficiency between two posterior samples.
}
\description{The function tests to determine adding a (set of) statistics is informative in ABC inference.}
\usage{
AS.test(grid = 10, x1, x2, supp=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid}{the number of bins into which to divide the posterior sample for the approximate sufficiency calculation.}
  \item{x1}{the posterior sample using the first set of summary statistics.}
  \item{x2}{the posterior sample using the second (alternative) set of summary statistics.}
%  \item{prior}{a single letter specifying the prior which generated the parameter values, e.g ."u" (uniform) or "n" (normal).}
  \item{supp}{the "support" of the prior (e.g. uniform bounds).}
}
\details{After dividing each posterior sample into a number of bins (specified by \code{grid}), the function computes the ratio of the posterior densities.  This is seen as a measure of information added (sufficiency) by using the alternative posterior sample instead of the first posterior sample.  If the ratio exceeds a particular threshold (a number of standard deviations away from the expected counts in each bin), then the alternative set of summaries is seen as being more informative.}
\value{
\item{extreme}{a boolean value indicating whether the alternative posterior sample is more informative than the first (i.e. the extra summary statistics add information.}
}
\references{
Blum, M. G. B, Nunes, M. A., Prangle, D. and Sisson, S. A. (2013) A
comparative review of dimension reduction methods in approximate
Bayesian computation. \emph{Stat. Sci.} \bold{28}, Issue 2, 189--208.\cr\cr
Joyce, P. and P. Marjoram (2008) Approximately sufficient statistics and Bayesian computation. \emph{Stat. Appl. Gen. Mol. Biol.} \bold{7} Article 26.
}
\author{Matt Nunes}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{AS.select}}}
\examples{

#create two fake posterior samples:

x1<-runif(10000)
x2<-rnorm(10000)

AS.test(x1=x1,x2=x2,supp=range(x2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
