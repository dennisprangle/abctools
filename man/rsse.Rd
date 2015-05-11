\name{rsse}
\alias{rsse}
\alias{sse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation error measures.}
\description{
Computes error measures between two sets of data.
}
\usage{
rsse(a, b, v = 1)
sse(a, b, v = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Dataset 1, of dimension (n x p)}
  \item{b}{Dataset 2: of dimension (1 x p) or a vector.}
  \item{v}{an optional factor to normalise the data before computation of the RSSE.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{The RSSE between dataset \code{a} and \code{b}.}
\author{Matt Nunes}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

a<-matrix(rnorm(1000),ncol=2)
b<-runif(2)

rsse(a,b)

}

\keyword{manip}
