\name{nn.ent}
\alias{nn.ent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Works out entropy of a sample.}
\description{
The function computes the k nearest neighbour sample entropy.}
\usage{
nn.ent(th, k=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{th}{The sample from which to compute the entropy.}
  \item{k}{The order (number of neighbours) of the sample entropy calculation.}
}
\details{The sample entropy gives a measure of information in a (posterior) sample, or lack of it.}
\value{ The k nearest neighbour entropy from the sample.}
\references{
Singh, H. et al. (2003) Nearest neighbor estimates of entropy. \emph{Am. J. Math. Man. Sci.},\bold{23}, 301--321.\cr\cr
Shannon, C. E. and Weaver, W. (1948) A mathematical theory of communication. \emph{Bell Syst. Tech. J.}, \bold{27}, 379--423.
}
\author{Matt Nunes}

\section{Warning}{
For high-dimensional posterior samples, the \code{nn.ent} calculation is quite computationally intensive.
}

\seealso{
	\code{\link{mincrit}}
}
\examples{

# create a dummy sample to calculate an entropy measure:

theta<-rnorm(10000)

nn.ent(theta)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
