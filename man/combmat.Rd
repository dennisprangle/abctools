\name{combmat}
\Rdversion{1.1}
\alias{combmat}
\alias{fillcomb}
\alias{combtable}
\title{table of combinations}
\description{This function creates a table of binary masks representing combinations of statistics.}
\usage{combmat(n,limit = NULL)}
\arguments{
  \item{n}{number of statistics}
%  \item{order}{logical value of whether to order the masks in order of subset size.}
  \item{limit}{an optional (integer) value indicating whether to limit the table to subsets up to a certain size.}
}
%\details{
%}
\value{
\item{m}{The matrix of binary masks.}
}
%\references{
%}
\author{
Matt Nunes
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#
# Find all binary masks of a set of statistics {C1,C2,C3,C4},
# listing all singlets, pairs, triples and then the whole set:

combmat(4,TRUE)

}

\keyword{manip}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
