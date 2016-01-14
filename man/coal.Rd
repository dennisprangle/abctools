\name{coal}
\alias{coal}
\alias{coalobs}
\docType{data}
\title{Examples of coalescent data}
\description{Data generated from a coalescent model for genetic variation.}
\usage{data(coal)}
\format{Matrices of parameters (first 2 columns) and summary statistics
  (remaining columns) from a coalescent model.}
\details{
"coal" contains 100,000 simulated datasets to use for an ABC
analysis. "coalobs" contains 100 which can each be used as observed
datasets. This allows methods to be tested on many different
observations.
  
The parameters are the scaled mutation rate (theta) and scaled
recombination rate (rho).
  
The summary statistics are the number of segregating sites (segsites),
independent Uniform[0,25] random noise (unif),
the pairwise mean number of nucleotidic differences (meandiff),
the mean R^2 across pairs separated by < 10\% of the simulated genomic regions (R2),
the number of distinct haplotypes (nhap),
the frequency of the most common haplotype (fhap),
and the number of singleton haplotypes (shap).

To simulate a dataset 5,001 basepair DNA sequences for 50 individuals
are generated from the coalescent model, with recombination, under the
infinite-sites mutation model, using the software ms (Hudson 2002)

Data of this form has been analysed in several ABC papers.  See the
references for more details.
}
\references{
Blum, M. G. B., M. A. Nunes, D. Prangle and S. A. Sisson `A comparative
review of dimension reduction methods in approximate Bayesian
computation' Statistical Applications in Genetics and
Molecular Biology 13 (2014)

Hudson, R. R. `Generating samples under a Wright-Fisher neutral
model of genetic variation' Bioinformatics 18 (2002)

Joyce, P. and P. Marjoram `Approximately sufficient statistics
and Bayesian computation' Statistical Applications in Genetics and
Molecular Biology 7 (2008)

Nunes, M. A. and D. J. Balding `On optimal selection of summary
statistics for approximate Bayesian computation' Statistical
Applications in Genetics and Molecular Biology 9 (2010)

Nunes, M. A. and Prangle, D. abctools: an R package for tuning
approximate Bayesian computation analyses. The R Journal 7 (2016)
}

\keyword{datasets}
