\name{recalibrationABC}
\alias{recalibrationABC}
\title{
  ABC inference with a recalibration adjustment
}
\description{
  This function post-processes ABC output with the aim of calibrating
  its credible intervals to have the correct probabilities. The approach
  can be thought of as extending coverage tests to correct deviations
  from the desired posterior. See the reference for details.
}
\usage{
  recalibrationABC(target, param, sumstat, eps, tol, method="rejection",
  multicore=FALSE, cores=NULL, abc.p.options=list(method="loclinear"), ...)
}

\arguments{
  \item{target}{
    A vector of the observed summary statistics.
  }
  \item{param}{
    A data frame of the simulated parameter values.
  }
  \item{sumstat}{
    A data frame of the simulated summary statistics.
  }
  \item{eps}{
    The ABC acceptance threshold i.e. max acceptable distance. This or
    \code{tol} must be specified (but not both).
  }
  \item{tol}{
    The ABC acceptance tolerance i.e. proportion of simulations to
    accept. This or \code{eps} must be specified (but not both).    
  }
  \item{method}{
    A character string indicating the type of ABC algorithm to be
    applied. Possible values are "rejection", "loclinear",
    "neuralnet" and "ridge".
  }

  \item{multicore}{
    Whether to use the \code{\link{parallel}} package to perform
    analyses of test datasets in parallel.
  }
  \item{cores}{
    Number of cores to use when \code{multicore==TRUE}.
  }
  \item{abc.p.options}{
    A list of further arguments to be supplied to the \code{abc} command
    when applying regression correction to the p-values based on the
    summary statistics.
  }
  \item{...}{
    Further arguments to be supplied to the \code{abc} command when
    performing inference for each simulated data set. Typically these will
    control regression correction. If omitted no correction is performed.
  }
}

\value{
  A list with the following components is returned.

  \code{sample.abc} is the ordinary ABC output sample (with any
  regression correction requested). The rows represent accepted samples
  and the columns represent the parameters.
  
  \code{sample.recal} is the ABC output sample following
  recalibration. The rows are Monte Carlo draws from an approximation to
  the posterior and the columns represent the parameters.

  \code{sample.regrecal} is the ABC output sample following coverage
  correction including regression correction of p-values. It has a
  similar form to \code{sample.recal}.

  \code{weights} are weights for the ABC output. These apply to all
  types of ABC.

  Each row of \code{pvalues} corresponds to a particular accepted
  dataset. It gives the p-values of the true parameters within the ABC
  sample generated from this data.

  Each row of \code{pvalues.reg} corresponds to \code{pvalues} after
  a regression-adjustment has been performed on them.

  \code{svalues} is a subset of the rows of \code{sumstat} corresponding
  to accepted datasets. These can be used in conjunction with
  \code{pvalues} to perform a recalibration correction manually.
}

\references{
  G. S. Rodrigues, D. Prangle and S. A. Sisson (2017)
  Recalibration: A post-processing method for approximate Bayesian
  computation. \emph{In submission}
}

\author{
  Dennis Prangle and Guilherme Rodrigues
}

\examples{
  data(musigma2)
  P <- data.frame(par.sim)
  S <- data.frame(stat.sim)
  out <- recalibrationABC(target=stat.obs, param=P, sumstat=S, tol=0.3)
  plot(rbind(out$sample.plain[1:500,], out$sample.recal[1:500,]),
  col=c(rep("red",500), rep("blue", 500)))
  ##Red shows plain ABC sample, blue shows recalibrated output
}
