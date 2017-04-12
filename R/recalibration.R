recalibrationABC <- function(target, param, sumstat, eps, tol, method="rejection", multicore=FALSE, cores=NULL, abc.p.options=list(method="loclinear"), ...) {
    ##Check which of tol and eps is used
    got_tol <- !missing(tol)
    got_eps <- !missing(eps)
    if (got_tol && got_eps) {
        stop("Only one of eps and tol should be supplied to recalibrationABC")
    } else if (!got_tol && !got_eps) {
        stop("One of eps or tol must be supplied to recalibrationABC")
    }

    ##Do ABC on observed data
    more.args <- list(...)
    abc.args <- list(target=target, param=param, sumstat=sumstat, tol=tol, method=method)
    if (!"transf" %in% names(more.args)) {
        more.args <- c(more.args, list(transf=rep("none", ncol(param)))) #Avoids abc outputing warning messages
    }
    abcout <- do.call(abc, c(abc.args, more.args))
    if (abc.args$method == "rejection") {
        sample.abc <- abcout$unadj.values
    } else {
        sample.abc <- abcout$adj.values
    }

    if (got_eps) {
        accepted <- which(abcout$dist<=eps)
        tol <- length(accepted)/nrow(param)
    } else {
        accepted <- order(abcout$dist)[1:ceiling(nrow(param)*tol)] #Matches what's used in abc package
    }

    sample.abc <- as.matrix(sample.abc) #(might previous have been a vector)
    ##Get p-values for accepted simulations
    cov.args <- list(param=param, sumstat=sumstat, testsets=accepted,
                     diagnostics=c(), multicore=multicore, cores=cores, type.cdf=2)
    if (got_eps) {
        cov.args <- c(cov.args, eps=eps)
    } else {
        cov.args <- c(cov.args, tol=tol)
    }
    covout.pi <- do.call(cov.pi, c(cov.args, more.args))
    pvalues <- as.matrix(covout.pi$raw[,-(1:3)])
    ##Apply regression correction to p-values
    LB <- matrix(nrow=ncol(param), ncol=2)
    LB[,1] <- 0
    LB[,2] <- 1
    abc.p.args <- list(target=target, param=pvalues, sumstat=sumstat[accepted,], tol=1, transf=rep("logit", ncol(param)), logit.bounds=LB)
    abc.p.args <- c(abc.p.args, abc.p.options)
    pvalues.reg <- do.call(abc, abc.p.args)$adj.values
    ##Store weights
    weights_used <- ("weights" %in% names(abcout))
    if (weights_used) {
        weights <- abcout$weights
        get_quantiles <- function(A, p) {
            Hmisc::wtd.quantile(x=A, weights=weights, probs=p)
        }
    } else {
        weights <- rep(1, nrow(sample.abc))
        get_quantiles <- function(A, p) {
            quantile(A, probs=p, type=8)
        }
    }
    ##Do coverage correction
    sample.recal <- sapply(1:ncol(pvalues), function(i) get_quantiles(sample.abc[, i], pvalues[, i]))
    sample.regrecal <- sapply(1:ncol(pvalues), function(i) get_quantiles(sample.abc[, i], pvalues.reg[, i]))
    colnames(sample.recal) <- colnames(sample.abc)
    colnames(sample.regrecal) <- colnames(sample.abc)
    ##Return results
    list(sample.abc=sample.abc, sample.recal=sample.recal, sample.regrecal=sample.regrecal,
         weights=weights,
         pvalues=pvalues, pvalues.reg=pvalues.reg,
         svalues=as.matrix(sumstat[accepted,,drop=FALSE])
         )
}
