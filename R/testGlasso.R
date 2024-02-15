
#' @name testGlasso
#' @title Test for for the equality of connectivity based on the Graphical lasso estimation.
#' @description This function utilizes Dynamic Connectivity Regression (DCR) algorithm proposed by Cribben el al. (2012) to test the equality of connectivity in two fMRI signals.
#' @param subY1 a sample of size length*dim
#' @param subY2 a sample of size length*dim
#' @param p  Gep(p) distribution controls the size of stationary bootstrap. The mean block length is 1/p
#' @param lambda two selections possible for optimal parameter of lambda. "bic" finds lambda from bic criteria, or user can directly input the penalty value.
#' @param nboot the number of bootstrap sample for pvalue. Default is 100.
#' @param n.cl number of cores in parallel computing. The default is (machine cores - 1)
#' @param bound bound of bic search in "bic" rule. Default is (.001, 1)
#' @param gridTF Utilize a grid search to optimize hyperparameters
#' @return \strong{pval} The empirical p-value for testing the equality of connectivity structure
#' @return \strong{rho} The sequence of penalty paramter based on the combined sample, subY1 and subY2.
#' @return \strong{fit0} Output of glasso for combied sample
#' @return \strong{fit1} Output of glasso for subY1
#' @return \strong{fit2} Output of glasso for subY2
#' @examples test1= testGlasso(testsim$X, testsim$Y, n.cl=1)
#' @export

testGlasso <- function(subY1, subY2, p, lambda = "bic", nboot = 100, n.cl, bound = c(.001, 1), gridTF=FALSE) {
  if (missing(n.cl)) {
    n.cl <- parallel::detectCores(logical = TRUE) - 1
  }
  if (missing(p)) {
    p <- .2
  }

  #library(glasso)
  #library(doParallel)

  ## Data dimension is Tt*q
  ### Auxilary functions
  glasso.bic2 <- function(y, khat, lambda, bound, debiasTF=TRUE) {
    st <- khat[1]
    end <- khat[2]
    id <- seq(st, end, by = 1)
    tj <- length(id)
    cy <- stats::cov(y[id, ])
    n <- nrow(y);

    
    if (lambda == "bic") {
      bic.rho <- function(rho, cy, tj) {
        B <- 0
        out <- glasso::glasso(cy, rho, maxit=400)
        ell1 <- 1 * (out$wi != 0)
        diag(ell1) <- rep(0, nrow(out$wi))
        ## tj is replacec by tj-1 by comparing result in the below.
        dwi <- sum(log(svd(out$wi)$d));
        B <- sum(diag((tj - 1) * cy %*% out$wi)) - tj * dwi + sum(ell1) * (log(tj))/2;
        return(B)
      }

      if(gridTF){
        ## Grid search?? Most stable..
        sqrho = seq(bound[1],  bound[2], length=20);
        est = lapply(sqrho, bic.rho, cy = cy, tj=tj);
        id = which.min(unlist(est)); rho = sqrho[id];
      } else{
        est <- stats::optim(par = .2, fn = bic.rho, cy = cy, tj = tj, method = "L-BFGS-B", lower = bound[1], upper = bound[2], control = list(maxit = 400, factr = 1e2, lmm = 4))
        rho <- est$par
      }

    } else {
      rho <- lambda
    } ## Fixed constant

    out <- glasso::glasso(cy, rho, maxit=400)
    idz <- which(out$wi == 0)

    ## If close to singular, do not debias
    if(det(out$w) < 1e-6){ debiasTF=FALSE; }

    ## Re-estimate precision matrix by imposing zero constraints
    ## It only valid when there is zero estimates in out()
    if ((length(idz) > 0)*debiasTF) {
      q <- dim(out$wi)[1]
      I <- matrix(rep(1:q, each = q), byrow = TRUE, ncol = q)
      J <- matrix(rep(1:q, each = q), ncol = q)
      zeroid <- cbind(I[idz], J[idz])
      out.zero <- glasso::glasso(cy, 0, zero = zeroid, trace = FALSE, maxit=10)
      ell1 <- 1 * (out.zero$wi != 0)
      diag(ell1) <- rep(0, nrow(out.zero$wi))
      det.value = sum(log(svd(out.zero$wi)$d));
      B = sum(diag((tj-1)*cy%*%out.zero$wi)) - tj*det.value + sum(ell1)*(log(n))/2
      out.zero$BIC <- B
    } else {
      ell1 <- 1 * (out$wi != 0)
      diag(ell1) <- rep(0, nrow(out$wi))
      dwi <- sum(log(svd(out$wi)$d));
      out$BIC <- sum(diag((tj - 1) * cy %*% out$wi)) - tj * dwi + sum(ell1) * (log(n))/2
      out.zero <- out
    }

    out.zero$debiasTF= debiasTF;
    out.zero$rho <- rho
    out.zero$wiold <- out$wi
    options(warn = 0)
    return(out.zero)
  }


  ### Speed up bootstrapping by using the same rho inside bootstrap sample

  boot.Cribben2 <- function(zz, tcp, p, rhos, n.cl, bound, nboot) {
    T <- dim(zz)[1]
    q <- dim(zz)[2]
    B <- nboot

    # building a stationary bootstrap sample
    if (missing(p)) {
      p <- .2
    }

    # similar but for Cribben et al

    #library(doParallel)
    cl <- parallel::makeCluster(n.cl)
    doParallel::registerDoParallel(cl)
    Bstat <- foreach::foreach(b = 1:B, .combine = c, .export = c("glasso.bic2"), .packages = c("glasso")) %dopar% {
      # building a stationary bootstrap sample
      size <- 0
      z_b <- NULL
      while (size < T) {
        b <- 1 + stats::rgeom(1, p)
        i <- sample(1:T, 1)
        if (i + b - 1 <= T) {
          z_b <- rbind(z_b, zz[i:(i + b - 1), ])
        } else {
          z_b <- rbind(z_b, zz[i:T, ], zz[(1:min(i - 1, i + b - T)), ])
        }
        size <- size + b
      }
      z_b <- z_b[(1:T), ]

      out <- glasso.bic2(z_b, khat = c(1, T), lambda = rhos[1], bound = bound)
      out1 <- glasso.bic2(z_b, khat = c(1, tcp), lambda = rhos[2], bound = bound)
      out2 <- glasso.bic2(z_b, khat = c(tcp + 1, T), lambda = rhos[3], bound = bound)

      dif <- out$BIC - out1$BIC - out2$BIC
    }
    parallel::stopCluster(cl)
    return(Bstat)
  }

  #################################
  ### Main function starts here
  # Combine sample
  #################################
  y <- rbind(subY1, subY2)
  n <- dim(y)[1]
  khat <- nrow(subY1)

  fit0 <- glasso.bic2(y, khat = c(1, n), lambda = lambda, bound = bound)
  fit1 <- glasso.bic2(y, khat = c(1, khat), lambda = lambda, bound = bound)
  fit2 <- glasso.bic2(y, khat = c(khat + 1, n), lambda = lambda, bound = bound)
  diff <- fit0$BIC - (fit1$BIC + fit2$BIC)
  rhos <- c(fit0$rho, fit1$rho, fit2$rho)
  bsample <- boot.Cribben2(y, tcp = khat, p, rhos = rhos, n.cl = n.cl, bound = bound, nboot = nboot)
  pval <- sum(bsample >= diff) / length(bsample)

  out <- list()
  out$bsample = bsample;
  out$fit0 <- fit0
  out$fit1 <- fit1
  out$fit2 <- fit2
  out$diff <- diff
  out$pval <- pval
  out$rho <- rhos

  return(out)
}
