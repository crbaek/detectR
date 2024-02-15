#' @name testMax
#' @title Max-type test for for the equality of connectivity
#' @description This function produces three test results based on max-type block boostrap (BMB), long-run variance block boostrapping with lagged-window estimator (LVBWR) and sum-type block bootstrap (BSUM). See Baek el al. (2019) for details.
#' @param subY1 a sample of size length*dim
#' @param subY2 a sample of size length*dim
#' @param diagTF include diagonal term of covariance matrix or not
#' @param nboot number of bootstrap sample, default is 2000
#' @param q methods in calculating long-run variance of the test statistic. Default is "fixed" = length^(1/3) or "andrews" implements data adaptive selection, or user specify the length
#' @param n.cl number of cores in parallel computing. The default is (machine cores - 1)
#' @return \strong{tstat} Test statisticlue for testing the equality of connectivity structure
#' @return \strong{pval} The pvalue for testing the equality of connectivity structure
#' @return \strong{q} The tuning parameter used in calulating long-run variance
#' @examples \donttest{test2 = testMax(testsim$X, testsim$Y, n.cl=1)}
#' @export

testMax <- function(subY1, subY2, diagTF = TRUE, nboot, q = "fixed", n.cl) {
  if (missing(n.cl)) {
    n.cl <- parallel::detectCores(logical = TRUE) - 1
  }
  if (missing(nboot)) {
    nboot <- 199
  }


  vecDF <- function(y1, y2, K = 0, diagTF, centerTF = TRUE) {
    if (missing(diagTF)) {
      diagTF <- FALSE
    } ## Include diagonal for lag 0? default is false
    # Data dimenskon of y1 & y2 is p*n

    ## centering the sequence
    if (centerTF) {
      y1 <- y1 - rowMeans(y1)
      y2 <- y2 - rowMeans(y2)
    }
    Tt <- ncol(y1)
    d <- nrow(y1)
    #  d1 = d+1*diagTF;  Zt = matrix(0, d1*(d1-1)/2, Tt-K)
    Zt <- NULL
    for (i in 1:(Tt - K)) {
      zz <- y1[, i] %*% t(y1[, i]) - y2[, i] %*% t(y2[, i])
      Z <- zz[upper.tri(zz, diag = diagTF)]

      if (K >= 1) {
        for (j in 1:K) {
          z2 <- y1[, (i + j)] %*% t(y1[, i]) - y2[, (i + j)] %*% t(y2[, i])
          Z <- c(Z, as.vector(z2))
        }
      }

      Zt <- cbind(Zt, Z)
    }

    return(Zt) # Output of Zt
  }


  #####################################
  # Block size selection with Andrews
  #####################################
  bsize.andrew <- function(Z) {
    p <- nrow(Z)
    q <- numeric(p)
    n <- ncol(Z)
    for (i in 1:p) {
      data <- Z[i, ]
      rhohat <- stats::ar.ols(data, FALSE, order = 1)
      rhohat <- rhohat$ar[1]
      q[i] <- 1.1447 * (n * 4 * rhohat^2 / (1 - rhohat^2)^2)^(1 / 3)
    }
    return(ceiling(2 * mean(q)))
  }

  ## Another LRV, using sliding method
  LRV.lag <- function(Z, m) {

    ## Lagged window estimation of long run variance
    var.internal <- function(z, m) {
      gmean <- mean(z)
      ff <- stats::filter(z, rep(1, m) / m, sides = 1, method = "convolution")
      f1 <- ff[-(1:m)]
      return(mean((f1 - gmean)^2))
    }

    Sig <- apply(Z, 1, var.internal, m = m)
    return(m * Sig)
  }

  sumLRV <- function(Z, M) {
    N <- ncol(Z)
    aa <- stats::acf(t(Z), M, plot = FALSE, type = ("covariance"))
    gammahat <- numeric(M + 1)
    for (i in 1:(M + 1)) {
      gammahat[i] <- sum(diag(aa$acf[i, , ]))
    }
    hh <- seq(from = 0, to = M, by = 1)
    wt <- N / (N - hh)
    k <- seq(from = 1, to = M, by = 1)
    bn <- 2 * (1 - (k - 1) / N)
    bn <- c(1, bn)
    sum(bn * wt * gammahat / N)
  }

  ##########################################################################
  ### Main function starts here #########################################

  Z <- vecDF(t(subY1), t(subY2), K = 0, diagTF = diagTF)
  p <- nrow(Z)
  N <- ncol(Z)
  Zbar <- rowMeans(Z)

  if (q == "fixed") {
    q <- floor(N^(1 / 3))
  } else if (q == "andrew") {
    q <- bsize.andrew(Z)
  } else {
    q <- q
  } # User definded
  bsize <- q

  blength <- floor(N / bsize) + 1
  # 1. BMB
  xihat <- sqrt(N) * max(abs(Zbar))

  # 2. LVBW / LVBWR
  Sig3 <- sqrt(LRV.lag(Z, bsize))
  xihat3 <- sqrt(N) * max(abs(Zbar) / Sig3)

  # 4. BSUM
  xihat5 <- sum(Zbar^2)

  if (n.cl > 1) {
    #library(doParallel)
    cl <- parallel::makeCluster(n.cl)
    doParallel::registerDoParallel(cl)
    bsample <- foreach::foreach(b = 1:nboot, .combine = rbind) %dopar% {
      weight <- rep(stats::rnorm(blength, 0, 1), each = bsize)
      weight <- weight[1:N]
      bZ <- 0 * Z
      bZ <- t(Z) * weight
      bZ <- t(bZ)
      bZbar <- rowMeans(bZ)
      bs <- sqrt(N) * max(abs(rowMeans(bZ)))

      ## Rademacher disribution
      wt.R <- rep(sample(c(1, -1), blength, replace = TRUE), each = bsize)
      wt.R <- wt.R[1:N]
      rZ <- 0 * Z
      rZ <- t(Z) * wt.R
      rZ <- t(rZ)
      rZbar <- rowMeans(rZ)

      bSig2 <- sqrt(LRV.lag(rZ, bsize))
      bs2 <- sqrt(N) * max(abs(rZbar) / bSig2)

      bs5 <- sum(bZbar^2)

      c(bs, bs2, bs5)
    }
    parallel::stopCluster(cl)

    pval1 <- sum(bsample[, 1] > abs(xihat)) / nboot # BMB
    pval2 <- sum(bsample[, 2] > abs(xihat3)) / nboot # LVBW
    pval3 <- sum(bsample[, 3] > abs(xihat5)) / nboot # BSUM

    pval <- c(pval1, pval2, pval3)
  } else {
    bsample <- matrix(0, nboot, 3)
    for (r in 1:nboot) {
      weight <- rep(stats::rnorm(blength, 0, 1), each = bsize)
      weight <- weight[1:N]
      bZ <- 0 * Z
      bZ <- t(Z) * weight
      bZ <- t(bZ)
      bZbar <- rowMeans(bZ)
      bsample[r, 1] <- sqrt(N) * max(abs(rowMeans(bZ)))

      wt.R <- rep(sample(c(1, -1), blength, replace = TRUE), each = bsize)
      wt.R <- wt.R[1:N]
      rZ <- 0 * Z
      rZ <- t(Z) * wt.R
      rZ <- t(rZ)
      rZbar <- rowMeans(rZ)

      bSig2 <- sqrt(LRV.lag(rZ, bsize))
      bsample[r, 2] <- sqrt(N) * max(abs(rZbar) / bSig2)
      bsample[r, 3] <- sum(bZbar^2)
    }

    pval1 <- sum(bsample[, 1] > abs(xihat)) / nboot # BMB
    pval2 <- sum(bsample[, 2] > abs(xihat3)) / nboot # LVBW
    pval3 <- sum(bsample[, 3] > abs(xihat5)) / nboot # BSUM
    pval <- c(pval1, pval2, pval3)
  }

  out <- list()
  out$bsample <- bsample
  out$tstat <- c(xihat, xihat3, xihat5)
  out$pval <- pval
  out$q <- bsize
  names(out$pval) <- names(out$tstat) <- c("BMB", "LVBWR", "BSUM")

  return(out)
}
