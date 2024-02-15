#' @name testPCA
#' @title PCA-based test for the equality of connectivity
#' @description  This function performs PCA-test for testing the equality of connectivity in two fMRI signals
#' @param subY1 a sample of size length*dim
#' @param subY2 a sample of size length*dim
#' @param L the number of factors
#' @param diagTF include diagonal term of covariance matrix or not
#' @param nlag is the number of ACF lag to be used in the test, default is 2, Default is nlag = floor(N^(1/3))
#' @return \strong{tstat} Test statistic
#' @return \strong{pval} Returns the p-value
#' @return \strong{df} The degree of freedom in PCA-best test
#' @return \strong{L} The number of factors used in the test
#' @return \strong{diagTF} If true, the diagonal entry of covarianc matrix is used in testing
#' @examples test3 = testPCA(testsim$X, testsim$Y, L=2)
#' @export


testPCA <- function(subY1, subY2, L = 2, nlag, diagTF = TRUE) {

  ######### Internal functions ###################
  # .........................................................
  # Lund et al. test
  # Inputs:
  #    y1,y2: two vector series of dimension qxN
  #    nlag: the number of ACF lags to consider in testing
  # .........................................................
  test.lund2 <- function(y1, y2, nlag, diagTF = TRUE) {
    q <- dim(y1)[1]
    N <- min(dim(y1)[2], ncol(y2))
    py1 <- stats::acf(t(y1), lag.max = nlag, type = "covariance", plot = FALSE)$acf
    py2 <- stats::acf(t(y2), lag.max = nlag, type = "covariance", plot = FALSE)$acf
    p1 <- py1 - py2
    py <- (py1 + py2) / 2
    tmp <- t(p1[1, , ])
    idx <- lower.tri(p1[1, , ], diag = TRUE)
    delGam <- tmp[idx]

    for (j in 1:nlag) {
      delGam <- c(delGam, as.vector(t(p1[(j + 1), , ])))
    }

    ## Generate index set for the covariance of delGam
    ## Index of row
    I <- t(matrix(rep(1:q, each = q), ncol = q, byrow = TRUE))[idx]
    I <- c(I, rep(rep(1:q, each = q), nlag))
    ## Index of column
    J <- matrix(rep(1:q, each = q), ncol = q, byrow = TRUE)[idx]
    J <- c(J, rep(rep(1:q, q), nlag))
    K <- I
    eL <- J
    H1 <- H2 <- c(rep(0, q * (q + 1) / 2), rep(1:nlag, each = q^2))

    ## Assuming Gaussian in (15) of paper
    M <- floor(N^(1 / 2))
    lmax <- min(M + nlag, N)

    py1 <- stats::acf(t(y1), lag.max = lmax, type = "covariance", plot = FALSE)$acf
    py2 <- stats::acf(t(y2), lag.max = lmax, type = "covariance", plot = FALSE)$acf
    py <- (py1 + py2) / 2
    D <- length(delGam)
    W1 <- matrix(0, D, D)

    ccfh2 <- function(py1, id1, id2, h) {
      if (h >= 0) {
        cf <- py1[h + 1, id1, id2]
      }
      if (h < 0) {
        cf <- py1[abs(h) + 1, id2, id1]
      }
      return(cf)
    }

    for (i in 1:D) {
      for (j in i:D) {
        S <- 0
        for (r in -M:M) {
          term <- ccfh2(py, I[i], K[j], r) * ccfh2(py, J[i], eL[j], r - H1[i] + H2[j])
          +ccfh2(py, I[i], eL[j], r + H2[j]) * ccfh2(py, J[i], K[j], r - H1[i])
          S <- S + term
        }
        W1[i, j] <- S
      }
    }

    if (diagTF) {
      W1 <- W1 + t(W1) - diag(diag(W1))
      tstat2 <- N / 2 * t(delGam) %*% solve(W1) %*% delGam
      pval2 <- 1 - stats::pchisq(tstat2, D)
    } else {
      ## Only include  off-diagonal term in gamma(0) and the rest of them are the full autocovariances
      Iq <- matrix(1, q, q)
      diag(Iq) <- rep(0, q)
      vC <- Iq[lower.tri(diag(q), diag = TRUE)]

      for (j in 1:nlag) {
        vC <- c(vC, matrix(1, q, q))
      }
      id.one <- which(vC == 1)
      C <- matrix(0, nrow = length(id.one), ncol = length(vC))
      for (i in 1:length(id.one)) {
        C[i, id.one[i]] <- 1
      }
      delGam3 <- C %*% delGam
      W3 <- C %*% W1 %*% t(C)
      tstat2 <- N / 2 * t(delGam3) %*% solve(W3) %*% delGam3
      pval2 <- 1 - stats::pchisq(tstat2, length(delGam3))
      delGam <- delGam3
    }

    out <- list()
    out$tstat <- tstat2
    out$pval <- pval2
    out$df <- length(delGam)
    return(out)
  }

  n <- min(nrow(subY1), nrow(subY2))
  if (missing(nlag)) {
    nlag <- floor(n^(1 / 3))
  }

  subY <- rbind(subY1, subY2) # Pooled sample
  vvz <- svd(stats::cov(subY))
  Lamz <- vvz$u[, 1:L]
  zL <- t(Lamz) %*% t(subY)
  n1 <- nrow(subY1)
  out <- test.lund2(zL[, (1:n1)], zL[, -(1:n1)], nlag = nlag, diagTF = diagTF)
  out$L <- L
  out$diagTF <- diagTF

  return(out)
}
