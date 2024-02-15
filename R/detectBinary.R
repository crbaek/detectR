#' @name detectBinary
#' @title Change point detection using PCA and binary segmentation
#' @description This function uses PCA-based method to find breaks. Simultaneous breaks are found from binary segmentation.
#' @param Y data: Y = length*dim
#' @param Del Delta away from the boundary restriction
#' @param L the number of factors
#' @param q methods in calculating long-run variance of the test statistic. Defaul is "fixed" = length^(1/3). Adaptive selection method is also available via "andrews", or user specify the length
#' @param alpha significance level of the test
#' @param nboot the number of bootstrap sample for pvalue. Defauls is 199.
#' @param n.cl number of cores in parallel computing. The default is (machine cores - 1)
#' @param bsize block size for the Block Wild Boostrapping. Default is log(length),  "sqrt" uses sqrt(length), "adaptive" deterines block size usign data dependent selection of Andrews
#' @param bootTF determine whether the threshold is calculated from bootstrap or asymptotic
#' @param scaleTF scale the variance into 1
#' @param diagTF include diagonal term of covariance matrix or not
#' @param plotTF Draw plot to see test statistic and threshold
#' @return \strong{tstathist} The complete history of test tsatistic
#' @return \strong{Brhist} The sequence of breakspoints found from binay splitting
#' @return \strong{L} The number of factors used in the procedure
#' @return \strong{q} The estimated vecorized autocovariance on each regime.
#' @return \strong{crit} The critical vlaue to identify change point
#' @return \strong{bsize} The block size of the bootstrap
#' @return \strong{diagTF} If TRUE, the diagonal entry of covariance matrix is used in detecting connectivity changes.
#' @return \strong{bootTF} If TRUE, boostrap is used to find critical value
#' @return \strong{scaleTF} If TRUE, the multivariate signal is studentized to have zero mean and unit variance.
#' @examples \donttest{out3= detectBinary(changesim, L=2, n.cl=1)}
#' @export


detectBinary <- function(Y, Del, L, q = "fixed", alpha = .05, nboot = 199, n.cl, bsize = "log", bootTF = TRUE, scaleTF = TRUE, diagTF = TRUE, plotTF = TRUE) {
  regpar <- par(no.readonly = TRUE)
  on.exit(par(regpar))
  if (scaleTF) {
    Y <- scale(Y, center = TRUE, scale = TRUE)
  }
  length <- n <- nrow(Y)
  if (missing(L)) {
    L <- 3
  }
  if (missing(n.cl)) {
    n.cl <- parallel::detectCores(logical = TRUE) - 1
  }

  # Extract factors
  vvz <- svd(stats::cov(Y))
  Lamz <- vvz$u[, 1:L]
  zL <- t(Lamz) %*% t(Y)
  zL <- scale(t(zL), center = TRUE, scale = FALSE)
  zL <- t(zL)
  ell <- nrow(zL)

  ## Now transform zL into covariance form; vech(ftf_t')
  if (!diagTF) {
    cL <- matrix(0, ell * (ell - 1) / 2, length)
    for (i in 1:length) {
      tp <- zL[, i] %*% t(zL[, i])
      cL[, i] <- tp[lower.tri(tp, diag = FALSE)]
    }
  } else {
    cL <- matrix(0, ell * (ell + 1) / 2, length)
    for (i in 1:length) {
      tp <- zL[, i] %*% t(zL[, i])
      cL[, i] <- tp[lower.tri(tp, diag = TRUE)]
    }
  }
  if (missing(Del)) {
    Del <- 40
  }

  Find.change <- function(cL, Del, q, n.cl) {
    # Find the change-point using sst, see Han and Inoue
    # Data is dim*length
    if (missing(Del)) {
      Del <- 40
    }
    dim <- nrow(cL)
    length <- n <- ncol(cL)
    csum.cL <- apply(cL, 1, cumsum) # it is length*dim
    gmean <- rowMeans(cL)
    ## Wald test; takes quite bit of time
    qq <- max(q, Del)
    inpi <- seq(from = qq + 1, to = length - qq, by = 1)
    W <- NULL

    if (n.cl > 1) {
      #require(doParallel)
      cl <- parallel::makeCluster(n.cl)
      doParallel::registerDoParallel(cl)
      W <- foreach::foreach(k = 1:length(inpi), .combine = c, .export = c("LRV.bartlett", "blocklength.andrew")) %dopar% {
        j <- inpi[k]
        S1 <- LRV.bartlett(cL[, (1:j)], q = q)
        S2 <- LRV.bartlett(cL[, -(1:j)], q = q)
        SS <- (j * S1$sn2 + S2$sn2 * (length - j)) / length
        invSS <- solve(diag(diag(SS)))
        v1 <- (csum.cL[j, ] - j * gmean) * sqrt(length) / (sqrt(j * (length - j)))
        t(v1) %*% invSS %*% v1
      }
      parallel::stopCluster(cl)
    } else {
      for (j in inpi) {
        S1 <- LRV.bartlett(cL[, (1:j)], q = q)
        S2 <- LRV.bartlett(cL[, -(1:j)], q = q)
        SS <- (j * S1$sn2 + S2$sn2 * (length - j)) / length
        invSS <- solve(diag(diag(SS)))
        v1 <- (csum.cL[j, ] - j * gmean) * sqrt(length) / (sqrt(j * (length - j)))
        W <- c(W, t(v1) %*% invSS %*% v1)
      }
    }

    k1 <- which.max(W)
    Wstat <- W[k1]
    khat1 <- k1 + qq
    out <- list()
    out$W <- W
    out$khat.W <- khat1
    out$xW <- inpi
    out$tstat.W <- Wstat
    return(out)
  }

  ### Binary segmentation
  Change.binary <- function(cL, Del, crit, n.cl, q, plotTF = TRUE) {
    if (missing(Del)) {
      Del <- 40
    }
    maxK <- 3
    # Data is dim*length
    n <- length <- ncol(cL)
    dim <- nrow(cL)
    out1 <- Find.change(t(scale(t(cL))), Del = Del, q = q, n.cl = n.cl)
    tstat <- out1$tstat.W
    khat <- out1$khat.W
    result <- list(tstathist = tstat, Brhist = khat)
    br <- c(0, n)
    nmax <- 1

    # Step1 : 1 break
    if (tstat > crit) {
      br <- c(0, khat, n)
      st <- sum(abs(tstat) > crit)

      # Step1 : 1 break
      spindex <- c(1, 1)

      while (st > 0 && nmax < maxK) {
        nmax <- nmax + 1
        lbr <- length(br)
        Ts <- NULL
        sst <- Br <- NULL
        brindex <- seq(from = 1, to = lbr - 1, by = 1)
        for (j in brindex) {
          if (spindex[j]) {
            id <- seq(from = br[j] + 1, to = br[j + 1], by = 1)
            dat <- cL[, id]
            if (ncol(dat) > 2 * Del + 1) { # set the minimum segment length as 2*Del
              out2 <- Find.change(t(scale(t(dat))), Del = Del, q = q, n.cl = n.cl)
              tstat2 <- out2$tstat.W
              khat2 <- out2$khat.W
              Br <- c(Br, br[j] + khat2)
              Ts <- c(Ts, tstat2)
              idf <- 1 * (tstat2 > crit)
              sst <- c(sst, rep(idf, idf + 1))
            } else {
              sst <- c(sst, 0)
            }
          } else {
            sst <- c(sst, 0)
          }
        }

        st <- sum(sst)
        if (!is.null(Ts)) {
          newbr <- (abs(Ts) > crit) * Br
          br <- unique(c(br, newbr))
        }
        br <- sort(br)
        spindex <- sst

        result$tstathist <- c(result$tstathist, Ts)
        result$Brhist <- c(result$Brhist, Br)
      }
    }

    if (plotTF) {
      graphics::par(mfrow = c(1, 1))
      graphics::par(mar = c(2, 2, 2, 1))
      graphics::par(cex = .7)
      graphics::plot(out1$xW, out1$W, type = "l", xlab = "", ylab = "Test statistic", main = "PCA binary")
      graphics::abline(v = br[br != n], col = "blue")
      graphics::abline(h = crit, col = "red")
    }

    result$br <- br
    result$iter <- length(br) - 1
    return(result)
  }

  ## For critical value calculation
  ## This is a critical value of Q
  critQ <- function(D, n, r) {
    if (missing(r)) {
      r <- 2000
    }
    if (missing(n)) {
      n <- 5000
    }

    dist <- rep(0, r)
    eps <- .18
    index <- seq(from = eps * n, to = (1 - eps) * n, by = 1)
    pp <- index / n
    wt <- pp * (1 - pp)
    for (k in 1:r) {
      BB <- matrix(0, D, length(index))
      for (i in 1:D) {
        e <- stats::rnorm(n)
        B <- cumsum(e)
        tmp <- B[index] - pp * B[n]
        BB[i, ] <- tmp / sqrt(n * wt)
      }
      b <- colSums(BB * BB)
      dist[k] <- max(b)
    }
    return(dist)
  }

  ##############################################

  ########################################################
  # Block bootstrap critical value calculation
  #########################################################
  crit.BWB <- function(cL, bsize = "log", Del = 40, alpha = .05, q, nboot, n.cl) {
    nB <- nboot
    dim <- p <- nrow(cL)
    length <- n <- ncol(cL)

    if (bsize == "log") {
      bsize <- log(n)
    } else if (bsize == "sqrt") {
      bsize <- sqrt(n)
    } else if (bsize == "adaptive") {
      K.d <- apply(cL, 1, blocklength.andrew)
      bsize <- min(1.5 * stats::median(K.d), log(n * p))
    } else {
      bsize <- bsize
    } ## User provide

    bsize <- floor(bsize)
    X1 <- scale(t(cL), center = TRUE, scale = FALSE)
    cL <- t(X1)

    BWB.cL <- function(cL, bsize) {
      ## Dimenison of X is dim*length
      dim <- p <- nrow(cL)
      length <- n <- ncol(cL)
      nk <- ifelse(n %% bsize == 0, floor(n / bsize), floor(n / bsize) + 1)
      wt <- rep(stats::rnorm(nk), each = bsize)[1:n]
      cL.WB <- cL * wt
      return(cL.WB)
    }

    #library(doParallel)
    cl <- parallel::makeCluster(n.cl)
    doParallel::registerDoParallel(cl)
    Bstat <- foreach::foreach(b = 1:nB, .combine = c, .export = c("Find.change", "blocklength.andrew", "LRV.bartlett")) %dopar% {
      resample <- BWB.cL(cL, bsize = bsize)
      out <- Find.change(resample, Del = Del, q = q, n.cl = 1)
      c(out$tstat.W)
    }
    parallel::stopCluster(cl)
    return(c(stats::quantile(Bstat, (1 - alpha)), bsize))
  }


  ##########################################################################
  ### Main function starts here #########################################

  if (q == "fixed") {
    q <- floor(n^(1 / 3))
  } else if (q == "andrew") {
    q <- apply(cL, 1, blocklength.andrew)
    q <- floor(stats::median(q))
    q <- min(q, floor(n / 3))
  } else {
    q <- q
  } # User definded

  if (bootTF) {
    bb <- crit.BWB(cL, Del = Del, bsize = bsize, alpha = alpha, q = q, n.cl = n.cl, nboot = nboot)
    crit <- bb[1]
  } else {
    crit <- stats::quantile(critQ(nrow(cL), n = 1000, r = 500), 1 - alpha)
  }

  out <- Change.binary(cL, Del = Del, crit, n.cl = n.cl, q = q, plotTF = plotTF)
  out$L <- L
  out$q <- q
  out$crit <- crit
  out$bsize <- ifelse(bootTF, bb[2], NA)
  out$diagTF <- diagTF
  out$bootTF <- bootTF
  out$scaleTF <- scaleTF

  return(out)
}

