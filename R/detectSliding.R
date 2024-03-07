#' @name detectSliding
#' @title Change point detection using PCA and sliding method
#' @importFrom foreach %dopar%
#' @param Y data: Y = length*dim
#' @param wd window size for sliding averages
#' @param Del Delta away from the boundary restriction
#' @param L the number of factors
#' @param q methods in calculating long-run variance of the test statistic. Default is "fixed" = length^(1/3) or "andrews" implements data adaptive selection, or user specify the length
#' @param alpha significance level of the test
#' @param nboot the number of bootstrap sample for p-value. Default is 199.
#' @param n.cl number of cores in parallel computing. The default is (machine cores - 1)
#' @param bsize block size for the Block Wild Bootstrapping. Default is log(length),  "sqrt" uses sqrt(length), "adaptive" determines block size using data dependent selection of Andrews
#' @param bootTF determine whether the threshold is calculated from bootstrap or asymptotic
#' @param scaleTF scale the variance into 1
#' @param diagTF include diagonal term of covariance matrix or not
#' @param plotTF Draw plot to see test statistic and threshold
#' @return \strong{sW} The test statistic
#' @return \strong{L} The number of factors used in the procedure
#' @return \strong{q} The estimated vectorized autocovariance on each regime.
#' @return \strong{crit} The critical value to identify change point
#' @return \strong{bsize} The block size of the bootstrap
#' @return \strong{diagTF} If TRUE, the diagonal entry of covariance matrix is used in detecting connectivity changes.
#' @return \strong{bootTF} If TRUE, bootstrap is used to find critical value
#' @return \strong{scaleTF} If TRUE, the multivariate signal is studentized to have zero mean and unit variance.
#' @examples \donttest{out4 = detectSliding(changesim, wd=40, L=2, n.cl=1)}
#' @export

detectSliding <- function(Y, wd = 40, L, Del, q = "fixed", alpha = .05, nboot = 199, n.cl, bsize = "log", bootTF = TRUE, scaleTF = TRUE, diagTF = TRUE, plotTF = TRUE) {
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


  Change.sliding <- function(cL, wd = 30, crit, q, n.cl, plotTF = TRUE) {
    # data cL is dim*length
    length <- ncol(cL)

    sliding.stat <- function(cL, Del, q) {
      length <- n <- dim(cL)[2]
      if (missing(Del)) {
        Del <- floor(length * .05)
      }
      cL <- scale(t(cL), center = TRUE, scale = FALSE) # centering
      cL <- t(cL)
      csum.cL <- apply(cL, 1, cumsum) # apply gives length*dim as output
      gmean <- rowMeans(cL)
      m <- floor(length / 2)
      S1 <- LRV.bartlett(cL[, (1:m)], q)
      S2 <- LRV.bartlett(cL[, -(1:m)], q)
      SS <- (S1$sn2 + S2$sn2) / 2
      invSS <- solve(SS)
      setid <- seq(from = Del + 1, to = n - Del, by = 1)
      W <- NULL

      for (j in setid) {
        v1 <- (csum.cL[j, ] - j * gmean) * sqrt(length) / (sqrt(j * (length - j)))
        W <- c(W, t(v1) %*% invSS %*% v1)
      }
      return(max(W))
    }

    time.seq <- seq(from = wd, to = length - wd)

    if (n.cl > 1) {
      cl <- parallel::makeCluster(n.cl)
      doParallel::registerDoParallel(cl)
      result <- foreach::foreach(k = 1:length(time.seq), .combine = c, .export = c("LRV.bartlett", "blocklength.andrew")) %dopar% {
        i <- time.seq[k]
        id1 <- seq(from = i - wd + 1, to = i + wd, by = 1)
        subY1 <- cL[, id1]
        out <- sliding.stat(subY1, q = q)
      }
      sW <- result
      parallel::stopCluster(cl)
    } else {
      sW <- NULL
      for (i in time.seq) {
        id1 <- seq(from = i - wd + 1, to = i + wd, by = 1)
        subY1 <- cL[, id1]
        out <- sliding.stat(subY1, q = q)
        sW <- c(sW, out)
      }
    }

    out <- list()
    out$time.seq <- time.seq
    out$sW <- sW

    if (!is.null(crit)) {
      khat.W <- locatecp(sW, crit = crit, Del = wd)$khat
      out$crit <- crit
      out$br <- khat.W

      if (plotTF) {
        graphics::par(mfrow = c(1, 1))
        graphics::par(mar = c(2, 2, 2, 1))
        graphics::par(cex = .7)
        graphics::plot(time.seq, sW, type = "l", ylab = "Test statistic", xlab = "")
        graphics::title("PCA sliding")
        graphics::abline(v = khat.W, col = "blue")
        graphics::abline(h = crit, col = "red")
      }
    }

    return(out)
  }


  locatecp <- function(ts1, crit, Del = Del) {
    khat <- which.max(ts1)
    tstat <- ts1[khat]
    khat <- khat + Del - 1
    ts1 <- c(rep(0, Del - 1), ts1, rep(0, Del))
    n <- length(ts1)
    br <- c(0, n)
    result <- list(tstathist = tstat, Brhist = khat)

    # Step1 : 1 break
    if (tstat > crit) {
      br <- c(0, khat, n)
      st <- sum(abs(tstat) > crit)

      # Step1 : 1 break
      spindex <- c(1, 1)
      while (st > 0) {
        lbr <- length(br)
        Ts <- NULL
        sst <- Br <- NULL
        brindex <- seq(from = 1, to = lbr - 1, by = 1)
        for (j in brindex) {
          if (spindex[j]) {
            id <- seq(from = br[j] + 1, to = br[j + 1], by = 1)
            dat <- ts1[id]
            if (length(dat) > 2 * Del + 1) { # set the minimum segment length as 2*Del+1
              nn <- length(dat)
              subdat <- dat[(Del + 1):(nn - Del)]
              khat2 <- Del + which.max(subdat)
              tstat2 <- max(subdat)
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
    ## Small refinement disregarding peaked breaks.
    if (length(br) > 2) {
      fbr <- 0
      for (i in 2:(length(br) - 1)) {
        id <- c(br[i] - 2, br[i] - 1, br[i], br[i] + 1, br[i] + 2) ## 5 consecutive to be a change-point
        if (sum(ts1[id] > crit) == 5) {
          fbr <- c(fbr, br[i])
        }
      }
      fbr <- c(fbr, n)
    } else {
      fbr <- br
    }

    result$khat <- fbr
    return(result)
  }

  ## For critical value calculation for sliding Q

  critQ.sliding <- function(del, D, n, r) {
    if (missing(r)) {
      r <- 2000
    }
    if (missing(n)) {
      n <- 5000
    }

    dist <- rep(0, r)
    eps <- floor((del / 2) * n)
    x <- floor(del * n)
    #  index = seq(from=eps*n, to = (1-eps)*n, by=1);
    index2 <- seq(from = 1, to = n - x - 1, by = 1)

    pp <- seq(from = 1, to = x - 1, by = 1) / x
    wt <- pp * (1 - pp)
    dd <- numeric(r)
    BB <- matrix(0, D, x - 1)

    #library(doParallel)
    n.cl <- parallel::detectCores(logical = TRUE) - 1
    cl <- parallel::makeCluster(n.cl)
    doParallel::registerDoParallel(cl)
    dist <- foreach::foreach(k = 1:r, .combine = c) %dopar% {
      e <- matrix(stats::rnorm(n * D), ncol = n) ## Generate D*n matrix
      cc <- NULL
      for (j in index2) {
        id <- seq(from = j, to = j + x - 1, by = 1)
        B <- apply(e[, id], 1, cumsum)
        for (i in 1:D) {
          tmp <- B[-x, i] - pp * B[x, i]
          BB[i, ] <- tmp / sqrt(x * wt)
        }
        b <- colSums(BB * BB)
        cc <- c(cc, max(b))
      }
      max(cc) ## maximum over index2
    }
    parallel::stopCluster(cl)
    return(dist)
  }

  ########################################################
  # Block bootstrap critical value calculation
  #########################################################
  crit.Sliding <- function(cL, wd = 40, bsize = "log", alpha = .05, q, n.cl, nboot) {
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
      bsize <- bsize ## User provide
    }

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
    Bstat <- foreach::foreach(b = 1:30, .combine = c, .export = c("Change.sliding", "blocklength.andrew", "LRV.bartlett")) %dopar% {
      resample <- BWB.cL(cL, bsize = bsize)
      out <- Change.sliding(resample, wd = wd, crit = NULL, q = q, n.cl = 1)
      c(stats::quantile(out$sW, 1 - alpha))
    }
    parallel::stopCluster(cl)
    return(c(mean(Bstat), bsize))
  }


  ##########################################################################
  ### Main function starts here #########################################
  nn <- 2 * wd
  if (q == "fixed") {
    q <- floor(nn^(1 / 3))
  } else if (q == "andrew") {
    q <- apply(cL, 1, blocklength.andrew)
    q <- floor(stats::median(q))
    q <- min(q, floor(nn / 3))
  } else {
    q <- q
  } # User definded

  if (bootTF) {
    bb <- crit.Sliding(cL, wd = wd, q = q, bsize = bsize, alpha = alpha, n.cl = n.cl, nboot = nboot)
    crit <- bb[1]
  } else {
    crit <- stats::quantile(critQ.sliding(nrow(cL), n = 1000, r = 500), 1 - alpha)
  }

  out <- Change.sliding(cL, wd = wd, crit = crit, q = q, plotTF = plotTF, n.cl = n.cl)
  out$L <- L
  out$q <- q
  out$crit <- crit
  out$bsize <- ifelse(bootTF, bb[2], NA)
  out$diagTF <- diagTF
  out$bootTF <- bootTF
  out$scaleTF <- scaleTF

  return(out)
}
