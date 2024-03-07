#' @name detectMaxChange
#' @title Change point detection using max-type statistic as in Jeong et. al (2016)
#' @param Y Input data matrix
#' @param m window sizes
#' @param margin margin
#' @param thre.localfdr threshold for local fdr
#' @param design.mat design matrix for analyzing task data
#' @param plotTF Draw plot to see test statistic and threshold
#' @param n.cl number of clusters for parallel computing
#' @return \strong{CLX} Test statistic corresponding to window size arranged in column
#' @return \strong{CLXLocalFDR} The Local FDR calculated for each time point
#' @return \strong{br} The final estimated break points
#' @examples \donttest{out2= detectMaxChange(changesim, m=c(30, 35, 40, 45, 50), n.cl=1)}
#' @export

detectMaxChange <- function(Y, m=c(30, 40, 50), margin = 30, thre.localfdr = 0.2, design.mat = NULL, plotTF = TRUE, n.cl) {

  ################## Internal Functions ###############################

  sp.mix.logconc <- function(z, doplot = FALSE, thre.localfdr,
                             tol = 1.0e-6, max.iter = 100) {
    #library(LogConcDEAD)

    ## Initial step
    q0 <- stats::quantile(z, probs = .75)
    p.0 <- mean(z <= q0)
    mu.0 <- mean(z[z <= q0])
    sig.0 <- stats::sd(z[z <= q0])

    mu.1 <- mean(z[z > q0])
    sig.1 <- stats::sd(z[z > q0])
    f1.tilde <- stats::dnorm(z, mean = mu.1, sd = sig.1)

    term1 <- p.0 * stats::dnorm(z, mu.0, sig.0)
    term2 <- term1 + (1 - p.0) * f1.tilde
    gam <- term1 / term2

    ## EM-step
    diff <- 100
    k <- 0
    while ((k < 3) | ((k < max.iter) & (diff > tol))) {
      k <- k + 1
      ## E-step
      term1 <- p.0 * stats::dnorm(z, mu.0, sig.0)
      term2 <- term1 + (1 - p.0) * f1.tilde
      new.gam <- term1 / term2

      ## M-step
      w.gam <- new.gam / sum(new.gam, na.rm = T)
      new.mu.0 <- sum(w.gam * z, na.rm = T)
      new.sig.0 <- sqrt(sum(w.gam * (z - new.mu.0)^2, na.rm = T))
      new.p.0 <- mean(new.gam, na.rm = T)

      which.z <- new.gam <= .95
      weight <- (1 - new.gam[which.z]) / sum(1 - new.gam[which.z])
      f1.tilde[which.z] <- exp(LogConcDEAD::mlelcd(z[which.z], w = weight)$logMLE)
      f1.tilde[!which.z] <- 0
      which.gam <- (new.gam <= 0.75) * (new.gam >= 0.01)
      diff <- max(abs(gam - new.gam)[which.gam])

      p.0 <- new.p.0
      mu.0 <- new.mu.0
      sig.0 <- new.sig.0
      gam <- new.gam
    }

    y0 <- p.0 * stats::dnorm(z, mu.0, sig.0)
    y1 <- (1 - p.0) * f1.tilde
    y <- y0 + y1
    Efdr <- sum(y0 * y1 / y) / sum(y1)

    which.z <- gam <= thre.localfdr
    thre <- min(z[which.z])

    res <- list(
      p.0 = p.0, mu.0 = mu.0, sig.0 = sig.0,
      localfdr = gam, thre = thre, Efdr = Efdr, iter = k
    )

    return(res)
  }

  CLX.test <- function(X, Y)
  # Reference: Cai, Liu and Xia (2013), JASA 108, 265-277.
  {
    X <- despike(X)
    Y <- despike(Y)
    n1 <- dim(X)[1]
    n2 <- dim(Y)[1]
    p <- dim(X)[2]
    W1 <- X - matrix(1 / n1, n1, n1) %*% X
    W2 <- Y - matrix(1 / n2, n2, n2) %*% Y
    S1 <- t(W1) %*% W1 / n1
    S2 <- t(W2) %*% W2 / n2
    Theta1 <- Theta2 <- matrix(0, p, p)
    for (i in 1:n1) {
      Theta1 <- Theta1 + (1 / n1) * (W1[i, ] %*% t(W1[i, ]) - S1)^2
    }
    for (i in 1:n2) {
      Theta2 <- Theta2 + (1 / n2) * (W2[i, ] %*% t(W2[i, ]) - S2)^2
    }
    W <- (S1 - S2) / sqrt(Theta1 / n1 + Theta2 / n2)
    M <- W^2
    M.n <- max(M)
    Test.stat <- M.n - 4 * log(p) + log(log(p))
    pvalue <- 1 - exp(-1 / sqrt(8 * pi) * exp(-Test.stat / 2))

    Support <- matrix(0, p, p)
    diag(Support) <- diag(M) >= 2 * log(p)
    off.diag <- as.logical(lower.tri(Support) + upper.tri(Support))
    q.alpha <- -log(8 * pi) - 2 * log(-log(1 - 0.2))
    Support[off.diag] <- (M[off.diag] >= (4 * log(p) - log(log(p)) + q.alpha))

    return(list(Test.stat = Test.stat, pvalue = pvalue, Support = Support))
  }

  Compute.pval <- function(X, m, margin, n.cl = 4) {
    #library(doParallel)

    cl <- parallel::makeCluster(n.cl)
    doParallel::registerDoParallel(cl)

    n <- nrow(X)
    Test.stat <- p.value <- n.edge <- z <- rep(NA, n)
    nonmargin <- (margin + 1):(n - margin)

    PAR <- foreach::foreach(tau = nonmargin, .export = c("CLX.test", "despike"), .combine = rbind) %dopar% {
      x1 <- X[max(1, tau - m):(tau - 1), ]
      x2 <- X[(tau + 1):min(tau + m, n), ]
      res <- CLX.test(x1, x2)
      data.frame(
        Test.stat = res$Test.stat,
        p.value <- res$pvalue,
        n.edge <- mean(res$Support),
        z <- min(max(-10, stats::qnorm(1 - res$pvalue)), 10)
      )
    }

    Test.stat[nonmargin] <- PAR$Test.stat
    p.value[nonmargin] <- PAR$p.value
    n.edge[nonmargin] <- PAR$n.edge
    z[nonmargin] <- PAR$z

    parallel::stopCluster(cl)

    return(data.frame(
      tau = 1:nrow(X),
      Test.stat = Test.stat,
      p.value = p.value,
      z = z,
      n.edge = n.edge
    ))
  }

  determineCHGPT <- function(localfdr, thre, ChunkSize = 10) {
    ind <- which(localfdr <= thre)
    if (is.na(ind[1])) {
      message("No change point has been detected.")
      return(0)
    } else {
      begining <- ind[c(1, which(diff(ind) > ChunkSize) + 1)]
      ending <- ind[c(which(diff(ind) > ChunkSize), length(ind))]
      ccc <- cbind(begining, ending)
      minat <- rep(0, nrow(ccc))
      for (j in 1:nrow(ccc)) {
        minat[j] <- (ccc[j, 1]:ccc[j, 2])[which.min(localfdr[ccc[j, 1]:ccc[j, 2]])]
      }

      return(minat)
    }
  }

  find.taskchange <- function(mat) {
    ans <- NULL
    for (i in 1:ncol(mat)) {
      ans <- c(ans, which(diff(mat[, i]) != 0))
    }
    ans <- sort(unique(ans))

    return(ans)
  }

  minmax.norm <- function(x) {
    M <- max(x, na.rm = TRUE)
    m <- min(x, na.rm = TRUE)

    return((x - m) / max(M - m, 1.0e-12))
  }

  root.sumsq <- function(x) sqrt(sum(x^2))

  L1.norm <- function(x) max(abs(x), na.rm = TRUE)


  plot.vec.ts <- function(X, title, xlab = "Time", ylab = "Node") {
    graphics::plot(1:nrow(X), minmax.norm(X[, 1]) - 1 / 2,
      type = "l",
      ylim = c(0, ncol(X) + 1),
      xlab = "Time", ylab = "Node",
      main = title, axes = FALSE
    )
    graphics::axis(1)
    graphics::axis(2)
    for (j in 2:ncol(X)) {
      graphics::lines(1:nrow(X), minmax.norm(X[, j]) + j - 1 / 2, col = j)
    }
  }

  despike <- function(X) {
    m <- apply(X, 2, stats::median, na.rm = TRUE)
    MAD <- apply(X, 2, mad, na.rm = TRUE)
    LB <- m - 2.5 * MAD
    UB <- m + 2.5 * MAD
    for (i in 1:ncol(X)) {
      X[, i] <- pmin(pmax(LB[i], X[, i]), UB[i])
      # X[(LB[i] > X[,i] | X[,i] > UB[i]),i] <- NA
    }
    # X <- na.omit(X)

    return(X)
  }

  Resid.matrix <- function(X) diag(nrow(X)) - X %*% solve(t(X) %*% X) %*% t(X)

  support.recovery <- function(X, chgpt) {
    n <- length(chgpt)
    X1 <- X[1:(chgpt[1] - 1), ]
    Support <- NULL
    S <- tmp1 <- stats::cov(X1)
    for (i in 1:(n - 1)) {
      X2 <- X[chgpt[i]:(chgpt[i + 1] - 1), ]
      Support <- rbind(Support, CLX.test(X1, X2)$Support)
      S <- rbind(S, stats::cov(X2))
      X1 <- X2
    }
    X2 <- X[chgpt[n]:nrow(X), ]
    Support <- rbind(Support, CLX.test(X1, X2)$Support)
    S <- rbind(S, stats::cov(X2))

    return(list(chgpt = chgpt, S = S, Support = Support))
  }


  #######################################################################
  ## Start of main function #############################################

  X <- Y
  n <- nrow(X)
  nonmargin <- (margin + 1):(n - margin)

  if (missing(n.cl)) {
    n.cl <- parallel::detectCores(logical = TRUE) - 1
  }

  if( length(m) < 2){ stop('More than two window sizes for m!')}

  # To compute CLX test statistics
  z <- test.stat <- N.edge <- matrix(NA, nrow(X), length(m))
  for (i in 1:length(m)) {
    res <- Compute.pval(X = X, m = m[i], margin = margin, n.cl = n.cl)
    test.stat[, i] <- res$Test.stat
    z[, i] <- res$z
    N.edge[, i] <- res$n.edge
  }

  # To compute the PC scores
  pc1 <- rep(NA, nrow(z))
  zs <- z[nonmargin, ]
  if (!is.null(design.mat)) {
    pc1[nonmargin] <- stats::princomp(zs, scores = TRUE)$scores[, 1]
    if (stats::cor(apply(zs, 1, mean), pc1[nonmargin]) < 0) {
      pc1[nonmargin] <- -pc1[nonmargin]
    }
  } else {
    pc1[nonmargin] <- apply(zs, 1, mean)
  }

  # Local FDR based on the semiparametric mixture model
  res <- sp.mix.logconc(pc1[nonmargin], thre.localfdr = thre.localfdr)
  localfdr <- rep(NA, nrow(z))
  localfdr[nonmargin] <- res$localfdr
  BF <- res$p.0 / (1 - res$p.0) * (1 / thre.localfdr - 1)

  # To determine the change points
  CHGPT <- determineCHGPT(localfdr = localfdr, thre = thre.localfdr)
  CHGPT <- c(0, CHGPT, n)
  if (!is.null(design.mat)) {
    task.change <- find.taskchange(design.mat)
    edge <- rep(NA, nrow(z))
    edge[nonmargin] <- apply(N.edge[nonmargin, ], 1, max)
  } else {
    task.change <- NULL
  }

  if (plotTF) {
  stats::plot.ts(localfdr, main = "Max", xlab = "Time", ylab = "Local FDR")
  graphics::abline(h = thre.localfdr, col = 2, lty = 2)
  graphics::abline(v = CHGPT, col = "blue")
}

out <- list()
  out$CLX = test.stat
  out$CLXLocalFDR = localfdr
  out$br = CHGPT
  return(out)
}
