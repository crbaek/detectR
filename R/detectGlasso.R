#' @name detectGlasso
#' @title Change point detection using Graphical lasso as in Cribben et al. (2012)
#' @description This function implements the Dynamic Connectivity Regression (DCR) algorithm proposed by Cribben el al. (2012) to locate changepoints.
#' @param Y Input data of dimension length*dim (T times d)
#' @param Del Delta away from the boundary restriction
#' @param p  Gep(p) distribution controls the size of stationary bootstrap. The mean block length is 1/p
#' @param lambda two selections possible for optimal parameter of lambda. "bic" finds lambda from bic criteria, or user can directly input the penalty value
#' @param nboot the number of bootstrap sample for pvalue. Default is 100.
#' @param n.cl number of cores in parallel computing. The default is (machine cores - 1)
#' @param bound bound of bic search in "bic" rule. Default is (.001, 1)
#' @param gridTF minimum bic is found by grid search. Default is FALSE
#' @param plotTF Draw plot to see test statistic
#' @return A list with component
#' @return \strong{br} The estimated breakpoints including boundary (0, T)
#' @return \strong{brhist} The sequence of breakspoints found from binay splitting
#' @return \strong{diffhist} The history of BIC reduction on each step
#' @return \strong{W} The estimated vecorized autocovariance on each regime.
#' @return \strong{WI} The estimated vecorized precision matrix on each regime.
#' @return \strong{lambda} The penalty parameter estimated on each regime.
#' @return \strong{pvalhist} The empirical p-values on each binary spltting.
#' @return \strong{fitzero} Detailed output at first stage. Useful in producing plot.
#' @examples \donttest{out1= detectGlasso(changesim, p=.2, n.cl=1)}
#' @export 

detectGlasso <- function(Y, Del, p, lambda = "bic", nboot = 100, n.cl, bound = c(.001, 1), gridTF=FALSE, plotTF=TRUE) {

  # y: Input data of dimension length*dim (T times d)
  # Del:  Delta away from the boundary restriction
  # p : Gep(p) distribution controls the size of stationary bootstrap. The mean block length
  #     is 1/p
  # lambda: two selections possible for optimal parameter of lambda. "bic" finds lambda
  # from bic criteria, or user can directly input the penalty value
  # nboot: the number of bootstrap sample for pvalue. Default is 100.
  # n.cl: number of cores in parallel computing. The default is (machine cores - 1)
  # bound: bound of bic search in "bic" rule. Default is (.001, 1)
  # gridTF: whether rho from BIC is found by grid search
  # debiasTF: wheter debiasing from graphical lasso estimates
  regpar <- par(no.readonly = TRUE)
  on.exit(par(regpar))
  
  if (missing(n.cl)) {
    n.cl <- parallel::detectCores(logical = TRUE) - 1
  }
  if (missing(p)) {
    p <- .2
  }
  if (missing(Del)) {
    Del <- 40
  }

  ## Data dimension is Tt*q
  ### Auxilary functions

  glasso.bic <- function(Y, khat, lambda, bound, debiasTF) {
    st <- khat[1]
    end <- khat[2]
    id <- seq(st, end, by = 1)
    tj <- length(id)
    cy <- stats::cov(Y[id, ])
    n <-  nrow(Y);

    
    if (lambda == "bic") {
      bic.rho <- function(rho, cy, tj) {
        B <- 0
        out <- glasso::glasso(cy, rho, maxit=400)
        ell1 <- 1 * (out$wi != 0)
        diag(ell1) <- rep(0, nrow(out$wi))
        ## tj is replacec by tj-1 by comparing result in the below.
        dwi <- sum(log(svd(out$wi)$d))
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

    if ((length(idz) > 0)*debiasTF) {
      ## Re-estimate precision matrix by imposing zero constraints
      ## It only valid when there is zero estimates in out()
      ## This takes quite a bit of time if dimension is high and close to singluar
      q <- dim(out$wi)[1]
      I <- matrix(rep(1:q, each = q), byrow = TRUE, ncol = q)
      J <- matrix(rep(1:q, each = q), ncol = q)
      zeroid <- cbind(I[idz], J[idz])
      out.zero <- glasso::glasso(cy, 0, zero = zeroid, trace = FALSE, maxit=100)
      ell1 <- 1 * (out.zero$wi != 0)
      diag(ell1) <- rep(0, nrow(out.zero$wi))
            det.value <- sum(log(svd(out.zero$wi)$d));
      out.zero$nlik <- sum(diag((tj-1)*cy%*%out.zero$wi)) - tj*det.value ;
      B <- out.zero$nlik + sum(ell1)*(log(n))/2;
      out.zero$BIC <- B;
      ## log(n) is used for the sum of two BICs
    } else {
      ell1 <- 1 * (out$wi != 0)
      diag(ell1) <- rep(0, nrow(out$wi))
      det.value <- sum(log(svd(out$wi)$d));
      out$nlik <- sum(diag((tj-1)*cy%*%out$wi)) - tj*det.value;
      out$BIC <- out$nlik + sum(ell1)*(log(n))/2
      out.zero <- out;
    }

    out.zero$debiasTF= debiasTF;
    out.zero$rho <- rho
    options(warn = 0)
    return(out.zero)
  }


  boot.Cribben2 <- function(zz, tcp, p, rhos, n.cl, bound, nboot) {
    T <- dim(zz)[1]
    q <- dim(zz)[2]

    # building a stationary bootstrap sample
    if (missing(p)) {
      p <- .2
    }

    # similar but for Cribben et al
    cl <- parallel::makeCluster(n.cl)
    doParallel::registerDoParallel(cl)
    Bstat <- foreach::foreach(b = 1:nboot, .combine = c, .export = c("glasso.bic"), .packages = "glasso") %dopar% {
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

      # Use the same rho to save time
      out <- glasso.bic(z_b, khat = c(1, T), lambda = rhos[1], bound = bound, debiasTF=TRUE)
      out1 <- glasso.bic(z_b, khat = c(1, tcp), lambda = rhos[2], bound = bound, debiasTF=TRUE)
      out2 <- glasso.bic(z_b, khat = c(tcp + 1, T), lambda = rhos[3], bound = bound, debiasTF=TRUE)

      dif <- out$BIC - out1$BIC - out2$BIC
    }
    parallel::stopCluster(cl)

    return(Bstat)
  }

  ########################

  glasso.break <- function(Y, Del, lambda, n.cl, bound) {

    ## Input data is Tt*q
    if (missing(Del)) {
      Del <- 40
    }

    nn <- dim(Y)[1]
    step0 <- glasso.bic(Y, khat = c(1, nn), lambda = lambda, bound = bound, debiasTF=TRUE)
    grid <- seq(Del + 1, nn - Del, by = 1)
    BIC = NULL;
    cl <- parallel::makeCluster(n.cl)
    doParallel::registerDoParallel(cl)

    BIC <- foreach::foreach(b = 1:length(grid), .combine = rbind, .export = c("glasso.bic"), .packages = "glasso") %dopar% {
      out1 <- glasso.bic(Y, khat = c(1, grid[b]), lambda = step0$rho, bound = bound, debiasTF=FALSE)
      out2 <- glasso.bic(Y, khat = c(grid[b] + 1, nn), lambda = step0$rho, bound = bound, debiasTF=FALSE)
      c(out1$nlik + out2$nlik, out1$BIC + out2$BIC)
    }
    parallel::stopCluster(cl)

    id = which.min(BIC[,1]); newkhat <- grid[id];
    step1 <- glasso.bic(Y, khat = c(1, newkhat), lambda = lambda, bound = bound, debiasTF=TRUE)
    step2 <- glasso.bic(Y, khat = c(newkhat+1, nn), lambda = lambda, bound = bound, debiasTF=TRUE)

    bicnew = step1$BIC + step2$BIC;
    return(list(BIC = BIC[,1], khat = newkhat, time.seq = grid, bic0 = step0$BIC, bicnew = bicnew))
  }

  ### Main function starts here

  n <- dim(Y)[1]
  out <- list()

  br <- newbr <- brhist <- c(0, n)
  st <- 1
  TF <- 1
  diffhist <- NULL
  nmax <- 1
  maxK = floor(log2(n/(2*Del)))+2;
  while (st > 0 && nmax < maxK) {
    newTF <- NULL
    for (i in 1:(length(br) - 1)) {
      if (TF[i] > 0 && br[i + 1] - br[i] > 2 * Del) {
        id <- seq(from = br[i] + 1, to = br[i + 1], by = 1)
        zz <- Y[id, ]
        fit <- glasso.break(zz, Del, lambda = lambda, n.cl = n.cl, bound = bound)
        if( nmax == 1 && i == 1){ fitzero = fit; }
        khat <- fit$khat
        bicdiff <- fit$bicnew - fit$bic0
        newk <- khat + br[i]

        if (bicdiff < 0 && khat > Del && (br[i + 1] - newk) > Del) {
          newbr <- c(newbr, newk)
          newTF <- c(newTF, 1, 1)
          brhist <- c(brhist, newk)
          diffhist <- rbind(diffhist, bicdiff)
        } else {
          newTF <- c(newTF, 0)
        }
      } else {
        newTF <- c(newTF, 0)
      }
    }
    br <- unique(newbr)
    br <- sort(br)
    TF <- newTF
    newbr <- br
    st <- sum(TF)
    nmax <- nmax + 1
  }

  ## Final tuning stage (Appendix 1  steps 7-9)
  cdiff <- 1
  pvalhist <- NULL
  while (cdiff > 0 & length(br) > 2) {
    cdiff <- 0
    finalbr <- NULL

    for (j in 1:(length(br) - 2)) {
      fit0 <- glasso.bic(Y, khat = c(br[j] + 1, br[j + 2]), lambda = lambda, bound = bound, debiasTF=TRUE)
      fit1 <- glasso.bic(Y, khat = c(br[j] + 1, br[j + 1]), lambda = lambda, bound = bound, debiasTF=TRUE)
      fit2 <- glasso.bic(Y, khat = c(br[j + 1] + 1, br[j + 2]), lambda = lambda, bound = bound, debiasTF=TRUE)
      diff <- fit0$BIC - (fit1$BIC + fit2$BIC)
      rhos <- c(fit0$rho, fit1$rho, fit2$rho)
      idx <- ((br[j] + 1):br[j + 2])

      zz <- Y[idx, ]
      bsample <- boot.Cribben2(zz, tcp = br[j + 1] - br[j], p, rhos = rhos, n.cl = n.cl, bound = bound, nboot = nboot)
      pval <- sum(bsample >= diff) / length(bsample)
      pvalhist <- c(pvalhist, pval)
      cdiff <- cdiff + 1 * (pval >= .05)
      if (pval <= 0.05) {
        finalbr <- c(finalbr, br[j + 1])
      }
    }

    br <- c(0, finalbr, n)
  }

  ## Get the estimates of final breaks
  W <- WI <- R <- NULL
  for (j in 1:(length(br) - 1)) {
    fit1 <- glasso.bic(Y, khat = c(br[j] + 1, br[j + 1]), lambda = lambda, bound = bound, debiasTF=TRUE)
    W <- rbind(W, as.vector(fit1$w))
    WI <- rbind(WI, as.vector(fit1$wi))
    R <- c(R, fit1$rho)
  }

  if(plotTF){
      graphics::par(cex = .7)
      plot(fitzero$time.seq, fitzero$BIC, type = "l", ylab = "Negative log-likelihood", xlab = "Time")
      graphics::title("DCR")
      graphics::abline(v = br, col = "blue")
  }

  out = list()
  out$br <- br
  out$brhist <- brhist
  out$diffhist <- diffhist
  out$W <- W
  out$WI <- WI
  out$lambda <- R
  out$pvalhist <- pvalhist
  out$fitzero = fitzero;

  return(out)
}
