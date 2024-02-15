
### Auxiliary functions#############

#' @name blocklength.andrew
#' @title Block length formula by Andrews
#' @keywords internal

blocklength.andrew <- function(data) {
  n <- length(data)
  rhohat <- stats::ar.ols(data, FALSE, order = 1)
  rhohat <- rhohat$ar[1]
  q <- 1.1447 * (n * 4 * rhohat^2 / (1 - rhohat^2)^2)^(1 / 3)
  return(round(q))
}

#' @name LRV.bartlett
#' @title (multivariate) Bartlett long-run variance calculation
#' @keywords internal
#' 
LRV.bartlett <- function(cL, q) {
  # input data is dim*length
  n <- dim(cL)[2]
  
  if (missing(q)) {
    q <- floor(n^(1 / 3))
    #    q = apply(cL, 1, blocklength.andrew);
    #    q = floor(stats::median(q));
    #    q = min(q, floor(n/3));
  } # Data dependent bw
  
  if (q == 0) {
    sn2 <- stats::cov(t(cL))
  } else {
    wq <- seq(from = q, to = 1, by = -1) / (q + 1)
    xcov <- stats::acf(t(cL), lag.max = q, type = "covariance", plot = FALSE, na.action = stats::na.fail)$acf
    sn2 <- xcov[1, , ]
    for (i in 2:(q + 1)) {
      sn2 <- sn2 + wq[i - 1] * (xcov[i, , ] + t(xcov[i, , ]))
    }
  }
  return(list(sn2 = sn2, q = q))
}
