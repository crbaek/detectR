#' @name preproc
#' @aliases preproc
#' @title Optional preprocessing step for regressing out noise variables.
#' @description This function saves the residuals after regressing the user-defined
#' noise variables from the user-defined "signal" variables. Noise variables can
#' be nuisance regressors (such as cebral spinal fluid in fMRI), time vectors,
#' or any other vector that the user wishes to regress out.
#' @keywords internal
#' @param Y data: Y = length*dim
#' @param noise string vector indicating which variables (columns) are regressors
#' @param signal string vector indicating which variables to regress out the noise variables from
#' @return \strong{signalmatrix}

preproc <- function(Y, noise, signal) {
  Y = as.data.frame(Y)
  networkframe <- data.frame(Y[, which(colnames(Y) %in% noise)], seq(nrow(Y)))
  nt1 <- as.matrix(networkframe, nrow = nrow(Y), byrow = TRUE)
  signalmatrix <- NULL
  for (i in signal) {
    signalmatrix <- as.data.frame(cbind(networkmatrix, stats::lm(Y[, i] ~ 1 + nt1)$residuals))
  }
  colnames(signalmatrix) <- signal
  return(signalmatrix)
}



