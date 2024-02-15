#' @name bandpass
#' @aliases bandpass
#' @title Optional bandpass filtering using the Butterworth filter.
#' @description Id
#'
#' @keywords internal 


bandpass <- function(data, butterfreq) {
  bf1 <- signal::butter(5, butterfreq, type = "pass", plane = "z")
  for (i in 1:ncol(data)) {
    data[, i] <- signal::filter(bf1, data[, i] - mean(data[, i]))
  }
  return(data)
}
