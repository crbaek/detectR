#' @name preprocess
#' @aliases preprocess
#' @title Data preparation for changepoint detection using functions in this package..
#' @description Id
#' @usage
#' preprocess(file = NULL,
#' header = NULL,
#' sep    = NULL,
#' signal = NULL,
#' noise = NULL,
#' butterfreq = NULL,
#' model = NULL)
#'
#' @param file a data matrix or file name with columns as variables and rows as observations across time.
#' @param header logical for whether or not there is a header in the data file.
#' @param sep The spacing of the data files.
#' "" indicates space-delimited,
#' "/t" indicates tab-delimited, "," indicates comma delimited. Only necessary
#' to specify if reading data in from physical directory.
#' @param signal (optional) a character vector containg the names of variables that contain signal i.e., which variables to use to detect change point. The default
#' (NULL) indicates all variables except those in 'noise' argument are considered signal. Example: signal = c("dDMN4", "vDMN5", "vDMN1",
# "antSal2", "antSal5", "postSal9", "postSal12") uses only these variables to detect change point.
#' @param noise (optional) a character vector containg the names of variables that contain noise. The signal variables will be regressed on these variables and residuals
#' used in change point detection. The deault (NULL) indicates there are no noise variables. Example: noise = c("White.Matter1", "CSF1")
#' @param butterfreq (optional) bandpass filter frequency ranges. Example: c(.04,.4)
#' @param model (optional) syntax indicating which variables belong to which networks for first pass of data reduction that is user-specified. If no header naming convention follows "V#". Notation should follow lavaan syntax style.

preprocess <- function(file = NULL,
                       header = NULL,
                       sep = NULL,
                       signal = NULL,
                       noise = NULL,
                       butterfreq = NULL,
                       model = NULL) {

  # read in/ assign data to Y
  if (is.data.frame(file)) {
    Y <- file
  } else {
    # Error check for header
    if (!is.data.frame(file) & is.null(header)) {
      stop(paste0(
        "detectR ERROR: a data file name is specified but a header argument is not. ",
        "Please specify a logical value for 'header' argument before continuing or provide data frame."
      ))
    }

    # Error check for sep
    if (!is.data.frame(file) & is.null(sep)) {
      stop(paste0(
        "detectR ERROR: a data file name is specified but a seperation argument is not. ",
        "Please specify 'sep' argument before continuing or provide data frame."
      ))
    }
    Y <- utils::read.table(file = file, sep = sep, header = header)
  }

  varnames <- c(signal, noise)
  Y <- (Y[names(Y)[names(Y) %in% varnames]])

  # preprocessing
  if (is.null(signal)) {
    signal <- colnames(Y)[!(colnames(Y) %in% noise)]
  }

  if (!is.null(noise)) {
    networkmatrix <- preproc(Y, noise, signal)
  } else {
    networkmatrix <- Y
  }

  # bandpass filtering
  if (!is.null(butterfreq)) {
    networkmatrix <- bandpass(networkmatrix, c(.04, .4))
  }

  if (!is.null(model)) {
    networkmatrix <- networkpca(model, networkmatrix)
  }

  out <- list(data.used = networkmatrix)
  return(out)
}
