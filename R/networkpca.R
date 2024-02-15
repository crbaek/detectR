#' If model is indicated, reduce data to components.
#'
#' @keywords internal

networkpca <- function(model, Y) {
  pcatable <- lavaan::lavParTable(Mymodel)
  pcatable <- pcatable[ which(pcatable$op == "=~"), ]
  netw_num <- length(unique(pcatable$lhs))
  netw_names <- unique(pcatable$lhs)

  # grab rhs for each network
  z <- list() # List because networks have different amounts of variables
  r <- 1 # how many components
  for (k in 1:netw_num) {
    id <- which(pcatable$lhs == netw_names[k])
    varnames <- pcatable[id, ]$rhs
    z[[k]] <- Y[, varnames]
  }


  new <- list() # new list created for compiling the sum of the Principle components
  newdf <- matrix(0, nrow(Y), length(netw_names)) # new data fram for concatanating the network sums

  for (i in 1:length(z)) {
    vvz <- svd(stats::cov(z[[i]]))
    Lamz <- vvz$u[, 1:1]
    zL <- t(Lamz) %*% t(z[[i]])
    newdf[, i] <- t(zL)
    newdf <- as.data.frame(newdf)
    colnames(newdf) <- netw_names
  }
  return(newdf)
}
