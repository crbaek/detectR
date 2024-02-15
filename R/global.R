#' @name global
#' @aliases global
#' @title Global Variables and functions
#' @description Defining ariables and functions used in the internal functions

utils::globalVariables(c("k","b","tau","mad","Mymodel", "networkmatrix"))
`%dopar` = foreach::`%dopar%`
`plot` = graphics::`plot`
`par` = graphics::`par`