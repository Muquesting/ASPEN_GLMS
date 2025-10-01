#' Logit transform with clipping
#'
#' @param p Probabilities in (0, 1).
#' @param eps Numerical guard to avoid infinities.
#' @return Transformed values on the logit scale.
logit <- function(p, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

#' Inverse-logit transform
#'
#' @param x Values on the logit scale.
#' @return Probabilities in (0, 1).
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

#' Convert logit coefficient to allelic ratio
#'
#' @param beta Value on the logit scale (e.g. model coefficient).
#' @return Allelic ratio between 0 and 1.
allelic_ratio_from_logit <- function(beta) {
  inv_logit(beta)
}

#' Convert allelic ratio to logit scale
#'
#' @param ar Allelic ratio between 0 and 1.
#' @return Logit-transformed value.
logit_from_allelic_ratio <- function(ar) {
  logit(ar)
}
