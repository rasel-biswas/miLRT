#' Returns the exact transformation of a chi-square random variable
#' to the standard normal random variable.
#'
#' @export

transform_to_std_norm <- function(chi_square_rv, df) {
  qnorm(p = pchisq(q = chi_square_rv, df = df) )
}
