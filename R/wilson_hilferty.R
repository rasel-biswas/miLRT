#' Returns the Wilson-Hilferty transformation of a chi-square random variable
#'
#' @references
#' Wilson, E.B., and Hilferty, M.M. (1931).
#' The distribution of chi-square.
#' \emph{Proceedings of the National Academy of Sciences of the United States of America},
#' 17, 684-688.
#'
#' @export

wilson_hilferty <- function(chi_square_rv, df) {
  ((chi_square_rv / df) ^ (1 / 3) - (1 - 2 / (9 * df))) / sqrt(2 / (9 * df))
}
