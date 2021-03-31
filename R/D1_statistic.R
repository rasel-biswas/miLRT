#' Compare two nested models using D1-statistic
#'
#' The D1 statistic is the multivariate Wald test.
#'
#' @references
#' Li, K. H., T. E. Raghunathan, and D. B. Rubin. 1991.
#' Large-Sample Significance Levels from Multiply Imputed Data Using
#' Moment-Based Statistics and an F Reference Distribution.
#' \emph{Journal of the American Statistical Association}, 86(416): 1065-73.
#'
#' @import mice
#' @export

D1_statistic <- function(model, null_model) {
  m <- length(model$analyses)
  moments <- get_completed_data_moments(model, null_model)

  theta_bar <- apply(moments$theta_hat, 1, mean)
  U_bar <- apply(moments$U_hat, c(1, 2), mean)

  B <- cov(t(moments$theta_hat))

  k <- length(theta_bar)
  rm <- (1 + m ^ (-1)) * sum(diag(B %*% solve(U_bar))) / k
  D1 <- t(theta_bar) %*% solve(U_bar) %*% theta_bar / (k * (1 + rm))
  D1
}
