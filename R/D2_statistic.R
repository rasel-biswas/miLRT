#' Compare two nested models using D2-statistic
#'
#' The D2 statistic pools test statistics from the repeated analyses.
#' The method is less powerful than the D1 and D3 statistics.
#'
#' @references
#' Li, K. H., X. L. Meng, T. E. Raghunathan, and D. B. Rubin. 1991.
#' Significance Levels from Repeated p-Values with Multiply-Imputed Data.
#' \emph{Statistica Sinica} 1 (1): 65-92.
#'
#' @import mice
#' @export

D2_statistic <- function(model, null_model) {
  m <- length(model$analyses)
  moments <- get_completed_data_moments(model, null_model)

  theta_hat <- moments$theta_hat
  U_hat <- moments$U_hat

  chi_square_statistics <- rep(NA, m)
  for (j in 1:m) {
    chi_square_statistics[j] <-
      t(theta_hat[, j]) %*% solve(U_hat[, , j]) %*% theta_hat[, j]
  }

  dW_bar <- mean(chi_square_statistics)

  r <- (1 + m ^ (-1)) * var(sqrt(chi_square_statistics))
  k <- dim(theta_hat)[1]

  D2 <- (dW_bar / k - (m + 1) / (m - 1) * r) / (1 + r)
  D2
}
