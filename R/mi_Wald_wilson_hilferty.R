#' Compare two nested models
#'
#' This method is our proposed method using
#' the Wilson-Hilferty transformation of the Wald statistics
#'
#'
#' @references
#' Wilson, E.B., and Hilferty, M.M. (1931).
#' The distribution of chi-square.
#' \emph{Proceedings of the National Academy of Sciences of the United States of America},
#' 17, 684-688.
#'
#' @import mice
#' @export

mi_Wald_wilson_hilferty <- function(model, null_model, df="Rubin1987") {

  m <- length(model$analyses)
  k <- length(model$analyses[[1]]$coef) - length(null_model$analyses[[1]]$coef)
  n <- nobs(object = model$analyses[[1]])

  moments <- get_completed_data_moments(model, null_model)

  theta_hat <- moments$theta_hat
  U_hat <- moments$U_hat

  Wald_statistics <- rep(NA, m)
  for (j in 1:m) {
    Wald_statistics[j] <-
      t(theta_hat[, j]) %*% solve(U_hat[, , j]) %*% theta_hat[, j]
  }  # vector of the 'm' Wald statistics

  transformed_Walds <- wilson_hilferty(Wald_statistics, k)
  # the transformed variable is a standardized normal variable
  # with mean zero and variance one

  z_bar <- mean(transformed_Walds)
  B <- var(transformed_Walds)
  T <- 1 + (1 + m ^ -1) * B
  t_statistic <- z_bar / sqrt(T)

  if(df=="BR1999") {
    lambda <- ((1 + m ^ -1) * B) / (1 + ((1 + m ^ -1) * B))
    lambda[lambda < 1e-04] <- 1e-04
    v_old <- (m-1) * lambda ^ -2
    v_complete <- n - k
    v_observed <- (v_complete + 1) / (v_complete + 3) * v_complete *
      (1 - ((1 + m ^ -1) * B) / T)
    df_BR <- (v_old * v_observed) / (v_old + v_observed)
    df <- df_BR
  } else {
    lambda <- ((1 + m ^ -1) * B) / (1 + ((1 + m ^ -1) * B))
    v_old <- (m-1) * lambda ^ -2
    #v_old <- (m - 1) * { 1 + 1 / ((1 + m ^ -1) * B) } ^ 2
    df <- v_old
  }
  p_value <- pt(q = t_statistic, df = df, lower.tail = FALSE)
  out <- data.frame(t.statistic = t_statistic,
                    df = df,
                    p.value = p_value)
  out
}
