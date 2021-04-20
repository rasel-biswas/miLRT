#' Compare two nested models
#'
#' This method is our proposed method using
#' the Wilson-Hilferty transformation of the LRT statistics
#'
#'
#' @references
#' #' Meng, X. L., and D. B. Rubin. 1992.
#' Performing Likelihood Ratio Tests with Multiply-Imputed Data Sets.
#' \emph{Biometrika}, 79 (1): 103-11.
#'
#' Wilson, E.B., and Hilferty, M.M. (1931).
#' The distribution of chi-square.
#' \emph{Proceedings of the National Academy of Sciences of the United States of America},
#' 17, 684-688.
#'
#' @import mice
#' @export

mi_lrt_new <- function(model, null_model, df="Rubin1987") {

  m <- length(model$analyses)
  k <- length(model$analyses[[1]]$coef) - length(null_model$analyses[[1]]$coef)
  n <- nobs(object = model$analyses[[1]])

  loglik1_m <- -2 * sapply(model$analyses, logLik)
  loglik0_m <- -2 * sapply(null_model$analyses, logLik)
  repeated_lrts <- loglik0_m - loglik1_m  # vector of the 'm' LRT statistics

  transformed_lrts <- transform_to_std_norm(repeated_lrts, k)
  # the transformed variable is a standardized normal variable
  # with mean zero and variance one

  z_bar <- mean(transformed_lrts)
  B <- var(transformed_lrts)
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
