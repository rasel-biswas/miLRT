#' Compare two nested models
#'
#' This method is our proposed method using
#' the Wilson-Hilferty transformation of the LRT statistics
#' This codes uses the `mice::pool.scalar()` function. This does the same thing as `miLRT::mi_lrt_wilson_hilferty()`.
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

mi_lrt_wilson_hilferty_another <- function(model, null_model, df="Rubin1987") {

  m <- length(model$analyses)
  k <- length(model$analyses[[1]]$coef) - length(null_model$analyses[[1]]$coef)
  n <- nobs(object = model$analyses[[1]])

  loglik1_m <- -2 * sapply(model$analyses, logLik)
  loglik0_m <- -2 * sapply(null_model$analyses, logLik)
  repeated_lrts <- loglik0_m - loglik1_m  # vector of the 'm' LRT statistics

  transformed_lrts <- wilson_hilferty(repeated_lrts, k)
  # the transformed variable is a standardized normal variable
  # with mean zero and variance one
  pool.scalar <- function (Q, U, n = Inf, k = 1) {
    m <- length(Q)
    qbar <- mean(Q)
    ubar <- mean(U)
    b <- var(Q)
    t <- ubar + (m + 1) * b/m
    barnard.rubin <- function(m, b, t, dfcom = 999999) {
      lambda <- (1 + 1 / m) * b / t
      lambda[lambda < 1e-04] <- 1e-04
      dfold <- (m - 1) / lambda^2
      dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
      dfold * dfobs / (dfold + dfobs)
    }
    df <- barnard.rubin(m, b, t, dfcom = n - k)
    r <- (1 + 1/m) * b/ubar
    fmi <- (r + 2/(df + 3))/(r + 1)
    list(m = m, qhat = Q, u = U, qbar = qbar, ubar = ubar, b = b,
         t = t, df = df, r = r, fmi = fmi)
  }

  pooled_rs <- pool.scalar(Q = transformed_lrts, U = 1, n = n, k = k)
  t_statistic <- pooled_rs$qbar / sqrt(pooled_rs$t)
  df <- pooled_rs$df

  p_value <- pt(q = t_statistic, df = df, lower.tail = FALSE)
  out <- data.frame(t.statistic = t_statistic,
                    df = df,
                    p.value = p_value)
  out
}
