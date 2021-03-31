#' The most straightforward Wald statistic
#'
#' The most straightforward Wald statistic which does not need the
#' Equal Fraction of Missing Information" assumption.
#'
#' @import mice
#' @export

D0_naive_test <- function(model, null_model) {
  ## This test do not need the EFMI assumption

  m <- length(model$analyses)
  moments <- get_completed_data_moments(model, null_model)

  theta_bar <- apply(moments$theta_hat, 1, mean)
  U_bar <- apply(moments$U_hat, c(1, 2), mean)

  B <- cov(t(moments$theta_hat))
  T_overall <- U_bar + (1+m^(-1))* B
  k <- dim(B)[1]
  D0 <- t(theta_bar) %*% solve(T_overall) %*% theta_bar
  p_value <- pchisq(q = D0, df = k, lower.tail = FALSE)
  out <- data.frame(Wald.statistic = D0,
                    df = k,
                    p.value = p_value)
  out
}
