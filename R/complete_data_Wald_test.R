#' The multivariate Wald test statistic for complete data.
#'
#' @export

complete_data_Wald_test <- function(model, null_model){

  names_1 <- names(coef(model))
  names_0 <- names(coef(null_model))
  theta_names <- setdiff(names_1, names_0)

  theta_hat <- coef(model)
  theta_hat <- theta_hat[theta_names]

  i <- which(names_1 %in% theta_names)
  k <- length(i)

  U_hat <- vcov(model)
  U_hat <- U_hat[i, i]

  Wald <- t(theta_hat) %*% solve(U_hat) %*% theta_hat
  p_value <- pchisq(q = Wald, df = k, lower.tail = FALSE)
  out <- data.frame(Wald.statistic = Wald,
                    df = k,
                    p.value = p_value)
  out
}
