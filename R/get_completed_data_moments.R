#' This function computes completed-data moments.
#'
#'  (This is an internal function).
#'
#' @return A list containing the completed-data estimates and
#' their variance-covariance matrices.
#'
#' @import mice
#' @export

get_completed_data_moments <- function(model, null_model, sandwich = FALSE){
  names_1 <- names(coef(model$analyses[[1]]))
  names_0 <- names(coef(null_model$analyses[[1]]))
  theta_names <- setdiff(names_1, names_0)

  theta_hat <- sapply(model$analyses, coef)
  theta_hat <- theta_hat[theta_names, ]

  p <- length(names_1)
  i <- which(names_1 %in% theta_names)
  U_hat <- vapply(model$analyses,
                  ifelse(sandwich == FALSE, vcov, sandwich::vcovHC),
                  FUN.VALUE = matrix(0, p, p))
  U_hat <- U_hat[i, i,]

  out <- list(theta_hat = theta_hat, U_hat = U_hat)
  out
}
