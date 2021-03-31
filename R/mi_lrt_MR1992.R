#' Compare two nested models using MI LRT
#'
#' The MI LRT (D3) test is a variant of likelihood ratio by Meng and Rubin (1992).
#'
#' @references
#' Meng, X. L., and D. B. Rubin. 1992.
#' Performing Likelihood Ratio Tests with Multiply-Imputed Data Sets.
#' \emph{Biometrika}, 79 (1): 103-11.
#'
#' @import mice
#' @export

mi_lrt_MR1992 <- function(model, null_model) {

  loglik1_m <- -2 * sapply(model$analyses, logLik)
  loglik0_m <- -2 * sapply(null_model$analyses, logLik)
  dbar_m <- mean(loglik0_m - loglik1_m)  # mean of the LR test statistics

  est1 <- mice::pool(model)          # using Rubin's rule
  est0 <- mice::pool(null_model)     # using Rubin's rule

  estimates_bar1 <- est1$pooled$estimate
  estimates_bar0 <- est0$pooled$estimate

  m <- length(model$analyses)
  k <- length(estimates_bar1) - length(estimates_bar0)

  # Next step: fitting full and null models with their parameters fixed to
  # \psi_1, and \psi_0, respectively. And calculate
  # the avg LR test statistics (under those constraints).

  fix_coef_gaussian <- function(model, psi) {
    beta <- psi$beta # final estimates of coefficients
    sig2 <- psi$sigma2 # final estimate of sigma^2
    ytrm <- attr(model$terms, "variables")[-1][attr(model$terms, "response")]
    y <- as.matrix(model$model[as.character(ytrm)]) # response vector, y
    X <- model.matrix(model$terms, model$model) # design matrix, X
    n <- length(y)
    L <- -(n / 2) * log(2 * pi * sig2) -
      (1 / (2 * sig2)) * t(y - X %*% beta) %*% (y - X %*% beta)
    as.numeric(L)
  }

  fix_coef_glm <- function(model, beta) {
    data <- model.frame(model)
    design_matrix <- model.matrix(formula(model), data = data)
    offset <- as.vector(design_matrix %*% beta)
    updated_model <- update(
      model,
      formula. = . ~ 1,
      offset = offset,
      data = cbind(data, offset = offset)
    )
    logLik(updated_model)
  }

  if (family(model$analyses[[1]])[1]$family == "gaussian") {
    sigma2_1 <- mean(sapply(model$analyses, resid) ^ 2)
    sigma2_0 <- mean(sapply(null_model$analyses, resid) ^ 2)
    psi_1 <- list(beta = estimates_bar1, sigma2 = sigma2_1)
    psi_0 <- list(beta = estimates_bar0, sigma2 = sigma2_0)
    loglik1_L <- sapply(model$analyses, fix_coef_gaussian, psi = psi_1)
    loglik0_L <- sapply(null_model$analyses, fix_coef_gaussian, psi = psi_0)
  } else {
    loglik1_L <- sapply(model$analyses, fix_coef_glm, beta = estimates_bar1)
    loglik0_L <- sapply(null_model$analyses, fix_coef_glm, beta = estimates_bar0)
  }

  dbar_L <- -2 * mean(loglik0_L - loglik1_L)

  rm <- ((m + 1) / (k * (m - 1))) * (dbar_m - dbar_L)
  Dm <- dbar_L / (k * (1 + rm))
  v <- k * (m - 1)
  if (v > 4) {
    w <- 4 + (v - 4) * ((1 + (1 - 2 / v) * (1 / rm)) ^ 2)
  } else {
    w <- v * (1 + 1 / k) * ((1 + 1 / rm) ^ 2) / 2
  }
  pvalue = ifelse(Dm>=0 & rm>=0, 1 - pf(Dm, k, w), NA)
  out <-
    data.frame(
      F.statistic = Dm,
      df1 = k,
      df2 = w,
      p.value = pvalue,
      RIV = rm
    )
  out
}
