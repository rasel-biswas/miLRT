#' Compare two nested models using D4 test
#'
#' This is the improved MI LRT by Chan and Meng (2019), of the original test by Meng and Rubin (1992).
#'
#'
#' @import mice
#' @export

mi_lrt_CM2019 <- function(model, null_model) {

  m <- length(model$analyses)

  all_data_rbinded <- model.frame(model$analyses[[1]])
  for (i in 2:m) {
    all_data_rbinded <-
      rbind(all_data_rbinded, model.frame(model$analyses[[i]]))
  }

  model_all <-
    glm(
      formula = formula(model$analyses[[1]]),
      data = all_data_rbinded,
      family = family(model$analyses[[1]])
    )
  model_all_null <-
    glm(
      formula = formula(null_model$analyses[[1]]),
      data = all_data_rbinded,
      family = family(null_model$analyses[[1]])
    )

  dhat_S <- 2 * as.numeric(logLik(model_all) - logLik(model_all_null)) / m


  k <- length(model_all$coefficients) - length(model_all_null$coefficients)
  h <- length(model_all$coefficients) + 1

  deltabar_S <- 2 * mean(sapply(model$analyses, logLik))
  deltahat_S <- 2 * as.numeric(logLik(model_all)) / m

  rm <- ((m + 1) / (h * (m - 1))) * (deltabar_S - deltahat_S)
  Dm <- dhat_S / (k * (1 + rm))


  fm <- rm / (1 + rm)
  df_hat <- h * (m - 1) / (fm^2)   ## from equation 2.17 of the paper. (robust)

  p_value <- 1 - pf(Dm, k, df_hat)
  out <-
    data.frame(
      Dm.statistic = Dm,
      df1 = k,
      df2 = df_hat,
      p.value = p_value,
      RIV = rm
    )
  out
}
