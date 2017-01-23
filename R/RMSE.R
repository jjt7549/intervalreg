#' The Root Mean-Square Error(RMSE)
#'
#' The lower bound root mean-square error and the upper bound root mean-square error proposed by Lima Neto and de Carvalho(2008) measure the differences between the predicted values and the observed values.
#' @param model \code{\link{imcmuni}} or \code{\link{imcmtn}} object, etc..
#'
#' @references Lima Neto, E.A. and De Carvalho, F.A.T(2010), Constrained linear regression models for symbolic interval-valued variables \emph{Computational Statistics and Data Analysis, 54}, 333-347
#'
#' @examples
#' set.seed(2017)
#' x1_L = rnorm(30, 3, 0.01) - rnorm(30, 0, 0.01)
#' x1_U = rnorm(30, 3, 0.01) + rnorm(30, 3, 0.01)
#' x2_L = runif(30, 1.5, 3) - runif(30, 0, 1)
#' x2_U = runif(30, 1.5, 3) + runif(30, 1, 2)
#' y_L = x1_L + x2_L
#' y_U = x1_U + x2_U
#' temp <- as.data.frame(cbind(y_L, y_U, x1_L, x1_U, x2_L, x2_U))
#' m1 <- imcmtn(cbind(y_L, y_U) ~ x1_L + x1_U + x2_L + x2_U, data = temp, b = 100)
#' RMSE(m1)
#' @export
RMSE <- function(model) {
  y = model$response
  yhat = model$fitted.values
  n = nrow(y)

  RMSE_L <- sqrt(sum((y[, 1]- yhat[, 1])^{2}) / n)
  RMSE_U <- sqrt(sum((y[, 2]- yhat[, 2])^{2}) / n)

  rmse <- cbind(RMSE_L, RMSE_U)
  return(rmse)
}
