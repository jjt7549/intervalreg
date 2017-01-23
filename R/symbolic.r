#' The Symbolic Correlation
#'
#' The symbolic correlation coefficient, \eqn{r}, proposed by Billard(2007, 2008) and applied to the regression problem by Xu(2010), measures the correlation between the predicted values and the observed values.
#' @param model \code{\link{imcmuni}} or \code{\link{imcmtn}} object, etc..
#'
#' @references Billard(2007, 2008). Sample covariance functions for complex quantitative data
#' @references Xu(2010), Symbolic Data Analysis: Interval-Valued Data Regression
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
#' symbolic.r(m1)
#' @export
symbolic.r <- function(model) {
  x = model$predictor
  y = model$response
  haty = model$fitted.values

  # sample deviation
  sample.dev <- function(lower, upper) {
    n = length(lower)
    s.cov = (sum(lower^2 + lower*upper + upper^2) / (3*n)) - (sum(lower + upper)^{2} / (4*n^{2}))

    return(s.cov)
  }

  n = nrow(x)
  meany = sum(y[, 1:2]) / (2*n)
  meanhaty = sum(haty[, 1:2]) / (2*n)

  # symbolic covariance
  sym.c = ( sum(2*(y[, 1] - meany)*(haty[, 1] - meanhaty) +
                  (y[, 1] - meany)*(haty[, 2] - meanhaty) +
                  (y[, 2] - meany)*(haty[, 1] - meanhaty) +
                  2*(y[, 2] - meany)*(haty[, 2] - meanhaty)) ) / ( 6*n )

  sym.r = sym.c / (sqrt(sample.dev(y[, 1], y[, 2])) * sqrt(sample.dev(haty[, 1], haty[, 2])))

  return(sym.r)
}
