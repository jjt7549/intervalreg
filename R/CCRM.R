#' Constrained Center and Range Method for Interval-Valued Data Regression
#'
#' \code{CCRM()} is used to fit a linear regression model based on the inequality constraints over the range variables.(Lima Neto and De Carvalho, 2010)
#' @param formula an object of class \code{\link[stats]{formula}}, a symbolic description of the model to be fitted.
#' @param data an data frame containing the variables in the model.
#'
#' @details
#' The disadvantages of the CRM model is that the predicted valu of the range model may be negative. Similar to Center and Range Method(CRM), but adds constraint condition that all estimative of the parameters of the range's model are positive.(based on inequality constraints)
#' There is no constraints over the parameters estimates for the center point regression equation.
#'
#' @note In dataset, a pair of the interval variables should always be composed in order from lower to upper bound. In order to apply this function, the data should be composed as follows:
#' \tabular{cccccc}{
#' \eqn{y_L}     \tab  \eqn{y_U}     \tab  \eqn{x1_L}    \tab  \eqn{x1_U}    \tab  \eqn{x2_L}    \tab  \eqn{x2_U}\cr
#' \eqn{y_L1}  \tab  \eqn{y_U1}  \tab  \eqn{x_L11} \tab  \eqn{x_U11} \tab  \eqn{x_L12} \tab  \eqn{x_U12}\cr
#' \eqn{y_L2}  \tab  \eqn{y_U2}  \tab  \eqn{x_L21} \tab  \eqn{x_U21} \tab  \eqn{x_L22} \tab  \eqn{x_U22}\cr
#' \eqn{y_L3}  \tab  \eqn{y_U3}  \tab  \eqn{x_L31} \tab  \eqn{x_U31} \tab  \eqn{x_L32} \tab  \eqn{x_U32}\cr
#' \eqn{y_L4}  \tab  \eqn{y_U4}  \tab  \eqn{x_L41} \tab  \eqn{x_U41} \tab  \eqn{x_L42} \tab  \eqn{x_U42}\cr
#' \eqn{y_L5}  \tab  \eqn{y_U5}  \tab  \eqn{x_L51} \tab  \eqn{x_U51} \tab  \eqn{x_L52} \tab  \eqn{x_U52}\cr
#' }
#' @note The upper limit value of the variable should be unconditionally greater than the lower limit value. Otherwise, it will be output as \code{NA} or \code{NAN}, and the value can not be generated.
#'
#' @return \item{coefficients.Center}{Coefficients for the Center variable.}
#' @return \item{coefficeints.Range}{Coefficients for the Range variable.}
#' @return \item{fitted.values}{The fitted values for the lower and upper interval bound.}
#' @return \item{residuals}{The residuals for the lower and upper interval bound.}
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
#' m1 <- CCRM(cbind(y_L, y_U) ~ x1_L + x1_U + x2_L + x2_U, data = temp)
#' m1
#' @seealso \code{\link{RMSE}} \code{\link{symbolic.r}}
#' @import stats
#' @importFrom mgcv pcls
#' @export
CCRM <- function(formula = formula, data = data) {

  # Non-negative Least Squares, Penalized Constrained Least Squares fitting
  requireNamespace("mgcv", quietly = TRUE)
  #use_package("mgcv")
  #library(mgcv)

  n = nrow(data)

  # input formula for regression
  formula1 = as.formula(formula)
  model_interval = model.frame(formula = formula1, data = data)

  # response variable interval value
  yinterval = model.response(model_interval)
  y = as.matrix(yinterval)
  y.C = (y[, 1] + y[, 2]) / 2   # center point of y
  y.R = y[, 2] - y[, 1]         # range point of y

  # predictor variable interval matrix (each left: Lower, each right: Upper)
  xinterval = model.matrix(attr(model_interval, "terms"), data = data)
  x = as.matrix(xinterval)

  p = ncol(x) - 1 # num. of predictor variable except intercept
  var.names = colnames(x)[(1:{p/2})*2]

  x.C <- matrix(NA, n, {p/2})   # center point of x
  colnames(x.C) <- var.names
  x.R <- matrix(NA, n, {p/2})   # range point of x
  colnames(x.R) <- var.names
  for (j in 1:{p/2}) {
    for (i in 1:n) {
      x.C[i, j] = (x[i, 2*j] + x[i, 2*j+1]) / 2
      x.R[i, j] = (x[i, 2*j+1] - x[i, 2*j])
    }
  }

  ### coefficient estimate STEP ###
  model.C = model.frame(y.C ~ x.C)
  y.C = model.response(model.C)
  x.C = model.matrix(attr(model.C, "terms"))
  model.R = model.frame(y.R ~ x.R)
  y.R = model.response(model.R)
  x.R = model.matrix(attr(model.R, "terms"))

  var.names = colnames(x)[(1:{p/2})*2]
  colnames(x.C) <- c("Intercept", var.names)
  colnames(x.R) <- c("Intercept", var.names)

  # Center point
  coef.C <- solve.qr(qr(x.C), y.C)
  coef.C <- as.matrix(coef.C)
  colnames(coef.C) <- c("coefficients")
  df.C = nrow(x.C) - ncol(x.C)  # degree of freedom

  # Range point (constraint valid)
  pcls.list <- list(y = y.R,
                    w = rep(1, n),
                    X = x.R,
                    C = matrix(0, 0, 0),
                    S = list(),
                    off = numeric(0),
                    sp = numeric(0),
                    p = rep(1, p - 1),
                    Ain = diag(rep(1, p - 1)),
                    bin = rep(0, p - 1))
  coef.R <- pcls(pcls.list)
  coef.R <- as.matrix(coef.R)
  rownames(coef.R) <- c("Intercept", var.names)
  colnames(coef.R) <- c("coefficients")
  df.R = nrow(x.R) - ncol(x.R)  # degree of freedom

  # fitted values (Lower & Upper)
  fitted_L <- (x.C %*% coef.C) - ((x.R %*% coef.R) / 2)
  fitted_U <- (x.C %*% coef.C) + ((x.R %*% coef.R) / 2)
  fitted_values <- cbind(fitted_L, fitted_U)
  colnames(fitted_values) <- c("fitted.Lower", "fitted.Upper")

  # residuals (Lower & Upper)
  residual_L <- (y.C - y.R / 2) - fitted_L
  residual_U <- (y.C + y.R / 2) - fitted_U
  residuals <- cbind(residual_L, residual_U)
  colnames(residuals) <- c("resid.Lower", "resid.Upper")

  ### final STEP ###
  result <- list(call = match.call(),
                 response = y,
                 predictor = x,
                 coefficients.Center = coef.C,
                 coefficeints.Range = coef.R,
                 fitted.values = fitted_values,
                 residuals = residuals)

  return(result)

}
