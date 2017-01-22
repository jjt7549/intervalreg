#' Symbolic Covariance Method
#'
#' \code{CRM()} is used to fit a linear regression model based on symbolic covariance matrix(Xu, 2010).
#' @param formula an object of class \code{\link[stats]{formula}}, a symbolic description of the model to be fitted.
#' @param data an data frame containing the variables in the model.
#'
#' @details
#' The SCM proposed by Xu(2010) is a method of estimating a regression coefficient using a symbolic covariance matrix. In SCM, the centralized linear regression model is used(model with centered variables).
#' The regression coefficient is estimated by least squares method, but it is used by symbolic covariance matrix. Because the process of calculating the symbolic sample covariance uses the lower and upper limits of each variable, the SCM reflects the variability of the interval.
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
#' @return \item{symbolic.covariance.Sxx}{Symbolic sample variance-covariance matrix between response variable Y and predictor variables X.}
#' @return \item{symbolic.covariance.Sxy}{Symbolic sample covariance vector between response variable Y and predictor variables X}
#' @return \item{coefficients}{regression coefficients}
#' @return \item{fitted.values}{The fitted values for the lower and upper interval bound.}
#' @return \item{residuals}{The residuals for the lower and upper interval bound.}
#'
#' @references Xu, W.(2010), Symbolic Data Analysis: Interval-Valued Data Regression
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
#' m1 <- SCM(cbind(y_L, y_U) ~ x1_L + x1_U + x2_L + x2_U, data = temp)
#' m1
#' @seealso \code{\link{RMSE}} \code{\link{symbolic.r}}
#' @import stats
SCM <- function(formula = formula, data = data) {

  n = nrow(data)

  # input formula for regression
  formula1 = as.formula(formula)
  model_interval = model.frame(formula = formula1, data = data)

  # response variable interval value
  yinterval = model.response(model_interval)
  y = as.matrix(yinterval)

  # predictor variable interval matrix (each left: Lower, each right: Upper)
  xinterval = model.matrix(attr(model_interval, "terms"), data = data)
  x = as.matrix(xinterval)

  p = ncol(x) - 1 # except num. of intercept

  ### SCM STEP ###

  # symbolic sample mean
  mean_x <- NULL
  for (j in 1:{p/2}) {
    mean_x[j] = sum(x[, {2*j}] + x[, {2*j+1}]) / (2*n)
  }

  mean_y <- NULL
  for (i in 1:n) {
    mean_y = sum(y[, 1] + y[, 2]) / (2*n)
  }

  # set lower & upper bounds
  x_L <- x_U <- matrix(NA, nrow = n, ncol = {p/2})
  for (j in {1:p/2}) {
    x_L[, j] = x[, 2*j]
    x_U[, j] = x[, 2*j+1]
  }
  y_L = y[, 1]
  y_U = y[, 2]

  # symbolic covariance
  symbolic.cov <- function(lower, upper) {

    n = nrow(lower)
    p = ncol(lower)

    x_mean = colSums(lower + upper) / (2*n)

    sym.c = matrix(NA, nrow = p, ncol = p)

    for (j in 1:p) {
      for (k in 1:p) {
        sym.c[j, k] = ( sum(2*(x_L[,j] - x_mean[j])*(x_L[,k] - x_mean[k]) +
                              (x_L[,j] - x_mean[j])*(x_U[,k] - x_mean[k]) +
                              (x_U[,j] - x_mean[j])*(x_L[,k] - x_mean[k]) +
                              2*(x_U[,j] - x_mean[j])*(x_U[,k] - x_mean[k]) ) / ( 6*n ) )
      }
    }

    return(sym.c)

  }

  # S_xx
  S_xx <- symbolic.cov(x_L, x_U)

  # S_xy
  S_xy <- matrix(NA, nrow = {p/2}, ncol = 1)
  for (j in 1:{p/2}) {
    S_xy[j, ] <- sum(2*(x_L[, j] - mean_x[j])*(y_L - mean_y) +
                       (x_L[, j] - mean_x[j])*(y_U - mean_y) +
                       (x_U[, j] - mean_x[j])*(y_L - mean_y) +
                       2*(x_U[, j] - mean_x[j])*(y_U - mean_y)) / (6*n)
  }

  ### coefficients estimate STEP ###
  beta <- t(solve(S_xx) %*% S_xy)
  intercept <- mean_y - sum(mean_x * beta)

  estcoef <- matrix(c(intercept, beta), ncol = 1)
  var.names = colnames(x)[(1:{p/2})*2]
  rownames(estcoef) = c("intercept", var.names)
  colnames(estcoef) = c("coefficients")

  # fitted values (Lower & Upper)
  fitted_L <- fitted_U <- NULL
  for (i in 1:n) {
    fitted_L[i] <- min(x[i, c(1, seq(2, p, by = 2))] %*% c(t(estcoef[, 1])), x[i, c(1, seq(3, p+1, by = 2))] %*% c(t(estcoef[, 1])))
    fitted_U[i] <- max(x[i, c(1, seq(2, p, by = 2))] %*% c(t(estcoef[, 1])), x[i, c(1, seq(3, p+1, by = 2))] %*% c(t(estcoef[, 1])))
  }
  fitted_values <- cbind(fitted_L, fitted_U)
  colnames(fitted_values) <- c("fitted.Lower", "fitted.Upper")

  # residuals (Lower & Upper)
  residual_L <- y[, 1] - fitted_L
  residual_U <- y[, 2] - fitted_U
  residuals <- cbind(residual_L, residual_U)
  colnames(residuals) <- c("resid.Lower", "resid.Upper")

  ### final STEP ###
  result <- list(call = match.call(),
                 response = y,
                 predictor = x,
                 symbolic.covariance.Sxx = S_xx,
                 symbolic.covariance.Sxy = S_xy,
                 coefficients = estcoef,
                 fitted.values = fitted_values,
                 residuals = residuals)

  return(result)

}
