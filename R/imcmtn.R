#' MCM for Interval-Valued Data Regression Using Truncated normal distribution
#'
#' \code{imcmtn()} is used to fit a linear regression model based on the Monte Carlo Method using truncated normal distribution.
#' @param formula an object of class \code{\link[stats]{formula}}, a symbolic description of the model to be fitted.
#' @param data an data frame containing the variables in the model.
#' @param b number of resampling (default : 100)
#'
#' @details Similar to \code{imcmuni()}, but uses a truncated normal distribution rather than a uniform distribution.
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
#' @return \item{resampling.coefficients}{B coefficient vectors obtained by resampling.}
#' @return \item{coefficients}{Average of B coefficients obtained by resampling, standard error and p-value.}
#' @return \item{fitted.values}{The fitted values for the lower and upper interval bound.}
#' @return \item{residuals}{The residuals for the lower and upper interval bound.}
#'
#' @references Ahn, J., Peng, M., Park, C., Jeon, Y.(2012), A Resampling Approach for Interval-Valued Data Regression. \emph{Statistical Analysis and Data Mining, 5, 336-348}
#' @references Lim, S., Kang, K.(2016), On Statistical Analysis of Interval-Valued Data
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
#' m1
#'
#' @seealso \code{\link[truncnorm]{rtruncnorm}} \code{\link{RMSE}} \code{\link{symbolic.r}}
#' @import stats
#' @importFrom truncnorm rtruncnorm
imcmtn <- function(formula = formula, data = data, b = 100) {

  requireNamespace("truncnorm", quietly = TRUE)
  #use_package("truncnorm")
  #library(truncnorm)

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

  p = ncol(x) - 1 # num. of predictor variable except intercept

  ### Resampling STEP (by using truncated normal distribution) ###
  # re-sampling response variable
  re_y <- matrix(NA, n, b)

  for (k in 1:b) {
    for (i in 1:n) {
      mean_y = ((y[i, 1] + y[i, 2]) / 2)
      sigma = (y[i, 2] - y[i, 1])^{2} / 12
      alpha = (y[i, 1] - mean_y) / sigma
      beta = (y[i, 2] - mean_y) / sigma
      re_y[i, k] = rtruncnorm(1, a = y[i, 1], b = y[i, 2],
                              mean = mean_y,
                              sd = sqrt(sigma^2 * (1 + (alpha*dnorm(alpha) - beta*dnorm(beta))/(pnorm(beta) - pnorm(alpha)) - ((dnorm(alpha) - dnorm(beta))/(pnorm(beta) - pnorm(alpha)))^{2})))
    }
  }

  # re-sampling predictor variables
  re_xx <- array(NA, c(n, {p/2}, b))
  re_x <- array(NA, c(n, {p/2}+1, b)) # include intercept term

  for (k in 1:b) {
    for (j in 1:{p/2}) {
      for (i in 1:n) {
        mean_x = (x[i, {2*j}] + x[i, {2*j}+1]) / 2
        sigma = (x[i, {2*j}+1] - x[i, {2*j}])^{2} / 12
        alpha = (x[i, {2*j}] - mean_x) / sigma
        beta = (x[i, {2*j}+1] - mean_x) / sigma
        re_xx[i, j, k] = rtruncnorm(1, a = x[i, {2*j}], b = x[i, {2*j}+1],
                                    mean = mean_x,
                                    sd = sqrt(sigma^2 * (1 + (alpha*dnorm(alpha) - beta*dnorm(beta))/(pnorm(beta) - pnorm(alpha)) - ((dnorm(alpha)-dnorm(beta))/(pnorm(beta) - pnorm(alpha)))^{2})))
      }
    }
    re_x[, , k] <- cbind(matrix(rep(1, n), ncol = 1), re_xx[, , k])
  }


  ### coefficient estimate STEP ###
  # MCM beta matrix
  re_coef <- matrix(NA, {p/2}+1, b)

  for (k in 1:b) {
    re_coef[, k] <- solve(t(re_x[, , k]) %*% re_x[, , k]) %*% t(re_x[, , k]) %*% re_y[, k]
  }

  estcoef <- apply(re_coef, 1, mean)

  # coefficient standard error
  coefse <- matrix(NA, nrow = {p/2}+1)

  for (j in 1:{p/2}) {
    for (k in 1:b) {
      coefse[j+1, ] = sqrt(sum((re_coef[j+1, k] - estcoef[j+1])^{2}) / (b-1))
    }
  }

  estcoef <- matrix(estcoef, ncol = 1)
  coef <- cbind(estcoef, coefse)

  var.names = colnames(x)[(1:{p/2})*2]
  rownames(coef) <- c("intercept", var.names)

  # p-value
  p.value <- NULL
  for (j in 1:{p/2}) {
    z = coef[j+1, 1] / coef[j+1, 2]
    p.value[j+1] =2 * (1 - pnorm(abs(z)))
  }

  coef <- cbind(coef, p.value)
  colnames(coef) <- c("Coeff.", "S.E.", "p.value")
  rownames(re_coef) <- c("intercept", var.names)

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


  ### final output STEP ###
  result <- list(call = match.call(),
                 response = y,
                 predictor = x,
                 resampling.coefficients = re_coef,
                 coefficients = coef,
                 fitted.values = fitted_values,
                 residuals = residuals)
  return(result)
}
