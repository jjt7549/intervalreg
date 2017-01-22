#' MCM for Interval-Valued Data Regression Using Uniform distribution
#'
#' \code{imcmuni()} is used to fit a linear regression model based on the Monte Carlo Method using uniform distribution.
#' @param formula an object of class \code{\link[stats]{formula}}, a symbolic description of the model to be fitted.
#' @param data an data frame containing the variables in the model.
#' @param b number of resampling (default : 100)
#'
#' @details
#' Ahn et al.(2012) propose a regression approach for interval-valued data based on resampling. So this method is to resample by randomly selecting a single-valued point within each observed intervals.
#' Then, fit a classical linear regression model on each single-valued points, and calculate the average of regression coefficients over the models. The use of the resampling approach method, called Monte Carlo method (MCM), has the advantage of estimating on sample distribution approximately, and statistical inference is possible using this.
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
#' @references Ahn, J., Peng, M., Park, C., Jeon, Y.(2012), A Resampling Approach for Interval-Valued Data Regression. \emph{Statistical Analysis and Data Mining, 5}, 336-348
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
#' m1 <- imcmuni(cbind(y_L, y_U) ~ x1_L + x1_U + x2_L + x2_U, data = temp, b = 100)
#' m1
#'
#' @seealso \code{\link{RMSE}} \code{\link{symbolic.r}}
#' @import stats
imcmuni <- function(formula = formula, data = data, b = 100) {

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

  ### Resampling STEP (by using uniform distribution) ###
  # re-sampling response variable
  re_y <- matrix(NA, n, b)

  for (k in 1:b) {
    for (i in 1:n) {
      re_y[i, k] = runif(1, min = y[i, 1], max = y[i, 2])
    }
  }

  # re-sampling predictor variables
  re_xx <- array(NA, c(n, {p/2}, b))
  re_x <- array(NA, c(n, {p/2}+1, b)) # include intercept term

  for (k in 1:b) {
    for (j in 1:{p/2}) {
      for (i in 1:n) {
        re_xx[i, j, k] = runif(1, min = x[i, {2*j}], max = x[i, {2*j}+1])
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


  ### final STEP ###
  result <- list(call = match.call(),
                 response = y,
                 predictor = x,
                 resampling.coefficients = re_coef,
                 coefficients = coef,
                 fitted.values = fitted_values,
                 residuals = residuals)
  return(result)
}
