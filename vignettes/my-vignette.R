## ---- echo = FALSE-------------------------------------------------------
library(truncnorm)

## ------------------------------------------------------------------------
set.seed(2017)
x1_L = rnorm(30, 3, 0.01) - rnorm(30, 0, 0.01)
x1_U = rnorm(30, 3, 0.01) + rnorm(30, 3, 0.01)
x2_L = runif(30, 1.5, 3) - runif(30, 0, 1)
x2_U = runif(30, 1.5, 3) + runif(30, 1, 2)
y_L = x1_L + x2_L + rnorm(30, 0, 0.03)
y_U = x1_U + x2_U + rnorm(30, 0, 0.03)

temp <- as.data.frame(cbind(y_L, y_U, x1_L, x1_U, x2_L, x2_U))

## ---- echo = FALSE, results == "asis"------------------------------------
knitr::kable(temp)

## ---- echo = FALSE-------------------------------------------------------
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

## ---- echo = FALSE-------------------------------------------------------
RMSE <- function(model) {
  y = model$response
  yhat = model$fitted.values
  n = nrow(y)

  RMSE_L <- sqrt(sum((y[, 1]- yhat[, 1])^{2}) / n)
  RMSE_U <- sqrt(sum((y[, 2]- yhat[, 2])^{2}) / n)

  rmse <- cbind(RMSE_L, RMSE_U)
  return(rmse)
}

## ---- echo = FALSE-------------------------------------------------------
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

## ------------------------------------------------------------------------
m1 <- imcmtn(cbind(y_L, y_U) ~ x1_L + x1_U + x2_L + x2_U, data = temp, b = 100)
m1
RMSE(m1)
symbolic.r(m1)

