require(nprobust)
require(mvtnorm)
require(CVXR)
require(hal9001)

###########################

#' Statistical inference for partially additive models
#'
#' Performs statistical inference for the mean regression function
#' in a partially additive model via the restricted score test.
#' Also computes confidence bands for the mean regression functions.
#' 
#' @param y Response variable, an n-dimensional vector.
#' @param x Covariate of interest, an n-dimensional vector.
#' @param w Adjustment variables, an (n * q)-dimensional matrix.  With w = NULL,
#'          no covariate adjustment is performed.
#' @param null.hypothesis The null hypothesis to be test, a univariate function
#' @param d Dimension of basis expension for mean regression function, an integer
#' @param norm Choice of norm for restricted score test.  
#'             Must be one of "sup" or "L2".
#' @param no.bs Number of bootstrap samples, an integer.
#' @param confidence.band Whether to construct a confidence band, a boolean.
#' @param alpha Significance level for confidence band (one minus the coverage rate),
#'              scalar between 0 and 1.
#' @param x0 Sequence of evaluation points for confidence band, a numerical vector.
#' @param n.lam Number of penalty parameters considered in fit of
#'              partially additive model, an integer.
#' @param lam.max Largest tuning parameter considered in fit of
#'                partially additive model, a positive scalar.
#' @param n.folds Number of folds in cross-validation for estimation of
#'                partially additive model
#' @param accept  Proportion of accepted number of Monte Carlo samples in
#'                calculation of L2 norm, a scalar between 0 and 1.
#'                Only relevant when norm = "L2"
#' @param B Number of Monte Carlo samples in calculation of L2 norm, an integer.
#' @param tol Tolerance for convergence of root-finding problem in calculation
#'            of supremum norm, a positive scalar.
#' @param rkhs.norm User-specified RKHS norm for mean regression function,
#'                  a positive scalar.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{x} Covariate of interest.
#'            \item{x0} Evaluation points for confidence band.
#'            \item{fitted} Fitted values for estimated mean regression function at x.
#'            \item{fitted.x0} Fitted values for estimated mean regression function at x0.
#'            \item{test.stat} Value of test statistic for the restricted score test.
#'            \item{bs.sample} Multiplier bootstrap sample for approximating the distribution
#'                             of the test statistic.
#'            \item{p.val} P-value for test against null hypothesis.
#'            \item{cb.lower} Lower confidence limits at evaluation points x0.
#'            \item{cb.upper} Upper confidence limits at evaluation points x0.
#'          }
#'
#' @examples 
#' # Example 1: Regression w/out Adjustments
#' theta0 <- function(x){sin(3 * pi * x/2) * (1 + 18*x^2 * (sign(x) + 1))^(-1)}
#' n <- 500
#' x <- runif(n, -1, 1)
#' y <- theta0(x) + rnorm(n, sd = 1)
#' out <-  PAM(y = y, x = x)
#' 
#' plot(x, y - mean(y))
#' lines(out$x0, out$cb.lower, col = "red")
#' lines(out$x0, out$cb.upper, col = "red")
#' lines(out$x0, out$fitted.x0, col = "green")
#' lines(out$x0, theta0(out$x0) - mean(theta0(x)), col = "blue")
#' 
#' # Example 2: Regression w/ Adjustments
#' theta0 <- function(x){sin(pi * x)}
#' mu.x <- function(w){w[1] + w[2]}
#' f0 <- function(w){-2 * w[1] - 2 * w[2]}
# 
#' n <- 500
#' w <- matrix(runif(2 * n, -1, 1), nrow = n, ncol = 2)
#' x <- apply(w, 1, mu.x) + rnorm(n)
#' x <- 2 * (x - min(x))/(max(x) - min(x)) - 1
#' f0.w <- apply(w, 1, f0)
#' y <- f0.w + theta0(x) + rnorm(n)
#' out <- PAM(y = y, x = x, w = w, d = 10)
#' 
#' plot(x, y - f0.w - mean(y - f0.w))
#' lines(out$x0, out$cb.lower, col = "red")
#' lines(out$x0, out$cb.upper, col = "red")
#' lines(out$x0, out$fitted.x0, col = "green")
#' lines(out$x0, theta0(out$x0) - mean(theta0(x)), col = "blue")
#' 
#' @export
PAM <- function(y, x, w = NULL, 
                null.hypothesis = function(x) 0 * x,
                d = 20, norm = "sup", no.bs = 5000,
                confidence.band = TRUE,  alpha = .05,
                x0 = seq(min(x), max(x), length.out = 50),
                n.lam = 50, lam.max = NULL, n.folds = 10,
                accept = .5, B = round(1000/accept), tol = 1e-04,
                rkhs.norm = NULL) {
   
  scale.y <- sd(y)
  y <- y/scale.y
  null.scale <- function(x){null.hypothesis(x)/scale.y}
  if(!is.null(rkhs.norm)) rkhs.norm <- rkhs.norm/(scale.y^2)
  
  x0 <- sort(x0)
  out.fit <- PAMFit(y = y, x = x, w = w, x0 = x0, d = d,
                     n.folds = n.folds, n.lam = n.lam, lam.max = lam.max)
  
  if(norm == "sup") {
    out.test <- PAMSupTest(x = x, y.w = out.fit$y.w, w = w,
                           H.w = out.fit$H.w, H = out.fit$H,
                           coef = out.fit$coef, pen = out.fit$pen, 
                           fitted.H.w = out.fit$fitted.H.w,
                           rkhs.norm = rkhs.norm,
                           null.hypothesis = null.scale, no.bs = no.bs)
    
    if(confidence.band) {
      out.cb <- PAMConfBand(x = x, y.w = out.fit$y.w, w = w,
                       H.w = out.fit$H.w, H = out.fit$H, coef = out.fit$coef, 
                       pen = out.fit$pen, fitted.H.w = out.fit$fitted.H.w, norm = "sup",
                       no.bs = no.bs, alpha = alpha, x0 = x0,
                       rkhs.norm = rkhs.norm,
                       accept = accept, B = B, tol = tol)
    } else {
      out.cb <- NULL
    }
  }
  if(norm == "L2") {
    out.test <- PAML2Test(x = x, y.w = out.fit$y.w, w = w,
                           H.w = out.fit$H.w, H = out.fit$H,
                           coef = out.fit$coef, pen = out.fit$pen, 
                           fitted.H.w = out.fit$fitted.H.w, null.hypothesis = null.scale,
                           B = B, rkhs.norm = rkhs.norm, no.bs = no.bs)
    if(confidence.band) {
      out.cb <- PAMConfBand(x = x, y.w = out.fit$y.w, w = w,
                       H.w = out.fit$H.w, H = out.fit$H, coef = out.fit$coef, 
                       pen = out.fit$pen, fitted.H.w = out.fit$fitted.H.w, norm = "L2",
                       no.bs = no.bs, alpha = alpha, x0 = x0,
                       rkhs.norm = rkhs.norm,
                       accept = accept, B = B, tol = tol)
    } else {
      out.cb <- NULL
    }
  }
  
  # return output
  if(confidence.band) {
    out <- list(x = x,
                x0 = x0,
                fitted = out.fit$fitted * scale.y,
                fitted.x0 = out.fit$fitted.x0 * scale.y,
                test.stat = out.test$test.stat,
                bs.sample = out.test$bs.samples,
                p.val = out.test$p.val,
                cb.lower = out.cb$cb.lower * scale.y,
                cb.upper = out.cb$cb.upper * scale.y)
  } else {
    out <- list(x = x,
                x0 = x0,
                fitted = out.fit$fitted * scale.y,
                fitted.x0 = out.fit$fitted.x0 * scale.y,
                test.stat = out.test$test.stat,
                bs.sample = out.test$bs.samples,
                p.val = out.test$p.val)
  }
  
  return(out)
}

#' Sobolev basis construction
#'
#' Evaluation of Sobolev basis functions at a fixed value.
#' 
#' @param x0 Evaluation point, a scalar value.
#' @param x Covariate of interest, an n-dimensional vector.
#' @param d Dimension of basis expension for mean regression function, an integer.
#'                  
#' @return A vector containing the evaluation of each basis function at x0. 
#'         
#' @export
SobBasis <- function(x0, x, d = 20) {
  
  n <- length(x)
  x.scale <- x - mean(x)
  x0.scale <- x0 - mean(x)
  
  x0.scale <- (x0.scale - ((n+1)/n) * min(x.scale))/
              (((n+1)/n) * (max(x.scale) - min(x.scale)))
  x.scale <- (x.scale - ((n+1)/n) * min(x.scale))/
             (((n+1)/n) * (max(x.scale) - min(x.scale)))
  
  phi <- numeric(d+1)
  phi[1] <- x0.scale
  phi[1 + seq(1, d, 2)] <- 2^(1/2) * cos(2 * (1:(d/2)) * pi * (x0.scale))
  phi[2 + seq(1, d, 2)] <- 2^(1/2) * sin(2 * (1:(d/2)) * pi * (x0.scale))
  
  H <- matrix(0, nrow = n, ncol = d + 1) #eigenbasis for H
  H[,1] <- x.scale
  for(i in 1:(d/2)) {
    H[,2*i] <- 2^(1/2) * cos(2 * i * pi * (x.scale))
    H[,2*i+1] <- 2^(1/2) * sin(2 * i * pi * (x.scale))
  }
  
  return(phi - apply(H, 2, mean))
}

#' Estimates mean regression function in partially additive model.
#'
#' @param y Response variable, an n-dimensional vector.
#' @param x Covariate of interest, an n-dimensional vector.
#' @param w Adjustment variables, an (n * q)-dimensional matrix.  With w = NULL,
#'          no covariate adjustment is performed.
#' @param d Dimension of basis expension for mean regression function, an integer.
#' @param x0 Sequence of evaluation points at which regression function is evaluated.
#' @param n.folds Number of folds in cross-validation for estimation of
#'                partially additive model 
#' @param n.lam Number of penalty parameters considered in fit of
#'              partially additive model, an integer.
#' @param lam.max Largest tuning parameter considered in fit of
#'                partially additive model, a positive scalar.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{coef} Vector of coefficients for regression function.
#'            \item{rkhs.norm} RKHS norm of estimate selected via cross-validation.
#'            \item{pen} Matrix with diagonals equal to inverse of
#'                       eigenvalues for the Soblev basis.
#'            \item{y.w} Difference between observed y and conditional mean given w
#'            \item{H} An (n * (d+1))-dimensional matrix, where the (i,j) element
#'                     is the evaluation of the j-th basis function for observation i.
#'            \item{H.w} Difference between H and conditional mean given w.
#'            \item{fitted.H.w} Product of H.w and coef.
#'            \item{x0} Points at which mean regression function is evaluated.
#'            \item{fitted.x0} Fitted values for estimated mean regression function at x0.
#'            \item{fitted} Fitted values for estimated mean regression function at x.
#'          }
#'
#' @export
PAMFit <- function(y, x, w = NULL, d = 20,
                   x0 = seq(min(x), max(x), length.out = 50),
                   n.folds = 10, n.lam = 50, lam.max = NULL) {
  
  n <- length(x)
  
  # scale x to be [0,1] variable
  x0.scale <- x0 - mean(x)
  x.scale <- x - mean(x)
  x0.scale <- (x0.scale - (n+1)/n * min(x.scale))/((n+1)/n * (max(x.scale) - min(x.scale)))
  x.scale <- (x.scale - (n+1)/n * min(x.scale))/((n+1)/n * (max(x.scale) - min(x.scale)))
  
  # construct eigenbasis for class of directions
  H <- matrix(0, nrow = n, ncol = d + 1) #eigenbasis for H
  H.eval <- matrix(0, nrow = length(x0), ncol = d + 1) #eigenbasis at evaluation points
  eig <- numeric(d) # eigenvalues
  H[,1] <- x.scale
  H.eval[,1] <- x0.scale
  for(i in 1:(d/2)) {
    H[,2*i] <- 2^(1/2) * cos(2 * i * pi * (x.scale))
    H[,2*i+1] <- 2^(1/2) * sin(2 * i * pi * (x.scale))
    H.eval[,2*i] <- 2^(1/2) * cos(2 * i * pi * (x0.scale))
    H.eval[,2*i+1] <- 2^(1/2) * sin(2 * i * pi * (x0.scale))
    eig[2*i-1] <- (pi * i)^(-4)
    eig[2*i] <- (pi * i)^(-4)
  }
  
  # center eigenbasis by conditional mean given w
  H.w <- matrix(0, nrow = n, ncol = d + 1)
  if(is.null(w)) {
    for(i in 1:(d+1)) {
      H.eval[,i] <- H.eval[,i] - mean(H[,i]) 
      H.w[,i] <- H[,i] - mean(H[,i])
      H[,i] <- H.w[,i]
    }
  } else {
    for(i in 1:(d+1)){
      H.eval[,i] <- H.eval[,i] - mean(H[,i]) 
      hal.fit <- fit_hal(X = w, Y = H[,i], yolo = FALSE)
      H.w[,i] <- H[,i] - predict(hal.fit, new_data = w)
      H[,i] <- H[,i] - mean(H[,i])
    }
  }
  
  # center y by conditional mean given w
  if(is.null(w)) {
    y.w <- y - mean(y)
  } else {
    hal.fit <- fit_hal(X = w, Y = y, yolo = FALSE)
    y.w <- y - predict(hal.fit, new_data = w)
  }
  
  # select tuning parameter via cross-validation
  folds <- sample(cut(1:n, breaks = n.folds, labels = FALSE))
  
  gram <- t(H.w) %*% H.w # gram matrix
  pen <- diag(c(0, 1/eig)) # penalty matrix
  
  if(is.null(lam.max)) lam.max <- sum(diag(gram))/(n/eig[1])
  lam.min <- .0005 * lam.max
  lam.seq <- exp(seq(log(lam.min), log(lam.max), length.out = n.lam)) # sequence of tuning parameters
  
  pred.error <- matrix(NA, nrow = n.lam, ncol = n.folds) # prediction error
  for(i in 1:n.folds) {
    for(j in 1:n.lam) {
      y.w.train <- y.w[folds != i]
      y.w.test <- y.w[folds == i]
      H.w.train <- H.w[folds != i,]
      H.w.test <- H.w[folds == i,]
      
      coef <- solve(t(H.w.train) %*% H.w.train + lam.seq[j] * pen) %*% 
              t(H.w.train) %*% y.w.train
      fit <- H.w.test %*% coef
      pred.error[j,i] <- mean((y.w.test - fit)^2)
    }
  }
  
  lam.opt <- lam.seq[which.min(apply(pred.error, 1, mean))] # optimal tuning parm
  
  mean.pred.error <- apply(pred.error, 1, mean)
  opt.pred.error <- min(mean.pred.error)
  se.pred.error <- sd(pred.error[which.min(mean.pred.error),])/sqrt(n.folds)
  lam.1se <- max(lam.seq[mean.pred.error < (opt.pred.error + se.pred.error)])
  lam.cv <- lam.seq[which.min(mean.pred.error)] 
  
  lam.opt <- lam.cv
  
  # obtain fit on full data
  coef <- solve(gram + lam.opt * pen) %*% t(H.w) %*% y.w
  fitted.H.w <- H.w %*% coef
  fitted <- H %*% coef
  fitted.x0 <- H.eval %*% coef
  rkhs.norm <- t(coef) %*% pen %*% coef
  
  # return output
  out <- list(coef = coef, rkhs.norm = rkhs.norm, pen = pen,
              fitted.H.w = fitted.H.w, y.w = y.w, H = H, H.w = H.w,
              x0 = x0, fitted.x0 = fitted.x0, fitted = fitted)
  return(out)
}

#' Restricted score test performed with the supremum norm.
#'
#' @param y.w Difference between response y and conditional mean given w
#' @param x Covariate of interest, an n-dimensional vector.
#' @param w Adjustment variables, an (n * q)-dimensional matrix.
#' @param H An (n * (d+1))-dimensional matrix, where the (i,j) element
#'          is the evaluation of the j-th basis function for observation i.
#' @param H.w Difference between H and conditional mean given w.
#' @param coef Vector of coefficients for regression function.
#' @param pen Matrix with diagonals equal to inverse of eigenvalues 
#'            for the Soblev basis.
#' @param fitted.H.w Product of H.w and coef.
#' @param null.hypothesis The null hypothesis to be test, a univariate function
#' @param no.bs Number of bootstrap samples, an integer.
#' @param tol Tolerance for convergence of root-finding problem in calculation
#'            of supremum norm, a positive scalar.
#' @param rkhs.norm User-specified RKHS norm for mean regression function,
#'                  a positive scalar.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat} Value of test statistic for the restricted score test.
#'            \item{bs.sample} Multiplier bootstrap sample for approximating the distribution
#'                             of the test statistic.
#'            \item{p.val} P-value for test against null hypothesis.
#'            \item{Pi} Matrix used for calculating the test statistic.
#'          }
#'
#' @export
PAMSupTest <- function(x, y.w, w = NULL, H, H.w, coef, pen, fitted.H.w,
                       null.hypothesis = function(x) {0 * x}, no.bs = 5000,
                       tol = 1e-04, rkhs.norm = NULL) {
  n <- length(y.w)
  
  # estimate conditional mean of theta given w
  theta.x <- null.hypothesis(x)
  theta.x <- theta.x - mean(theta.x)
  if(!is.null(w) & (var(theta.x) > 0)) {
    hal.fit <- fit_hal(X = w, Y = theta.x, yolo = FALSE)
    theta.xw <- theta.x - predict(hal.fit, new_data = w)
  } else {
    theta.xw <- theta.x - mean(theta.x)
  }
  
  theta.coef <- solve(t(H) %*% H) %*% t(H) %*% theta.x
  
  # calculate smoothness parm
  V <- t(H.w) %*% diag((y.w - c(fitted.H.w))^2) %*% H.w/n
  if(is.null(rkhs.norm)) {
    rkhs.norm <- c(t(coef - theta.coef) %*% pen %*% (coef - theta.coef))
    g <- c(t(coef - theta.coef) %*% pen %*% (coef - theta.coef))/
         c(t(coef - theta.coef) %*% V %*% (coef - theta.coef))
  } else {
    g <- c(rkhs.norm)/c(t(coef - theta.coef) %*% V %*% (coef - theta.coef))
  }
  
  # lam <- 1/g
  lam.old <- 1/g
  lam.converge <- FALSE
  step <- 1e-06
  max.iter <- 1000
  iter <- 0

  while(!lam.converge) {
    iter <- iter + 1
    a <- solve(V + lam.old * pen) %*% t(H.w) %*% (y.w - theta.xw)/n
    a.1 <- solve(V + (lam.old + step) * pen) %*% t(H.w) %*% (y.w - theta.xw)/n
    a.2 <- solve(V + (lam.old - step) * pen) %*% t(H.w) %*% (y.w - theta.xw)/n

    deriv.approx <-  (((t(a.1) %*% pen %*% a.1)/(t(a.1) %*% V %*% a.1)) -
                      ((t(a.2) %*% pen %*% a.2)/(t(a.2) %*% V %*% a.2)))/
                     (2 * step)
    lam <- lam.old - c((t(a) %*% pen %*% a)/(t(a) %*% V %*% a) - g)/c(deriv.approx)
    if(lam <= 0) {lam <- lam.old * .5}

    lam.converge <- abs(log(lam) - log(lam.old)) < tol
    lam.old <- lam
    if(iter > 1000) break
  }
  
  # compute test statistic
  Pi <- H.w %*% solve(V + lam * pen) %*% t(H.w)/n
  test.stat <- c(t(y.w - theta.xw) %*% Pi %*% (y.w - theta.xw))
  
  # calculate bootstraped p-values
  bs.samples <- numeric(no.bs)
  
  for(i in 1:no.bs) {
    xi <- rnorm(n)
    mult.bs <- (xi - mean(xi)) * (y.w - fitted.H.w)
    bs.samples[i] <- t(mult.bs) %*% Pi %*% mult.bs
  }
  
  # calculate p-value
  p.val <- mean(c(test.stat) < bs.samples)
  
  # return output
  out <- list(test.stat = test.stat,
              bs.samples = bs.samples,
              Pi = Pi,
              p.val = p.val)
  
  return(out)
}

#' Restricted score test performed with the L2 norm.
#'
#' @param y.w Difference between response y and conditional mean given w
#' @param x Covariate of interest, an n-dimensional vector.
#' @param w Adjustment variables, an (n * q)-dimensional matrix.
#' @param H An (n * (d+1))-dimensional matrix, where the (i,j) element
#'          is the evaluation of the j-th basis function for observation i.
#' @param H.w Difference between H and conditional mean given w.
#' @param coef Vector of coefficients for regression function.
#' @param pen Matrix with diagonals equal to inverse of eigenvalues 
#'            for the Soblev basis.
#' @param fitted.H.w Product of H.w and coef.
#' @param null.hypothesis The null hypothesis to be test, a univariate function
#' @param no.bs Number of bootstrap samples, an integer.
#' @param accept  Proportion of accepted number of Monte Carlo samples in
#'                calculation of L2 norm, a scalar between 0 and 1.
#' @param B Number of Monte Carlo samples in calculation of L2 norm, an integer.
#' @param rkhs.norm User-specified RKHS norm for mean regression function,
#'                  a positive scalar.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat} Value of test statistic for the restricted score test.
#'            \item{bs.sample} Multiplier bootstrap sample for approximating the distribution
#'                             of the test statistic.
#'            \item{p.val} P-value for test against null hypothesis.
#'            \item{Pi} Matrix used for calculating the test statistic.
#'          }
#'
#' @export
PAML2Test <- function(x, y.w, w = NULL, H.w, H, coef, pen, fitted.H.w,
                       null.hypothesis = function(x) {0 * x},
                       no.bs = 5000, accept = .5, B = round(1000/accept),
                       rkhs.norm = NULL) {
  n <- length(y.w)
  
  # estimate conditional mean of theta given w
  theta.x <- null.hypothesis(x)
  theta.x <- theta.x - mean(theta.x)
  if(!is.null(w) & (var(theta.x) > 0)) {
    hal.fit <- fit_hal(X = w, Y = theta.x, yolo = FALSE)
    theta.xw <- theta.x - predict(hal.fit, new_data = w)
  } else {
    theta.xw <- theta.x - mean(theta.x)
  }
  
  theta.coef <- solve(t(H) %*% H) %*% t(H) %*% theta.x
  
  V <- t(H.w) %*% diag((y.w - c(fitted.H.w))^2) %*% H.w/n
  
  # calculate smoothness parm
  V <- t(H.w) %*% diag((y.w - c(fitted.H.w))^2) %*% H.w/n
  if(is.null(rkhs.norm)) {
    rkhs.norm <- c(t(coef - theta.coef) %*% pen %*% (coef - theta.coef))
    g <- c(t(coef - theta.coef) %*% pen %*% (coef - theta.coef))/
         c(t(coef - theta.coef) %*% V %*% (coef - theta.coef))
  } else {
    g <- c(rkhs.norm)/c(t(coef - theta.coef) %*% V %*% (coef - theta.coef))
  }
  
  # generate Monte Carlo sample
  d <- ncol(H.w) - 1
  lam <- 1/g
  good.mc <- FALSE
  density.mc <- numeric(B)
  var.mc <- numeric(B)
  
  while(!good.mc) {
    U <- rmvnorm(B, mean = rep(0, d+1), sigma = solve(V + lam * pen))
    smooth.mc <- numeric(B)
    for(b in 1:B){
      U[b,] <- U[b,]/c((t(U[b,]) %*%  V %*% U[b,])^.5)
      smooth.mc[b] <- (t(U[b,]) %*% pen %*% U[b,])/(t(U[b,]) %*% V %*% U[b,])
      var.mc[b] <- t(U[b,]) %*% V %*% U[b,]
    }
    if(mean(smooth.mc < g) >= accept) {
        good.mc <- TRUE
     } else {
        lam <- lam * 1.1
     }
  }
  
  U <- U[smooth.mc <= g,]
  var.mc <- var.mc[smooth.mc <= g]
  smooth.mc <- smooth.mc[smooth.mc <= g]
  
  # approximate density of smooth.mc
  density.f <- approxfun(density(smooth.mc))
  density.mc <- density.f(smooth.mc)
  density.mc <- density.mc/max(density.mc)

  # compute test statistic
  # Pi.1 <- t(U) %*% diag(1/(density.mc * var.mc)) %*% U/B
  Pi <- H.w %*% t(U) %*% diag(1/(density.mc * var.mc)) %*% U %*% t(H.w)/B
  Pi <- Pi/sum(diag(Pi))
  # Pi <- Pi/max(eigen(Pi)$values)
  test.stat <- c(t(y.w - theta.xw) %*% Pi %*% (y.w - theta.xw))
   
  # calculate bootstrap samples
  bs.samples <- numeric(no.bs)
  
  for(i in 1:no.bs) {
    xi <- rnorm(n)
    mult.bs <- (xi - mean(xi)) * (y.w - fitted.H.w)
    bs.samples[i] <- t(mult.bs) %*% Pi %*% mult.bs
  }
  
  # calculate p-value
  p.val <- mean(test.stat < bs.samples)
  
  # return output
  out <- list(test.stat = test.stat,
              bs.samples = bs.samples,
              p.val = p.val,
              Pi = Pi)
  return(out)
}

#' Confidence band construction for the mean regression function in
#' a partially additive model.
#'
#' @param y.w Difference between response y and conditional mean given w
#' @param x Covariate of interest, an n-dimensional vector.
#' @param w Adjustment variables, an (n * q)-dimensional matrix.
#' @param H An (n * (d+1))-dimensional matrix, where the (i,j) element
#'          is the evaluation of the j-th basis function for observation i.
#' @param H.w Difference between H and conditional mean given w.
#' @param coef Vector of coefficients for regression function.
#' @param pen Matrix with diagonals equal to inverse of eigenvalues 
#'            for the Soblev basis.
#' @param fitted.H.w Product of H.w and coef.
#' @param no.bs Number of bootstrap samples, an integer.
#' @param alpha Significance level for confidence band (one minus the coverage rate),
#'              scalar between 0 and 1.
#' @param x0 Sequence of evaluation points for confidence band, a numerical vector.
#' @param accept  Proportion of accepted number of Monte Carlo samples in
#'                calculation of L2 norm, a scalar between 0 and 1.
#'                Only relevant when norm = "L2".
#' @param B Number of Monte Carlo samples in calculation of L2 norm, an integer.
#' @param rkhs.norm User-specified RKHS norm for mean regression function,
#'                  a positive scalar.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{test.stat} Value of test statistic for the restricted score test.
#'            \item{bs.sample} Multiplier bootstrap sample for approximating the distribution
#'                             of the test statistic.
#'            \item{p.val} P-value for test against null hypothesis.
#'            \item{Pi} Matrix used for calculating the test statistic.
#'          }
#'
#' @export
PAMConfBand <- function(x, y.w, w = NULL, H.w, H, coef, pen, fitted.H.w, norm = "sup",
                   no.bs = 5000, alpha = .05,
                   x0 = seq(min(x), max(x), length.out = 50),
                   rkhs.norm = NULL,
                   accept = .5, B = 10000, tol = 1e-04) {
  
  d <- ncol(H.w) - 1
  n <- length(y.w)
  
  cb.lower <- numeric(length(x0))
  cb.upper <- numeric(length(x0))
  
  if(is.null(rkhs.norm)) rkhs.norm <- c(t(coef) %*% pen %*% coef)
  
  # perform restricted score test at theta = 0
  if(norm == "sup") {
    test.zero <- PAMSupTest(x = x, y.w = y.w, w = w, H.w = H.w, H = H,
                            coef = coef, pen = pen, fitted.H.w = fitted.H.w,
                            null.hypothesis = function(x) {0 * x}, no.bs = no.bs,
                            tol = tol, rkhs.norm = rkhs.norm)
  }
  if(norm == "L2") {
    test.zero <- PAML2Test(x = x, y.w = y.w, w = w, H.w = H.w, H = H,
                           coef = coef, pen = pen, fitted.H.w = fitted.H.w,
                           null.hypothesis  = function(x) {0 * x}, no.bs = no.bs,
                           B = B, accept = accept, rkhs.norm = rkhs.norm)
  }
  
  Pi <- test.zero$Pi
  bs.samples <- test.zero$bs.samples
  t.star <- quantile(bs.samples, 1-alpha)
  
  # scale Pi to make CVXR less finnicky
  Pi <- Pi * rkhs.norm/t.star
  t.star <- rkhs.norm
  
  # calculate bands
  A <- t(H.w) %*% Pi %*% H.w # constraints for opt problems
  B <- 2 * t(y.w) %*% Pi %*% H.w
  C <- t(y.w) %*% Pi %*% y.w
  for(i in 1:length(x0)) {
    phi <- SobBasis(x0 = x0[i], x = x, d = ncol(H.w) - 1)
    
    # upper band
    a <- Variable(d + 1)
    objective <- Minimize(-t(phi) %*% a)
    prob <- Problem(objective,
                    list(quad_form(a, pen) - rkhs.norm <= 0,
                         quad_form(a, A) - B %*% a + C - t.star <= 0))
    CVXR.result <- solve(prob)
    cb.upper[i] <- ifelse(CVXR.result$status == "solver_error",
                          NA, -CVXR.result$value)
    
    # lower band
    a <- Variable(d + 1)
    objective <- Minimize(t(phi) %*% a)
    prob <- Problem(objective,
                    list(quad_form(a, pen) - rkhs.norm <= 0,
                         quad_form(a, A) - B %*% a + C - t.star <= 0))
    CVXR.result <- solve(prob)
    cb.lower[i] <- ifelse(CVXR.result$status == "solver_error",
                          NA, CVXR.result$value)
  }
  
  # return output
  out <- list(x0 = x0,
              cb.lower = cb.lower,
              cb.upper = cb.upper)
  return(out)
  
}

lprobust.uniform <- function(x, y, alpha = .05, no.bs = 5000,
                             null.function = function(x){0 * x},
                             x0 = seq(min(x), max(x), length.out = 50)) {
  
  # obtain de-biased kernel regression estimate
  fit <- lprobust(y - mean(y), x, eval = x0, covgrid = TRUE)
  cov.fit <- diag(1/fit$Estimate[,8]) %*% fit$cov.rb %*% diag(1/fit$Estimate[,8])
  
  # bootstrap sample max of gaussians with covariance sigma
  bs.samples <- replicate(no.bs,
                          max(abs(rmvnorm(1, mean = rep(0, length(x0)), sigma = cov.fit))))
  
  # calculate bands
  cb.upper <- fit$Estimate[,6] + fit$Estimate[,8] * quantile(bs.samples, 1-alpha)
  cb.lower <- fit$Estimate[,6] - fit$Estimate[,8] * quantile(bs.samples, 1-alpha)
  
  # test null hypothesis
  test.stat <- max(abs(fit$Estimate[,6] - null.function(x0))/fit$Estimate[,8])
  p.val <- mean(test.stat < bs.samples)
  
  # return output
  out <- list(x0 = x0,
              fitted = fit$Estimate[,6],
              cb.lower = cb.lower,
              cb.upper = cb.upper,
              test.stat = test.stat,
              bs.samples = bs.samples,
              p.val = p.val)
  return(out)
}

# Example 1: Regression w/out Adjustments
# set.seed(206)
# theta0 <- function(x){sin(3 * pi * x/2) * (1 + 18*x^2 * (sign(x) + 1))^(-1)}
# n <- 500
# x <- runif(n, -1, 1)
# y <- theta0(x) + rnorm(n, sd = 1)
# out <-  PAM(y = y, x = x, norm = "sup", d = 20)
#
# plot(x, y - mean(y))
# lines(out$x0, out$cb.lower, col = "red")
# lines(out$x0, out$cb.upper, col = "red")
# lines(out$x0, out$fitted.x0, col = "green")
# lines(out$x0, theta0(out$x0) - mean(theta0(x)), col = "blue")

# Example 2: Regression w/ Adjustments
# theta0 <- function(x){sin(pi * x)}
# mu.x <- function(w){w[1] + w[2]}
# f0 <- function(w){-2 * w[1] - 2 * w[2]}
# 
# n <- 500
# w <- matrix(runif(2 * n, -1, 1), nrow = n, ncol = 2)
# x <- apply(w, 1, mu.x) + rnorm(n)
# x <- 2 * (x - min(x))/(max(x) - min(x)) - 1
# f0.w <- apply(w, 1, f0)
# y <- f0.w + theta0(x) + rnorm(n)
# out <- PAM(y = y, x = x, w = w, d = 10, null.hypothesis = theta0)
# 
# plot(x, y - f0.w - mean(y - f0.w))
# lines(out$x0, out$cb.lower, col = "red")
# lines(out$x0, out$cb.upper, col = "red")
# lines(out$x0, out$fitted.x0, col = "green")
# lines(out$x0, theta0(out$x0) - mean(theta0(x)), col = "blue")
