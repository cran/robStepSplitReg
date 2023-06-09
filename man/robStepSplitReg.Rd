% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robStepSplitReg.R
\name{robStepSplitReg}
\alias{robStepSplitReg}
\title{Robust Stepwise Split Regularized Regression}
\usage{
robStepSplitReg(
  x,
  y,
  n_models = 1,
  model_saturation = c("fixed", "p-value")[1],
  alpha = 0.05,
  model_size = NULL,
  robust = TRUE,
  compute_coef = FALSE,
  en_alpha = 1/4
)
}
\arguments{
\item{x}{Design matrix.}

\item{y}{Response vector.}

\item{n_models}{Number of models into which the variables are split.}

\item{model_saturation}{Criterion to determine if a model is saturated. Must be one of "fixed" (default) or "p-value".}

\item{alpha}{P-value used to determine when the model is saturated}

\item{model_size}{Size of the models in the ensemble.}

\item{robust}{Argument to determine if robust measures of location, scale and correlation are used. Default is TRUE.}

\item{compute_coef}{Argument to determine if coefficients are computed (via adaptive PENSE) for each model. Default is FALSE.}

\item{en_alpha}{Elastic net mixing parmeter for parameters shrinkage. Default is 1/4.}
}
\value{
An object of class robStepSplitReg.
}
\description{
\code{robStepSplitReg} performs robust stepwise split regularized regression.
}
\examples{
# Required library
library(mvnfast)

# Simulation parameters
n <- 50
p <- 500
rho <- 0.5
p.active <- 100
snr <- 1
contamination.prop <- 0.2

# Setting the seed
set.seed(0)

# Simulation of beta vector
true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))

# Simulation of uncontaminated data 
sigma.mat <- matrix(0, nrow = p, ncol = p)
sigma.mat[1:p.active, 1:p.active] <- rho
diag(sigma.mat) <- 1
x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
sigma <- as.numeric(sqrt(t(true.beta) \%*\% sigma.mat \%*\% true.beta)/sqrt(snr))
y <- x \%*\% true.beta + rnorm(n, 0, sigma)

# Contamination of data 
contamination_indices <- 1:floor(n*contamination.prop)
k_lev <- 2
k_slo <- 100
x_train <- x
y_train <- y
beta_cont <- true.beta
beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
for(cont_id in contamination_indices){
  
  a <- runif(p, min = -1, max = 1)
  a <- a - as.numeric((1/p)*t(a) \%*\% rep(1, p))
  x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
    k_lev * a / as.numeric(sqrt(t(a) \%*\% solve(sigma.mat) \%*\% a))
  y_train[cont_id] <- t(x_train[cont_id,]) \%*\% beta_cont
}

# Ensemble models
ensemble_fit <- robStepSplitReg(x_train, y_train,
                                n_models = 5,
                                model_saturation = c("fixed", "p-value")[1],
                                alpha = 0.05, model_size = n - 1,
                                robust = TRUE,
                                compute_coef = TRUE,
                                en_alpha = 1/4)

# Ensemble coefficients
ensemble_coefs <- coef(ensemble_fit, group_index = 1:ensemble_fit$n_models)
sens_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/p.active
spec_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/sum(ensemble_coefs[-1]!=0)

# Simulation of test data
m <- 2e3
x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
y_test <- x_test \%*\% true.beta + rnorm(m, 0, sigma)

# Prediction of test samples
ensemble_preds <- predict(ensemble_fit, newx = x_test, 
                          group_index = 1:ensemble_fit$n_models,
                          dynamic = FALSE)
mspe_ensemble <- mean((y_test - ensemble_preds)^2)/sigma^2

}
\seealso{
\code{\link{coef.robStepSplitReg}}, \code{\link{predict.robStepSplitReg}}
}
\author{
Anthony-Alexander Christidis, \email{anthony.christidis@stat.ubc.ca}
}
