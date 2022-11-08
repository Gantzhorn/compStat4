library(tidyverse)
ggplot2::theme_set(theme_bw())
library(microbenchmark)
library(profvis)
library(Rcpp)
library(xtable)
library(numDeriv)
library(CSwR)

## Key-functions ##

#Density
densY <- function(x, alpha, beta, gamma, rho){
  gamma + (rho-gamma)/(1+exp(beta*x-alpha)) 
}

# Nonlinear least squares

#Gradient
grad_calc <- function(x, y, alpha, beta, gamma, rho){
  N <- length(x)
  grad_alpha <- 2/N*sum(
    (y-gamma+(gamma-rho)/(1+exp(beta*x-alpha)))*
      ((gamma-rho)*exp(beta*x-alpha))/((exp(beta*x-alpha)+1)^2)
  )
  
  grad_beta <- -2/N*sum(
    (y-gamma+(gamma-rho)/(1+exp(beta*x-alpha)))*
      ((gamma-rho)*x*exp(beta*x-alpha))/((exp(beta*x-alpha)+1)^2)
  )
  
  grad_gamma <- 2/N*sum(
    (y-gamma+(gamma-rho)/(1+exp(beta*x-alpha)))*
      (1/(1+exp(beta*x-alpha))-1)
  )
  
  grad_rho <- -2/N*sum(
    (y-gamma+(gamma-rho)/(1+exp(beta*x-alpha)))*
      (1/(1+exp(beta*x-alpha)))
  )
  c(grad_alpha, grad_beta, grad_gamma, grad_rho)
}

#True parameters
alpha0 <- 1
beta0 <- 2
gamma0 <- -1
rho0 <- 0.5


trueval_tibble <- tibble(Parameter = c("alpha", "beta", "gamma", "rho"), Value = c(alpha0, beta0, gamma0, rho0))

# Sim parameters
sigma <- 0.5

omega <- 5

# Grad-test
test_points <- c(1,1)


obj_func <- function(x, y, param){
  1/length(x)*(sum((y-densY(x, param[1], param[2], param[3], param[4]))^2))
}
obj_func_test <- function(param, x_inp, y_inp){
  1/length(x_inp)*(sum((y_inp-densY(x_inp, param[1], param[2], param[3], param[4]))^2))
}



grad_calc(test_points[1], test_points[2], alpha0, beta0, gamma0, rho0)

numDeriv::grad(obj_func_test, x = c(alpha0, beta0, gamma0, rho0), x_inp = test_points[1], y_inp = test_points[2])
