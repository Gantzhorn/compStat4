library(tidyverse)
ggplot2::theme_set(theme_bw())
library(microbenchmark)
library(profvis)
library(Rcpp)
library(splines)
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
beta0 <- 0.5
gamma0 <- -0.5
rho0 <- -1


trueval_tibble <- tibble(Parameter = c("alpha", "beta", "gamma", "rho"), Value = c(alpha0, beta0, gamma0, rho0))

# Sim parameters
sigma <- 0.1

omega <- 5

# Grad-test
test_points <- c(1,1)


obj_func <- function(x, y, param){
  1/length(x)*(sum((y-densY(x, param[1], param[2], param[3], param[4]))^2))
}
obj_func_test <- function(param){
  1/length(test_points[1])*(sum((test_points[2]-densY(test_points[1], param[1], param[2], param[3], param[4]))^2))
}



grad_calc(test_points[1], test_points[2], alpha0, beta0, gamma0, rho0)

numDeriv::grad(obj_func_test, x = c(alpha0, beta0, gamma0, rho0))