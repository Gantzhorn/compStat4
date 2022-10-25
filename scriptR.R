library(tidyverse)
ggplot2::theme_set(theme_bw())
library(microbenchmark)
library(profvis)
library(Rcpp)
library(splines)
library(xtable)
library(Matrix)

## Key-functions ##

#Density
densY <- function(x, alpha, beta, gamma, rho){
  gamma + (rho-gamma)/(1+exp(beta*log(x)-alpha)) #Note that x is log-normal
}

#Gradient
grad <- function(x, alpha, beta, gamma, rho){
  N <- length(x)
  grad_alpha <- 2/N*sum(Y-gamma+(gamma-rho)/(1+exp(beta*log(x)-alpha))*
                          ((gamma-rho)*exp(beta*log(x)-alpha))/(exp(beta*log(x)-alpha)+1)^2)
  
  grad_beta <- -2/N*sum(Y-gamma+(gamma-rho)/(1+exp(beta*log(x)-alpha))*
                         ((gamma-rho)*log(x)*exp(beta*log(x)-alpha))/(exp(beta*log(x)-alpha)+1)^2)
  
  grad_gamma <- 2/N*sum(Y-gamma+(gamma-rho)/(1+exp(beta*log(x)-alpha))*
                           (1/(1+exp(beta*log(x)-alpha))-1))
  
  grad_rho <- -2/N*sum(Y-gamma+(gamma-rho)/(1+exp(beta*log(x)-alpha))*
                         1/(1+exp(beta*log(x)-alpha)))
  c(grad_alpha, grad_beta, grad_gamma, grad_rho)
}


## Simuliation ##

#Number of observations
N <- 5000

sigma <- 0.2
epsilon <- rnorm(N, mean = 0, sd = sigma)

omega <- 0.5

X <- rlnorm(N, mean = 0, sd = omega)

#True parameters
alpha0 <- 0.5
beta0 <- 0.3
gamma0 <- 0.1
rho0 <- 1.5

# Y-values
Y <- densY(X, alpha0, beta0, gamma0, rho0) + epsilon

grad(X, alpha0, beta0, gamma0, rho0)

#Horsedata
horse_data <- readr::read_csv("4_Horses.csv", col_types = cols(dead = col_logical()))
horse_missing <- horse_data %>%
  mutate(missing = ifelse(is.na(Temperature) == TRUE, TRUE,FALSE)) %>%
  group_by(dead,missing) %>%
  summarise(n = n())



