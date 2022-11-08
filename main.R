library(tidyverse)
ggplot2::theme_set(theme_bw())
library(microbenchmark)
library(profvis)
library(Rcpp)
library(RcppArmadillo)
library(xtable)
library(dqrng)
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

source("~/Desktop/Skole/Compstat/Assignment4/compStat4/online_SGD.R")

source("~/Desktop/Skole/Compstat/Assignment4/compStat4/batch_SGD.R")

source("~/Desktop/Skole/Compstat/Assignment4/compStat4/miniBatch_SGD.R")

source("~/Desktop/Skole/Compstat/Assignment4/compStat4/deterministic_gradient_descent.R")


# Comparison of SGD's vs gradient descent
par_init <- c(100, 1000, 5000, 10000, 25000, 50000, 75000, 100000)
errorpath_batch <- vector("list", length(par_init))
errorpath_mini <- vector("list", length(par_init))
errorpath_deter <- vector("list", length(par_init))
set.seed(1)

for(i in seq_along(par_init)){
  N <- par_init[i]
  
  X <- rnorm(N, mean = 0, sd = omega)
  
  epsilon <- rnorm(N, mean = 0, sd = sigma) 
  
  Y <- densY(X, alpha0, beta0, gamma0, rho0) + epsilon
  
  deter_tracer <- CSwR::tracer(c("par"), N = 0)
  mini_SG_tracer <- CSwR::tracer("par", N = 0)
  batch_SG_tracer <- CSwR::tracer("par", N = 0)
  
  SG_mini(initpar, N = N, gamma = rate, cb = mini_SG_tracer$tracer, maxiter = 200)
  errorpath_mini[[i]] <- summary(mini_SG_tracer) %>% tibble(., N = N, type = "minibatch")
  mini_SG_tracer$clear()
  
  SG(initpar, N = N, gamma = rate, cb = batch_SG_tracer$tracer, maxiter = 200)
  errorpath_batch[[i]] <- summary(batch_SG_tracer) %>% tibble(., N = N, type = "basicbatch")
  batch_SG_tracer$clear()
  
  gradient_descent(initpar, gamma = rate, cb = deter_tracer$tracer, maxiter = 200)
  errorpath_deter[[i]] <- summary(deter_tracer) %>% tibble(., N = N, type = "Deterministic")
  deter_tracer$clear()
}

bind_rows(errorpath_batch,
          errorpath_deter,
          errorpath_mini) %>% 
  mutate(squarederror = (par.1-alpha0)^2+(par.2-beta0)^2+(par.3-gamma0)^2+(par.4-rho0)^2,
         time = .time) %>% 
  ggplot(aes(x = log(time), y = log(squarederror), col = type)) + geom_line(size = .9) +
  facet_wrap(~N,ncol = 4, labeller = ) + theme(axis.title = element_text(face = "bold", size = 14),
                                  axis.text = element_text(face = "bold", size = 12),
                                  legend.title = element_text(face = "bold", size = 16),
                                  legend.text = element_text(face = "bold", size =12),
                                  strip.text = element_text(face = "bold", size = 12))

# Profiling
profvis::profvis(source("~/Desktop/Skole/Compstat/Assignment4/compStat4/batch_SGD.R"))

