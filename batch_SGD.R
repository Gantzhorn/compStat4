# Batch stochastic gradient algorithm

# Simulation

N <- 5000

X <- rnorm(N, mean = 0, sd = omega)

epsilon <- rnorm(N, mean = 0, sd = sigma) 

Y <- densY(X, alpha0, beta0, gamma0, rho0) + epsilon

SG <- function(
    par,
    N,
    gamma,
    maxiter = 100,
    sampler = sample,
    cb = NULL,
    ...
){
  gamma <- if(is.function(gamma))gamma(1:maxiter) else rep(gamma, maxiter)
  for (k in 1:maxiter){
    if(!is.null(cb)) cb()
    samp <- sampler(N)
    for (j in 1:N){
      i <- samp[j]
      par <- par - gamma[k] * grad_calc(X[i], Y[i], par[1], par[2], par[3], par[4])
      
    }
  }
  par
}

# Tracing

loglogis_SG_tracer <- tracer("par", N = 0)

res <- SG(initpar,
          N = N,
          gamma = rate,
          cb = loglogis_SG_tracer$tracer)

# Plot estimates against true values
summary(loglogis_SG_tracer) %>% as_tibble() %>% rename("alpha" = par.1, "beta" = par.2, "gamma" = par.3, "rho" = par.4, "time" = .time) %>% 
  pivot_longer(cols = -time, names_to = "Parameter", values_to = "Estimate") %>% 
  ggplot(aes(x = time, y = Estimate, col = Parameter)) + geom_line(size = 1) +
  geom_hline(data = trueval_tibble, aes(yintercept = Value, col = Parameter), linetype = "dashed", size = 1)

# Objective function
