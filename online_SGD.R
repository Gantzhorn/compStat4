decay_scheduler <- function(gamma1 = 1, a = 1, K = 1, gamma2, n1) {
  force(a)
  if (!missing(gamma2) && !missing(n1))
    K <- n1^a * gamma2 / (gamma0 - gamma2)
  b <- gamma1 * K
  function(n) b / (K + n^a)
}

rate <- decay_scheduler(gamma1 = 0.5, K = 100)

# Online learning

initpar <- c(alpha0 + runif(1,-1,1), beta0 + runif(1,-1,1), gamma0 + runif(1,-1,1), rho0 + runif(1,-1,1))

online_SG <- function(
    par0,
    gamma,
    maxiter = 10000,
    cb = NULL,
    ...
){
  gamma <- if(is.function(gamma))gamma(1:maxiter) else rep(gamma, (maxiter))
  par <- par0
  for (i in 1:maxiter){
    if(!is.null(cb)) cb()
    # Simulate online data 
    X <- rnorm(1, mean = 0, sd = omega)
    
    epsilon <- rnorm(1, mean = 0, sd = sigma) 
    
    Y <- densY(X, alpha0, beta0, gamma0, rho0) + epsilon
    
    par <- par - gamma[i] * grad_calc(X, Y, par[1],par[2], par[3], par[4])
  }
  par
}


online_loglogis_SG_tracer <- tracer("par", N = 0)
res <- online_SG(initpar, rate, cb = online_loglogis_SG_tracer$tracer)

summary(online_loglogis_SG_tracer) %>% as_tibble() %>% rename("alpha" = par.1, "beta" = par.2, "gamma" = par.3, "rho" = par.4, "time" = .time) %>% 
  pivot_longer(cols = -time, names_to = "Parameter", values_to = "Estimate") %>% 
  ggplot(aes(x = time, y = Estimate, col = Parameter)) + geom_line(size = 1) +
  geom_hline(data = trueval_tibble, aes(yintercept = Value, col = Parameter), linetype = "dashed", size = 1)
  

