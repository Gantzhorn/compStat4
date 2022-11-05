SG_mini <- function(
    par, 
    N,                 
    gamma,             
    epoch = batch,     
    ...,               
    maxiter = 100,    
    sampler = sample,  
    cb = NULL
) {
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter) 
  for(k in 1:maxiter) {
    if(!is.null(cb)) cb()
    samp <- sampler(N)
    par <- epoch(par, samp, gamma[k], ...)
  }
  par
}

batch <- function(
    par, 
    samp,
    gamma,            
    m = 50,            # Mini-batch size
    ...
) {
  M <- floor(length(samp) / m)
  for(j in 0:(M - 1)) {
    i <- samp[(j * m + 1):(j * m + m)]
    par <- par - gamma * grad_calc(X[i], Y[i], par[1], par[2], par[3], par[4])
  }
  par
}

SG_tracer <- tracer("par", N = 0)
SG_mini(initpar,
        N = length(X),
        gamma = 5e-2,
        maxiter = 200,
        cb = SG_tracer$tracer
        )

summary(SG_tracer) %>% as_tibble() %>% rename("alpha" = par.1, "beta" = par.2, "gamma" = par.3, "rho" = par.4, "time" = .time) %>% 
  pivot_longer(cols = -time, names_to = "Parameter", values_to = "Estimate") %>% 
  ggplot(aes(x = time, y = Estimate, col = Parameter)) + geom_line(size = 1) +
  geom_hline(data = trueval_tibble, aes(yintercept = Value, col = Parameter), linetype = "dashed", size = 1)
