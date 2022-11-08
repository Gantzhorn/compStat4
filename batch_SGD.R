# Batch stochastic gradient algorithm

# Simulation

N <- 10000

X <- rnorm(N, mean = 0, sd = omega)

epsilon <- rnorm(N, mean = 0, sd = sigma) 

Y <- densY(X, alpha0, beta0, gamma0, rho0) + epsilon

SG <- function(
    par,
    N,
    gamma,
    gradient = grad_calc,
    maxiter = 500,
    sampler = sample,
    cb = NULL,
    epsilon = 1e-6,
    ...
){
  gamma <- if(is.function(gamma))gamma(1:maxiter) else rep(gamma, maxiter)
  for (k in 1:maxiter){
    if(!is.null(cb)) cb()
    samp <- sampler(N)
    for (j in 1:N){
      par0 <- par
      i <- samp[j]
      par <- par0 - gamma[k] * gradient(X[i], Y[i], par0[1], par0[2], par0[3], par0[4])
      if(sum((par-par0)^2)<epsilon) break
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




# Plot estimates against true values and initial values
# initpar_tibble <- tibble(initpar = initpar, Parameter = c("alpha", "beta", "gamma", "rho"))
# 
# 
# summary(loglogis_SG_tracer) %>% as_tibble() %>% rename("alpha" = par.1, "beta" = par.2, "gamma" = par.3, "rho" = par.4, "time" = .time) %>% 
#   pivot_longer(cols = -time, names_to = "Parameter", values_to = "Estimate") %>% 
#   ggplot(aes(x = time, y = Estimate, col = Parameter)) + geom_line(size = 1) +
#   geom_hline(data = trueval_tibble, aes(yintercept = Value, col = Parameter), linetype = "dashed", size = 1) + 
#   theme(axis.title = element_text(face = "bold", size = 14),
#         axis.text = element_text(face = "bold", size = 12),
#         legend.title = element_text(face = "bold", size = 16),
#         legend.text = element_text(face = "bold", size =12),
#         strip.text.x = element_text(face = "bold", size = 12))

# # Convergence for different starting values of alpha, beta
# loglogis_SG_tracer <- tracer("par", N = 0)
# 
# 
# par_init <- expand.grid(
#   alpha = c(-4,-2, 0, 1, 1.5, 2, 5, 8),
#   beta = c(-4, -2, -0.5, 0.5, 2, 4, 8)
# ) %>% as.matrix()
# 
# 
# set.seed(1)
# 
# paths <- vector("list", nrow(par_init))
# 
# for(i in 1:nrow(par_init)){
#   SG(c(par_init[i,1], par_init[i,2], initpar[3], initpar[4]),
#      N = N,
#      gamma = rate,
#      cb = loglogis_SG_tracer$tracer)
#   paths[[i]] <- summary(loglogis_SG_tracer) %>% tibble(., alphainit = par_init[i,1], betainit = par_init[i,2], .name_repair = "unique")
#   loglogis_SG_tracer$clear()
# }
# 
# paths_tibble <- bind_rows(paths) %>% select(-c("par....3", "par....4")) %>%
#   rename("alpha" = par.alpha, "beta" = par.beta, "time" = .time) %>% 
#   mutate(alphainit = as.factor(alphainit), betainit = as.factor(betainit)) %>%
#   pivot_longer(-c(time,alphainit, betainit), names_to = "Parameter")
# 
# 
# 
# paths_tibble %>% ggplot(aes(x = time, y = value, col = Parameter)) + geom_line(size = 1) + 
#   facet_grid(betainit~alphainit) +
#   geom_hline(data = trueval_tibble[1:2,], aes(yintercept = Value, col = Parameter), linetype = "dashed") +
#   facet_grid(betainit~alphainit) + 
#   theme(axis.title = element_text(face = "bold", size = 14),
#         axis.text = element_text(face = "bold", size = 12),
#         legend.title = element_text(face = "bold", size = 16),
#         legend.text = element_text(face = "bold", size =12),
#         strip.text = element_text(face = "bold", size = 12)) +
#   scale_y_continuous(breaks = seq(-3,6, by = 3))


