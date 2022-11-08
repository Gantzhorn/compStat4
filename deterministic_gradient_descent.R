gradient_descent <- function(par,
                             gamma = 0.5,
                             maxiter = 400,
                             epsilon = 1e-06,
                             cb = NULL){
  
  gamma <- if(is.function(gamma))gamma(1:maxiter) else rep(gamma, (maxiter))
  
  for (i in 1:maxiter){
    par0 <- par
    par <- par0 - gamma[i]*grad_calc(X, Y, par[1], par[2], par[3], par[4])
    if(!is.null(cb)) cb()
    if(sum((par0-par)^2) <= epsilon) break
  }
  if(i == maxiter)(base::warning("Maximum number of iterations reached"))
  else{print(glue::glue("Convergence reached after ", i, " steps."))}
  par
}

deter_tracer <- CSwR::tracer(c("par"))

gradient_descent(par = c(4.5,3,0.5,2), gamma = rate, cb = deter_tracer$tracer)

summary(deter_tracer) %>% as_tibble() %>% rename("alpha" = par.1, "beta" = par.2, "gamma" = par.3, "rho" = par.4, "time" = .time) %>% 
  pivot_longer(cols = -time, names_to = "Parameter", values_to = "Estimate") %>% 
  ggplot(aes(x = time, y = Estimate, col = Parameter)) + geom_line(size = 1) +
  geom_hline(data = trueval_tibble, aes(yintercept = Value, col = Parameter), linetype = "dashed", size = 1)
