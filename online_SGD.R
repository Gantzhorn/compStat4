decay_scheduler <- function(gamma1 = 1, a = 1, K = 1, gamma2, n1) {
  force(a)
  if (!missing(gamma2) && !missing(n1))
    K <- n1^a * gamma2 / (gamma0 - gamma2)
  b <- gamma1 * K
  function(n) b / (K + n^a)
}

rate <- decay_scheduler(gamma1 = 0.5, K = 100)

# Online learning
omega <- 4

sigma <- 1

max_iter <- 20000

par <- numeric(max_iter*4)

dim(par) <- c(max_iter, 4)

initpar <- c(alpha0+runif(1,-1,1), beta0+runif(1,-1,1), gamma0+runif(1,-1,1), rho0+runif(1,-1,1))

par[1, ] <- initpar

for (i in 2:max_iter){
  # Simulate online data 
  X <- rnorm(1, mean = 0, sd = omega)
  
  epsilon <- rnorm(1, mean = 0, sd = sigma) 
  
  Y <- densY(X, alpha0, beta0, gamma0, rho0) + epsilon
  
  par[i, ] <- par[(i-1), ] - rate(i) * grad_calc(X, Y, par[(i-1), 1],
                                                 par[(i-1), 2],
                                                 par[(i-1), 3],
                                                 par[(i-1), 4])
}

par

ggplot(tibble(x = seq_along(par[, 1]), y = par[, 1]), aes(x = x, y = y)) + geom_line() + geom_hline(yintercept = alpha0, col = "red")
ggplot(tibble(x = seq_along(par[, 1]), y = par[, 2]), aes(x = x, y = y)) + geom_line() + geom_hline(yintercept = beta0, col = "red")
ggplot(tibble(x = seq_along(par[, 1]), y = par[, 3]), aes(x = x, y = y)) + geom_line() + geom_hline(yintercept = gamma0, col = "red")
ggplot(tibble(x = seq_along(par[, 1]), y = par[, 4]), aes(x = x, y = y)) + geom_line() + geom_hline(yintercept = rho0, col = "red")
  

