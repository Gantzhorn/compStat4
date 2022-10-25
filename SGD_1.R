decay_scheduler <- function(gamma1 = 1, a = 1, K = 1, gamma2, n1) {
  force(a)
  if (!missing(gamma2) && !missing(n1))
    K <- n1^a * gamma2 / (gamma0 - gamma2)
  b <- gamma1 * K
  function(n) b / (K + n^a)
}

rate <- decay_scheduler(gamma1 = 0.0004, K = 100)

beta <- vector("list", N)

initpar <- c(0.5,0.3, 0.1, 1.5)

beta[[1]] <- initpar
for (i in 2:N){
  beta[[i]] <- beta[[i-1]] - rate(i) * grad(X[i], alpha0, beta0, gamma0, rho0)
}

