simuHMM <- function(T, mu, sigma, theta, p){
    res <- rep(NA, T)
    state <- 1 + (runif(1) > p[1])
    res[1] <- rnorm(1,mu[state],sigma[state])
    for(t in 2:T){
        state <- 1 + (runif(1) > theta[state,1])
        res[t] <- rnorm(1,mu[state],sigma[state])
    }
    return(res)
}
