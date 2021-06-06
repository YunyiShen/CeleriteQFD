library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

t <- seq(0,100,1)
f <- sin(t/(2*pi))
plot(t,f)
y <- 5*f + rnorm(length(t),0,1)
plot(t,y)


test_GP_data <- list(N=length(t), x = t, y = y, foo = c(0,1))

fit <- stan(file = './Stan/GP_latent.stan', data = test_GP_data, 
                iter = 2000, control = list(adapt_delta = .8))
stan_array <- as.array(fit)
summ_fit <- summary(fit)
plot(summ_fit[[1]][4:105,1], type = "l")
lines(10*f, col = "red")
points(test_GP_data$y, col = "blue")
