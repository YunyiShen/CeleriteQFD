library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = F)

N <- 600
tt <- seq(0,10,length.out = N)
f <- 5*(sin(0.5*tt)+1)*sin(3 * tt)
plot(tt,f)
y <- f + rnorm(length(tt),0,1)
plot(tt,y)
lines(tt,f)
celerite_data <- list(N=N, t = tt, y = y,
                     B_prior = c(-10,0),
                     L_prior = c(1.5,5),
                     P_prior = c(-3,5),
                     C_prior = c(-5,5),
                     err_prior = c(0.01,0.01),
                     diag = 1e-6+0*tt)

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celeriteQP.stan', 
            model_name = "celeritQP", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fit <- sampling(modelcelerite, data = celerite_data,control = list(adapt_delta = 0.9,max_treedepth=10), iter = 2000,chains = 2)
summ_fit <- summary(fit)
plot(summ_fit[[1]][1:N + (N+9),1])
plot(y-summ_fit[[1]][1:N + (N+9),1])

fit_array <- as.array(fit)
plot(fit_array[,1,N+7])
par(mfrow = c(1,2))

plot(t,10*f, type = "l")
lines(t,summ_fit[[1]][1:101 + 111,1],col = "red")
