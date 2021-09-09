library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = F)

N <- 300
tt <- seq(0,5,length.out = N)
f <- 5*sin(3 * tt)
#plot(t,f)
y <- f + rnorm(length(tt),0,1)
plot(tt,y)
lines(tt,f)
celerite_data <- list(N=N, t = tt, y = y,
                     S0_prior = c(-10,10),
                     w0_prior = c(-10,10),
                     Q_prior = c(-10,10),
                     err_prior = c(0.01,0.01),
                     diag = 1e-6+0*tt)

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celeriteSHO.stan', 
            model_name = "celeritSHO", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))
fit2 <- vb(modelcelerite, data = celerite_data,tol_rel_obj=0.001)
summ_celerite2 <- summary(fit2)

fitopt <- optimizing(modelcelerite, data = celerite_data)
para <- summ_celerite2[[1]][1:23 + N,1]

lomega <- seq(-5,5,0.01 )
omega <- 10^(lomega) 
ps1 <- sqrt(2/pi) *  (para["S1"] * para["w1"]^4)/((omega^2-para["w1"]^2)^2+para["w1"]^2*(omega^2)/(para["Q1"]^2))
ps2 <- sqrt(2/pi) *  (para["S2"] * para["w2"]^4)/((omega^2-para["w2"]^2)^2+para["w2"]^2*(omega^2)/(para["Q2"]^2))

plot(lomega,log10(ps1), type = "l")
lines(lomega,log10(ps2), col = "red")


summ_celerite2[[1]][1:22 + N,1]
plot(summ_celerite2[[1]][1:N + (N+23),1], type = "l")
points(y,col = "red")

fit <- sampling(modelcelerite, data = celerite_data,control = list(adapt_delta = 0.9,max_treedepth=10), iter = 4000)
summ_fit <- summary(fit)
plot(summ_fit[[1]][1:N + (N+9),1])
plot(y-summ_fit[[1]][1:N + (N+9),1])

fit_array <- as.array(fit)
plot(fit_array[,1,N+6])
par(mfrow = c(1,2))

plot(t,10*f, type = "l")
lines(t,summ_fit[[1]][1:101 + 111,1],col = "red")
