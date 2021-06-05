stan_model(file = "./Stan/celerite.stan", model_name = "celerit", 
            allow_undefined = TRUE,
           includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))



mc <- '
functions { real logLikSHO(vector t, vector y,real S0, real w0, real Q, real eps, vector diag){
    real fw;
    fw = sqrt(S0);
    return(0);
} }
model {} // use the fib() function somehow
'

writeLines(stanc(model_code = mc, allow_undefined = TRUE)$cppcode,"./minimum.hpp")
