functions {
   real logLikSHO(vector t, vector y, vector jitter,real S0, real w0, real Q, real eps){
       return(w0);
   }

   real logLikRotation(vector t, vector y, vector jitter,real sigma, real period, real Q0, real dQ, real eps){
       return(0);
   }
}
