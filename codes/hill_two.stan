data {
  int <lower=0> len;          // number of dose points 
  real <lower=0> d[len];      // dose levels  
  real <lower=0> y[len];      // observed responses
  real p_b[2];                // prior for ec50
  real p_n[2];                // prior for hill coefficient
  real p_sig[2];              // prior for sigma
} 
parameters {
  real <lower=0> b;
  real <lower=0> n;
  real <lower=0> sigma;
}
model {
  vector[len] yp;
  b ~ uniform_lpdf (p_b[1], p_b[2]);
  n ~ uniform (p_n[1], p_n[2]);
  sigma ~ normal (p_sig[1], p_sig[2]);
  for (j in 1:len)
    yp[j] = 1 / (1 + d[j] / b )^n;
    y ~ student_t(5, yp, sigma);  
}
