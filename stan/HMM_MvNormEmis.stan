functions {
  
  /**  generate a N * K matrix with
  multivariate Gaussian emissions for t=1,..,N and j=1,..,K.
  */
  matrix mvnormEmissions(int N, int K, int D, vector[] y, vector[] mu, 
                          vector[] tau, matrix[] L_Omega) 
  {
    matrix[N, K] allprobs; 
    for (k in 1:K) {
      matrix[D, D] Omega;
      matrix[D, D] Sigma;
      Omega = multiply_lower_tri_self_transpose(L_Omega[k]);
      Sigma = quad_form_diag(Omega, tau[k]);
      for (n in 1:N) {
        allprobs[n, k] = exp(multi_normal_lpdf(y[n] |
                             mu[k], Sigma));
      }
    }
    return allprobs; 
  }
  
    /** perform forward algorithm to 
  generate alpha_t (j) for t=1,..,N and j=1,.., M = sum(m)
  via forward dynamic programming - p.243 (Zucchini et al.).
  */
  real forwardMessages(int N, int K, matrix emissions, matrix gamma_mat)
  {
    vector[K] foo;
    real sumfoo;
    real lscale;
    vector[K] gamma_init; 
    // alpha_1
    for (k in 1:K) {
      foo[k] = emissions[1, k];
    }
    sumfoo = sum(foo);
    lscale = log(sumfoo);
    foo = foo/sumfoo;
    // alpha_t, t = 2, ..., N
    for (i in 2:N) {
      foo = (foo'*gamma_mat .* emissions[i, :])';
      sumfoo = sum(foo);
      lscale = lscale + log(sumfoo);
      foo = foo/sumfoo;
    }
    return lscale;
  }
  
  /** convert simplex gamma to a K*K matrix. 
  */
  matrix gammaToMat(int K, vector[] gamma)
  {
    matrix[K, K] gamma_mat;
    // Creating Gamma Matrix for Forward
    for(i in 1:K)
      for(j in 1:K)
        gamma_mat[i, j] = gamma[i][j];
    return gamma_mat;
  }

  /** compute likelihood p(y_1, .., y_T  | )
  */
  real llk_lp(int N, int K, int D, vector[] y, vector[] mu,
              vector[] tau, matrix[] L_Omega, vector[] gamma)
  {
    matrix[N, K] emissions = mvnormEmissions(N, K, D, y, mu, tau, L_Omega);
    matrix[K, K] gamma_mat = gammaToMat(K, gamma);
    real llk =  forwardMessages(N, K, emissions, gamma_mat);
    return llk;
  }
}

data {
  int<lower=0> N; // length time series
  int<lower=1> D; // dimensionality // lower bound is 1 becuase of the deomposition of mu
  int<lower=0> K; // number of states
  vector[D] y[N]; // data

  // hyperparms
  real mu_loc; 
  real<lower=0> mu_scale; 
  real tau_loc; 
  real<lower=0> tau_scale;
  real<lower=0> Omega_shape;
  vector<lower=0>[K] alpha_0[K]; // prior dirichlet probs
}

parameters {
  simplex[K] gamma[K]; // transition probs
  ordered[K] mu_first; // mean gauss emission - first dimension 
  vector[D - 1] mu_rest[K]; // mean gauss emission - rest of the dimension 
  cholesky_factor_corr[D] L_Omega[K]; // corr gauss emission
  vector<lower=0>[D] tau[K];  // scales gauss emission
}


transformed parameters { 
  vector[D] mu[K];
  for(k in 1:K){
    mu[k] = append_row(mu_first[k], mu_rest[k]);
  }
}


model {
  // priors
  for (k in 1:K) {
    target += lkj_corr_cholesky_lpdf(L_Omega[k] | Omega_shape);
    target += student_t_lpdf(tau[k] | 4, tau_loc, tau_scale);
    target += dirichlet_lpdf(gamma[k] | alpha_0[k]);
    target += student_t_lpdf(to_vector(mu_rest[k]) | 4,  mu_loc, mu_scale);
  }
  target += student_t_lpdf(to_vector(mu_first) | 4,  mu_loc, mu_scale);
  // likelihood
  target += llk_lp(N, K, D, y, mu, tau, L_Omega, gamma); 
}
