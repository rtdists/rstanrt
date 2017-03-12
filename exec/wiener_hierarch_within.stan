// general model
data {
  int<lower=1> N; // number of participants
  int<lower=1> n_par[4]; // number of regression parameters per diffusion parameter
  int<lower=1> n_intercept[4];
  int<lower=0> n_other_pars[4];
  int pos_intercept[4, max(n_intercept)];
  int pos_other_pars[4, max(n_other_pars)];
  
  int<lower=1> K_u; // number of trials reaching upper boundary
  vector[K_u] Y_u; // response times reaching upper boundary
  matrix[K_u,n_par[1]] mmatrix_alpha_u; // model matrix for upper boundary
  matrix[K_u,n_par[2]] mmatrix_tau_u; // model matrix for upper boundary
  matrix[K_u,n_par[3]] mmatrix_beta_u; // model matrix for upper boundary
  matrix[K_u,n_par[4]] mmatrix_delta_u; // model matrix for upper boundary
  int<lower=1> p_boundaries_u[N, 2];  // boundaries of each participant in Y_u

  int<lower=1> K_l; // number of trials reaching lower boundary
  vector[K_l] Y_l; // response times reaching lower boundary
  matrix[K_l,n_par[1]] mmatrix_alpha_l; // model matrix for upper boundary
  matrix[K_l,n_par[2]] mmatrix_tau_l; // model matrix for upper boundary
  matrix[K_l,n_par[3]] mmatrix_beta_l; // model matrix for upper boundary
  matrix[K_l,n_par[4]] mmatrix_delta_l; // model matrix for upper boundary
  int<lower=1> p_boundaries_l[N, 2];  // boundaries of each participant in Y_l

}
transformed data{
int nparams;
nparams = sum(n_par);
}
parameters {
  
  // correlation among parameters and individual effects;
  matrix[nparams,N] deltahat_tilde; 
  cholesky_factor_corr[nparams] L_Omega; 
  vector<lower=0>[nparams] sigma; 
  
  // hyper parameters
  vector[n_par[1]] alpha;
  vector[n_par[2]] tau;
  vector[n_par[3]] beta;
  vector[n_par[4]] delta;
}
transformed parameters{
  matrix[nparams,N] deltahat;
  
  vector<lower=0>[K_u] alpha_u;
  vector<lower=0>[K_u] tau_u;
  vector<lower=0, upper=1>[K_u] beta_u;
  vector[K_u] delta_u;
  
  vector<lower=0>[K_l] alpha_l;
  vector<lower=0>[K_l] tau_l;
  vector<lower=0, upper=1>[K_l] beta_l;
  vector[K_l] delta_l;
  
  deltahat = diag_pre_multiply(sigma, L_Omega) * deltahat_tilde; 
  
  for (i in 1:N) {
    alpha_u[ p_boundaries_u[i, 1]:p_boundaries_u[i, 2] ] =  mmatrix_alpha_u[p_boundaries_u[i, 1]:p_boundaries_u[i, 2],] * (alpha + deltahat[1:n_par[1],i]);
    alpha_l[ p_boundaries_l[i, 1]:p_boundaries_l[i, 2] ] =  mmatrix_alpha_l[p_boundaries_l[i, 1]:p_boundaries_l[i, 2],] * (alpha + deltahat[1:n_par[1],i]);
    
    tau_u[ p_boundaries_u[i, 1]:p_boundaries_u[i, 2] ] = mmatrix_tau_u[p_boundaries_u[i, 1]:p_boundaries_u[i, 2],] * (tau + deltahat[(n_par[1]+1):(sum(n_par[1:2])),i]);
    tau_l[ p_boundaries_l[i, 1]:p_boundaries_l[i, 2] ] = mmatrix_tau_l[p_boundaries_l[i, 1]:p_boundaries_l[i, 2],] * (tau + deltahat[(n_par[1]+1):(sum(n_par[1:2])),i]);
    
    beta_u[ p_boundaries_u[i, 1]:p_boundaries_u[i, 2] ] = mmatrix_beta_u[p_boundaries_u[i, 1]:p_boundaries_u[i, 2],] * (beta + deltahat[(sum(n_par[1:2])+1):(sum(n_par[1:3])),i]);
    beta_l[ p_boundaries_l[i, 1]:p_boundaries_l[i, 2] ] = mmatrix_beta_l[p_boundaries_l[i, 1]:p_boundaries_l[i, 2],] * (beta + deltahat[(sum(n_par[1:2])+1):(sum(n_par[1:3])),i]);
    
    delta_u[ p_boundaries_u[i, 1]:p_boundaries_u[i, 2] ] = mmatrix_delta_u[p_boundaries_u[i, 1]:p_boundaries_u[i, 2],] * (delta + deltahat[(sum(n_par[1:3])+1):(sum(n_par[1:4])),i]);
    delta_l[ p_boundaries_l[i, 1]:p_boundaries_l[i, 2] ] = mmatrix_delta_l[p_boundaries_l[i, 1]:p_boundaries_l[i, 2],] * (delta + deltahat[(sum(n_par[1:3])+1):(sum(n_par[1:4])),i]);
  }
}
model {
  L_Omega ~ lkj_corr_cholesky(1); 
  sigma ~ cauchy(0, 4); 
  to_vector(deltahat_tilde) ~ normal(0, 1);
  
  alpha[pos_intercept[1,1:n_intercept[1]]] ~ cauchy(1, 2);
  if (n_other_pars[1] > 0)
    alpha[pos_other_pars[1,1:n_other_pars[1]]] ~ cauchy(0, 1);
  
  tau[pos_intercept[2,1:n_intercept[2]]] ~ normal(0.2, 0.2);
  if (n_other_pars[2] > 0)
    tau[pos_other_pars[2,1:n_other_pars[2]]] ~ normal(0, 0.2);
  
  beta[pos_intercept[3,1:n_intercept[3]]] ~ normal(0.5, 0.5);
  if (n_other_pars[3] > 0)
    beta[pos_other_pars[3,1:n_other_pars[3]]] ~ normal(0, 0.5);
  
  delta[pos_intercept[4,1:n_intercept[4]]] ~ cauchy(1, 2);
  if (n_other_pars[4] > 0)
    delta[pos_other_pars[4,1:n_other_pars[4]]] ~ cauchy(0, 2);
  
  Y_u ~ wiener(alpha_u, tau_u, beta_u, delta_u);
  Y_l ~ wiener(alpha_l, tau_l, (1-beta_l), -(delta_l));

}
generated quantities {
  corr_matrix[nparams] Omega;
  //real log_lik[K_u + K_l];
  
  //log_lik[1:K_u] = wiener_lpdf(Y_u| alpha_u, tau_u, beta_u, delta_u + drift_crit_u);
  //log_lik[(K_u+1):(K_u+K_l)] = wiener_lpdf(Y_l| alpha_l, tau_l, (1-beta_l), -(delta_l + drift_crit_l));

  // Post-Processing Means, Standard Deviations, Correlations
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}

