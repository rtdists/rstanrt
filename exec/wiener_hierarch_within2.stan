// general model
data {
  int<lower=1> N; // number of participants
  int<lower=1> n_par[4]; // number of regression parameters per diffusion parameter
  int<lower=1> n_intercept[4];
  int<lower=0> n_other_pars[4];
  int pos_intercept[4, max(n_intercept)];
  int pos_other_pars[4, max(n_other_pars)];
  
  int<lower=1> K_u; // number of trials
  int<lower=1> K_l; // number of trials
  matrix[K_u,3] data_u; // [id, response, rt]
  matrix[K_l,3] data_l;
  int map_u[K_u,4]; // one column per parameter
  int map_l[K_l,4];
  int<lower=1> p_boundaries_u[N, 2];  // boundaries of each participant upper (data_map_u)
  int<lower=1> p_boundaries_l[N, 2];  // boundaries of each participant lower
  
  int<lower=0> rows_mmatrix[4];
  matrix[rows_mmatrix[1],n_par[1]] mmatrix_alpha; // model matrix fixed effects
  matrix[rows_mmatrix[2],n_par[2]] mmatrix_tau; // model matrix 
  matrix[rows_mmatrix[3],n_par[3]] mmatrix_beta; // model matrix 
  matrix[rows_mmatrix[4],n_par[4]] mmatrix_delta; // model matrix 
  
  int<lower=0> rows_mmatrix_indiv[4];
  int mmatrix_alpha_indiv[rows_mmatrix_indiv[1],2]; // mappings for individual model matrices (random effects)
  int mmatrix_tau_indiv[rows_mmatrix_indiv[2],2]; 
  int mmatrix_beta_indiv[rows_mmatrix_indiv[3],2]; 
  int mmatrix_delta_indiv[rows_mmatrix_indiv[4],2]; 
  int<lower=1> p_boundaries_2[N, 8];  // boundaries of each participant in random effects model matrix
  
}
transformed data{
  int nparams;
  int param_bounds_tau[2];
  int param_bounds_beta[2];
  int param_bounds_delta[2];
  nparams = sum(n_par);
  param_bounds_tau[1] = n_par[1]+1;
  param_bounds_tau[2] = sum(n_par[1:2]);
  param_bounds_beta[1] = sum(n_par[1:2])+1;
  param_bounds_beta[2] = sum(n_par[1:3]);
  param_bounds_delta[1] = sum(n_par[1:3])+1;
  param_bounds_delta[2] = sum(n_par[1:4]);
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
  matrix[nparams,N] deltahat; // individual level effects
  
  // individual level predictions (upper bound)
  vector<lower=0>[rows_mmatrix_indiv[1]] alpha_pred_indiv;
  vector<lower=0>[rows_mmatrix_indiv[2]] tau_pred_indiv;
  vector<lower=0, upper=1>[rows_mmatrix_indiv[3]] beta_pred_indiv;
  vector[rows_mmatrix_indiv[4]] delta_pred_indiv;
  
  // individual level predictions (lower bound, where they differ)
  vector<lower=0, upper=1>[rows_mmatrix_indiv[3]] beta_pred_indiv_lower;
  vector[rows_mmatrix_indiv[4]] delta_pred_indiv_lower;
  
  deltahat = diag_pre_multiply(sigma, L_Omega) * deltahat_tilde; // individual level effects
  
  // calculate individual predicted parameters (fixed + random):
  for (i in 1:N) {
    alpha_pred_indiv[ p_boundaries_2[i, 1]:p_boundaries_2[i, 2] ] =  
      mmatrix_alpha[mmatrix_alpha_indiv[p_boundaries_2[i, 1]:p_boundaries_2[i, 2],2],] *
      (alpha + deltahat[1:n_par[1],i]);
    tau_pred_indiv[ p_boundaries_2[i, 3]:p_boundaries_2[i, 4] ] = 
      mmatrix_tau[mmatrix_tau_indiv[p_boundaries_2[i, 3]:p_boundaries_2[i, 4],2],] *
      (tau + deltahat[param_bounds_tau[1]:param_bounds_tau[2],i]);
    beta_pred_indiv[ p_boundaries_2[i, 5]:p_boundaries_2[i, 6] ] =
      mmatrix_beta[mmatrix_beta_indiv[p_boundaries_2[i, 5]:p_boundaries_2[i, 6],2],] *
      (beta + deltahat[param_bounds_beta[1]:param_bounds_beta[2],i]);
    delta_pred_indiv[ p_boundaries_2[i, 7]:p_boundaries_2[i, 8] ] = 
      mmatrix_delta[mmatrix_delta_indiv[p_boundaries_2[i, 7]:p_boundaries_2[i, 8],2],] *
      (delta + deltahat[param_bounds_delta[1]:param_bounds_delta[2],i]);
  }
  beta_pred_indiv_lower = (1-beta_pred_indiv);
  delta_pred_indiv_lower = -delta_pred_indiv;
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
  
  data_u[,3] ~ wiener(
    alpha_pred_indiv[map_u[,1]], 
    tau_pred_indiv[map_u[,2]], 
    beta_pred_indiv[map_u[,3]],
    delta_pred_indiv[map_u[,4]]);
  data_l[,3] ~ wiener(
    alpha_pred_indiv[map_l[,1]], 
    tau_pred_indiv[map_l[,2]], 
    beta_pred_indiv_lower[map_l[,3]],
    delta_pred_indiv_lower[map_l[,4]]);
}
generated quantities {
  corr_matrix[nparams] Omega;
  //real log_lik[K_u + K_l];
  
  //log_lik[1:K_u] = wiener_lpdf(Y_u| alpha_u, tau_u, beta_u, delta_u + drift_crit_u);
  //log_lik[(K_u+1):(K_u+K_l)] = wiener_lpdf(Y_l| alpha_l, tau_l, (1-beta_l), -(delta_l + drift_crit_l));

  // Post-Processing Means, Standard Deviations, Correlations
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}

