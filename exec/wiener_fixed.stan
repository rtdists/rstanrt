// general model
data {
  //int<lower=1> N; // number of participants
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
  //int<lower=1> p_boundaries_u[N, 2];  // boundaries of each participant in Y_u

  int<lower=1> K_l; // number of trials reaching lower boundary
  vector[K_l] Y_l; // response times reaching lower boundary
  matrix[K_l,n_par[1]] mmatrix_alpha_l; // model matrix for upper boundary
  matrix[K_l,n_par[2]] mmatrix_tau_l; // model matrix for upper boundary
  matrix[K_l,n_par[3]] mmatrix_beta_l; // model matrix for upper boundary
  matrix[K_l,n_par[4]] mmatrix_delta_l; // model matrix for upper boundary
  //int<lower=1> p_boundaries_l[N, 2];  // boundaries of each participant in Y_l

}
transformed data{
int nparams;
nparams = sum(n_par);
}
parameters {
  
  // correlation among parameters and individual effects;
  //matrix[nparams,N] deltahat_tilde; 
  //cholesky_factor_corr[nparams] L_Omega; 
  //vector<lower=0>[nparams] sigma; 
  
  // hyper parameters
  vector[n_par[1]] alpha;
  vector[n_par[2]] tau;
  vector[n_par[3]] beta;
  vector[n_par[4]] delta;
}
transformed parameters{
  //matrix[nparams,N] deltahat;
  
  vector<lower=0>[K_u] alpha_u;
  vector<lower=0>[K_u] tau_u;
  vector<lower=0, upper=1>[K_u] beta_u;
  vector[K_u] delta_u;
  
  vector<lower=0>[K_l] alpha_l;
  vector<lower=0>[K_l] tau_l;
  vector<lower=0, upper=1>[K_l] beta_l;
  vector[K_l] delta_l;
  
  //deltahat = diag_pre_multiply(sigma, L_Omega) * deltahat_tilde; 
  
  alpha_u =  mmatrix_alpha_u * alpha;
  alpha_l =  mmatrix_alpha_l * alpha;
    
  tau_u = mmatrix_tau_u * tau;
  tau_l = mmatrix_tau_l * tau;
    
  beta_u = mmatrix_beta_u * beta;
  beta_l = mmatrix_beta_l * beta;
    
  delta_u = mmatrix_delta_u * delta;
  delta_l = mmatrix_delta_l * delta;
}
model {
  //L_Omega ~ lkj_corr_cholesky(1); 
  //sigma ~ cauchy(0, 4); 
  //to_vector(deltahat_tilde) ~ normal(0, 1);
  
  alpha[pos_intercept[1,1:n_intercept[1]]] ~ cauchy(1, 2);
  if (n_other_pars[1] > 0)
    alpha[pos_other_pars[1,1:n_other_pars[1]]] ~ cauchy(0, 1);
  
  tau[pos_intercept[2,1:n_intercept[2]]] ~ normal(0.2, 0.2);
  if (n_other_pars[2] > 0)
    tau[pos_other_pars[2,1:n_other_pars[2]]] ~ normal(0, 0.2);
  
  beta[pos_intercept[3,1:n_intercept[3]]] ~ normal(0.5, 0.5);
  if (n_other_pars[3] > 0)
    beta[pos_other_pars[3,1:n_other_pars[3]]] ~ normal(0, 0.5);
  
  delta[pos_intercept[4,1:n_intercept[4]]] ~ cauchy(0, 3);
  if (n_other_pars[4] > 0)
    delta[pos_other_pars[4,1:n_other_pars[4]]] ~ cauchy(0, 2);
  
  Y_u ~ wiener(alpha_u, tau_u, beta_u, delta_u);
  Y_l ~ wiener(alpha_l, tau_l, (1-beta_l), -(delta_l));

}
