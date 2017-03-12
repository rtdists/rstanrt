
create_data_wiener_within <- function(alpha, tau, beta, delta, 
                                      data, rt, response, id,
                                      pos_intercepts) {
  N <- length(unique(data[,id]))
  if (!is.numeric(data[,id])) data[,id] <- as.numeric(data[,id])
  data_ordered <- data[order(data[,id]),] ## make sure data is ordered along id
  d_split <- split(data_ordered, f = data_ordered[,response])
  
  ## model matrices for each parameter
  mmatrix_alpha <- split(as.data.frame(model.matrix(alpha, data = data_ordered)), f = data_ordered[,response])
  mmatrix_tau <- split(as.data.frame(model.matrix(tau, data = data_ordered)), f = data_ordered[,response])
  mmatrix_beta <- split(as.data.frame(model.matrix(beta, data = data_ordered)), f = data_ordered[,response])
  mmatrix_delta <- split(as.data.frame(model.matrix(delta, data = data_ordered)), f = data_ordered[,response])
  
  # rts per boundary
  K_u <- nrow(d_split[[2]])
  K_l <- nrow(d_split[[1]])
  #browser()
  ## mapping of participants to trials
  rle_p_upper <- cumsum(rle(d_split[[2]][,id])$lengths)
  rle_p_lower <- cumsum(rle(d_split[[1]][,id])$lengths)
  
  # count parameters and create intercept-other mapping:
  n_par <- vapply(list(mmatrix_alpha, mmatrix_tau, mmatrix_beta, mmatrix_delta), function(x) ncol(x[[1]]), 1)
  if (pos_intercepts=="all") {
    pos_intercepts <- lapply(n_par, function(x) seq(1, x))
  } else if (is.null(pos_intercepts)) {
    pos_intercepts <- as.list(rep(1, 5))
  } 
  n_intercept <- vapply(pos_intercepts, length, 1)
  n_other_pars <- n_par - n_intercept
  pos_intercept <- matrix(0, nrow = 4, ncol = max(n_intercept))
  pos_other_pars <- matrix(0, nrow = 4, ncol = max(n_other_pars))
  for (i in seq_along(n_par)) {
    pos_intercept[i, 1:length(pos_intercepts[[i]])]  <- pos_intercepts[[i]]
    pos_other_pars[i, seq_len(n_par[i] - length(pos_intercepts[[i]]))]  <- seq_len(n_par[i])[-pos_intercept[i, 1:length(pos_intercepts[[i]])]]
  }
  list(N=N,
       K_u=K_u,
       K_l=K_l,
       Y_u=d_split[[2]][,rt],
       Y_l=d_split[[1]][,rt],
       p_boundaries_u = cbind(c(1,rle_p_upper[-length(rle_p_upper)]+1),rle_p_upper),
       p_boundaries_l = cbind(c(1,rle_p_lower[-length(rle_p_lower)]+1),rle_p_lower),
       mmatrix_alpha_u = as.matrix(mmatrix_alpha[[2]]),
       mmatrix_alpha_l = as.matrix(mmatrix_alpha[[1]]),
       mmatrix_tau_u = as.matrix(mmatrix_tau[[2]]),
       mmatrix_tau_l = as.matrix(mmatrix_tau[[1]]),
       mmatrix_beta_u = as.matrix(mmatrix_beta[[2]]),
       mmatrix_beta_l = as.matrix(mmatrix_beta[[1]]),
       mmatrix_delta_u = as.matrix(mmatrix_delta[[2]]),
       mmatrix_delta_l = as.matrix(mmatrix_delta[[1]]),
       n_par = n_par,
       n_intercept = n_intercept,
       n_other_pars = n_other_pars,
       pos_intercept = pos_intercept,
       pos_other_pars = pos_other_pars
  )
}



wiener_within_init <- function(data) {
  N <- data$N
  nparams <- sum(data$n_par)
  tmp_out <- list(deltahat_tilde=matrix(rnorm(N * nparams, 0, 0.01), nparams, N), 
                  L_Omega=diag(nparams),
                  sigma = runif(nparams)
  )
  tmp_out$alpha <- rep(0, data$n_par[1])
  tmp_out$alpha[data$pos_intercept[1,seq_len(data$n_intercept[1])]] <- runif(data$n_intercept[1], 0.5, 1.5)
  tmp_out$alpha[data$pos_other_pars[1,seq_len(data$n_other_pars[1])]] <- rnorm(data$n_other_pars[1], 0, 0.1)
  
  tmp_out$tau <- rep(0, data$n_par[2])
  tmp_out$tau[data$pos_intercept[2,seq_len(data$n_intercept[2])]] <- runif(data$n_intercept[2], 0.0, min(c(data$Y_u, data$Y_l)-0.05))
  tmp_out$tau[data$pos_other_pars[2,seq_len(data$n_other_pars[2])]] <- rnorm(data$n_other_pars[2], 0, 0.01)
  
  
  tmp_out$beta <- rep(0, data$n_par[3])
  tmp_out$beta[data$pos_intercept[3,seq_len(data$n_intercept[3])]] <- runif(data$n_intercept[3], 0.3, 0.7)
  tmp_out$beta[data$pos_other_pars[3,seq_len(data$n_other_pars[3])]] <- rnorm(data$n_other_pars[3], 0, 0.01)
  
  
  tmp_out$delta <- rep(0, data$n_par[4])
  tmp_out$delta[data$pos_intercept[4,seq_len(data$n_intercept[4])]] <- runif(data$n_intercept[4], 0.5, 4)
  tmp_out$delta[data$pos_other_pars[4,seq_len(data$n_other_pars[4])]] <- rnorm(data$n_other_pars[4], 0, 1)
  
  for (i in seq_along(tmp_out)) {
    if (is.null(dim(tmp_out[[i]]))) {
      tmp_out[[i]] <- as.array(tmp_out[[i]])
    }
  }
  
  return(tmp_out)
}