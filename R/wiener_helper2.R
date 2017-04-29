
create_data_wiener_within_2 <- function(..., 
                                      data, rt, response, id,
                                      pos_intercepts) {
  N <- length(unique(data[,id]))
  if (!is.numeric(data[,id])) data[,id] <- as.numeric(data[,id])
  data_ordered <- data[order(data[,id]),] ## make sure data is ordered along id
  
  dots <- list(...) # formulas for parameters
  dots <- dots[c("alpha", "tau", "beta", "delta")] # make sure dots has correct ordering

  # mapping matrix that maps row of identical model matrices to participants, responses, and rts
  dat_pass <- data_ordered[,c(id, response, rt), drop=FALSE]
  dat_pass[,id] <- as.numeric(dat_pass[,id])
  dat_pass[,response] <- as.numeric(dat_pass[,response])
  
  par_maps <- matrix(-1, nrow = nrow(data_ordered), ncol = length(dots)*2, 
                     dimnames = list(NULL, c(paste0(names(dots), "_f"), paste0(names(dots), "_i"))))
  
  # lists of model matrices
  mms_full <- vector("list", length(dots))   # full
  mms_unique <- vector("list", length(dots)) # unique
  names(mms_unique) <- paste0("mmatrix_", names(dots))
  rows_mmatrix <- vector("numeric", length(dots))
  mms_unique_indiv <- vector("list", length(dots)) # unique
  names(mms_unique_indiv) <- paste0("mmatrix_", names(dots), "_indiv")
  rows_mmatrix_indiv <- vector("numeric", length(dots))
  
  #head(par_maps)
  ## create model matrices:
  # 1. full model matrix
  # 2. unique elements for fixed effects
  # 3. map unique elements to rows per parameter for fixed effects
  # 4. unique elements for each participant
  # 5. map unique elements to rows per parameter for each participant and parameter
  for (i in seq_along(dots)) {
    mms_full[[i]] <- model.matrix(dots[[i]], data = data_ordered) # 1.
    mms_unique[[i]] <- unique(mms_full[[i]]) # 2.
    rows_mmatrix[i] <- nrow(mms_unique[[i]])
    par_maps[,i] <- match(apply(mms_full[[i]], 1, paste, collapse ='/t'), 
                           apply(mms_unique[[i]], 1, paste, collapse ='/t'))  # 3.
    
    # step 4:
    tmp_indiv <- vector("list", N)
    for (j in seq_len(N)) {
      tmp_indiv[[j]] <- cbind(id = j, mms = sort(unique(par_maps[dat_pass[,id]==j,i])))
    }
    mms_unique_indiv[[i]] <- do.call("rbind", tmp_indiv)
    rows_mmatrix_indiv[i] <- nrow(mms_unique_indiv[[i]])
    ## step 5:
    par_maps[,i+length(dots)] <- match(apply(cbind(dat_pass[,id], par_maps[,i]), 1, paste, collapse ='/t'), 
                                        apply(mms_unique_indiv[[i]], 1, paste, collapse ='/t'))
    #mms_unique_indiv[[i]] <- mms_unique_indiv[[i]][,-1] # remove id column for cleaner model matrix 
  }
  
  # Split map_mat into upper and lower part
  dat_pass_split <- split(as.data.frame(dat_pass), f = dat_pass[,response])
  data_lower <- as.matrix(dat_pass_split[[1]])
  data_upper <- as.matrix(dat_pass_split[[2]])
  
  map_split <- split(as.data.frame(par_maps), f = dat_pass[,response])
  map_lower <- as.matrix(map_split[[1]])
  map_upper <- as.matrix(map_split[[2]])
  
  # numbers of observations
  K_l <- nrow(data_lower)
  K_u <- nrow(data_upper)
  
  ## mapping of participants to data/model matrices
  rle_p_lower <- cumsum(rle(data_lower[,1])$lengths)
  rle_p_upper <- cumsum(rle(data_upper[,1])$lengths)
  rle_p2 <- matrix(-1, nrow = N, ncol = length(dots)*2)
  for (i in seq_along(dots)) {
    tmp_rle <- cumsum(rle(mms_unique_indiv[[i]][,1])$lengths)
    rle_p2[,((i-1)*2+1):((i-1)*2+2)] <- cbind(c(1,tmp_rle[-length(tmp_rle)]+1),tmp_rle)
  }
  # count parameters and create intercept-other mapping:
  n_par <- vapply(mms_unique, ncol, 1)
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
    pos_other_pars[i, seq_len(n_par[i] - length(pos_intercepts[[i]]))]  <- 
      seq_len(n_par[i])[-pos_intercept[i, 1:length(pos_intercepts[[i]])]]
  }
  #browser()
  #str(
    c(N=list(N),
      K_u=list(K_u),
      K_l=list(K_l),
      p_boundaries_u = list(cbind(c(1,rle_p_upper[-length(rle_p_upper)]+1),rle_p_upper)),
      p_boundaries_l = list(cbind(c(1,rle_p_lower[-length(rle_p_lower)]+1),rle_p_lower)),
      p_boundaries_2 = list(rle_p2),
      data_u = list(data_upper),
      data_l = list(data_lower),
      map_u = list(map_upper[,5:8]),
      map_l = list(map_lower[,5:8]),
      mms_unique,
      mms_unique_indiv,
      rows_mmatrix = list(rows_mmatrix),
      rows_mmatrix_indiv = list(rows_mmatrix_indiv),
      n_par = list(n_par),
      n_intercept = list(n_intercept),
      n_other_pars = list(n_other_pars),
      pos_intercept = list(pos_intercept),
      pos_other_pars = list(pos_other_pars)
  )
  #,1)
}


wiener_within_init2 <- function(data) {
  N <- data$N
  nparams <- sum(data$n_par)
  tmp_out <- list(deltahat_tilde=matrix(rnorm(N * nparams, 0, 0.001), nparams, N), 
                  L_Omega=diag(nparams),
                  sigma = runif(nparams, 0.1, 0.5)
  )
  tmp_out$alpha <- rep(0, data$n_par[1])
  tmp_out$alpha[data$pos_intercept[1,seq_len(data$n_intercept[1])]] <- runif(data$n_intercept[1], 0.5, 1.5)
  tmp_out$alpha[data$pos_other_pars[1,seq_len(data$n_other_pars[1])]] <- rnorm(data$n_other_pars[1], 0, 0.1)
  
  tmp_out$tau <- rep(0, data$n_par[2])
  tmp_out$tau[data$pos_intercept[2,seq_len(data$n_intercept[2])]] <- pmin(runif(data$n_intercept[2], 0.0, min(c(data$data_u[,3], data$data_u[,3])-0.05)), 0.05)
  tmp_out$tau[data$pos_other_pars[2,seq_len(data$n_other_pars[2])]] <- rnorm(data$n_other_pars[2], 0, 0.01)
  
  
  tmp_out$beta <- rep(0, data$n_par[3])
  tmp_out$beta[data$pos_intercept[3,seq_len(data$n_intercept[3])]] <- runif(data$n_intercept[3], 0.3, 0.7)
  tmp_out$beta[data$pos_other_pars[3,seq_len(data$n_other_pars[3])]] <- rnorm(data$n_other_pars[3], 0, 0.01)
  
  
  tmp_out$delta <- rep(0, data$n_par[4])
  tmp_out$delta[data$pos_intercept[4,seq_len(data$n_intercept[4])]] <- rnorm(data$n_intercept[4], 0, 2)
  tmp_out$delta[data$pos_other_pars[4,seq_len(data$n_other_pars[4])]] <- rnorm(data$n_other_pars[4], 0, 1)
  
  for (i in seq_along(tmp_out)) {
    if (is.null(dim(tmp_out[[i]]))) {
      tmp_out[[i]] <- as.array(tmp_out[[i]])
    }
  }
  
  return(tmp_out)
}