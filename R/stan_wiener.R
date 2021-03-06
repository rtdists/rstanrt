#' Fitting a completely within-subjects Wiener model
#' 
#' For all functions implemented here, all parameters are assumed to vary within subjects and the full variance-covariance matrix is calculated (the correlation parameter is \code{Omega}). \code{stan_wiener_within2} uses a slightly different model implementation that could be faster (at least 10%), but is not tested very thoroughly.
#' 
#' @import rstan
#' @import rstantools
#' @import Rcpp
#' @useDynLib rstanrt, .registration = TRUE
#' 
#' @example examples/examples_wiener_within.R
#' 
#' @rdname stan_wiener_within
#' @export
stan_wiener_within_fixed_beta <- function(alpha, tau, delta, 
                               data, rt, response, id,
                               pos_intercepts = "all",
                               chains = 4, warmup = 1000, iter = 2000,
                               thin = 1,
                               control = list(adapt_delta = 0.99, max_treedepth = 15), 
                               pars = c("alpha", "tau", "delta", "sigma", "deltahat", "Omega"),
                               seed = sample.int(.Machine$integer.max, 1)
                               ) 
  {
  
  data_prep <- create_data_wiener_within(alpha = alpha, tau = tau, beta = ~1, delta = delta,
                                         data = data, rt = rt, response = response, id = id,
                                         pos_intercepts = pos_intercepts)
  
  out <- sampling(stanmodels$wiener_hierarch_within_fixed_beta, 
                  data = data_prep,
                  init=lapply(seq_len(chains), function(x) wiener_within_init(data = data_prep)),
                  pars = pars,
                  warmup = warmup,
                  iter = iter,
                  chains = chains,
                  thin = thin,
                  control = control,
                  seed=seed)
  return(out)
}

#' @rdname stan_wiener_within
#' @export
stan_wiener_within <- function(alpha, tau, beta, delta, 
                               data, rt, response, id,
                               pos_intercepts = "all",
                               chains = 4, warmup = 1000, iter = 2000,
                               thin = 1,
                               control = list(adapt_delta = 0.99, max_treedepth = 15), 
                               pars = c("alpha", "tau", "beta", "delta", "sigma", "deltahat", "Omega"),
                               seed = sample.int(.Machine$integer.max, 1)
                               ) 
  {
  
  data_prep <- create_data_wiener_within(alpha = alpha, tau = tau, beta = beta, delta = delta,
                                         data = data, rt = rt, response = response, id = id,
                                         pos_intercepts = pos_intercepts)
  
  out <- sampling(stanmodels$wiener_hierarch_within, 
                  data = data_prep,
                  init=lapply(seq_len(chains), function(x) wiener_within_init(data = data_prep)),
                  pars = pars,
                  warmup = warmup,
                  iter = iter,
                  chains = chains,
                  thin = thin,
                  control = control,
                  seed=seed)
  return(out)
}

#' @export
#' @rdname stan_wiener_within
stan_wiener_within2 <- function(alpha, tau, beta, delta, 
                               data, rt, response, id,
                               pos_intercepts = "all",
                               chains = 4, warmup = 1000, iter = 2000,
                               thin = 1,
                               control = list(adapt_delta = 0.99, max_treedepth = 15), 
                               pars = c("alpha", "tau", "beta", "delta", "sigma", "deltahat", "Omega"),
                               seed = sample.int(.Machine$integer.max, 1)
                               ) 
  {
  
  data_prep <- create_data_wiener_within_2(alpha = alpha, tau = tau, beta = beta, delta = delta,
                                         data = data, rt = rt, response = response, id = id,
                                         pos_intercepts = pos_intercepts)
  
  #ppp <- wiener_within_init2(data = data_prep)
  out <- sampling(stanmodels$wiener_hierarch_within2, 
  #out <- stan("exec/wiener_hierarch_within2.stan", 
                  data = data_prep,
                  init=lapply(seq_len(chains), function(x) wiener_within_init2(data = data_prep)),
                  pars = pars,
                  warmup = warmup,
                  iter = iter,
                  chains = chains,
                  thin = thin,
                  control = control,
                  seed=seed)
  return(out)
}