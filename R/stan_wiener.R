#' Fitting a completely within-subjects Wiener model
#' 
#' @import rstan
#' @import rstantools
#' @import Rcpp
#' @useDynLib rstanrt
#' 
#' @example examples/examples_wiener_within.R
#' 
#' @export
stan_wiener_within <- function(alpha, tau, beta, delta, 
                               data, rt, response, id,
                               pos_intercepts = "all",
                               chains = 4, warmup = 1000, iter = 2000,
                               thin = 1,
                               control = list(adapt_delta = 0.99, max_treedepth = 15), 
                               pars = c("alpha", "tau", "beta", "delta", "sigma", "deltahat", "Omega")
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
                  control = control)
  return(out)
}
