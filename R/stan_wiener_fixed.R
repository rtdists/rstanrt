#' Fitting a fixed effects Wiener model
#' 
#' 
#' 
#' @example examples/examples_wiener_within.R
#' 
#' @rdname stan_wiener_fixed
#' @export
stan_wiener_fixed <- function(alpha, tau, beta, delta, 
                               data, rt, response, id,
                               pos_intercepts = "all",
                               chains = 4, warmup = 1000, iter = 2000,
                               thin = 1,
                               control = list(adapt_delta = 0.99, max_treedepth = 15), 
                               pars = c("alpha", "tau", "beta","delta"),
                               seed = sample.int(.Machine$integer.max, 1)
                               ) 
  {
  
  data_prep <- create_data_wiener_within(alpha = alpha, tau = tau, beta = beta, delta = delta,
                                         data = data, rt = rt, response = response, id = id,
                                         pos_intercepts = pos_intercepts)
  
  out <- sampling(stanmodels$wiener_fixed, 
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
