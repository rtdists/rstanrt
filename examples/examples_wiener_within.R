
\dontrun{
# Speed Accuracy data from Wagenmakers et al. (2008):
data(speed_acc, package = "rtdists")
speed_acc <- droplevels(speed_acc[!speed_acc$censor,])
str(speed_acc)

## use only a subset of the data:
sa_red <- droplevels(speed_acc[speed_acc$frequency %in% c("high", "nw_high"),])
sa_red <- droplevels(sa_red[sa_red$id %in% as.character(1:10),])
sa_red <- droplevels(sa_red[ sample(seq_len(nrow(sa_red)), ceiling(nrow(sa_red)/20)),])
str(sa_red)
table(sa_red$response)

options(mc.cores = parallel::detectCores()) # package need to be installed.
fit_rr <- stan_wiener_within(alpha = ~0+condition, 
                             tau = ~1, 
                             beta = ~1, 
                             delta = ~1, 
                             data = sa_red, 
                             rt = "rt", response = "response", id = "id",
                             chains = 2, warmup = 250, iter = 500)
# shows some divergent transitions, but should work better with more data.

hyper_pars <- c("alpha", "tau", "beta", "delta")
hyper_struct_pars <- c("sigma", "Omega")

# run on windows to avoid warnings:
windowsFonts(Arial=windowsFont("TT Arial"))

rstan::stan_trace(fit_rr, pars = hyper_pars)
rstan::stan_trace(fit_rr, pars = hyper_struct_pars)

print(fit_rr, pars = c(hyper_pars, hyper_struct_pars))

# use bayesplot for more diagnostic plots
library("bayesplot")
posterior <- as.array(fit_rr)
np2 <- nuts_params(fit_rr)
mcmc_trace(posterior, regex_pars = "sigma", divergences = np2)
mcmc_trace(posterior, regex_pars = "alpha\\[", divergences = np2)


}
