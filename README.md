# rstanrt
Bayesian Estimation of Response Time Models

## Installation

* From CRAN: In the future...

* Development version from Github:  
For installing the development version package `devtools` is needed which may require some additional software (see [here](http://r-pkgs.had.co.nz/intro.html) section "Getting started")
`devtools::install_github("rtdists/rstanrt")`

## Functionality

The package has currently only one exported function that allows to fit completely within-subject models `stan_wiener_within`. Currently it does not allow custom priors or other setting, but will in the future.

```R
#### Can we recover parameters with small N.
## The following small simulation should:
## Recover means quite well, 
## Recover variances okay
## Do not recover the correlations (too few data).

require(RWiener) # for random generation.
N_participants <- 10
N_trials <- 20
mus <- c(alpha = 1, tau = 0.3, beta = 0.5, delta = 1.5)  # hyperparameter means
vars <- c(alpha = 0.2, tau = 0.05, beta = 0.1, delta = 0.5) # SD of individual level parameters
cors <- matrix(0, 4, 4)  ## create correlation matrix of parameters
cors[upper.tri(cors)] <- c(0, 0.2, 0.3, 0.7, 0.5, 0)
cors <- cors + t(cors)
diag(cors) <- 1
Sigma <- MBESS::cor2cov(cors, vars)  # package MBESS required
pars <- MASS::mvrnorm(N_participants,mu = mus, Sigma = Sigma)

d_prep <- vector("list", N_participants)
for (i in seq_len(N_participants)) {
d_prep[[i]] <-cbind( id = factor(i), rwiener(N_trials, alpha = pars[i,"alpha"], tau = pars[i,"tau"], beta = pars[i,"beta"], delta = pars[i,"delta"]))
}
d_prep <- do.call("rbind", d_prep)
str(d_prep)

options(mc.cores = parallel::detectCores()) # package need to be installed.
recov_pars <- stan_wiener_within(alpha = ~1, 
                             tau = ~1, 
                             beta = ~1, 
                             delta = ~1, 
                             data = d_prep, 
                             rt = "q", response = "resp", id = "id",
                             chains = 4, warmup = 250, iter = 750)
# very few divergent iterations.

hyper_pars <- c("alpha", "tau", "beta", "delta")
hyper_struct_pars <- c("sigma", "Omega")

# run on windows to avoid warnings:
windowsFonts(Arial=windowsFont("TT Arial"))

rstan::stan_trace(recov_pars, pars = hyper_pars)
rstan::stan_trace(recov_pars, pars = hyper_struct_pars)

print(recov_pars, pars = c(hyper_pars, hyper_struct_pars))
# Inference for Stan model: wiener_hierarch_within.
# 4 chains, each with iter=750; warmup=250; thin=1; 
# post-warmup draws per chain=500, total post-warmup draws=2000.
# 
#             mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
# alpha[1]    0.93    0.00 0.11  0.72  0.86  0.92  0.99  1.15   653 1.01
# tau[1]      0.29    0.00 0.02  0.25  0.28  0.29  0.31  0.34   589 1.00
# beta[1]     0.49    0.00 0.07  0.36  0.45  0.49  0.54  0.63   938 1.00
# delta[1]   -1.49    0.01 0.31 -2.12 -1.67 -1.48 -1.29 -0.86  1361 1.00
# sigma[1]    0.29    0.00 0.10  0.15  0.22  0.27  0.34  0.55   748 1.01
# sigma[2]    0.07    0.00 0.02  0.04  0.06  0.07  0.08  0.12   648 1.00
# sigma[3]    0.19    0.00 0.07  0.08  0.14  0.18  0.22  0.33   477 1.01
# sigma[4]    0.53    0.02 0.43  0.02  0.21  0.44  0.75  1.60   632 1.00
# Omega[1,1]  1.00    0.00 0.00  1.00  1.00  1.00  1.00  1.00  2000  NaN
# Omega[1,2] -0.31    0.01 0.30 -0.79 -0.53 -0.34 -0.11  0.34   622 1.01
# Omega[1,3] -0.27    0.01 0.32 -0.81 -0.51 -0.30 -0.05  0.41   881 1.00
# Omega[1,4]  0.05    0.01 0.41 -0.73 -0.25  0.05  0.37  0.80  1630 1.00
# Omega[2,1] -0.31    0.01 0.30 -0.79 -0.53 -0.34 -0.11  0.34   622 1.01
# Omega[2,2]  1.00    0.00 0.00  1.00  1.00  1.00  1.00  1.00  1993 1.00
# Omega[2,3] -0.29    0.01 0.29 -0.77 -0.50 -0.31 -0.08  0.33  1296 1.00
# Omega[2,4] -0.14    0.01 0.40 -0.83 -0.45 -0.17  0.14  0.68  2000 1.00
# Omega[3,1] -0.27    0.01 0.32 -0.81 -0.51 -0.30 -0.05  0.41   881 1.00
# Omega[3,2] -0.29    0.01 0.29 -0.77 -0.50 -0.31 -0.08  0.33  1296 1.00
# Omega[3,3]  1.00    0.00 0.00  1.00  1.00  1.00  1.00  1.00  1917 1.00
# Omega[3,4]  0.02    0.01 0.43 -0.76 -0.31  0.01  0.35  0.81  2000 1.00
# Omega[4,1]  0.05    0.01 0.41 -0.73 -0.25  0.05  0.37  0.80  1630 1.00
# Omega[4,2] -0.14    0.01 0.40 -0.83 -0.45 -0.17  0.14  0.68  2000 1.00
# Omega[4,3]  0.02    0.01 0.43 -0.76 -0.31  0.01  0.35  0.81  2000 1.00
# Omega[4,4]  1.00    0.00 0.00  1.00  1.00  1.00  1.00  1.00  1323 1.00
# 
# Samples were drawn using NUTS(diag_e) at Mon Mar 13 17:16:41 2017.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
# convergence, Rhat=1).

# use bayesplot for more diagnostic plots
library("bayesplot")
posterior_recov <- as.array(recov_pars)
np_recov <- nuts_params(recov_pars)
mcmc_trace(posterior_recov, regex_pars = "sigma", divergences = np_recov)
mcmc_trace(posterior_recov, regex_pars = c("alpha\\[", "beta\\[", "tau\\[", "delta\\["), divergences = np_recov)


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
```


