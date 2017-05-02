
require(RWiener) # for random generation.
N_participants <- 10
N_trials <- 20
mus <- c(alpha = 1, tau = 0.3, beta = 0.5, delta = 1.5)
vars <- c(alpha = 0.2, tau = 0.05, beta = 0.1, delta = 0.5)
cors <- matrix(0, 4, 4)
cors[upper.tri(cors)] <- c(0, 0.2, 0.3, 0.7, 0.5, 0)
cors <- cors + t(cors)
diag(cors) <- 1
Sigma <- MBESS::cor2cov(cors, vars)
set.seed(5)
pars <- MASS::mvrnorm(N_participants,mu = mus, Sigma = Sigma)

d_prep <- vector("list", N_participants)
for (i in seq_len(N_participants)) {
d_prep[[i]] <-cbind( id = factor(i), rwiener(N_trials, alpha = pars[i,"alpha"], tau = pars[i,"tau"], beta = pars[i,"beta"], delta = pars[i,"delta"]))
}
d_prep <- do.call("rbind", d_prep)

set.seed(666)
system.time(
  recov_pars <- stan_wiener_within(alpha = ~1, 
                             tau = ~1, 
                             beta = ~1, 
                             delta = ~1, 
                             data = d_prep, 
                             rt = "q", response = "resp", id = "id",
                             chains = 1, warmup = 250, iter = 750, seed=666)

)
#    user  system elapsed 
# 701.907   0.462 703.052 

#require(rstan)
#model2 <- stan_model("exec/wiener_hierarch_within2.stan")
set.seed(666)
system.time(
  recov_pars2 <- stan_wiener_within2(alpha = ~1, 
                             tau = ~1, 
                             beta = ~1, 
                             delta = ~1, 
                             data = d_prep, 
                             rt = "q", response = "resp", id = "id",
                             chains = 1, warmup = 250, iter = 750, seed=666)
)
                              chains = 2, warmup = 250, iter = 500)
#    user  system elapsed 
# 662.459   0.568 663.352 

fit_rr <- stan_wiener_within2(alpha = ~0+condition, 
                              tau = ~1, 
                              beta = ~1, 
                              delta = ~1, 
                              data = sa_red, 
                              rt = "rt", response = "response", id = "id",


data(speed_acc, package = "rtdists")
speed_acc <- droplevels(speed_acc[!speed_acc$censor,])
str(speed_acc)

## use only a subset of the data:
sa_red <- droplevels(speed_acc[speed_acc$frequency %in% c("high", "nw_high"),])
#sa_red <- droplevels(sa_red[sa_red$id %in% as.character(1:10),])
#sa_red <- droplevels(sa_red[ sample(seq_len(nrow(sa_red)), ceiling(nrow(sa_red)/20)),])


require(rstan)
data_prep_old <- create_data_wiener_within(alpha = ~0+condition, 
                             tau = ~0+condition, 
                             beta = ~0+condition, 
                             delta = ~0+condition, 
                             data = sa_red, 
                             rt = "rt", response = "response", id = "id",
                             pos_intercepts = "all")

data_prep <- create_data_wiener_within_2(alpha = ~0+condition, 
                             tau = ~0+condition, 
                             beta = ~0+condition, 
                             delta = ~0+condition, 
                             data = sa_red, 
                             rt = "rt", response = "response", id = "id",
                             pos_intercepts = "all")
tmp_pars <- wiener_within_init2(data = data_prep)

model_old <- sampling(stanmodels$wiener_hierarch_within,
                      data = data_prep_old, 
                      chains = 0)

model_new <- sampling(stanmodels$wiener_hierarch_within2,
                      data = data_prep, 
                      chains = 0)
library(microbenchmark)
microbenchmark(
log_prob(model_old, upars = unconstrain_pars(model_old, pars = tmp_pars))
,
log_prob(model_new, upars = unconstrain_pars(model_new, pars = tmp_pars))
, times = 100)
