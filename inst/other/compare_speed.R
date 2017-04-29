
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
#     user   system  elapsed 
# 1256.319    2.445 1273.813 

require(rstan)
model2 <- stan_model("exec/wiener_hierarch_within2.stan")
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

fit_rr <- stan_wiener_within2(alpha = ~0+condition, 
                              tau = ~1, 
                              beta = ~1, 
                              delta = ~1, 
                              data = sa_red, 
                              rt = "rt", response = "response", id = "id",
                              chains = 2, warmup = 250, iter = 500)

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
