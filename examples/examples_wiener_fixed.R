

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
fit_rr_fixed <- stan_wiener_fixed(alpha = ~0+condition, 
                             tau = ~1, 
                             beta = ~1, 
                             delta = ~1, 
                             data = sa_red, 
                             rt = "rt", response = "response", id = "id",
                             chains = 1, warmup = 250, iter = 500)
fit_rr_fixed
}