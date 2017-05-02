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

fit_rr <- stan_wiener_within_fixed_beta(alpha = ~0+condition, 
                             tau = ~1, 
                             delta = ~1, 
                             data = sa_red, 
                             rt = "rt", response = "response", id = "id",
                             chains = 1, warmup = 250, iter = 500)

pars_alphs <- extract_tbl(fit_rr, "alpha", levels = levels(sa_red$condition))

require(latticeExtra)
lattice.options(default.theme = standard.theme(color = FALSE))
lattice.options(default.args = list(as.table = TRUE))

bwplot(condition ~ value, pars_alphs, panel = bayesian_panelfunction, level = 0.05)
bwplot(value~condition, pars_alphs, panel = bayesian_panelfunction, level = 0.05)


print_signif_par(fit_rr) # no correlation ssignificant

print_signif_par(fit_rr, pars = "sigma") # shows other parameters that are != 0.
}
