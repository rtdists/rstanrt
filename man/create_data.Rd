\name{create_data}
\alias{create_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create_data(alpha, tau, beta, delta, drift_crit, data, rt, response, boundary, id, pos_intercepts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{drift_crit}{
%%     ~~Describe \code{drift_crit} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{rt}{
%%     ~~Describe \code{rt} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
  \item{boundary}{
%%     ~~Describe \code{boundary} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{pos_intercepts}{
%%     ~~Describe \code{pos_intercepts} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (alpha, tau, beta, delta, drift_crit, data, rt, response, 
    boundary, id, pos_intercepts) 
{
    N <- length(unique(data[, id]))
    data_ordered <- data[order(data[, id]), ]
    d_split <- split(data_ordered, f = data_ordered[, response])
    mmatrix_alpha <- split(as.data.frame(model.matrix(alpha, 
        data = data_ordered)), f = data_ordered[, response])
    mmatrix_tau <- split(as.data.frame(model.matrix(tau, data = data_ordered)), 
        f = data_ordered[, response])
    mmatrix_beta <- split(as.data.frame(model.matrix(beta, data = data_ordered)), 
        f = data_ordered[, response])
    mmatrix_delta <- split(as.data.frame(model.matrix(delta, 
        data = data_ordered)), f = data_ordered[, response])
    mmatrix_drift_crit <- split(as.data.frame(model.matrix(drift_crit, 
        data = data_ordered)), f = data_ordered[, response])
    K_u <- nrow(d_split[[2]])
    K_l <- nrow(d_split[[1]])
    rle_p_upper <- cumsum(rle(d_split[[2]][, id])$lengths)
    rle_p_lower <- cumsum(rle(d_split[[1]][, id])$lengths)
    if (missing(pos_intercepts)) {
        pos_intercepts <- as.list(rep(1, 5))
    }
    n_par <- vapply(list(mmatrix_alpha, mmatrix_tau, mmatrix_beta, 
        mmatrix_delta, mmatrix_drift_crit), function(x) ncol(x[[1]]), 
        1)
    n_intercept <- vapply(pos_intercepts, length, 1)
    n_other_pars <- n_par - n_intercept
    pos_intercept <- matrix(0, nrow = 5, ncol = max(n_intercept))
    pos_other_pars <- matrix(0, nrow = 5, ncol = max(n_other_pars))
    for (i in seq_along(n_par)) {
        pos_intercept[i, 1:length(pos_intercepts[[i]])] <- pos_intercepts[[i]]
        pos_other_pars[i, seq_len(n_par[i] - length(pos_intercepts[[i]]))] <- seq_len(n_par[i])[-pos_intercept[i, 
            1:length(pos_intercepts[[i]])]]
    }
    list(N = N, K_u = K_u, K_l = K_l, Y_u = d_split[[2]][, rt], 
        Y_l = d_split[[1]][, rt], p_boundaries_u = cbind(c(1, 
            rle_p_upper[-length(rle_p_upper)] + 1), rle_p_upper), 
        p_boundaries_l = cbind(c(1, rle_p_lower[-length(rle_p_lower)] + 
            1), rle_p_lower), mmatrix_alpha_u = as.matrix(mmatrix_alpha[[2]]), 
        mmatrix_alpha_l = as.matrix(mmatrix_alpha[[1]]), mmatrix_tau_u = as.matrix(mmatrix_tau[[2]]), 
        mmatrix_tau_l = as.matrix(mmatrix_tau[[1]]), mmatrix_beta_u = as.matrix(mmatrix_beta[[2]]), 
        mmatrix_beta_l = as.matrix(mmatrix_beta[[1]]), mmatrix_delta_u = as.matrix(mmatrix_delta[[2]]), 
        mmatrix_delta_l = as.matrix(mmatrix_delta[[1]]), mmatrix_drift_crit_u = as.matrix(mmatrix_drift_crit[[2]]), 
        mmatrix_drift_crit_l = as.matrix(mmatrix_drift_crit[[1]]), 
        boundary_u = (as.numeric(factor(d_split[[2]][, boundary])) - 
            1.5) * 2, boundary_l = (as.numeric(factor(d_split[[1]][, 
            boundary])) - 1.5) * 2, n_par = n_par, n_intercept = n_intercept, 
        n_other_pars = n_other_pars, pos_intercept = pos_intercept, 
        pos_other_pars = pos_other_pars)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
