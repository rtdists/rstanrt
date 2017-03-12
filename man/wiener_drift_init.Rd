\name{wiener_drift_init}
\alias{wiener_drift_init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wiener_drift_init(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (data) 
{
    N <- data$N
    nparams <- sum(data$n_par)
    tmp_out <- list(deltahat_tilde = matrix(rnorm(N * nparams, 
        0, 0.01), nparams, N), L_Omega = diag(nparams), sigma = runif(nparams))
    tmp_out$alpha <- rep(0, data$n_par[1])
    tmp_out$alpha[data$pos_intercept[1, seq_len(data$n_intercept[1])]] <- runif(data$n_intercept[1], 
        0.5, 1.5)
    tmp_out$alpha[data$pos_other_pars[1, seq_len(data$n_other_pars[1])]] <- rnorm(data$n_other_pars[1], 
        0, 0.1)
    tmp_out$tau <- rep(0, data$n_par[2])
    tmp_out$tau[data$pos_intercept[2, seq_len(data$n_intercept[2])]] <- runif(data$n_intercept[2], 
        0.2, 0.25)
    tmp_out$tau[data$pos_other_pars[2, seq_len(data$n_other_pars[2])]] <- rnorm(data$n_other_pars[2], 
        0, 0.01)
    tmp_out$beta <- rep(0, data$n_par[3])
    tmp_out$beta[data$pos_intercept[3, seq_len(data$n_intercept[3])]] <- runif(data$n_intercept[3], 
        0.4, 0.6)
    tmp_out$beta[data$pos_other_pars[3, seq_len(data$n_other_pars[3])]] <- rnorm(data$n_other_pars[3], 
        0, 0.01)
    tmp_out$delta <- rep(0, data$n_par[4])
    tmp_out$delta[data$pos_intercept[4, seq_len(data$n_intercept[4])]] <- runif(data$n_intercept[4], 
        1, 2)
    tmp_out$delta[data$pos_other_pars[4, seq_len(data$n_other_pars[4])]] <- rnorm(data$n_other_pars[4], 
        0, 1)
    tmp_out$drift_crit <- rep(0, data$n_par[5])
    tmp_out$drift_crit[data$pos_intercept[5, seq_len(data$n_intercept[5])]] <- rnorm(data$n_intercept[5], 
        0, 1)
    tmp_out$drift_crit[data$pos_other_pars[5, seq_len(data$n_other_pars[5])]] <- rnorm(data$n_other_pars[5], 
        0, 0.1)
    return(tmp_out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
