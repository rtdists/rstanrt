#' Analysis Functions for Stan Models
#' 
#' @example examples/examples_analysis.R
#' 
#' @rdname analysis_functions
#' @export
extract_tbl <- function(model, pars, levels) {
  #browser()
  if (!requireNamespace("dplyr")) stop("dplyr required for exporting a tbble.")
  parl <- rstan::extract(model, pars = pars)
  #str(parl)
  # dim1 <- dim(parl[[1]])
  # for (i in seq_along(parl)) {
  #   if (!all(dim(parl[[i]]) == dim1)) stop("dimensions of all parameters need to be identical!")
  #   if (!missing(dimensions)) {
  #     names(dimnames(parl[[i]]))[-1] <- names(dimensions)
  #     dimnames(parl[[i]])[-1] <- dimensions 
  #   }
  # }
  parl2 <- lapply(parl, reshape2::melt)
  for (i in seq_along(parl)) {
    parl2[[i]]$parameter <- factor(names(parl2)[i])
  }
  par_df <- dplyr::tbl_df(do.call("rbind", parl2))
  colnames(par_df)[2] <- "condition"
  par_df$condition <- factor(par_df$condition, labels = levels)
  par_df
}

#' @rdname analysis_functions
#' @export
bayesian_panelfunction <- function(x, y, ..., box.ratio, level = 0.05, abline) {
  if (!requireNamespace("multcompView")) stop("package multcompView required.")
  dots <- list(...)
  if (is.numeric(x)) {
    ad <- tapply(x, INDEX = list(y), FUN = get_tripple, interval = (1-level)*100)
    lattice::panel.violin(x, y, col = "transparent", lty = 2, border = "grey",
                          varwidth = FALSE, box.ratio = box.ratio)
    for (i in seq_len(nrow(ad))) {
      lattice::panel.arrows(x0 = ad[[i]][2], x1 = ad[[i]][3], y0 = i, y1=i, ends="both", angle = 90, length = 0.1)
      lattice::panel.points(x = ad[[i]][1], y = i, pch =19, cex = 2.5)
    }
    
    tmp_m <- tapply(x, INDEX = list(y), FUN = identity, simplify = FALSE)
    attributes(tmp_m) <- NULL
    names(tmp_m) <- levels(y)
    tmp_m <- as.matrix(as.data.frame(tmp_m))
    p_matrix <- get_p_matrix(tmp_m)
    lp_matrix <- (p_matrix < (level/2) | p_matrix > (1-(level/2)))
    cld <- multcompView::multcompLetters(lp_matrix)$Letters
    limits <- lattice::current.panel.limits()
    lattice::panel.text(x = limits$xlim[2] - diff(limits$xlim)/10, 
                        y = seq_along(levels(y))+0.1, 
                        labels = cld)
  }
  else {
    ad <- tapply(y, INDEX = list(x), FUN = get_tripple, interval = (1-level)*100)
    lattice::panel.violin(x, y, col = "transparent", lty = 2, border = "grey",
                          varwidth = FALSE, box.ratio = box.ratio, horizontal = FALSE)
    for (i in seq_len(nrow(ad))) {
      lattice::panel.arrows(y0 = ad[[i]][2], y1 = ad[[i]][3], x0 = i, x1=i, ends="both", angle = 90, length = 0.1)
      lattice::panel.points(y = ad[[i]][1], x = i, pch =19, cex = 2.5)
    }
    
    tmp_m <- tapply(y, INDEX = list(x), FUN = identity, simplify = FALSE)
    attributes(tmp_m) <- NULL
    names(tmp_m) <- levels(y)
    tmp_m <- as.matrix(as.data.frame(tmp_m))
    p_matrix <- get_p_matrix(tmp_m)
    lp_matrix <- (p_matrix < (level/2) | p_matrix > (1-(level/2)))
    cld <- multcompView::multcompLetters(lp_matrix)$Letters
    limits <- lattice::current.panel.limits()
    lattice::panel.text(y = limits$ylim[2] - diff(limits$ylim)/10, 
                        x = seq_along(levels(x))+0.1, 
                        labels = cld)
  }
  
 
}

get_tripple <- function(x, interval = 95) {
  if (!requireNamespace("hdrcde")) stop("package hdrcde required.")
  #c(mean(x), TeachingDemos::emp.hpd(x, conf = interval))
  tmp <- hdrcde::hdr(x, prob = interval)
  c(tmp$mode, tmp$hdr)
}

#' @rdname analysis_functions
#' @export
print_signif_par <- function(model, pars = c("Omega"), level = 0.05, use_knitr = TRUE) {
  if (!requireNamespace("stringr")) stop("package stringr required.")
  if (use_knitr) 
    if (!requireNamespace("knitr"))  {
      use_knitr <- FALSE
      warning("knitr not available for nice printing.")
    }
  for (name in pars) {
    corrs1 <- summary(model, pars = name, probs = c(level/2, 1-level/2))
    corrs1 <- corrs1$summary[corrs1$summary[,4] > 0 | corrs1$summary[,5] < 0,]

    ns <- stringr::str_extract_all(rownames(corrs1), "\\d+")
    if (all(vapply(ns, length, 1) == 2)) 
      corrs1 <- round(corrs1[!vapply(ns, function(x) as.numeric(x[1]) >= as.numeric(x[2]), NA),, drop = FALSE], 2)
    if (use_knitr) print(knitr::kable(corrs1, caption = name))
    else print(corrs1)
  }
}

get_p_matrix <- function(matrix, only_low = TRUE) {
  out <- matrix(-1, nrow = ncol(matrix), ncol = ncol(matrix), dimnames = list(colnames(matrix), colnames(matrix)))
  #str(df)
  for (i in seq_len(ncol(matrix))) {
    for (j in seq_len(ncol(matrix))) {
      out[i, j] <- mean(matrix[,i] < matrix[,j])
    }
  }
  if (only_low) out[out > .5] <- 1- out[out > .5] 
  out
}