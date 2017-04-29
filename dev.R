require(devtools)
load_all()
load_all(recompile = TRUE, reset = TRUE)

#require(roxygen2)
roxygen2::roxygenize(clean = TRUE)
