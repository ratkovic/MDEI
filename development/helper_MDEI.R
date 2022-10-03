Rcpp::compileAttributes('~/Dropbox/Github/MDEI')           # this updates the Rcpp layer from C++ to R
roxygen2::roxygenize('~/Dropbox/Github/MDEI', roclets="namespace", clean=TRUE)  # this updates the documentation based on roxygen comments
roxygen2::roxygenize('~/Dropbox/Github/MDEI')  # this updates the documentation based on roxygen comments

Rcpp::compileAttributes('~/Downloads/MDEI')           # this updates the Rcpp layer from C++ to R
roxygen2::roxygenize('~/Downloads/MDEI', roclets="namespace", clean=TRUE)  # this updates the documentation based on roxygen comments
roxygen2::roxygenize('~/Downloads/MDEI')  # this updates the documentation based on roxygen comments

