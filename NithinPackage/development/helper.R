#"C:\Users\think\OneDrive\Documents\NithinPackage"

Rcpp::compileAttributes('C:\\Users\\think\\OneDrive\\Documents\\NithinPackage')           # this updates the Rcpp layer from C++ to R
roxygen2::roxygenize('C:\\Users\\think\\OneDrive\\Documents\\NithinPackage', roclets="namespace", clean=TRUE)  # this updates the documentation based on roxygen comments

Rcpp::compileAttributes('~/Dropbox/Github/Warmup/NithinPackage')           # this updates the Rcpp layer from C++ to R
roxygen2::roxygenize('~/Dropbox/Github/Warmup/NithinPackage', roclets="namespace", clean=TRUE)  # this updates the documentation based on roxygen comments