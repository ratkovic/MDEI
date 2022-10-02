<!-- badges: start -->
[![R build status](https://github.com/ratkovic/MDEI/workflows/R-CMD-check/badge.svg)](https://github.com/ratkovic/MDEI/actions)
<!-- badges: end -->

# MDEI
This R package implements the Method of Direct Information and Inference from the manuscript ``Estimation and Inference on Nonlinear and Heterogeneous Effects'' by Ratkovic and Tingley (2023, accepted at the Journal of Politics). The method calculates takes an outcome, variable of theoretical interest (treatment), and set of covariates and returns a partial derivative (marginal effect) of the treatment variable at each point along with uncertainty intervals.

The approach offers two advances. First, a split-sample approach is used as a guard against overfitting. Second, the method uses a data-driven interval derived from conformal inference, rather than relying on a normality assumption on the error terms.

A help file illustrates the method's use and the associated manuscript describes the approach and how to apply it in real world data.

# Installation #

```r
# To Install
devtools::install_github("ratkovic/MDEI")
```
