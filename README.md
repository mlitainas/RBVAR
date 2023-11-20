
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RBVAR

<!-- badges: start -->
<!-- badges: end -->

RBVAR is a set of functions for the estimation of BVAR models.

## Installation

You can install the development version of RBVAR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mlitainas/RBVAR")
```
In the following example I estimate a SVAR using the Inverse-Wishart prior implemented using dummy observations ala Banbura 2007. The identification is recursive. IRFs show the shock for the first variable. 
``` r
data %>% 
  BVAR_estimation_NIW(lags = 4,reps = 2000,burn = 1000, lamda = 1, tau = 4, epsilon = 0.1 ) %>%  
  BVAR_irf_chol(shock = 1)


```