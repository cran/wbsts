# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cusum <- function(x) {
    .Call(`_wbsts_cusum`, x)
}

finner_prod_maxp <- function(x, p) {
    .Call(`_wbsts_finner_prod_maxp`, x, p)
}

across_fip <- function(X, tau, p, epp, p1, Ts) {
    .Call(`_wbsts_across_fip`, X, tau, p, epp, p1, Ts)
}

multi_across_fip <- function(X, M, min_draw, tau, p, epp, Ts) {
    .Call(`_wbsts_multi_across_fip`, X, M, min_draw, tau, p, epp, Ts)
}

