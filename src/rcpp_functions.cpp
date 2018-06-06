// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
arma::mat cunfold(arma::cube Q) {
    int n = Q.slice(1).n_rows;
    int d = Q.slice(1).n_cols;
    int p = Q.n_slices;
    arma::mat m(n*d, p, arma::fill::none);
    m = reshape(Q,n*d,p,1);
    return m;
}