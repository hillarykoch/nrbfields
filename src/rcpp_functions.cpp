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


// [[Rcpp::export]]
arma::cube copula_field(int gridlen, int nbasis, int l, arma::mat b, bool rev){
    arma::cube ops = arma::zeros<arma::cube>(gridlen, gridlen, pow(nbasis,2));
    arma::mat subb(gridlen, gridlen, arma::fill::none);
    int count = 0;
    if (rev){
        for(int i = 0; i < nbasis; ++i){
            for(int j = 0; j < nbasis; ++j){
                subb.set_size(gridlen-l*i, gridlen-l*j);
                subb = b.submat(0, 0, gridlen-l*i-1, gridlen-l*j-1);
                ops.slice(count).submat(l*i, l*j, gridlen-1, gridlen-1) = subb;
                // to rotate, transpose, reverse cols, transpose, reverse cols
                ops.slice(count) = reverse(reverse(ops.slice(count).t(),0).t(),0);
                ++count;
            }
        }
    } else{
        for(int i = 0; i < nbasis; ++i){
            for(int j = 0; j < nbasis; ++j){
                subb.set_size(gridlen-l*i, gridlen-l*j);
                subb = b.submat(0, 0, gridlen-l*i-1, gridlen-l*j-1);
                ops.slice(count).submat(l*i, l*j, gridlen-1, gridlen-1) = subb;            
                ++count;
            }
        }
    }
    return ops;
}