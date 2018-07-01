// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// find corners
// [[Rcpp::export]]
arma::mat ccornercall(arma::mat cannymat){
    int x = cannymat.n_rows;
    int y = cannymat.n_cols;
    arma::mat keepmat(x, y, arma::fill::zeros);
    for(int i = 1; i < x - 1; ++i){
        for(int j = i+1; j < y - 1; ++j){
            if(cannymat(i,j) + cannymat(i,j-1) + cannymat(i+1,j) == 3){
                keepmat(i,j) = 1;
            }
        }
    }
    return keepmat;
}
