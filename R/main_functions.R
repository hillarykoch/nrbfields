# Simulate data, given a field
rfield <- function(gridlen, coefs, shape=NULL, rate=NULL, scale=NULL,
                   nbasis, type = c("gamma", "weibull"), rev = FALSE,
                   errtype = c("gaussian"), errvar = 1){
    if(type == "gamma" & any(is.null(c(shape, rate)))){
        stop("for gamma, shape and rate must be specified.")
    }
    if(type == "weibull" & any(is.null(c(shape, scale)))){
        stop("for weibull, shape and scale must be specified.")
    }
    if(nbasis^2 != length(coefs)-1){
        stop("number of coefficients must equal nbasis^2 + 1.")
    }
    if(type == "gamma"){
        fieldlist <- create.gamma.field(gridlen, shape, rate, nbasis, rev=rev)
    } else if(type == "weibull"){
        fieldlist <- create.weibull.field(gridlen, shape, scale, nbasis, rev=rev)
    }
    
    unfld <- cbind(1, mode3unfold(fieldlist))
    sclfield <- unfld %*% coefs
    refold <- matrix(sclfield, nrow = gridlen)
    if(errtype == "gaussian"){
        errs <- matrix(rnorm(gridlen^2, 0, sqrt(errvar)), nrow = gridlen)
    }
    refold + errs
}