# Simulate data, given a field
rfield <- function(gridlen, coefs,
                   shape=NULL, rate=NULL, scale=NULL, meanlog = NULL, sdlog = NULL,
                   nbasis, type = c("gamma", "weibull", "lognormal"), rev = FALSE,
                   copula = FALSE, copulaType = NULL, param = NULL,
                   errtype = c("gaussian"), errvar = 1){
    if(type == "gamma" & any(is.null(c(shape, rate)))){
        stop("for gamma, shape and rate must be specified.")
    }
    if(type == "weibull" & any(is.null(c(shape, scale)))){
        stop("for weibull, shape and scale must be specified.")
    }
    if(type == "lognormal" & any(is.null(c(meanlog, sdlog)))){
        stop("for lognormal, meanlog and sdlog must be specified.")
    }
    if(nbasis^2 != length(coefs)-1){
        stop("number of coefficients must equal nbasis^2 + 1.")
    }
    if(type == "gamma"){
        fieldlist <- create.gamma.field(gridlen, shape, rate, nbasis, rev=rev,
                                        copula=copula, copulaType=copulaType, param=param)
    }
    if(type == "weibull"){
        fieldlist <- create.weibull.field(gridlen, shape, scale, nbasis, rev=rev,
                                          copula=copula, copulaType=copulaType, param=param)
    }
    if(type == "lognormal"){
        fieldlist <- create.lognormal.field(gridlen, meanlog, sdlog, nbasis, rev=rev,
                                            copula=copula, copulaType=copulaType, param=param)
    }
    if(copula & any(is.null(c(copulaType, param)))){
        stop("copulaType and dependence parameter must be specified when copula = TRUE.")
    }
    
    # Unfold tensor into n^2 x nbasis^2 matrix, where n is the grid length
    unfld <- cbind(1, cunfold(simplify2array(fieldlist)))
    sclfield <- unfld %*% coefs
    refold <- matrix(sclfield, nrow = gridlen)
    if(errtype == "gaussian"){
        errs <- matrix(rnorm(gridlen^2, 0, sqrt(errvar)), nrow = gridlen)
    }
    refold + errs
}
