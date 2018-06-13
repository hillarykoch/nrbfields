create.weibull.basis <- function(gridlen, shape, scale, nbasis,
                                 rev = FALSE, copula = FALSE,
                                 copulaType = NULL, param = NULL){
    
    if(copula){
        cop <- archmCopula(copulaType, param = param)
        if(rev){
            indep <- outer(dweibull(seq(gridlen), shape, scale) %>% rev,
                           dweibull(seq(gridlen), shape, scale) %>% rev) # If x and y were treated as independent
            u <- pweibull(seq(gridlen), shape, scale) %>% rev # Compute distribution of data to get uniform RVs
        } else {
            indep <- outer(dweibull(seq(gridlen), shape, scale),
                           dweibull(seq(gridlen), shape, scale)) # If x and y were treated as independent
            u <- pweibull(seq(gridlen), shape, scale) # Compute distribution of data to get uniform RVs   
        }
        dep <- matrix(rep(NA, gridlen^2), nrow = gridlen) # Compute the copula
        for(i in seq(gridlen)){
            dep[i,] <- dCopula(cbind(u, u[i]), cop)
        }
        dep*indep # multiply independent density by dependent part
    } else{
        d <- dweibull(seq(gridlen), shape, scale)
        l <- round(gridlen/nbasis) # Compute how to stagger the basis
        S <- matrix(rep(0, nbasis*gridlen), ncol = nbasis)
        for(i in (seq(nbasis)-1)){
            S[((l*i)+1):gridlen,i+1] <- d[1:(length((l*i):gridlen)-1)]
        }
        if(rev){
            S <- apply(S, 2, rev)
        }
        S
    }
}

create.gamma.basis <- function(gridlen, shape, rate, nbasis,
                               rev = FALSE, copula = FALSE,
                               copulaType = NULL, param = NULL){
    if(shape <= 1){
        warning("Gamma distribution with shape <= 1 will result in very spiky basis.")
    }
    
    if(copula){
        cop <- archmCopula(copulaType, param = param)
        if(rev){
            indep <- outer(dgamma(seq(gridlen), shape=shape, rate=rate) %>% rev,
                           dgamma(seq(gridlen), shape=shape, rate=rate) %>% rev)
            u <- pgamma(seq(gridlen), shape=shape, rate=rate) %>% rev
        } else {
            indep <- outer(dgamma(seq(gridlen), shape=shape, rate=rate),
                           dgamma(seq(gridlen), shape=shape, rate=rate))
            u <- pgamma(seq(gridlen), shape=shape, rate=rate)
        }
        dep <- matrix(rep(NA, gridlen^2), nrow = gridlen)
        for(i in seq(gridlen)){
            dep[i,] <- dCopula(cbind(u, u[i]), cop)
        }
        dep*indep
    } else{
        d <- dgamma(seq(gridlen), shape=shape, rate=rate)
        l <- round(gridlen/nbasis)
        S <- matrix(rep(0, nbasis*gridlen), ncol = nbasis)
        for(i in (seq(nbasis)-1)){
            S[((l*i)+1):gridlen,i+1] <- d[1:(length((l*i):gridlen)-1)]
        }
        if(rev){
            S <- apply(S, 2, rev)
        }
        S   
    }
}

create.lognormal.basis <- function(gridlen, meanlog, sdlog, nbasis,
                                   rev = FALSE, copula = FALSE,
                                   copulaType = NULL, param = NULL){
    if(copula){
        cop <- archmCopula(copulaType, param = param)
        if(rev){
            indep <- outer(dlnorm(seq(gridlen), meanlog, sdlog) %>% rev,
                           dlnorm(seq(gridlen), meanlog, sdlog) %>% rev)
            u <- plnorm(seq(gridlen), meanlog, sdlog) %>% rev
        } else {
            indep <- outer(dlnorm(seq(gridlen), meanlog, sdlog),
                           dlnorm(seq(gridlen), meanlog, sdlog))
            u <- plnorm(seq(gridlen), meanlog, sdlog)
        }
        dep <- matrix(rep(NA, gridlen^2), nrow = gridlen)
        for(i in seq(gridlen)){
            dep[i,] <- dCopula(cbind(u, u[i]), cop)
        }
        dep*indep
    } else{
        d <- dlnorm(seq(gridlen), meanlog=meanlog, sdlog=sdlog)
        l <- round(gridlen/nbasis)
        S <- matrix(rep(0, nbasis*gridlen), ncol = nbasis)
        for(i in (seq(nbasis)-1)){
            S[((l*i)+1):gridlen,i+1] <- d[1:(length((l*i):gridlen)-1)]
        }
        if(rev){
            S <- apply(S, 2, rev)
        }
        S   
    }
}

create.weibull.field <- function(gridlen, shape, scale, nbasis,
                                 rev=FALSE, reduce = FALSE, copula = FALSE,
                                 copulaType = NULL, param = NULL){
    if(copula & is.null(copulaType)){
        stop("A copula type and copula dependence parameter must be specified.")
    }
    
    if(copula){
        l <- round(gridlen/nbasis) # Compute how to stagger the basis
        ops <- list()
        b <- create.weibull.basis(gridlen, shape, scale, nbasis, rev=rev,
                                      copula=TRUE, copulaType=copulaType,
                                      param=param)
        count <- 1
        if(rev){
            for(i in rev(seq(nbasis)-1)){
               for(j in rev(seq(nbasis)-1)){
                  ops[[count]] <- matrix(rep(0, gridlen^2), nrow = gridlen)
                  subb <- b[gridlen:(l*i+1), gridlen:(l*j+1)]
                  ops[[count]][(gridlen-l*i):1, (gridlen-l*j):1] <- subb
                  count <- count+1
               } 
            }
        } else{
            for(i in (seq(nbasis)-1)){
                for(j in (seq(nbasis)-1)){
                    ops[[count]] <- matrix(rep(0, gridlen^2), nrow = gridlen)
                    subb <- b[1:(gridlen-l*i), 1:(gridlen-l*j)]
                    ops[[count]][((l*i)+1):gridlen, ((l*j)+1):gridlen] <- subb
                    count <- count+1
                }
            }    
        }
        
    } else{
        S <- create.weibull.basis(gridlen, shape, scale, nbasis, rev=rev)
        combos <- combn(nbasis,2)
        combos <- cbind(combos, rbind(combos[2,], combos[1,]), rbind(1:nbasis, 1:nbasis))
        
        # outer product of independent bases
        ops <- lapply(seq(ncol(combos)), function(X) outer(S[,combos[1,X]], S[,combos[2,X]]))
    }
    
    if(reduce){
        ops <- Reduce('+', ops)
    }
    ops   
}

create.gamma.field <- function(gridlen, shape, rate, nbasis,
                               rev = FALSE, reduce = FALSE, copula = FALSE,
                               copulaType = NULL, param = NULL){
    if(copula & is.null(copulaType)){
        stop("A copula type and copula dependence parameter must be specified.")
    }
    
    if(copula){
        l <- round(gridlen/nbasis) # Compute how to stagger the basis
        ops <- list()
        b <- create.gamma.basis(gridlen, shape, rate, nbasis, rev=rev,
                                  copula=TRUE, copulaType=copulaType,
                                  param=param)
        count <- 1
        if(rev){
            for(i in rev(seq(nbasis)-1)){
                for(j in rev(seq(nbasis)-1)){
                    ops[[count]] <- matrix(rep(0, gridlen^2), nrow = gridlen)
                    subb <- b[gridlen:(l*i+1), gridlen:(l*j+1)]
                    ops[[count]][(gridlen-l*i):1, (gridlen-l*j):1] <- subb
                    count <- count+1
                } 
            }
        } else{
            for(i in (seq(nbasis)-1)){
                for(j in (seq(nbasis)-1)){
                    ops[[count]] <- matrix(rep(0, gridlen^2), nrow = gridlen)
                    subb <- b[1:(gridlen-l*i), 1:(gridlen-l*j)]
                    ops[[count]][((l*i)+1):gridlen, ((l*j)+1):gridlen] <- subb
                    count <- count+1
                }
            }    
        }
        
    } else{
        S <- create.gamma.basis(gridlen, shape, rate, nbasis, rev=rev)
        combos <- combn(nbasis,2)
        combos <- cbind(combos, rbind(combos[2,], combos[1,]), rbind(1:nbasis, 1:nbasis))
        
        # outer product of independent bases
        ops <- lapply(seq(ncol(combos)), function(X) outer(S[,combos[1,X]], S[,combos[2,X]]))   
    }
    
    if(reduce){
        ops <- Reduce('+', ops)
    }
    ops
}

create.lognormal.field <- function(gridlen, meanlog, sdlog, nbasis,
                                   rev = FALSE, reduce = FALSE, copula = FALSE,
                                   copulaType = NULL, param = NULL){
    if(copula & is.null(copulaType)){
        stop("A copula type and copula dependence parameter must be specified.")
    }
    
    if(copula){
        l <- round(gridlen/nbasis) # Compute how to stagger the basis
        ops <- list()
        b <- create.lognormal.basis(gridlen, meanlog, sdlog, nbasis, rev=rev,
                                copula=TRUE, copulaType=copulaType,
                                param=param)
        count <- 1
        if(rev){
            for(i in rev(seq(nbasis)-1)){
                for(j in rev(seq(nbasis)-1)){
                    ops[[count]] <- matrix(rep(0, gridlen^2), nrow = gridlen)
                    subb <- b[gridlen:(l*i+1), gridlen:(l*j+1)]
                    ops[[count]][(gridlen-l*i):1, (gridlen-l*j):1] <- subb
                    count <- count+1
                } 
            }
        } else{
            for(i in (seq(nbasis)-1)){
                for(j in (seq(nbasis)-1)){
                    ops[[count]] <- matrix(rep(0, gridlen^2), nrow = gridlen)
                    subb <- b[1:(gridlen-l*i), 1:(gridlen-l*j)]
                    ops[[count]][((l*i)+1):gridlen, ((l*j)+1):gridlen] <- subb
                    count <- count+1
                }
            }    
        }
    } else{
        S <- create.lognormal.basis(gridlen, meanlog, sdlog, nbasis, rev=rev)
        combos <- combn(nbasis,2)
        combos <- cbind(combos, rbind(combos[2,], combos[1,]), rbind(1:nbasis, 1:nbasis))
        
        # outer product of independent bases
        ops <- lapply(seq(ncol(combos)), function(X) outer(S[,combos[1,X]], S[,combos[2,X]]))   
    }
    if(reduce){
        ops <- Reduce('+', ops)
    }
    ops
}

