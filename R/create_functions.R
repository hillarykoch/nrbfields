create.weibull.basis <- function(gridlen, shape, scale, nbasis, rev = FALSE){
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

create.gamma.basis <- function(gridlen, shape, rate, nbasis, rev = FALSE){
    if(shape <= 1){
        warning("Gamma distribution with shape <= 1 will result in very spiky basis.")
    }
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

create.weibull.field <- function(gridlen, shape, scale, nbasis, rev=FALSE, reduce = FALSE){
    S <- create.weibull.basis(gridlen, shape, scale, nbasis, rev=rev)
    combos <- combn(nbasis,2)
    combos <- cbind(combos, rbind(combos[2,], combos[1,]), rbind(1:nbasis, 1:nbasis))
    
    # outer product of independent bases
    ops <- lapply(seq(ncol(combos)), function(X) outer(S[,combos[1,X]], S[,combos[2,X]]))
    
    if(reduce){
        ops <- Reduce('+', ops)
    }
    ops
}

create.gamma.field <- function(gridlen, shape, rate, nbasis, rev = FALSE, reduce = FALSE){
    S <- create.gamma.basis(gridlen, shape, rate, nbasis, rev=rev)
    combos <- combn(nbasis,2)
    combos <- cbind(combos, rbind(combos[2,], combos[1,]), rbind(1:nbasis, 1:nbasis))
    
    # outer product of independent bases
    ops <- lapply(seq(ncol(combos)), function(X) outer(S[,combos[1,X]], S[,combos[2,X]]))
    
    if(reduce){
        ops <- Reduce('+', ops)
    }
    ops
}
