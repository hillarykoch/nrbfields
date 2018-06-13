library(copula)

# bivariate lognormal density for independent lognormal distributions
n <- 100
x <- y <- dlnorm(seq(.001,10, length.out = n), mean = 3, sdlog = 2)
#contour(outer(x,y))
persp(outer(x,y), theta = 95, phi = 15)

# uniform data on [0,1] (quantile function applied to lognormal data)
u <- v <- seq(1/(n+1),n/(n+1),length.out = n)
unif <- outer(u,v)
contour(unif)

# Apply generator for given theta (theta is the copula dependence parameter)
generator <- function(u, theta){
    psi <- log((1-theta+theta*u)/u)
    psi[psi < 0] <- 0
    psi
}
theta <- 0.5
psiu <- psiv <- generator(u,theta)

# Apply pseudo-inverse generator
# The inv[psi < 0] part isn't necessary here I don't think because the generator is always non-negative
generator_inv <- function(psi, theta){
    inv <- (1-theta)/(exp(psi)-theta)
    inv[psi < 0] <- 0
    inv
}

psiinvu <- psiinvv <- generator_inv(psiu, theta)

# Apply Clayton copula
clayton <- function(u,v,theta, density = FALSE){
    if(!density){
        clay <- u^-theta + v^-theta-1
        clay[clay < 0] <- 0
    } else{
        clay <- (theta+1)*(u*v)^(-theta-1)*(u^-theta + v^-theta-1)^((-2*theta-1)/theta)
    }
    clay
}


theta <- -4
cop <- archmCopula("frank", param = theta)
sim <- rCopula(1000, cop)
plot(sim)

formula <- "dlnorm(seq(.001,10, length.out = n), mean = 0, sdlog = 2)"
indep <- cv <- outer(formula,
               dlnorm(seq(.001,10, length.out = n), mean = 0, sdlog = 2))
#indep <- cv <- outer(dgamma(seq(.001,10, length.out = n), shape = 3, rate = 2),
#                     dgamma(seq(.001,10, length.out = n), shape = 3, rate = 2))
# indep <- cv <- outer(dweibull(seq(.001,10, length.out = n), shape = 1.5, scale = 2),
#                      dweibull(seq(.001,10, length.out = n), shape = 1.5, scale = 2))
contour(indep)
plot(indep[1,])

u <- plnorm(seq(.001,10, length.out = n), mean = 0, sdlog = 1.5)
#u <- pgamma(seq(.001,10, length.out = n), shape = 3, rate = 2)
#u <- pweibull(seq(.001,10, length.out = n), shape = 2, scale = 3)

for(i in seq(n)){
    cv[i,] <- dCopula(cbind(u,u[i]), cop)
}

contour(indep*cv)
persp(seq(0,10,length.out = n), seq(0,10,length.out=n), indep*cv, 
      theta = 80, phi = 15,
      xlab = "", ylab = "", zlab = "")

# testfun <- function(param){
#     t1 <- param[1]
#     t2 <- param[2]
#     
#     u <- plnorm(t1)
#     v <- plnorm(t2)
#     cv <- dCopula(c(u,v),cop)
#     indep <- dlnorm(t1)*dlnorm(t2)
#     
#     indep*cv
# }
# 
# testint <- hcubature(testfun, c(0,0), c(100,100))
