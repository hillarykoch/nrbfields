library(glmnet)

nb <- 50

# Simulate weighted basis functions in lognormal copula field
coefs <- VGAM::rgpd(nb^2+1, location = 0, scale = 5, shape = 0.65)
coefs[sample(1:(nb^2+1), round(.7*(nb^2+1)), replace = F)] <- 0
sim <- rfield(gridlen = 100,
              coefs = coefs,
              meanlog = .5, sdlog = 3,
              nbasis = nb,
              type = "lognormal",
              rev = TRUE,
              copula = TRUE, copulaType = "frank", param = -4,
              errvar = 0.5)
sim <- abs(sim) ## WHY IS SIM EVER NEGATIVE??

# Create unweighted basis representation for the field
# Unfold tensor into n^2 x nbasis^2 matrix, where n is the grid length
b <- create.lognormal.field(gridlen = 100,
                            meanlog = .5, sdlog = 2,
                            nbasis = nb,
                            rev=TRUE, reduce=FALSE,
                            copula=TRUE, copulaType = "frank", param = -4)
unfld <- cbind(1, cunfold(simplify2array(b)))
fit <- glmnet(x = unfld, y = as.vector(sim),
              family = "poisson",
              alpha = 1,
              lambda = 10)
