test_unfold <- function(){
    fieldlist <- create.gamma.field(gridlen=100, shape=3, rate=0.5, nbasis=10, rev = FALSE, reduce = FALSE)
    checkEqualsNumeric(length(fieldlist), 10^2)
    unfld <- mode3unfold(fieldlist)
    checkEqualsNumeric(ncol(unfld), 100)
    checkEqualsNumeric(nrow(unfld), 100*10^2)
}
