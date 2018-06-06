# Unfold tensor into n^2 x nbasis^2 matrix, where n is the grid length
mode3unfold <- function(fieldlist){
    cunfold(simplify2array(fieldlist))
}
