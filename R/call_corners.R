call_corners <- function(x, alpha, pow=1, t=1, sigma=2){
    # Smooth the data matrix
    sm_x <- (x + abs(min(x))) %^% t
    
    # Edge filter with Canny algorithm
    # (adjust alpha to be more or less strict in filtering)
    can <- cannyEdges(as.cimg(sm_x^pow), 
                      alpha = alpha,
                      sigma = sigma)
    
    # Call corners with a naive algorithm
    ccornercall(as.matrix(can))
}
