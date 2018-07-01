# Example of extracting likely corner points to reduce computational burden of nrbfields
library(readr)
library(imager)
library(expm)
library(Matrix)
library(dplyr)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("~/Box Sync/School/research - Qunhua/3C/nrbfields/src/toy_rcpp.cpp")

# read in data
resolution <- 50000
rep1d <- list()
rep1d[[1]] <- read_table2("~/Box Sync/School/research - Qunhua/3C/nrbfields/data/mES_r1_chr5_50kb.txt", col_name = FALSE)

# Turn into full matrix
fulld <- list()
for(i in seq_along(rep1d)){
    b1 <- diff((range(rep1d[[i]][,1:2])/resolution)) + 1
    correction1 <- min(rep1d[[i]][,1:2])
    
    # Convert 3 column format into sparse matrix representation
    # sparseMatrix is a command found in the Matrix package
    sd1 <- sparseMatrix(i=(rep1d[[i]]$X1-correction1)/resolution,
                        j=(rep1d[[i]]$X2-correction1)/resolution,
                        x = rep1d[[i]]$X3,
                        dims = c(b1,b1),
                        index1 = FALSE)
    utri1 <- as.matrix(sd1)
    full1 <- utri1 + t(utri1)
    diag(full1) <- diag(full1)/2
    
    fulld[[i]] <- full1
}

# Smooth the data matrix
d <- fulld[[1]]
sm_d <- d %^% 3

# Edge filter with Canny algorithm
# (adjust alpha to be more or less strict in filtering)
can <- cannyEdges(as.cimg(sm_d^.2), alpha = 0.6)

# Hough transform to find likely lines
# ntheta = 5 guarantees only horizontal and vertical lines
hl <- hough_line(can, ntheta = 5, data.frame = TRUE)

# Filter out duplicate lines (going in opposing directions)
hor <- dplyr::filter(hl, theta %in% c(0, pi, 2*pi))
hor <- hor[match(unique(hor$rho), hor$rho),]
ver <- dplyr::filter(hl, theta %in% c(pi/2, 3*pi/2))
ver <- ver[match(unique(ver$rho), ver$rho),]

hl <- data.frame(rbind(hor, ver))

# Visualize smoothed data,
# filtered input data,
# and lines found on them
par(mfrow = c(2,2))

plot(as.cimg(sm_d^.2))
plot(as.cimg(sm_d^.2))
with(subset(hl,score > quantile(score,.99)), nfline(theta, rho, col="red", lwd = 0.5))
plot(can)
plot(can)
with(subset(hl, score > quantile(score,.99)), nfline(theta, rho, col="red", lwd = 0.5))


# Call corners with a naive algorithm
keepmat <- ccornercall(as.matrix(can))
keepmat_pad

# Where are the corners (only one half is computed by ccornercall)
plot(as.cimg(keepmat))

