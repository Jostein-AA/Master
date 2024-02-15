rm(list = ls())

library(bigDM)
library(MASS)
library(RColorBrewer)
library(sf)
library(splines)
library(tmap)

data("Carto_SpainMUN")
carto <- Carto_SpainMUN

############################################################
## Simulate risks from a spatial two-dimensional P-spline ##
############################################################

#Function that takes two kronecker products and ??? 
Rten <- function(X1, X2) {
  #What in the hell?
  one1 <- matrix(1, 1, ncol(X1))      # Two vectors containing only ones!
  one2 <- matrix(1, 1, ncol(X2))
  kronecker(X1, one2) * kronecker(one1, X2)   #???
}

#What is bdeg? Guess it is B-spline degree (ie. cubic)
bdeg <- 3

#Get coordinates of each centriod
locs <- st_coordinates(st_centroid(carto))

#Make the coordinate values between 0 and 1!
locs <- apply(locs, 2, function(x) (x - min(x)) / (max(x) - min(x)))

#Extract the coordinates
x1 <- locs[, 1]
x2 <- locs[, 2]

#Specify number of knots in x1-direction
ndx1 <- 40

#Find 40 equi-spaced knots in x1-direction
dis1 <- (max(x1) - min(x1)) / ndx1
x1l <- min(x1) - dis1 * 0.05
x1r <- max(x1) + dis1 * 0.05
dx1 <- (x1r - x1l) / ndx1
knots1 <- seq(x1l - bdeg * dx1, x1r + bdeg * dx1, by = dx1)

#Specify number of knots in x2-direction
ndx2 <- 40

#Find 40 equi-spaced knots in x2-direction
dis2 <- (max(x2) - min(x2)) / ndx2
x2l <- min(x2) - dis2 * 0.05
x2r <- max(x2) + dis2 * 0.05
dx2 <- (x2r - x2l) / ndx2
knots2 <- seq(x2l - bdeg * dx2, x2r + bdeg * dx2, by = dx2)

#Create design matrices for marginal B-splines
B1 <- spline.des(knots1, x1, bdeg + 1, 0 * x1)$design  ##Base marginal longitude
B2 <- spline.des(knots2, x2, bdeg + 1, 0 * x2)$design  ##Base marginal latitude

#Take some form of Kronecker product involving both B2 and B1
Bs <- Rten(B2, B1)

#what are these??
k1 <- dim(B1)[2]
k2 <- dim(B2)[2]
k <- dim(Bs)[2]

# Precision matrix 
D1 <- diff(diag(k1), differences = 2)
P1 <- t(D1) %*% D1 #RW2 precision matrix

D2 <- diff(diag(k2), differences = 2)
P2 <- t(D2) %*% D2 #RW2 precision matrix

#Identity matrices k1 x k1 and k2 x k2
I1 <- diag(k1)
I2 <- diag(k2)

#Combined penalty matrix?
Ps <- kronecker(P2, I1) + kronecker(I2, P1)

## Simulate risks from a multivariate normal distribution ##
set.seed(19890502)

alpha <- 0             #alpha is intercept
sig2 <- 0.03           #sig2 is variance of parameters
#Sigma is covariance matrix (take generalized inverse of penalty matrix)
Sigma <- sig2*ginv(Ps)  #scale by sig2

#Draw risk as a loglink of 0 + 2-dim P-spline
risk <- exp(as.vector(alpha+Bs%*%mvrnorm(1,rep(0,k),Sigma)))
summary(risk)

carto$Real_risk <- risk

tmap_mode("plot")
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.67,0.77,0.83,1,1.20,1.30,1.50,Inf)

Map <- tm_shape(carto) +
  tm_polygons(col="Real_risk", palette=paleta,id="ID", title="Real risk", legend.show=T, border.col="transparent",
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              popup.vars=c("ID","name","Real_risk")) +
  tm_layout(main.title="Real risk", main.title.position="center",
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            outer.margins=c(0.02,0.01,0.02,0.01))
print(Map)
