rm(list = ls())

library(bigDM)
library(MASS)
library(RColorBrewer)
library(sf)
library(splines)
library(tmap)

data("Carto_SpainMUN")
carto <- Carto_SpainMUN
carto <- carto[carto$region %in% c("Navarra","Pais Vasco"),]

############################################################

#Weird Kronecker product
Rten <- function(X1, X2) {
  one1 <- matrix(1, 1, ncol(X1))
  one2 <- matrix(1, 1, ncol(X2))
  kronecker(X1, one2) * kronecker(one1, X2)
}

# B-spline basis functions degree
bdeg <- 3

#Get coordinates of each centroid
locs <- st_coordinates(st_centroid(carto))

#Scale coordinates
locs <- apply(locs, 2, function(x) (x - min(x)) / (max(x) - min(x)))

#Extract coordinates lat/lon
x1 <- locs[, 1]
x2 <- locs[, 2]

#Specify number of knots in x1 direction and their locations
ndx1 <- 20
dis1 <- (max(x1) - min(x1)) / ndx1
x1l <- min(x1) - dis1 * 0.05
x1r <- max(x1) + dis1 * 0.05
dx1 <- (x1r - x1l) / ndx1
knots1 <- seq(x1l - bdeg * dx1, x1r + bdeg * dx1, by = dx1)

#Specify number of knots in x2 direction and their locations
ndx2 <- 20
dis2 <- (max(x2) - min(x2)) / ndx2
x2l <- min(x2) - dis2 * 0.05
x2r <- max(x2) + dis2 * 0.05
dx2 <- (x2r - x2l) / ndx2
knots2 <- seq(x2l - bdeg * dx2, x2r + bdeg * dx2, by = dx2)

#Make marginal design matrices for B-splines in x1 and x2 direction 
B1 <- spline.des(knots1, x1, bdeg + 1, 0 * x1)$design  ##Base marginal longitude
B2 <- spline.des(knots2, x2, bdeg + 1, 0 * x2)$design  ##Base marginal latitude

#Do that weird Kronecker product
Bs <- Rten(B2, B1)

#Find dimensions
k1 <- dim(B1)[2]
k2 <- dim(B2)[2]
k <- dim(Bs)[2]

#Define spatial domain (1, 25)
xt <- 1:25
xt <- (xt - min(xt)) / (max(xt) - min(xt))
dist <- (max(xt) - min(xt)) / 6
xtl <- min(xt) - dist * 0.05
xtr <- max(xt) + dist * 0.05
dxt <- (xtr - xtl) / 6
knots <- seq(xtl - bdeg * dxt, xtr + bdeg * dxt, by = dxt)
Bt <- spline.des(knots, xt, bdeg + 1)$design
kt <- ncol(Bt)

#Define the 3-dimensional P-spline basis function
Bst <- kronecker(Bt, Bs)

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

#Precision matrix of parameters for temporal P-spline
Dt <- diff(diag(kt), differences = 2)
Pt <- t(Dt) %*% Dt #RW2 precision matrix

#Kronecker between identity w. kt dim, identity with k2 dim and P1
RR1 <- kronecker(diag(kt), kronecker(diag(k2), P1)) 
RR2 <- kronecker(diag(kt), kronecker(P2, diag(k1))) #...
RR3 <- kronecker(Pt, kronecker(diag(k2), diag(k1))) #...

#Combined penalization matrix
Pst <- RR1 + RR2 + RR3

## Simulate risks from a multivariate normal distribution ##
set.seed(07101983)
n <- 1 # number of simulations
alpha <- 0
sig2st <- 0.05
Sigma_st <- sig2st * ginv(Pst)
risk <- exp(as.vector(alpha + Bst %*% mvrnorm(1, rep(0, kt*k), Sigma_st)))
summary(risk)

risk <- matrix(risk, nrow = nrow(carto), ncol = length(xt), byrow = F)
result <- cbind(carto, risk)

tmap_mode("plot")
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.67,0.77,0.83,1,1.20,1.30,1.50,Inf)

Map <- tm_shape(result) +
  tm_polygons(col=c("X1","X5","X10","X15","X20","X25"), palette=paleta,id="ID", title="Real risk", legend.show=T, border.col="transparent",
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              popup.vars=c("ID","name","x1")) +
  tm_layout(main.title="Real risk", main.title.position="center",
            legend.outside=T, legend.outside.position="right", legend.frame=F, legend.outside.size=0.2,
            outer.margins=c(0.02,0.01,0.02,0.01))
print(Map)


