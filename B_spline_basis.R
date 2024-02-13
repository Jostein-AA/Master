#Clear environment
rm(list = ls())

library(splines)
library(ggfortify)

x <- seq(-0.5, 1.5, by=0.001)

spl <- bs(x, knots = c(-0.25, 0, 0.25, 0.5, 0.75, 1, 1.25), Boundary.knots = c(-0.5, 1.5), degree = 3) 

#Save as 7.5 by 3.5
plot(spl[,4]~x, 
     ylim=c(0,max(spl)), 
     type='l', 
     lwd=2, col=1, 
     xlab="Cubic B-spline basis",
     ylab="")
  lines(spl[, 5]~x, lwd = 2, col = 1)
  lines(spl[, 6]~x, lwd = 2, col = 1)


#autoplot(spl)
  

