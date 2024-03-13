#Clear environment
rm(list = ls())

library(splines)
library(ggfortify)

x = seq(0, 1, length.out = 100)


bdeg = 1

## Specify number of intervals in time (knots - 1)
n_intervals = 2

## Find distance between knots
dx <- (max(x) - min(x)) / n_intervals

## Define knot locations in time
knots <- seq(min(x) - bdeg * dx, max(x) + bdeg * dx, by = dx)

## Define temporal marginal B-spline basis design matrix
B <- spline.des(knots, x, bdeg + 1)$design

# Save to pdf as 7.5 by 3.5, name: B_spline_basis_func_deg_1
plot(B[, 2] ~ x, type = "l", ylim = c(0, 1), 
     ylab = "B(x)", xlab = "x", xaxt='n')
axis(1, at= c(0, 0.5, 1))


bdeg = 3

## Specify number of intervals in time (knots - 1)
n_intervals = 4

## Find distance between knots
dx <- (max(x) - min(x)) / n_intervals

## Define knot locations in time
knots <- seq(min(x) - bdeg * dx, max(x) + bdeg * dx, by = dx)

## Define temporal marginal B-spline basis design matrix
B <- spline.des(knots, x, bdeg + 1)$design


# Save to pdf as 7.5 by 3.5, name: B_spline_basis_func_deg_3
plot(B[, 4] ~ x, type = "l", ylim = c(0, 1), 
     ylab = "B(x)", xlab = "x", xaxt='n')
axis(1, at=knots[knots >= 0 & knots <= 1])




x = seq(1, 10, length.out = 100)


bdeg = 3

## Specify number of intervals in time (knots - 1)
n_intervals = 8

## Find distance between knots
dx <- (max(x) - min(x)) / n_intervals

## Define knot locations in time
knots <- seq(min(x) - bdeg * dx, max(x) + bdeg * dx, by = dx)

## Define temporal marginal B-spline basis design matrix
B <- spline.des(knots, x, bdeg + 1)$design


# Save to pdf as 7.5 by 3.5, name: B_spline_basis_deg_3
plot(B[, 4] ~ x, type = "l", ylim = c(0, 1), 
     ylab = "B(x)", xlab = "x", xaxt='n')
for(i in 4:(dim(B)[2] - 3)){lines(B[, i] ~ x)}
axis(1, at=knots)








