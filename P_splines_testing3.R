#Clear environment
rm(list = ls())

#Load in data/functions and do necessary preliminaries
source('Preliminaries.R')
source("Utilities.R")

################################################################################
#Make a grid over Germany and make a B-spline basis

## Get a polygon_grid over Germany
polygon_grid = st_as_sf(st_make_grid(germany_border, 
                                     n = c(75, 75), 
                                     square = FALSE),
                        crs = crs_$srid)

## Add unique polygon_id
polygon_grid$polygon_id = 1:nrow(polygon_grid)

## Remove polygons outside of Germany
polygon_grid2 <- sf::st_intersection(polygon_grid, 
                                     st_make_valid(germany_border))

## Get the centroids of each polygon in the grid
centroids_grid <- st_centroid(polygon_grid)

## Get the coordinates of the centroids
gpts = sf::st_coordinates(centroids_grid)

## standardize x and y
locs <- apply(gpts, 2, function(x) (x - min(x)) / (max(x) - min(x)))

## Extract coordinates lat/lon
x <- locs[, 1]
y <- locs[, 2]

## Define number of knots in 1) the spatial directions and 2) temporal direction 
n_knots_s = 10; n_knots_t = 5

## define B-spline basis degree
bdeg = 3

# Specify the knot locations and create the B-spline basis design matrices
## Find the distance between each knot in x- and y-directions
disx <- (max(x) - min(x)) / n_knots_s; disy <- (max(y) - min(y)) / n_knots_s

## Find the left (l) and right (r) end and a bitmore spline not 0 at ends!
xl <- min(x) - disx * 0.05; yl <- min(y) - disy * 0.05
xr <- max(x) + disx * 0.05; yr <- max(y) + disy * 0.05

## This is distance between each knot?
dx <- (xr - xl) / n_knots_s; dy <- (yr - yl) / n_knots_s

## Specify the knot locations in x1 and x2 direction
knotsx <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
knotsy <- seq(yl - bdeg * dy, yr + bdeg * dy, by = dy)

## Make marginal design matrices for B-splines in x and y direction 
Bx <- spline.des(knotsx, x, bdeg + 1, 0 * x)$design  ##Base marginal x
By <- spline.des(knotsy, y, bdeg + 1, 0 * y)$design  ##Base marginal y

## Perform row-wise Kronecker product to get the spatial B-spline basis
Bs = row_wise_Kronecker(Bx, By)

## Find dimensions (spatial)
kx <- dim(Bx)[2]; ky <- dim(By)[2]; ks <- dim(Bs)[2]

## Define spatial domain (which discretely is: 1,...,13, and cont [0, 13])
t_axis <- seq(0.0001, 13, length.out = 57); xt_ = seq(0, 13, length.out = 57)

## Scale the temporal 'coordinates'
xt <- (t_axis - min(t_axis)) / (max(t_axis) - min(t_axis))

## Find distance 
dist <- (max(xt) - min(xt)) / n_knots_t
xtl <- min(xt) - dist * 0.05; xtr <- max(xt) + dist * 0.05
dxt <- (xtr - xtl) / n_knots_t

## Define knot locations in time
knots <- seq(xtl - bdeg * dxt, xtr + bdeg * dxt, by = dxt)

## Define temporal marginal B-spline basis design matrix
Bt <- spline.des(knots, xt, bdeg + 1)$design
kt <- ncol(Bt)

## Define the 3-dimensional P-spline basis function
Bst <- kronecker(Bt, Bs)

################################################################################
# Define the penalization matrices and get the risk field
## Define the marginal penalization matrices (RW1 precision matrices)
Px <- INLA:::inla.rw(n = kx, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Py <- INLA:::inla.rw(n = ky, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

Pt <- INLA:::inla.rw(n = kt, order = 1, 
                     scale.model = FALSE, 
                     sparse = TRUE)

## Identity matrices k1 x k1, k2 x k2, and kt x kt
Ix <- diag(kx); Iy <- diag(ky); It <- diag(kt)


## Combined spatial penalty matrix
tau_s = 25 # Isotropic spatial smoothness param
tau_t = 150 # Temporal smoothness parameter

Pst <- tau_s *(It %x% Py %x% Ix + It %x% Iy %x% Px) + 
  tau_t * (Pt %x% Iy %x% Ix)

## Dont know if needed, but I'll keep it (make it proper so that inla.qsample works)
Pst <- Pst + diag(x = 0.0001, nrow = ncol(Pst))

## 
n <- 1 # number of simulations
intercept = -7.6 #Intercept

## Sampled space-time P-spline parameters
set.seed(07101984)
sampled_theta_st = inla.qsample(n = 1, Q = Pst)

## Sum-to-zero (since IGMRF approx)
sampled_theta_st = sampled_theta_st - mean(sampled_theta_st)
rownames(sampled_theta_st) <- 1:nrow(sampled_theta_st)

## Get the space-time interaction
interactions = Bst %*% sampled_theta_st[, 1]

## Make a linearly increasing temporal effect
beta_t = 0.0142
tmp_ = xt_ * beta_t; temporal_effect = rep(0, length(interactions))
for(t in 1:length(xt)){
  temporal_effect[((t - 1) * (dim(Bs)[1]) + 1):(t * dim(Bs)[1])] = tmp_[t]
}

## Get the risk-field
Lambda_st <- exp(as.vector(intercept + temporal_effect + interactions))
################################################################################
#Display the risk-field

## Make a yxt-grid
### First make xt-grid (not with y), then add y
yxt_grid = expand.grid(x = gpts[, 1], t = t_axis)
y_ = rep(gpts[, 2], length(t_axis))
yxt_grid$y = y_

## add the risk-values to yxt_grid
yxt_grid$values = Lambda_st

## Add polygon_id values
yxt_grid$polygon_id = rep(centroids_grid$polygon_id, length(t_axis))

## Make geometric points for each grid point and add values
risk_surface.list = st_as_sf(yxt_grid, 
                             coords = c("x", "y"),
                             crs = crs_$srid)

risk_surface.list$x = yxt_grid$x; risk_surface.list$y = yxt_grid$y 

## join risk_surface.list to polygon_grid2 based on polygon_id in order to plot
test = data.frame(values = risk_surface.list$values, 
                  t = risk_surface.list$t,
                  polygon_id = risk_surface.list$polygon_id)

test_2 = data.frame(polygon_id = polygon_grid2$polygon_id,
                    geometry = polygon_grid2$x)

test_3 <- merge(test, test_2)
tmp_ = st_set_geometry(test_3[, c("t", "values", "polygon_id")], 
                         test_3$geometry)

# Make it so that each heatmap is plotted on similar color scale 
scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 50 
scale = scale_col[c(3, 7, 10, 15, 20, 24, 27, 30)] #Select color scale to be more red
risk.min = min(tmp_$values); risk.max = max(tmp_$values) 
hardcoded_bins =  round(seq(risk.min, risk.max, length.out = 8), 6)


for(t in t_axis){
  if(t == t_axis[1]){files = c(); i = 1}
  
  #Extract the values for time = t
  tmp2_ = tmp_[tmp_$t == t, ]
  
  #Create a filename (and save filename) for where current time risk-surface is stored
  filename = paste("./Plots/P_spline_tests/", toString(i), ".png", sep = "")
  print(filename); files[i] = filename
  
  p <- ggplot(data = tmp2_) +  
              geom_sf(aes(fill = values), 
              alpha = 1,
              color="lightgrey") + ggtitle(round(t, 1)) + 
    theme(plot.title = element_text(size = 15),
          axis.title.x = element_blank(), #Remove axis and background grid
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.margin =  unit(c(0, 0, 0, 0), "inches"),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
          panel.spacing = unit(1, 'lines')) +
    guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
    binned_scale( #Scaling the color
      aesthetics = "fill",
      scale_name = "gradientn",
      palette = function(x) c(scale),
      labels = function(x){x},
      breaks = hardcoded_bins,
      guide = "colorscale")
  
  ## Save the plot to filename
  ggsave(filename = filename, plot=p, 
         width=4,height=4,units="in",scale=1)
  
  #update counter i
  i = i + 1
}

### Make a GIF (for fun)
img = c(image_read(files[1]))
for(name in files[2:length(files)]){
  img = c(img, image_read(name))
}

image_append(image_scale(img, "x200"))

my.animation <- image_animate(image_scale(img, 
                                          "400x400"),
                              fps = 1,
                              dispose = "previous")

image_write(my.animation, "test.gif")

################################################################################

## Drop points outside Germany
#risk_surface.list = sf::st_intersection(risk_surface.list, 
#                                       st_make_valid(germany_border))

