library(quantreg)
library(boot)
library(sp)
library(raster)
library(spatialEco)
library(rgeos)
library(png)
library(lqmm)

options(stringsAsFactors = FALSE)

make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

# growth functions

a_func <- function(r) {
  pi * r^2
}

r_func <- function(a) {
  sqrt(a / pi)
}

circularity <- function(area, perimeter) {
  (4 * pi * area)/(perimeter^ 2)
}

g_func <- function(x, y, qr0, p1, p2, p3) {
  # Actual radial growth
  g <- a_func(r_func(10^x) + qr0)
  # Partial mortality
  pm <- inv.logit(dnorm(logit(1-(10^y / g)), p1 + x * p2, p3))
  return(1 * pm)
}
