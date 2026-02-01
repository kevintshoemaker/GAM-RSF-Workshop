# Poisson point process for presence only data featuring Joshua Trees

rm(list=ls())

# load packages  ------
pkgs <- c("here", "readr", "janitor", "mgcv", "ppgam", "gratia", "dplyr",
          "ggplot2", "mvnfast", "patchwork", "tidyr")
vapply(pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)


# load data -------

# data - from Godsoe et al 2009 New Phytologist
# https://doi.org/10.1111/j.1469-8137.2009.02942.x
# download from dryad https://doi.org/10.5061/dryad.6s67t

joshua <- read_tsv(here("data/joshua-tree.txt"),
                   col_types = "ddc")


# plot/explore -------

ggplot(joshua, aes(x = longitude, y = latitude)) +
  geom_point() +
  coord_quickmap()


library(sf)
joshua_sp = st_as_sf(joshua,coords = c("longitude","latitude"),remove=F,crs=4326)

# plot(joshua_sp$geometry)
  
library(leaflet)
leaflet(joshua_sp) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addCircleMarkers(
    radius = 5,
    color = "red",
    fillOpacity = 0.7,
    stroke = TRUE,
    weight = 1
  )

# [note: there is nothing really random about the sampling process here. look closely!]


#  Use poisson point process model ---------

# fit the Poisson point process  [note: probably not an appropriate model here]
# about 5 seconds

m <- ppgam(~ s(longitude, latitude, k = 100), data = joshua)

summary(m)
gratia::draw(m, rug = FALSE)

# library(DHARMa)    # DHARMa tests don't work on this kind of model
# r=simulateResiduals(m)

     # try duchon spline - m parameter allows more control over the smooth
        # here we penalize the first derivative (gradient). 
             # second term in 'm' is harder to understand (penalty on null space)
m_ds <- ppgam(~ s(longitude, latitude, k = 100, bs = "ds", m = c(1,0.5)),
              data = joshua)

# summary(m_ds)
gratia::draw(m_ds,rug = FALSE)
  
     # try gaussian process smoother
m_gp <- ppgam(~ s(longitude, latitude, k = 100, bs = "gp"),
              data = joshua)

summary(m_gp)
gratia::draw(m_gp,rug = FALSE)   # looks a lot different!


# lattice method- break up into grid and count points  -----

library(sf)
sf::st_coordinates(joshua_sp)

library(terra)

# Define study area from buffered points
# Buffer distance in meters (e.g., 10000 = 10 km)
# First transform to a projected CRS for accurate buffering
joshua_sp_proj <- st_transform(joshua_sp, crs = 32611)  # UTM Zone 11N for California/Nevada

# Buffer the points
joshua_buffered <- st_buffer(joshua_sp_proj, dist = 25000)  # 25 km buffer

# Dissolve all buffered points into a single polygon
study_area <- st_union(joshua_buffered)

plot(study_area)

# Simplify if needed to reduce complexity
# study_area <- st_simplify(study_area, dTolerance = 100)

# Transform back to WGS84 for raster creation
study_area_wgs84 <- st_transform(study_area, crs = 4326)

# Visualize the study area
plot(st_geometry(study_area_wgs84), col = "lightblue", border = "blue")
plot(st_geometry(joshua_sp), add = TRUE, pch = 16, cex = 0.3, col = "red")


# Now create raster within the study area
# Convert study area to SpatVector for terra
study_area_vect <- vect(study_area_wgs84)

# Create raster from study area extent
joshua_raster <- rast(study_area_vect, res = 0.1)

# Count observations in each cell
joshua_vect <- vect(joshua_sp)
joshua_count <- rasterize(joshua_vect, joshua_raster, fun = "length", background=0)

joshua_count = mask(joshua_count,study_area_vect)

plot(joshua_count)

library(ape)

# Convert raster to data frame with x, y coordinates
joshua_df <- as.data.frame(joshua_count, xy = TRUE,cells=T,na.rm=T)

# View the result
head(joshua_df)
nrow(joshua_df)


# Rename the count column if desired
colnames(joshua_df) <- c("cell", "x", "y", "count")

# Visualize the result
ggplot(joshua_df, aes(x = x, y = y, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_quickmap() +
  labs(title = "Joshua Tree Observation Density",
       fill = "Count") +
  theme_minimal()

hist(joshua_df$count)  # looks zero inflated

# look at Moran's I statistic (spatial autocorrelation)

?Moran.I
d = 1/as.matrix(dist(joshua_df[,c("x","y")]))
dim(d)
diag(d) = 0

names(joshua_df)
ape::Moran.I(joshua_df$count,d)  # test for spatial correlation


## run count regression on the gridded data ------

josh1 <- gam(count ~ s(x, y, k = 100),
             family=ziP(),  data = joshua_df, method = 'REML')

library(gratia)
gratia::draw(josh1)


# try zi poisson with smoother on both terms

   # ziplss- uses location scale shape framework- allows for smoother on zi and p models
josh2 <- gam(list(count ~ s(x, y, k = 100),~1),
              family=ziplss(),  data = joshua_df, method = 'REML')
gratia::draw(josh2)


# try using sdmTMB package --------

library(sdmTMB)

names(joshua_df)

joshua_df = sdmTMB::add_utm_columns(joshua_df,c("x","y"))

mesh <- make_mesh(joshua_df, xy_cols = c("X", "Y"), cutoff = 10)
mesh$mesh$n   # 839 mesh nodes vs 5765 grid cells

plot(mesh)

library(inlabru)
ggplot() +
  inlabru::gg(mesh$mesh) +
  geom_point(data=joshua_df, aes(X,Y,col=count),size=0.4)


##  fit sdmTMB model ----

pc <- pc_matern(range_gt = 10, sigma_lt = 5)

fit <- sdmTMB(
  count ~ 1,
  data = joshua_df,
  mesh = mesh,
  priors = sdmTMBpriors(matern_s = pc),   # note the prior to help with convergence
  family = delta_truncated_nbinom2(),   # hurdle model (like a zero inflation model)
  spatial = "on"
)

summary(fit)
sanity(fit)

tidy(fit, model=1,effects = "ran_pars",conf.int = T)  # what if we want the count model params
tidy(fit, model=2,effects = "ran_pars",conf.int = T)

x_range <- range(joshua_df$X)
y_range <- range(joshua_df$Y)

p <- predict(fit, type="link")  

# p <- predict(fit)

   # zero inflation
names(p)
hist(p$est1); summary(p$est1)
hist(p$est2)
hist(p$x)
hist(p$y)

ggplot(p, aes(x, y, z = est1)) + geom_contour_filled() 

   # count model
ggplot(p, aes(x, y, z = est2)) + geom_contour_filled()
   

# note: could easily plot out the spatial field [omega], or the spatial field components, or the fixed effects, etc.

   # DHARMa tests 

s = simulate(fit, nsim = 500, type = "mle-mvn")
res = dharma_residuals(s, fit, return_DHARMa = T)  

library(DHARMa)

testResiduals(res)

testSpatialAutocorrelation(res,joshua_df$X, joshua_df$Y)


   # still not great fit... 


# Try just fitting the non-zero grid cells (truncated) -------

joshua_df2 = subset(joshua_df,count>0)
mesh2 <- make_mesh(joshua_df2, xy_cols = c("X", "Y"), cutoff = 10)

plot(mesh2)

pc <- pc_matern(range_gt = 10, sigma_lt = 5)
fit2 <- sdmTMB(
  count ~ 1,
  data = joshua_df2,
  mesh = mesh2,
  family = truncated_nbinom2(),   # the non-zero component of the hurdle model
  spatial = "on"
)

summary(fit2)
sanity(fit2)

p <- predict(fit2, newdata=joshua_df, type="response")

# p <- predict(fit)

# expected count intensity
ggplot(p, aes(x, y, z = est)) + geom_contour_filled()


# DHARMa tests 

s = simulate(fit, nsim = 500, type = "mle-mvn")
dharma_residuals(s, fit, test_uniformity = TRUE)  # okay?

res = dharma_residuals(s, fit, return_DHARMa = T)

testResiduals(res)

testSpatialAutocorrelation(res,joshua_df$X, joshua_df$Y)


# test without spatial random field

fit_null = update(fit2,spatial="off")

AIC(fit2,fit_null)  # spatial model is better based on AIC



# END SCRIPT






