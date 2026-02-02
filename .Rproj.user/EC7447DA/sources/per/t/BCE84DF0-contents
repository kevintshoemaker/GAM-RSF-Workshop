# Generate Additional Synthetic Raster Layers

#   this script generates "true" underlying static habitat suitability and avoidance layers 
     #  with varying degrees of spatial autocorrelation. 

# To do this we will use sdmTMB to generate spatial random fields, 

rm(list = ls())
gc()

# Load packages ----

library(terra)
library(sf)
library(tidyverse)
library(gstat)      # For spatial interpolation/kriging
# library(RandomFields) # For generating spatial random fields
# If RandomFields doesn't work, we'll use alternative methods

library(sdmTMB)

# Load existing data ----

data_dir <- "RSF_SIM/mojave_spatialdata_forsim"

# Load the covariate stack
env_stack <- rast(file.path(data_dir, "covariate_stack.tif"))
study_area <- st_read(file.path(data_dir, "study_area.shp"), quiet = TRUE)

random_points <- st_sample(study_area, size = 200)
random_points <- st_coordinates(random_points) 

env_stack
env_stack2 <- aggregate(env_stack,4)
# res(env_stack2)

# Extract resolution and extent info
template1 <- env_stack[[1]]
template2 <- env_stack2[[1]]
res_x <- res(template1)[1]
res_y <- res(template1)[2]
n_cells <- ncell(template1)


# PREPARE DATA FOR SDMTMB -------

    # convert rasters to data frame
# predictor_dat = terra::extract(env_stack,vect(random_points),xy=T,cells=T)

predictor_dat = as.data.frame(env_stack2, xy=T,cells=T)

predictor_dat[,c("x2","y2")] =  predictor_dat[,c("x","y")]/1000   # convert to km

range(predictor_dat$x2); range(predictor_dat$y2)

names(predictor_dat)
mesh = make_mesh(predictor_dat,xy_cols=c("x2","y2"),cutoff=1)
mesh$mesh$n
plot(mesh)

# CREATE 'TRUE' HABITAT SUITABILITY LAYER -----------

# Create a weighted combination of existing variables
# This represents a hypothetical "resource availability" index

# Normalize key variables to 0-1 scale
normalize <- function(x) {
  x_min <- global(x, "min", na.rm = TRUE)[1,1]
  x_max <- global(x, "max", na.rm = TRUE)[1,1]
  (x - x_min) / (x_max - x_min)
}

elev_norm <- normalize(env_stack$elevation)
slope_norm <- normalize(env_stack$slope)
shrub_norm <- normalize(env_stack$SHRUB)
pfg_norm <- normalize(env_stack$PFG)
slopeasp_norm <- normalize(env_stack$slope*(-env_stack$aspect_cos))

# Create composite index with ecological reasoning
# Higher at mid-elevations, moderate slopes, good vegetation
resource_index <- (
  0.3 * (1 - abs(elev_norm - 0.5) * 2) +  # Prefer mid-elevations
  0.2 * (1-abs(slope_norm - 0.7) * 2) +      # Prefer steeper slopes
  0.25 * shrub_norm +                      # More shrubs = more resources
  0.2 * pfg_norm +                          # More perennial grass = more resources
  0.15 * (1 - abs(elev_norm - 0.4) * 2)   # Prefer moderate south facing slopes - avoid very steep slopes
)

## generate a spatial smooth (create additional spatial autocorrelation)

resource_index <- focal(resource_index, w = 5, fun = mean, na.rm = TRUE)

plot(resource_index)   # store this for comparison with model results

# CREATE 'TRUE' AVOIDANCE LAYER --------------

names(env_stack)

# Start with a moderate correlation to existing variables
rough_norm <- normalize(env_stack$roughness)  # lack of rougnness
rough_norm2 = template1
rough_norm2[] = 0
rough_norm2[rough_norm<0.04] = 1
plot(rough_norm2)
tree_norm <- 1-normalize(env_stack$TREE)     # lack of trees
plot(tree_norm)
herb_norm <- 1-normalize(env_stack$AFG+env_stack$PFG)     # lack of herbaceous veg
plot(herb_norm)
trshr_norm <- 1-normalize(env_stack$SHRUB+env_stack$TREE)  # lack of tree and shrub cover
plot(trshr_norm)

avoidance_index <- (
  0.3 * herb_norm + 
  0.2 * trshr_norm + 
  0.2 * tree_norm +
  0.3 * rough_norm2
)

avoidance_index <- focal(avoidance_index, w = 7, fun = mean, na.rm = TRUE)
plot(avoidance_index)


# LAYER 1: SRF WITH NO CORRELATION WITH EXISTING VARIABLES -------------


# ?sdmTMB_simulate
sim_dat1 = sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link="identity"),
  range = 5,   # try changing this (spatial range)
  sigma_O = 0.2,  # try changing this (spatial sd)
  phi = 0.01,  # not used
  B=0   # fixed effect coefficient values
)

class(sim_dat1)

ggplot(sim_dat1, aes(x2,y2,fill=mu)) +
  geom_raster() +
  scale_fill_gradient2()

## Try playing around with the range and sigma_O parameters and plotting the results!

nrow(predictor_dat)
nrow(sim_dat1)
template2
srf1_nocovars = template2
srf1_nocovars[predictor_dat$cell] = sim_dat1$mu
srf1_nocovars = terra::resample(srf1_nocovars,template1)

plot(srf1_nocovars)

## FOR FUN: play with spatiotemporal random fields -------
     # ignore if we already did this earlier
     
# strf = sdmTMB_simulate(
#   formula = ~ 1,
#   data = predictor_dat,
#   mesh = mesh,
#   family = gaussian(link="identity"),
#   range = 5,   # try changing this (spatial range)
#   sigma_O = 0.2,  # try changing this (spatial sd)
#   phi = 0.01,  # not used
#   B=0   # fixed effect coefficient values
# )
# 
# class(sim_dat1)
# 
# ggplot(sim_dat1, aes(x2,y2,fill=mu)) +
#   geom_raster() +
#   scale_fill_gradient2()


# LAYER 2: SRF WITH HIGH CORRELATION WITH HABITAT SUITABILITY-----------

temp = terra::resample(resource_index,template2)
# temp = as.data.frame(temp,xy=T,cells=T)

predictor_dat$res_ndx = temp[predictor_dat$cell][,1]

temp = terra::resample(avoidance_index,template2)
predictor_dat$avd_ndx = temp[predictor_dat$cell][,1]

# ?sdmTMB_simulate
sim_dat2 = sdmTMB_simulate(
  formula = ~ res_ndx + avd_ndx,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link="identity"),
  range = 2,   # try changing this (spatial range)
  sigma_O = 0.05,  # try changing this (spatial sd)
  phi = 0.01,  # not used
  B=c(0,0.4,-0.2)   # fixed effect coefficient values
)

class(sim_dat2)

ggplot(sim_dat2, aes(x2,y2,fill=mu)) +
  geom_raster() +
  scale_fill_gradient2()

names(sim_dat2)
ggplot(sim_dat2, aes(x2,y2,fill=omega_s)) +
  geom_raster() +
  scale_fill_gradient2()

srf2_hs1 = template2
srf2_hs1[predictor_dat$cell] = sim_dat2$mu
srf2_hs1 = terra::resample(srf2_hs1,template1)

plot(srf2_hs1)

# LAYER 3: SRF WITH LOWER CORRELATION WITH HABITAT SUITABILITY LAYERS ------------

# ?sdmTMB_simulate
sim_dat3 = sdmTMB_simulate(
  formula = ~ res_ndx + avd_ndx,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link="identity"),
  range = 4,   # try changing this (spatial range)
  sigma_O = 0.1,  # try changing this (spatial sd)
  phi = 0.01,  # not used
  B=c(0,0.3,-0.3)   # fixed effect coefficient values
)

class(sim_dat3)

names(sim_dat3)    # plot spatial random field
ggplot(sim_dat3, aes(x2,y2,fill=omega_s)) +
  geom_raster() +
  scale_fill_gradient2()

sim_dat3$fixef = as.matrix(sim_dat3[,c("(Intercept)","res_ndx","avd_ndx")]) %*% matrix(c(0,0.3,-0.3),ncol=1)
names(sim_dat3)    # plot fixed effects
ggplot(sim_dat3, aes(x2,y2,fill=fixef)) +
  geom_raster() +
  scale_fill_gradient2()

ggplot(sim_dat3, aes(x2,y2,fill=mu)) +
  geom_raster() +
  scale_fill_gradient2()

srf3_hs2 = template2
srf3_hs2[predictor_dat$cell] = sim_dat3$mu
srf3_hs2 = terra::resample(srf3_hs2,template1)

plot(srf3_hs2)

# LAYER 4: SRF WITH CORRELATION ONLY WITH ELEVATION ----------

temp = terra::resample(elev_norm,template2)
# temp = as.data.frame(temp,xy=T,cells=T)

predictor_dat$elev_ndx = temp[predictor_dat$cell][,1]

# ?sdmTMB_simulate
sim_dat4 = sdmTMB_simulate(
  formula = ~ elev_ndx + I(elev_ndx^2),
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link="identity"),
  range = 4,   # try changing this (spatial range)
  sigma_O = 0.01,  # try changing this (spatial sd)
  phi = 0.01,  # not used
  B=c(0,0.3,-0.3)   # fixed effect coefficient values
)

names(sim_dat4)    # plot spatial random field
ggplot(sim_dat4, aes(x2,y2,fill=omega_s)) +
  geom_raster() +
  scale_fill_gradient2()

sim_dat4$fixef = as.matrix(sim_dat4[,c("(Intercept)","elev_ndx","I(elev_ndx^2)")]) %*% matrix(c(0,0.3,-0.3),ncol=1)
# names(sim_dat4)    # plot fixed effects
ggplot(sim_dat4, aes(x2,y2,fill=fixef)) +
  geom_raster() +
  scale_fill_gradient2()

ggplot(sim_dat4, aes(x2,y2,fill=mu)) +
  geom_raster() +
  scale_fill_gradient2()

srf4_elev1 = template2
srf4_elev1[predictor_dat$cell] = sim_dat4$mu*5
srf4_elev1 = terra::resample(srf4_elev1,template1)

plot(srf4_elev1)


# COMBINE AND ANALYZE ------------

# Create a stack of new layers
new_layers <- c(resource_index,avoidance_index, srf1_nocovars, srf2_hs1, srf3_hs2,srf4_elev1)
names(new_layers) = c("resource_index","avoidance_index", "srf1_nocovars", "srf2_hs1", "srf3_hs2","srf4_elev1")

# Calculate correlations with existing layers

# Sample points to calculate correlation (use all non-NA cells)
sample_df <- as.data.frame(c(env_stack, new_layers), na.rm = TRUE)

# Calculate correlations for each new layer
for(new_layer_name in names(new_layers)) {
  cat(new_layer_name, ":\n")
  
  # Select a few key existing variables
  key_vars <- c("elevation", "slope", "SHRUB", "TREE", "AFG", "PFG")
  
  for(var in key_vars) {
    if(var %in% names(sample_df)) {
      cor_val <- cor(sample_df[[var]], sample_df[[new_layer_name]], 
                     use = "complete.obs")
      cat(sprintf("  vs %-12s: r = %6.3f\n", var, cor_val))
    }
  }
  cat("\n")
}

# SAVE OUTPUTS -------------------

# Save individual layers

writeRaster(resource_index, 
            file.path(data_dir, "resource_index.tif"), 
            overwrite = TRUE)
writeRaster(avoidance_index, 
            file.path(data_dir, "avoidance_index.tif"), 
            overwrite = TRUE)
writeRaster(srf1_nocovars, 
            file.path(data_dir, "srf1_nocovars.tif"), 
            overwrite = TRUE)
writeRaster(srf2_hs1, 
            file.path(data_dir, "srf2_hs1.tif"), 
            overwrite = TRUE)
writeRaster(srf3_hs2, 
            file.path(data_dir, "srf3_hs2.tif"), 
            overwrite = TRUE)
writeRaster(srf4_elev1, 
            file.path(data_dir, "srf_elev1.tif"), 
            overwrite = TRUE)

# Save as multi-layer stack
writeRaster(new_layers, 
            file.path(data_dir, "unseen_layers.tif"), 
            overwrite = TRUE)



# END SCRIPT  -------