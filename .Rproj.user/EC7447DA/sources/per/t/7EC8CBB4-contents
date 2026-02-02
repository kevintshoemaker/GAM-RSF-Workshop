# SCRIPT: prep real mule deer data from Mojave National Preserve

#   Example will use central Mojave National Preserve

# Clear workspace
rm(list = ls())
gc()

# Load required packages ----
library(terra)        # For raster operations
library(sf)           # For spatial vector operations
library(elevatr)      # For downloading elevation data
library(tidyverse)    # For data manipulation
library(leaflet)      # For interactive mapping
library(geodata)      # Alternative for elevation data
library(rapr)         # Access RAP products

# Global variables --------

## Set coordinate reference systems ----
crs_proj <- 32611     # NAD 1983 / UTM Zone 11N (Mojave is in Zone 11)
crs_geog <- 4269      # NAD 1983 (Geographic, degrees)
crs_wgs84 <- 4326     # WGS84 for leaflet and elevatr


# Define study area ----

study_area <- st_read("RSF_REAL/Mojave_originalfiles/Midhills.shp")

study_area_proj <- st_transform(study_area,crs_proj)
study_area_proj <- st_buffer(study_area_proj,2500)

study_area_g <- st_transform(study_area_proj,crs_wgs84)

# Create leaflet map with topographic relief layer
leaflet(study_area_g) %>%
  # Add OpenTopoMap (great topographic relief)
  addProviderTiles(providers$OpenTopoMap) %>%
  # Add your shapefile polygons
  addPolygons(
    fillColor = "blue",
    fillOpacity = 0.3,
    color = "darkblue",
    weight = 2,
    popup = ~paste("Feature ID:", row.names(.))
  ) %>%
  # Fit map bounds to your data
  fitBounds(
    lng1 = ~min(st_bbox(study_area_g)[1]),
    lat1 = ~min(st_bbox(study_area_g)[2]),
    lng2 = ~max(st_bbox(study_area_g)[3]),
    lat2 = ~max(st_bbox(study_area_g)[4])
  )

# Verify area
actual_area_km2 <- as.numeric(st_area(study_area_proj)) / 1000000
cat("Actual study area: ", round(actual_area_km2, 2), "km²\n\n")

# Visualize with leaflet ----

leaflet(study_area_g) %>%
  addProviderTiles("Esri.WorldTopoMap") %>%
  addPolygons(data = study_area_g, 
              color = "red", 
              weight = 3, 
              fill = TRUE,
              fillOpacity = 0.2) %>%
  fitBounds(
    lng1 = ~min(st_bbox(study_area_g)[1]),
    lat1 = ~min(st_bbox(study_area_g)[2]),
    lng2 = ~max(st_bbox(study_area_g)[3]),
    lat2 = ~max(st_bbox(study_area_g)[4])
  )


# PART 1: TOPOGRAPHIC DATA -------------------

## Download elevation data ----

# Download elevation data at ~30m resolution (z=12)
elevation <- get_elev_raster(study_area_g, z = 12, src = "aws")

# Convert to terra SpatRaster
elev_terra <- rast(elevation)

# Reproject to UTM
elev_utm <- terra::project(elev_terra, paste0("EPSG:", crs_proj), method = "bilinear")

# Crop to exact study area
elev_utm <- crop(elev_utm, vect(study_area_proj))
elev_utm <- mask(elev_utm, vect(study_area_proj))

# Set proper layer name
names(elev_utm) <- "elevation"

## Calculate topographic derivatives ----

# Slope (degrees)
slope <- terrain(elev_utm, v = "slope", unit = "degrees")
names(slope) <- "slope"

# Aspect (degrees)
aspect <- terrain(elev_utm, v = "aspect", unit = "degrees")
names(aspect) <- "aspect"

# Convert aspect to cosine and sine components
aspect_rad <- aspect * pi / 180
aspect_cos <- cos(aspect_rad)
aspect_sin <- sin(aspect_rad)
names(aspect_cos) <- "aspect_cos"
names(aspect_sin) <- "aspect_sin"

# Roughness (TRI - Terrain Ruggedness Index)
tri <- terrain(elev_utm, v = "TRI")
names(tri) <- "roughness"

# Topographic Position Index
tpi <- terrain(elev_utm, v = "TPI")
names(tpi) <- "tpi"

## Create a multi-layer raster stack ----
topo_stack <- c(elev_utm, slope, aspect, aspect_cos, aspect_sin, tri, tpi)
names(topo_stack) <- c("elevation", "slope", "aspect", "aspect_cos", 
                       "aspect_sin", "roughness", "tpi")

# PART 2: VEGETATION DATA -------------------

# Set years
years_to_download <- 2008:2016
cat("Years:", paste(years_to_download, collapse = ", "), "\n\n")

## Get RAP data --------

   # Takes a minute
rap <- get_rap(study_area_g, product = "vegetation-cover", years = years_to_download, verbose = FALSE ) 

temp=gsub("vegetation-cover_v3_","",names(rap))

afg_ndx = grepl("annual_forb_and_grass",temp)
pfg_ndx = grepl("perennial_forb_and_grass",temp)
shrub_ndx = grepl("shrub",temp)
tree_ndx = grepl("tree",temp)

rapdf = data.frame(
  year = readr::parse_number(temp)
)

rapdf$layer = NA
rapdf$layer[afg_ndx] = "AFG"
rapdf$layer[pfg_ndx] = "PFG"
rapdf$layer[shrub_ndx] = "SHRUB"
rapdf$layer[tree_ndx] = "TREE"

rapdf$fullname = names(rap)

rapdf <- na.omit(rapdf)

cover_types = c("AFG","PFG","SHRUB","TREE")
veg_stack = list()
veg_stack_yr = list()
veg_stack_all = list()  # this time I'll make a separate list element for each year

c=1
for(c in 1:length(cover_types)){
  thiscov=cover_types[c]
  thisdf = subset(rapdf,layer==thiscov)
  veg_stack_yr[[thiscov]] = rap[[thisdf$fullname]]
  names(veg_stack_yr[[thiscov]]) = paste0(thiscov,years_to_download)
  
  veg_stack[[thiscov]] = terra::app(veg_stack_yr[[thiscov]],mean)
}
veg_stack = rast(veg_stack)

veg_stack = terra::project(veg_stack,topo_stack)
plot(veg_stack)
c=1
for(c in 1:length(cover_types)){
  thiscov=cover_types[c]
  veg_stack_yr[[thiscov]] = terra::project(veg_stack_yr[[thiscov]],topo_stack)
  
}

veg_stack_yr$AFG

veg_stack_resampled <- terra::resample(veg_stack, topo_stack, method = "bilinear")

c=1
for(c in 1:length(cover_types)){
  thiscov=cover_types[c]
  veg_stack_yr[[thiscov]] = terra::resample(veg_stack_yr[[thiscov]],topo_stack,method = "bilinear")
}

# Put veg and topo variables together --------

# Create master environmental stack
env_stack <- c(topo_stack, veg_stack_resampled)

# Change to 90m resolution
if(all(res(env_stack)<30)){
  env_stack <- aggregate(env_stack,6)
  c=1
  for(c in 1:length(cover_types)){
    thiscov=cover_types[c]
    veg_stack_yr[[thiscov]] = terra::aggregate(veg_stack_yr[[thiscov]],6)
    
  }
} 
res(env_stack)
res(veg_stack_yr$AFG)

## crop to study area --------

env_stack <- crop(env_stack,study_area_proj)
env_stack <- mask(env_stack,study_area_proj)

plot(env_stack$roughness)

# organize by year - full veg covars each year
y=1
for(y in 1:length(years_to_download)){
  thisy = as.character(years_to_download[y])
  temp <- list()
  c=1
  for(c in 1:length(cover_types)){
    thiscov=cover_types[c]
    temp[[thiscov]] = veg_stack_yr[[thiscov]][[sprintf("%s%s",thiscov,thisy)]]
  }
  veg_stack_all[[thisy]] = rast(temp)
  veg_stack_all[[thisy]] = crop(veg_stack_all[[thisy]],vect(study_area_proj))
  veg_stack_all[[thisy]] = mask(veg_stack_all[[thisy]],vect(study_area_proj))
}

length(veg_stack_all)
names(veg_stack_all) = years_to_download

plot(veg_stack_all$`2010`)   # time varying covariate

# PART 3: WATER ------

peren_water <- st_read("RSF_REAL/Mojave_originalfiles/Perennial_water.shp")
peren_water_proj <- st_transform(peren_water,crs_proj)

peren_water_g <- st_transform(peren_water,crs_wgs84)

# Create leaflet map with topographic relief layer
leaflet(study_area_g) %>%
  # Add OpenTopoMap (great topographic relief)
  addProviderTiles(providers$OpenTopoMap) %>%
  # Add your shapefile polygons
  addCircleMarkers(
    data = peren_water_g,
    radius = 8,
    color = "black",
    fillColor = "yellow",
    fillOpacity = 0.7,
    weight = 2
  ) %>%
  addPolygons(
    fillColor = "blue",
    fillOpacity = 0.3,
    color = "darkblue",
    weight = 2,
    popup = ~paste("Feature ID:", row.names(.))
  ) |> 
  # Fit map bounds to your data
  fitBounds(
    lng1 = ~min(st_bbox(study_area_g)[1]),
    lat1 = ~min(st_bbox(study_area_g)[2]),
    lng2 = ~max(st_bbox(study_area_g)[3]),
    lat2 = ~max(st_bbox(study_area_g)[4])
  )


## summarize density of nearby water sources -----------

water_vect <- vect(peren_water_proj)

# Create a raster template
# Define extent and resolution based on your study area
r <- env_stack$elevation
values(r) = 1
r <- crop(r,vect(study_area_proj))
r <- mask(r,vect(study_area_proj))
plot(r)

# This applies a Gaussian kernel with distance weighting
density_raster <- rasterize(
  water_vect, 
  r, 
  fun = "length",  # count points
  background = 0
)

plot(density_raster )

# c(env_stack,density_raster)

# Apply distance-weighted kernel density
# sigma controls the bandwidth (distance decay)
density_weighted <- focal(
  density_raster,
  w = focalMat(r, d = 1500, type = "Gauss"),
  fun = sum,
  na.rm = TRUE
)

density_weighted = mask(density_weighted,vect(study_area_proj))
plot(density_weighted, main = "Distance-Weighted Water Source Density")

names(density_weighted) = "water"

names(env_stack)
env_stack = c(env_stack,density_weighted*1000)


# summary stats --------

cat("Summary statistics for all layers:\n")
for(i in 1:nlyr(env_stack)) {
  layer_name <- names(env_stack)[i]
  layer_stats <- global(env_stack[[i]], c("min", "max", "mean", "sd"), na.rm = TRUE)
  cat(sprintf("%-12s: min=%7.2f  max=%7.2f  mean=%7.2f  sd=%6.2f\n",
              layer_name,
              layer_stats$min,
              layer_stats$max,
              layer_stats$mean,
              layer_stats$sd))
}
cat("\n")

# Plot the layers ----

# Set up plotting parameters
# par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))

# # Plot each layer
# plot(elev_utm, main = "Elevation (m)", col = terrain.colors(100))
# plot(slope, main = "Slope (degrees)", col = heat.colors(100))
# plot(aspect, main = "Aspect (degrees)", col = rainbow(100))
# plot(aspect_cos, main = "Aspect Cosine", col = terrain.colors(100))
# plot(aspect_sin, main = "Aspect Sine", col = terrain.colors(100))
# plot(tri, main = "Roughness (TRI)", col = heat.colors(100))
# plot(tpi, main = "TPI", col = terrain.colors(100))

# Add study area boundary to one plot
# plot(elev_utm, main = "Study Area Boundary", col = terrain.colors(100))
# plot(vect(study_area_proj), add = TRUE, border = "red", lwd = 3)

# par(mfrow = c(1, 1))

# Save rasters ----

# Create output directory
output_dir <- "RSF_REAL/mojave_spatialdata"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

y=1
for(y in 1:length(years_to_download)){
  thisy = as.character(years_to_download[y])
  this = veg_stack_all[[thisy]]
  writeRaster(this, 
              file.path(output_dir, sprintf("veg_%s.tif",thisy)), 
              overwrite = TRUE)
}

# Save as multi-layer stack
writeRaster(env_stack, 
            file.path(output_dir, "covariate_stack.tif"), 
            overwrite = TRUE)

# Save study area boundary
st_write(study_area_proj, 
         file.path(output_dir, "study_area.shp"), 
         delete_dsn = TRUE)

cat("\nData saved to:", output_dir, "\n")


# END SCRIPT --------





