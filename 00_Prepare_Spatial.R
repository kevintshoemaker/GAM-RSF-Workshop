# Mojave National Preserve - Topographic Data Extraction
# Define 50 km² study area and extract topographic layers

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

cat("========================================\n")
cat("INTEGRATED ENVIRONMENTAL DATA EXTRACTION\n")
cat("Mojave National Preserve\n")
cat("========================================\n\n")

# Set coordinate reference systems ----
crs_proj <- 32611     # NAD 1983 / UTM Zone 11N (Mojave is in Zone 11)
crs_geog <- 4269      # NAD 1983 (Geographic, degrees)
crs_wgs84 <- 4326     # WGS84 for leaflet and elevatr

# Define study area parameters ----
# 50 km² = 50,000,000 m²
# Square: 7071 m × 7071 m (approximately 7.071 km × 7.071 km)

area_size_km2 <- 50
side_length_m <- sqrt(area_size_km2 * 1000000)  # 7071.068 meters

cat("Study area: ", area_size_km2, "km²\n")
cat("Side length: ", round(side_length_m), "meters\n\n")

# Define center point in Mojave National Preserve ----
# Using a location in the central part of the preserve
# Coordinates: approximately 35.2°N, -115.5°W (near Kelso Dunes)

center_lat <- 34.996
center_lon <- -115.5232

# Create center point as sf object ----
center_point <- st_point(c(center_lon, center_lat)) %>%
  st_sfc(crs = crs_wgs84) %>%
  st_sf(id = 1, geometry = .)

# Transform to projected coordinates
center_point_proj <- st_transform(center_point, crs_proj)
center_coords <- st_coordinates(center_point_proj)

# Create square study area ----
half_side <- side_length_m / 2

# Define corners of the square
xmin <- center_coords[1] - half_side
xmax <- center_coords[1] + half_side
ymin <- center_coords[2] - half_side
ymax <- center_coords[2] + half_side

# Create polygon
study_area_proj <- st_polygon(list(
  matrix(c(xmin, ymin,
           xmax, ymin,
           xmax, ymax,
           xmin, ymax,
           xmin, ymin), 
         ncol = 2, byrow = TRUE)
)) %>%
  st_sfc(crs = crs_proj) %>%
  st_sf(id = 1, geometry = .)

# Verify area
actual_area_km2 <- as.numeric(st_area(study_area_proj)) / 1000000
cat("Actual study area: ", round(actual_area_km2, 2), "km²\n\n")

# Transform to WGS84 for elevation download
study_area_wgs84 <- st_transform(study_area_proj, crs_wgs84)

# Visualize with leaflet ----
cat("Creating interactive map...\n")

leaflet() %>%
  addProviderTiles("Esri.WorldTopoMap") %>%
  addPolygons(data = study_area_wgs84, 
              color = "red", 
              weight = 3, 
              fill = TRUE,
              fillOpacity = 0.2,
              popup = paste("Study Area: ", round(actual_area_km2, 2), "km²")) %>%
  addMarkers(data = center_point, 
             popup = "Center Point") %>%
  setView(lng = center_lon, lat = center_lat, zoom = 12)


# ========================================
# PART 1: TOPOGRAPHIC DATA
# ========================================

cat("\n========================================\n")
cat("DOWNLOADING TOPOGRAPHIC DATA\n")
cat("========================================\n\n")

# Download elevation data ----
cat("\nDownloading elevation data...\n")
cat("(This may take a minute)\n\n")

# Download elevation data at ~30m resolution (z=12)
elevation <- get_elev_raster(study_area_wgs84, z = 12, src = "aws")

# Convert to terra SpatRaster
elev_terra <- rast(elevation)

# Reproject to UTM
cat("Reprojecting to UTM Zone 11N...\n")
elev_utm <- project(elev_terra, paste0("EPSG:", crs_proj), method = "bilinear")

# Crop to exact study area
elev_utm <- crop(elev_utm, vect(study_area_proj))
elev_utm <- mask(elev_utm, vect(study_area_proj))

# Set proper layer name
names(elev_utm) <- "elevation"

cat("\nElevation statistics:\n")
cat("  Min: ", round(minmax(elev_utm)[1]), "m\n")
cat("  Max: ", round(minmax(elev_utm)[2]), "m\n")
cat("  Mean: ", round(global(elev_utm, "mean", na.rm = TRUE)[1,1]), "m\n")
cat("  Resolution: ", round(res(elev_utm)[1]), "m ×", round(res(elev_utm)[2]), "m\n")
cat("  Dimensions: ", nrow(elev_utm), "rows ×", ncol(elev_utm), "cols\n\n")

# Calculate topographic derivatives ----
cat("Calculating topographic variables...\n\n")

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

# Create a multi-layer raster stack ----
topo_stack <- c(elev_utm, slope, aspect, aspect_cos, aspect_sin, tri, tpi)
names(topo_stack) <- c("elevation", "slope", "aspect", "aspect_cos", 
                       "aspect_sin", "roughness", "tpi")

cat("Topographic layers created:\n")
print(names(topo_stack))
cat("\n")


# ========================================
# PART 2: VEGETATION DATA
# ========================================

cat("========================================\n")
cat("DOWNLOADING VEGETATION DATA (RAP)\n")
cat("========================================\n\n")

# Set years
years_to_download <- 2018:2024
cat("Years:", paste(years_to_download, collapse = ", "), "\n\n")

# Get RAP data --------

rap <- get_rap(study_area_wgs84, product = "vegetation-cover", years = years_to_download, verbose = FALSE ) 

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
c=1
for(c in 1:length(cover_types)){
  thiscov=cover_types[c]
  thisdf = subset(rapdf,layer==thiscov)
  tempstack = rap[[thisdf$fullname]]
  veg_stack[[thiscov]] = terra::app(tempstack,mean)
}
veg_stack = rast(veg_stack)

veg_stack = terra::project(veg_stack,topo_stack)
plot(veg_stack)


# Put veg and topo variables together --------

cat("Resampling vegetation data to match topography...\n")
veg_stack_resampled <- terra::resample(veg_stack, topo_stack, method = "bilinear")

# Create master environmental stack
env_stack <- c(topo_stack, veg_stack_resampled)


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
cat("Creating plots...\n")

# Set up plotting parameters
par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))

# Plot each layer
plot(elev_utm, main = "Elevation (m)", col = terrain.colors(100))
plot(slope, main = "Slope (degrees)", col = heat.colors(100))
plot(aspect, main = "Aspect (degrees)", col = rainbow(100))
plot(aspect_cos, main = "Aspect Cosine", col = terrain.colors(100))
plot(aspect_sin, main = "Aspect Sine", col = terrain.colors(100))
plot(tri, main = "Roughness (TRI)", col = heat.colors(100))
plot(tpi, main = "TPI", col = terrain.colors(100))

# Add study area boundary to one plot
plot(elev_utm, main = "Study Area Boundary", col = terrain.colors(100))
plot(vect(study_area_proj), add = TRUE, border = "red", lwd = 3)

par(mfrow = c(1, 1))

# Save rasters ----
cat("Saving raster data...\n")

# Create output directory
output_dir <- "mojave_topographic_data"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Save individual layers
writeRaster(elev_utm, 
            file.path(output_dir, "elevation.tif"), 
            overwrite = TRUE)
writeRaster(slope, 
            file.path(output_dir, "slope.tif"), 
            overwrite = TRUE)
writeRaster(aspect, 
            file.path(output_dir, "aspect.tif"), 
            overwrite = TRUE)
writeRaster(aspect_cos, 
            file.path(output_dir, "aspect_cos.tif"), 
            overwrite = TRUE)
writeRaster(aspect_sin, 
            file.path(output_dir, "aspect_sin.tif"), 
            overwrite = TRUE)
writeRaster(tri, 
            file.path(output_dir, "roughness.tif"), 
            overwrite = TRUE)
writeRaster(tpi, 
            file.path(output_dir, "tpi.tif"), 
            overwrite = TRUE)

c=1
for(c in 1:length(cover_types)){
  this = env_stack[[cover_types[c]]]
  writeRaster(this, 
              file.path(output_dir, sprintf("%s.tif",cover_types[c])), 
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

# Create summary dataframe ----
summary_df <- data.frame(
  layer = names(topo_stack),
  min = sapply(1:nlyr(topo_stack), function(i) {
    global(topo_stack[[i]], "min", na.rm = TRUE)[1,1]
  }),
  max = sapply(1:nlyr(topo_stack), function(i) {
    global(topo_stack[[i]], "max", na.rm = TRUE)[1,1]
  }),
  mean = sapply(1:nlyr(topo_stack), function(i) {
    global(topo_stack[[i]], "mean", na.rm = TRUE)[1,1]
  }),
  sd = sapply(1:nlyr(topo_stack), function(i) {
    global(topo_stack[[i]], "sd", na.rm = TRUE)[1,1]
  })
)

write.csv(summary_df, 
          file.path(output_dir, "topographic_summary.csv"), 
          row.names = FALSE)

cat("\nSummary statistics saved to:", 
    file.path(output_dir, "topographic_summary.csv"), "\n")

# Print summary ----
cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Study area: 50 km² in Mojave National Preserve\n")
cat("Center: ", center_lat, "°N, ", center_lon, "°W\n")
cat("CRS: UTM Zone 11N (EPSG:", crs_proj, ")\n")
cat("Resolution: ~", round(res(elev_utm)[1]), "m\n")
cat("Layers created: ", nlyr(topo_stack), "\n")
cat("Output directory: ", output_dir, "\n")
cat("========================================\n")


# END SCRIPT --------