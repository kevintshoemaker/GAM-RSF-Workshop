
# TELEMETRY SIMULATION 

## ideas for future

##  try simulating a 'pure' point process from the habitat suitability layer. Compare with more realistic correlated movement model

rm(list = ls())
gc()

library(terra)
library(sf)
library(tidyverse)
library(lubridate)
library(circular)

# LOAD DATA----------------

data_dir <- "mojave_spatialdata_forsim"
study_area <- st_read(file.path(data_dir, "study_area.shp"), quiet = TRUE)
unseen_layers <- rast(file.path(data_dir, "unseen_layers.tif"))

# Prepare habitat layers
normalize_raster <- function(r) {
  r_min <- global(r, "min", na.rm = TRUE)[1,1]
  r_max <- global(r, "max", na.rm = TRUE)[1,1]
  (r - r_min) / (r_max - r_min)
}

habitat <- normalize_raster(unseen_layers$resource_index - 0.5 * unseen_layers$avoidance_index)


# PARAMETERS ---------------

crs_proj <- 32611     # NAD 1983 / UTM Zone 11N (Mojave is in Zone 11)

# Population
n_animals <- 25
years <- 2015:2025

# Tracking
min_years <- 1
max_years <- 5
fixes_per_day <- 2  # Simplified: all animals same rate

# Movement (simplified for speed)
step_mean <- 100
step_sd <- 50
turn_concentration <- 0.5
habitat_beta <- 3      # Selection coefficient
hr_penalty <- -0.003    # Distance penalty

set.seed(123)

cat("  Parameters set\n\n")

# SETUP INDIVIDUALS --------------

cat("Generating individuals...\n")

individuals <- tibble(
  id = sprintf("DEER%03d", 1:n_animals),
  n_years = sample(min_years:max_years, n_animals, replace = TRUE),
  start_year = sample(years[1]:(years[length(years)] - max_years + 1), 
                     n_animals, replace = TRUE)
) %>%
  mutate(end_year = start_year + n_years - 1)


# Home range centers (vectorized selection)
hr_prob <- exp(5 * habitat)
hr_prob <- hr_prob / global(hr_prob, "sum", na.rm = TRUE)[1,1]
hr_cells <- sample(ncell(habitat), n_animals, 
                  prob = values(hr_prob, mat = FALSE), replace = TRUE)
hr_xy <- xyFromCell(habitat, hr_cells)

individuals$hr_x <- hr_xy[, 1]
individuals$hr_y <- hr_xy[, 2]

cat("  ", n_animals, "individuals created\n\n")

# Get coordinates of home range centers
hr_coords <- xyFromCell(habitat, hr_cells)

plot(habitat)
points(hr_coords)


# MOVEMENT SIMULATION --------------

start_time <- Sys.time()

# Pre-calculate total steps needed
individual_years <- individuals %>%
  rowwise() %>%
  do({
    tibble(
      id = .$id,
      year = .$start_year:.$end_year,
      hr_x = .$hr_x,
      hr_y = .$hr_y
    )
  }) %>%
  ungroup()

# Calculate steps per individual-year
individual_years <- individual_years %>%
  mutate(
    start_date = if_else(year == min(year), 
                        as.Date(paste0(year, "-03-01")),
                        as.Date(paste0(year, "-01-01"))),
    end_date = if_else(year == max(year),
                      as.Date(paste0(year, "-10-31")),
                      as.Date(paste0(year, "-12-31"))),
    n_days = as.numeric(end_date - start_date + 1),
    n_steps = round(n_days * fixes_per_day),
    id_year = paste(id, year, sep = "_")
  )

total_steps <- sum(individual_years$n_steps)
cat("Total steps to simulate:", total_steps, "\n\n")

# Simulate all paths
all_paths <- vector("list", nrow(individual_years))

pb <- txtProgressBar(min = 0, max = nrow(individual_years), style = 3)

i=1
for(i in 1:nrow(individual_years)) {
  
  ind_yr <- individual_years[i, ]
  
  # Starting location
  if(ind_yr$year == individuals$start_year[individuals$id == ind_yr$id]) {
    start_x <- ind_yr$hr_x + rnorm(1, 0, 500)
    start_y <- ind_yr$hr_y + rnorm(1, 0, 500)
  } else {
    # Continue from previous year
    prev_path <- all_paths[[i - 1]]
    start_x <- tail(prev_path$x, 1)
    start_y <- tail(prev_path$y, 1)
  }
  
  n <- ind_yr$n_steps
  
  # Pre-allocate
  x <- numeric(n)
  y <- numeric(n)
  x[1] <- start_x
  y[1] <- start_y
  
  # Pre-generate movement parameters (vectorized!)
  steps <- rgamma(n - 1, 
                 shape = (step_mean/step_sd)^2,
                 scale = step_sd^2/step_mean)
  
  angles <- rwrappedcauchy(n - 1, mu = circular(0), rho = turn_concentration)
  
  current_angle <- runif(1, 0, 2*pi)
  
  # Simplified selection: only 5 candidates per step
  for(j in 2:n) {
    
    # Generate candidates
    candidate_angles <- current_angle + angles[j-1] + 
                       seq(-pi/4, pi/4, length.out = 5)
    
    sl <- steps[j-1]
    cand_x <- x[j-1] + sl * cos(candidate_angles)
    cand_y <- y[j-1] + sl * sin(candidate_angles)
    
    # Batch extract habitat (faster than one-by-one)
    hab_vals <- extract(habitat, cbind(cand_x, cand_y))[, 1]
    hab_vals[is.na(hab_vals)] <- -100
    
    # Distance from HR
    dist <- sqrt((cand_x - ind_yr$hr_x)^2 + (cand_y - ind_yr$hr_y)^2)
    
    # Weights
    weights <- exp(habitat_beta * hab_vals + hr_penalty * dist)
    weights <- weights / sum(weights)
    
    # Select
    sel <- sample(1:5, 1, prob = weights)
    x[j] <- cand_x[sel]
    y[j] <- cand_y[sel]
    current_angle <- candidate_angles[sel]
  }
  
  # Create timestamps
  timestamps <- seq(as.POSIXct(ind_yr$start_date),
                   as.POSIXct(ind_yr$end_date + 1),
                   length.out = n)
  
  # Store
  all_paths[[i]] <- tibble(
    animal_id = ind_yr$id,
    timestamp = timestamps,
    x = x,
    y = y,
    year = ind_yr$year
  )
  
  setTxtProgressBar(pb, i)
}

close(pb)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\n\nSimulation complete in", round(runtime, 1), "minutes!\n\n")

# COMBINE AND PROCESS -----------

telemetry <- bind_rows(all_paths)

# Add variables
telemetry <- telemetry %>%
  mutate(
    date = as.Date(timestamp),
    month = month(timestamp),
    yday = yday(timestamp)
  )

## filter to only keep individuals within the study area 

telemetry_sp <- st_as_sf(telemetry,coords=c("x","y"),crs=crs_proj,remove=FALSE)

study_area <- st_read(file.path(data_dir, "study_area.shp"), quiet = TRUE)

sf_use_s2(FALSE)
study_area_buf = st_buffer(study_area,dist=-1000)

# plot(habitat)
# # plot(study_area,alpha=0.2,add=T)
# plot(study_area_buf,add=T,col="yellow",alpha=.9)

telemetry_crop <- st_intersection(telemetry_sp,study_area_buf)

telemetry <- st_drop_geometry(telemetry_crop)

cat("  Processing complete\n\n")

# VISUALIZE ----------------

cat("Creating plots...\n")

sample_ids <- sample(unique(telemetry$animal_id), 5)
sample_data <- telemetry %>% filter(animal_id %in% sample_ids)

p1 <- ggplot() +
  geom_raster(data = as.data.frame(habitat, xy = TRUE),
              aes(x = x, y = y, fill = resource_index)) +
  scale_fill_viridis_c() +
  geom_path(data = sample_data, aes(x = x, y = y, color = animal_id),
            alpha = 0.6) +
  labs(title = "Sample Paths") +
  theme_minimal() +
  coord_equal()

print(p1)

allids=sort(unique(telemetry$animal_id))

# SAVE ----------------

cat("\nSaving outputs...\n")

output_dir <- "telemetry_data"
if(!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE)

write_csv(telemetry, file.path(output_dir, "sim_telemetry_data.csv"))
write_csv(individuals, file.path(output_dir, "sim_individual_info.csv"))

telemetry_sf <- st_as_sf(telemetry, coords = c("x", "y"),
                        crs = st_crs(study_area), remove = FALSE)
st_write(telemetry_sf, file.path(output_dir, "telemetry_data.gpkg"),
        delete_dsn = TRUE, quiet = TRUE)

# ggsave(file.path(output_dir, "paths.png"), p1, width = 10, height = 8)
# ggsave(file.path(output_dir, "selection.png"), p2, width = 8, height = 6)

cat("  Outputs saved to:", output_dir, "\n\n")

cat("Total runtime:", round(runtime, 1), "minutes\n")
cat("Ready for analysis!\n\n")

# END SCRIPT
