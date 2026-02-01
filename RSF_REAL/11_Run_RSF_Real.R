
# RSF Analysis (real data from Mojave National Preserve)

# clear workspace ------------

rm(list = ls())

# install packages----------------------------------

library(Matrix)
# library(adehabitatHR)
library(terra)
library(DHARMa)
library(sf)
library(mgcv)
# library(INLA)
library(gratia)
# library(parallel)
# library(VSURF)   #new RFE package
library(amt)
# library(ctmm)
# library(mgcViz)
library(stringr)
# library(ranger)
library(tidyverse)

# set global vars--------------------------

## Set coordinate reference systems ----
crs_proj <- 32611     # NAD 1983 / UTM Zone 11N (Mojave is in Zone 11)
crs_geog <- 4269      # NAD 1983 (Geographic, degrees)
crs_wgs84 <- 4326     # WGS84 for leaflet and elevatr

resolution <- c(250, 250)   # for now: 0.25 km x 0.25 km square pixels 

allyears <- 2009:2016 

timevar_covars_list <- c("AFG","PFG","SHRUB","TREE")


# Load data  -------------

data_dir <- "RSF_REAL/mojave_spatialdata"
study_area <- st_read(file.path(data_dir, "study_area.shp"), quiet = TRUE)
study_area_extent <- ext(vect(study_area))

cov_layers <- rast(file.path(data_dir, "covariate_stack.tif"))

covars_list <- names(cov_layers)

## rasterize study area  -------

temp = rast(study_area_extent, resolution = resolution, crs=crs(vect(study_area)) )
study_area$Id = 1
template_rast <- rasterize(vect(study_area), temp, field = "Id")

# plot(template_rast)

ncell(template_rast)    # 12,648 cells

## Read in GPS collar data ------------------

pts_dir = "RSF_REAL"
realpoints <- readxl::read_excel(file.path(pts_dir, "real_telemetry_data.xlsx"))
names(realpoints)
ndx = complete.cases( realpoints[,c("UTM_X","UTM_Y")])

realpoints <- st_as_sf(realpoints[ndx,],coords=c("UTM_X","UTM_Y"),remove=F,crs=crs_proj)

## remove points outside of study area

library(sf)
realpoints <- st_intersection(realpoints,study_area)

plot(st_geometry(study_area),lwd=2)
plot(st_geometry(realpoints),col="darkgreen",add=T,pch=20,cex=.2)

## convert to data frame--------------------------------

xy = st_coordinates(realpoints)   # extract UTMs
proj = st_crs(realpoints)   # extract CRS

dat <- st_drop_geometry(realpoints)   # ~200k points

length(which(!complete.cases(dat)))  # lots of rows with one or more NA
temp = dat[which(!complete.cases(dat)),]  # datetime is the column with missing data

class(dat$Date)
dat$YEAR = year(dat$Date)
dat = subset(dat,YEAR %in% allyears)  # remove years outside the specified range

unique(dat$Area)
dat = subset(dat,Area=="MIDH")  # retain only the mid hills study area

dat$YDAY = yday(dat$Date)
dat$MONTH = month(dat$Date)
dat$WEEK = week(dat$Date)
dat$HOUR = hour(dat$Time)
dat$MINUTE = minute(dat$Time)

dat <- dat %>%
  mutate(
    DATE = ymd(paste0(YEAR, "-01-01")) + days(YDAY - 1),
    # Add hours and minutes
    DATETIME = DATE + hours(HOUR) + minutes(MINUTE)
  )

names(dat)
dat <- dat |> 
  select(all_of(c("ID", "UTM_Y","UTM_X","AREA", "YEAR","WEEK","MINUTE",    
                "YDAY","MONTH","HOUR","DATE","DATETIME")))

length(which(!complete.cases(dat)))
dat <- na.omit(dat)   # no missing data...


# process data -----------

## make home ranges ----------

## Nest by individual ID

all_ids = sort(unique(dat$ID))
dat$ID = factor(dat$ID,levels=all_ids)  # make ID factor

names(dat)
trk_grouped <- dat %>%
  nest(data = -ID) %>%  # group by individual ID
  mutate(track = map(data, ~ make_track(.x, UTM_X, UTM_Y, DATETIME, crs = crs_proj)))

## Calculate home ranges for each individual --------------------------------------------

hr <- trk_grouped %>%
  mutate(hr = map(track, ~ hr_kde(.x, levels = 0.9, h=c(1000,1000), keep.data=T )))  

i=10
hr$ID[i]
# terra::plot(hr_north$hr[[i]]$ud)
plot(hr$hr[[i]])
class(hr$hr[[i]])
r1 = random_points(hr$hr[[i]], n = 500, presence = hr$hr[[i]]$data)
# table(r1$case_)
# a=hr_north$hr[[i]]$data;View(a)
plot(r1)

#  Convert to sf for plotting

# subset(dat,ID==hr_north$ID[164])

hr_sf <- hr %>%
  mutate(hr_poly = map(hr, hr_isopleths)) %>%
  dplyr::select(ID, hr_poly) %>%
  unnest(cols = hr_poly) %>%
  st_as_sf()  # ensure it is a proper sf object

# # Add a buffer (1500 meters)    
hr_sf <- hr_sf %>%
  mutate(geometry = st_buffer(geometry, dist = 1500))  # distance in meters (if CRS is UTM)
plot(hr_sf$geometry)

# #clip to the study area- there was only small sections out of the study area
hr_clipped <- st_intersection(hr_sf, study_area)

## reduce dataset -------------

dat = st_as_sf(dat,coords = c("UTM_X","UTM_Y"),crs=crs_proj,remove=F)

## only keep points in study area
# temp = st_intersection(dat,study_area)  # takes a long time to run

temp = terra::extract(template_rast,vect(dat))
tokeep = !is.na(temp$Id)

if(length(tokeep)>0) dat = dat[tokeep,]

# names(dat)
dat$WDAY = wday(dat$DATE)   # add numeric day of week

dat <- dat |>      # thin to one or a few obs per ID per week each year
  group_by(ID, YEAR, YDAY) |> 
  # filter(WDAY%in%c(2,5)) |> 
  slice_sample(n=1) |> 
  ungroup()

dat <- dat |> 
  group_by(ID, YEAR, WEEK) |> 
  filter(WDAY %in% c(2,3,5))  |>    # select observations spaced apart
  ungroup()

dat$IDYEAR = paste0(as.character(dat$ID),"_",dat$YEAR)

all_idyrs = sort(unique(dat$IDYEAR))
idyrlist = lapply(all_idyrs, function(t) subset(dat,IDYEAR==t) )
nobs = sapply(idyrlist,nrow)
names(nobs)= all_idyrs
sort(nobs,decreasing = T)    # KTS: some ID/YEARs have less than 20 observations- let's remove some of the ones wiht very little information. I'll try using 200 as a minimum for now...

days <- sapply(idyrlist,function(t) interval(min(t$DATE),max(t$DATE))/days() )
names(days) = all_idyrs
sort(round(days),decreasing = T)   # a few that are only a few days long... let's only keep those with more than 50 days of observations...

# remove individuals with too few data points
tokeep <- intersect(names(days[days>=25]),names(nobs[nobs>=25]))

dat <- subset(dat,IDYEAR%in%tokeep)    # KTS: filter data to only include those with enough data...
nrow(dat)  # 8903 observation remaining

allids = sort(unique(dat$ID))
hr_clipped <- subset(hr_clipped, ID%in%allids)   # remove ids from spatial dataset as well

# table(dat$ID)  # check

all_ids = sort(unique(dat$ID))      # make ID factor
dat$ID = factor(dat$ID,levels=all_ids)
hr_clipped$ID = factor(hr_clipped$ID,levels=all_ids)

ggplot() +
  geom_sf(data = hr_clipped, aes(fill = ID), alpha = 0.5) +
  geom_sf(data = study_area, fill = NA, color = "black", size = 1) +
  geom_sf(data = dat, color = "red", size = 0.8, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Home Ranges")+
  theme(legend.position="none")

## remove points outside the MCP-------------

# Create a list of individual IDs

# Create list of unique IDs
library(sf)

class(dat)
# dat <- st_as_sf(dat, coords = c("X", "Y"), crs = crs_proj)  # or use your correct CRS
# dat <- st_transform(dat, st_crs(hr_clipped))

all_ids <- unique(dat$ID)   
length(all_ids)   

# table(dat$ID)

# Loop over each ID and retain used points within their MCP
usedpoints_hr <- map_dfr(all_ids, function(ID) {
  # Subset MCP and used points for this individual
  kde_i <- hr_clipped %>% filter(ID == !!ID)     # the !! is the 'bang-bang' operator, which 'unquotes' a variable name to refer to an object rather than the object's name
  used_i <- dat %>% filter(ID == !!ID)
  
  # Keep only used points inside the MCP
  used_in_hr <- used_i[st_within(used_i, kde_i, sparse = FALSE)[, 1], ]
  
  return(used_in_hr)
})

# sort(table(usedpoints_hr$ID))
# sort(table(usedpoints_hr$ID))

ggplot() +
  geom_sf(data = hr_clipped, aes(fill = as.factor(ID)), alpha = 0.3, color = "black") +
  geom_sf(data = usedpoints_hr, color = "darkgreen", size = 1, shape = 21, fill = "green") +
  theme_minimal() +
  labs(
    title = "Used Points Within Each Individual's hr",
    fill = "ID"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  ) +
  theme(legend.position = "none")

# rename dat with only points in mcp 

dat <- usedpoints_hr
nrow(dat)   # 8882 used points

all_ids <- sort(unique(dat$ID)) 
length(all_ids)  # 67 individuals

all_yrs <- sort(unique(dat$YEAR))
length(all_yrs)   # 8 years

ggplot() +
  geom_sf(data = hr_clipped, aes(fill = ID), alpha = 0.5) +
  geom_sf(data = study_area, fill = NA, color = "black", size = 1) +
  geom_sf(data = usedpoints_hr, color = "red", size = 0.8, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Home Ranges") +
  theme(legend.position = 'none')

##  process covariate rasters ---------

plot(cov_layers$elevation)
covars_list

### reproject all covariate rasters to the template raster dimension ----------------0

res(cov_layers)
res(template_rast)

thisstk = terra::project(cov_layers[[covars_list]],template_rast)   
res(thisstk)
covars_df = as.data.frame(thisstk, xy=T, cells = TRUE, na.rm = TRUE)

# 7661 observations...  with years will be 61,288 observations. Hopefully won't take too long to run

# data_dir
# allyears
dat_byyear = lapply(allyears,function(t) subset(dat,YEAR==t) )
names(dat_byyear) = allyears

veg_byyear = lapply(allyears, function(t) rast(sprintf("%s/veg_%s.tif",data_dir,t)) )
names(veg_byyear) = allyears

veg_byyear = lapply(veg_byyear, function(t) terra::project(t,template_rast) )
sapply(veg_byyear,res)

# make master data frame for analysis ------------

hr_clipped$hectares <- hr_clipped$area/10000

## make raster of point density (response) --------------

# terra::plot(template_rast)

## Convert sf points to SpatVector
points_vect <- lapply(allyears, function(t)  vect(dat_byyear[[as.character(t)]]) )
names(points_vect) = allyears

lapply(points_vect,plot)  # plot out points for each year

## Count telemetry points per cell  ---------

point_counts <- lapply(allyears, function(t){
  (template_rast-1) + rasterize(points_vect[[as.character(t)]], template_rast, 
                                fun = "count",background=0)
} )
names(point_counts) = allyears

terra::plot(point_counts$`2009`)   # plot count raster (response variable)

lapply(point_counts, plot)  # plot out count raster for all years

## Create master data frame for analysis --------------

master_df=NULL

if(any(timevar_covars_list %in% names(covars_df))){
  covars_df = covars_df |>    # remove time varying covars from covars_df
    select(!all_of(timevar_covars_list))
}

names(covars_df)

y=1
for(y in 1:length(allyears)){
  thisy = as.character(allyears[y])
  thisdf = as.data.frame(point_counts[[thisy]], cells = TRUE, na.rm = TRUE, xy=TRUE)
  names(thisdf) = c("cell","x","y","n_points")
  thisdf$year = allyears[y]
  thisdf = left_join(thisdf,covars_df,by=c("cell","x","y"))
  thisveg = veg_byyear[[thisy]] 
  thisvegdf = as.data.frame(thisveg,cells = TRUE, na.rm = TRUE, xy=TRUE)
  thisdf = left_join(thisdf,thisvegdf,by=c("cell","x","y"))
  master_df = rbind(master_df,thisdf)
}

nrow(master_df)  # about 60k observations...

## make raster of survey effort (exposure) -----------

poly_vect = vect(hr_clipped)
names(poly_vect)
# poly_vect$ID = as.numeric(as.character(poly_vect$ID))

# Extract: get which cells each polygon intersects
    # this returns a data frame with one row for each grid cell intersected by each hr polygon

extracted <- terra::extract(template_rast, poly_vect, cells = TRUE, ID = TRUE)

# Reverse grouping: for each cell, list all associated polygon IDs
cells_with_polygons <- extracted %>%
  group_by(cell) %>%
  summarise(
    polygon_ids = list(ID),  # List of all polygon IDs
    n_polygons = n()          # Count of polygons
  )

# exposure: sum across individuals overlapping this cell: 
#  specifically: for each cell, the total weeks of exposure per square km of home range  
#    takes a minute to run...

#    since we want to run a spatiotemporal model, we adjust effort based on the year. Some cells have more 'exposure' in one year vs another


all_idyrs
idyrlist = lapply(all_idyrs, function(t) subset(dat,IDYEAR==t) )
names(idyrlist) = all_idyrs
sapply(idyrlist,nrow)
all_idyrs = all_idyrs[sapply(idyrlist,nrow)>0]
idyrlist = idyrlist[all_idyrs]
days <- sapply(idyrlist,function(t) interval(min(t$DATE),max(t$DATE))/days() )
names(days) = all_idyrs

ids_byyear = lapply(allyears,function(t){temp = subset(dat,YEAR==t) ; levels(droplevels(temp$ID)) }  )
names(ids_byyear) = allyears
ids_byyear$`2009`

master_df$exposure = 0.0    # initialize exposure variables (total exposure and number of ids)
master_df$nids = 0.0       # not currently used- authors of sdmTMB are working on implementing a 'dispformula' and this would be used for that purpose... 
r=57
for(r in 1:nrow(cells_with_polygons)){
  cl = cells_with_polygons$cell[r]    # cell id
  pids = cells_with_polygons$polygon_ids[r][[1]]
  p=1
  for(p in 1:length(pids)){
    thispid = pids[p]
    thisid = as.character(hr_clipped$ID[thispid])  # individual id with range overlap
    thishr = hr_clipped[thispid,]$geometry
    thisarea_sqkm = as.numeric(hr_clipped[thispid,]$area)/1000000 
    y=2
    for(y in 1:length(allyears)){
      thisy = as.character(allyears[y])
      if(thisid %in% ids_byyear[[thisy]]){
        thisr = which(master_df$cell==cl & master_df$year == allyears[y])
        master_df$nids[thisr] = master_df$nids[thisr] + 1
        thisidyr = paste0(thisid,"_",thisy)
        thistime_weeks = as.numeric(days[thisidyr])/7
        # plot(thishr)
        master_df$exposure[thisr] <- master_df$exposure[thisr] + thistime_weeks/thisarea_sqkm
      }
    }
  }
}    # note: takes a couple minutes to run- I'm pretty sure I could make this more efficient but I don't have the time

## visualize the response variable ------

hist(master_df$n_points)   # lots of zeros- zero inflated model may be appropriate...

## save master df version 1 ----------

hist(master_df$exposure)

master_df <- subset(master_df, exposure>0.000)   # only keep cells 'exposed' to at least one individual for one year
     # about 30k observations total

master_df$offset = log(master_df$exposure)     # offset terms for count model

thisfile = "RSF_REAL/master_count_df_250m.csv"
write_csv(master_df,thisfile)

# process and scale covariates ----------

#  note: can start here, if you first read in master count df.

master_df <- read.csv(thisfile)

names(master_df)
varstouse = covars_list

## a priori variable thinning --------

varstouse  # remove sin aspect and raw aspect
varstouse = setdiff(varstouse,c("aspect","aspect_sin"))

## correlation matrix -------

library(corrplot)
library(ggcorrplot)
library(RColorBrewer) 

# Compute the correlation matrix
names(master_df)
corr_matrix <- cor(master_df[,varstouse]) 

corrplot(corr_matrix, 
         method = "circle",      # Use circles (other options: "square", "number", "pie", "shade", "color")
         type = "upper",         # Display only the upper triangle
         order = "hclust",       # Reorder variables using hierarchical clustering
         addrect = 2,            # Add rectangles around hierarchical clusters
         col = brewer.pal(n = 8, name = "RdBu"), # Use a diverging color palette
         tl.col = "black",       # Color of text labels 
         tl.srt = 45,            # Rotate text labels by 45 degrees
         p.mat = cor_pmat(master_df[,varstouse]), # Add p-values (requires ggcorrplot package for this function)
         sig.level = 0.01,       # Mark insignificant correlations (e.g., with an 'X')
         insig = "blank"         # Hide insignificant correlations
)

# alternative version
corrplot(corr_matrix, 
         method = "color",
         type = "upper",           # Show only upper triangle
         order = "hclust",         # Cluster similar variables
         addCoef.col = "black",    # Add correlation coefficients
         tl.col = "black",         # Text label color
         tl.srt = 45,              # Text label rotation
         diag = FALSE)             # Hide diagonal


    # remove roughness- too correlated with slope
varstouse = setdiff(varstouse,c("roughness"))

varstouse2 = c("ELEV_s","SLOPE_s","ASPECT_s","TPI_s","AFG_s","PFG_s","SHRUB_s","TREE_s","WATER_s")

##  transform, scale and center covariates  

transforms <- rep("none",length(varstouse))  # store transformations (for backtransformations)
names(transforms) <- varstouse2

offsets <- rep(0,length(varstouse))   # store offsets (for backtransformations)
names(offsets) <- varstouse2

hist(master_df$elevation)
master_df$ELEV_s <- scale(master_df$elevation)[,1]
hist(master_df$ELEV_s)   # okay untransformed

hist(master_df$tpi)
master_df$TPI_s <- scale(master_df$tpi)[,1]
master_df$TPI_s[master_df$TPI_s>4] = 4
master_df$TPI_s[master_df$TPI_s<(-4)] = -4
hist(master_df$TPI_s)

hist(master_df$slope)
master_df$SLOPE_s = scale(log(master_df$slope+1))[,1]
hist(master_df$SLOPE_s)
transforms["SLOPE_s"] <- "log"
offsets["SLOPE_s"] <- 1

hist(master_df$TREE)
master_df$TREE_s = scale(log(master_df$TREE+.5))[,1]
hist(master_df$TREE_s)
transforms["TREE_s"] <- "log"
offsets["TREE_s"] <- 0.5

hist(master_df$SHRUB)
master_df$SHRUB_s = scale(master_df$SHRUB)[,1]
hist(master_df$SHRUB_s)

hist(master_df$PFG)
master_df$PFG_s = scale(master_df$PFG)[,1]
hist(master_df$PFG_s)

hist(master_df$AFG)
master_df$AFG_s = scale(master_df$AFG)[,1]
hist(master_df$AFG_s)

hist(master_df$aspect_cos)
master_df$ASPECT_s = scale(master_df$aspect_cos)[,1]
hist(master_df$ASPECT_s)

hist(master_df$water)
master_df$WATER_s = scale(log(master_df$water+.1))[,1]
hist(master_df$WATER_s)
transforms["WATER_s"] <- "log"
offsets["WATER_s"] <- 0.1

### make data frame for backtransformations etc

means <- sapply(1:length(varstouse), function(t) {
  if(transforms[t]=="log"){
    mean(log(master_df[[varstouse[t]]]+ offsets[t]) ,na.rm=T)
  }else{
    mean(master_df[[varstouse[t]]],na.rm=T)
  }
})
names(means) =varstouse2

sds <- sapply(1:length(varstouse), function(t) {
  if(transforms[t]=="log"){
    sd(log(master_df[[varstouse[t]]]+ offsets[t]) ,na.rm=T)
  }else{
    sd(master_df[[varstouse[t]]],na.rm=T)
  }
})
names(sds) =varstouse2

transform_df = data.frame(
  origvar = varstouse,
  modvar = varstouse2,
  transform = transforms,
  offset = offsets,
  mean = means,
  sd = sds,
  min = sapply(varstouse,function(t) min(master_df[[t]],na.rm=T) ),   # NOTE: this is on original scale...
  max = sapply(varstouse,function(t) max(master_df[[t]],na.rm=T) ),
  min2 = sapply(varstouse2,function(t) min(master_df[[t]],na.rm=T) ),   # NOTE: this is on modified scale...
  max2 = sapply(varstouse2,function(t) max(master_df[[t]],na.rm=T) )
  
)


## save data frames -----------

write_csv(master_df,"RSF_REAL/master_count_df_250m.csv")
write_csv(transform_df,"RSF_REAL/transform_df_count_250m.csv")


# fit models -----------

## fit with gam/mgcv (nonspatial)  ---------

length(which(!complete.cases(master_df)))  # 0 rows with one or more NA
# master_df <- na.omit(master_df)   # ~ 30,000 cells

names(master_df)

fm = formula(n_points ~ 
             WATER_s + ASPECT_s + TREE_s + SLOPE_s + PFG_s
             + s(ELEV_s, bs = "cr",k=5) 
             + offset(offset)  ) 

#

mod_gam = gam(
  fm,
  data = master_df,
  family = nb(link = "log"),
  method = "REML"
)

summary(mod_gam)   # nonspatial GAM

# plot(mod_gam)
summary(mod_gam)
library(gratia)
# draw(mod_gam)  # takes too long

## model checking --------

gam.check(mod_gam)  # may need more knots for elevation

gratia::appraise(mod_gam)   # okay

library(DHARMa)
r = DHARMa::simulateResiduals(mod_gam)   # check with DHARMa
DHARMa::testResiduals(r)   # looks surprisingly good
DHARMa::testZeroInflation(r)  # no zero inflation detected

hist(master_df$n_points)

library(ggeffects)

pp = ggpredict(mod_gam, c("SLOPE_s [all]") )
plot(pp)

pp = ggpredict(mod_gam, c("WATER_s [all]") )
plot(pp)

### assess spatial correlation of residuals -------

library(spdep)   # note: might want to install spDataLarge as well 
library(sf)

# Assuming you have:
# gam_model - your fitted GAM
# spatial_data - sf object with coordinates

# Extract residuals
residuals <- residuals(mod_gam, type = "response")

master_df$res <- residuals

ggplot(master_df,aes(x,y,col=res)) +   # plot residuals spatially
  geom_point()

master_df_sp  = st_as_sf(master_df,coords = c("x","y"),remove = F )

# Get coordinates
coords <- st_coordinates(master_df_sp)

# Create spatial weights (k-nearest neighbors)
knn <- knearneigh(coords, k = 8)
nb <- knn2nb(knn)
weights <- nb2listw(nb, style = "W")

# Moran's I test
moran_test <- moran.test(residuals, weights)
print(moran_test)     # lots of spatial autocorrelation!

# Monte Carlo simulation for significance
# moran_mc <- moran.mc(residuals, weights, nsim = 999)
# print(moran_mc)   # tests are signicant. Morans I of 0.37 indicates substantial autocorrelation

moran_plot <- moran.plot(residuals, weights, 
                         main = "Moran's I Plot for GAM Residuals")


## Calculate correlogram (Moran's I at different lags) (takes some time)
correlogram <- sp.correlogram(nb, residuals, 
                              order = 8,          # Number of lags
                              method = "I",        # Moran's I
                              zero.policy = TRUE)

# Base R plot
plot(correlogram, main = "Spatial Correlogram")

## run GAM with spatial smoother to reduce autocorrelation  --------

names(master_df)

# tokeep = master_df$cell %% 2 == 0
# 
# master_df2 = master_df[tokeep,]

mod_spgam = update(mod_gam, .~.+s(x,y,k=100))   # takes a couple minutes...

gratia::appraise(mod_spgam)   # okay

summary(mod_spgam)


# plot(mod_spgam)
# draw(mod_spgam)


### model checking --------
gam.check(mod_spgam)

k.check(mod_spgam)

r = DHARMa::simulateResiduals(mod_spgam)   # check with DHARMa
DHARMa::testResiduals(r)   # worse fit with spatial smooth?
DHARMa::testZeroInflation(r)  # not zero inflated

### assess spatial correlation of residuals -------

# Extract residuals
residuals <- residuals(mod_spgam, type = "pearson")
master_df$res <- residuals

ggplot(master_df,aes(x,y,col=log(abs(res+0.01)) )) +   # plot residuals spatially
  geom_point()    # better

# # Moran's I test
# moran_test <- moran.test(residuals, weights)  # still nearly as high as before
# print(moran_test)
# 
# # try plotting out residuals?
# 
# # Monte Carlo simulation for significance
# # moran_mc <- moran.mc(residuals, weights, nsim = 999)
# # print(moran_mc)   # tests are significant. Morans I of 0.37 indicates substantial autocorrelation
# 
# moran_plot <- moran.plot(residuals, weights, 
#                          main = "Moran's I Plot for GAM Residuals")
# 
# 
# ## Calculate correlogram (Moran's I at different lags)
# correlogram <- sp.correlogram(nb, residuals, 
#                               order = 10,          # Number of lags
#                               method = "I",        # Moran's I
#                               zero.policy = TRUE)
# 
# plot(correlogram, main = "Spatial Correlogram")  # ### looks better -------

## make raster surface for predicted use intensity  ----------

# pick a year 
thisyr = 2010

nrow(covars_df)
names(covars_df)

thisveg = veg_byyear[[as.character(thisyr)]]
thisvegdf = as.data.frame(thisveg,xy=T,cells=T)
nrow(thisvegdf)
covars_df2 = left_join(covars_df,thisvegdf,by = c("x","y","cell"))

t = 6
for(t in 1:nrow(transform_df)){
  thisvar=transform_df$origvar[t]
  thisvar2=transform_df$modvar[t]
  if(is.null(covars_df2[[thisvar2]])){
    if(transform_df$transform[t]=="log"){
      covars_df2[[thisvar2]] = log(covars_df2[[thisvar]]+transform_df$offset[t])
      covars_df2[[thisvar2]] = (covars_df2[[thisvar2]] - transform_df$mean[t])/transform_df$sd[t]
    } else{
      covars_df2[[thisvar2]] = (covars_df2[[thisvar]] - transform_df$mean[t])/transform_df$sd[t]
    }
  }
  covars_df2[[thisvar2]] = pmin(pmax(covars_df2[[thisvar2]],transform_df$min2[t]),transform_df$max2[t])
  
}

thismod= mod_spgam

covars_df2$offset=0

pred = predict(thismod,newdata = covars_df2, type="response") 
preddf = as.data.frame(cbind(cell=covars_df$cell,pred=pred))

pred_rast = template_rast
pred_rast[preddf$cell] = preddf$pred

# hist(preddf$pred)

terra::plot(pred_rast)

my_breaks <- c(-Inf,0.15,0.2, 0.3, 0.5, 0.75, 1, Inf) # 5 bins
my_colors_classified <- viridis::viridis(n = 5, option = "D")

r = classify(pred_rast,my_breaks )
terra::plot(r , col = my_colors_classified, main = "Classified Raster Plot")

# plot univariate partial effects plots

pred_elev <- ggpredict(mod_spgam, terms = "ELEV_s [all]")
plot(pred_elev) +
  labs(title = "Effect of elev on intensity",
       x = "Elevation", y = "Predicted Count")

pred_slope <- ggpredict(mod_spgam, terms = "SLOPE_s [all]")
plot(pred_slope) +
  labs(title = "Effect of slope on intensity",
       x = "slope", y = "Predicted Count")

pred <- ggpredict(mod_spgam, terms = "WATER_s [all]")
plot(pred) +
  labs(title = "Effect of water on intensity",
       x = "water", y = "Predicted Count")

pred <- ggpredict(mod_spgam, terms = "TREE_s [all]")
plot(pred) +
  labs(title = "Effect of tree on intensity",
       x = "tree", y = "Predicted Count")

# pred_pfg <- ggpredict(mod_spgam, terms = "PFG_s [all]")
# plot(pred_pfg) +
#   labs(title = "Effect of pfg on intensity",
#        x = "rd", y = "Predicted Count")

summary(mod_spgam)

## use sdmTMB to build more sophisticated spatial regression models -------

library(sdmTMB)

# Scale coordinates (recommended for GP models)
master_df$x_s <- master_df$x/1000
master_df$y_s <- master_df$y/1000

mesh <- make_mesh(master_df, xy_cols = c("x_s", "y_s"), cutoff = 1)
mesh$mesh$n   # compare mesh nodes to number of observations (nodes should be quite a bit smaller!)

names(master_df)
library(inlabru)
ggplot() +
  inlabru::gg(mesh$mesh) 
  # geom_point(data=master_df, aes(x_s,y_s,col=n_points),size=0.4)


plot(mesh)  # simple mesh plot

varstouse

# takes about a minute to fit (which is pretty lightning fast actually)
fit <- sdmTMB(
  n_points ~ poly(ELEV_s,2) + WATER_s + SLOPE_s + TREE_s + SHRUB_s + PFG_s ,
  data = master_df,
  mesh = mesh,
  offset = "offset",
  family = nbinom2(),   
  time = "year"   # include spatiotemporal and spatial random field
)

sanity(fit)

summary(fit)

tidy(fit, effects = "fixed",conf.int = T)
tidy(fit, effects = "ran_pars",conf.int = T)

x_range <- range(master_df$x_s)
y_range <- range(master_df$y_s)


p <- predict(fit, type="link")  

p$x = p$x_s
p$y = p$y_s

# p <- predict(fit)

names(p)
ggplot(p, aes(x, y, z = omega_s)) + geom_contour_filled() + facet_wrap(~year)
ggplot(p, aes(x, y, z = epsilon_st)) + geom_contour_filled() + facet_wrap(~year)
ggplot(p, aes(x, y, z = est_non_rf)) + geom_contour_filled() + facet_wrap(~year)
ggplot(p, aes(x, y, z = est_rf)) + geom_contour_filled() + facet_wrap(~year)
ggplot(p, aes(x, y, z = est)) + geom_contour_filled() + facet_wrap(~year)
# scale_fill_viridis_c(trans = "sqrt")


# basic residual tests 

p$resids <- residuals(fit) # randomized quantile residuals
hist(p$resids)

qqnorm(p$resids);abline(a = 0, b = 1)    # looks good

  # check spatial correlation in residuals

ggplot(p, aes(x, y, col = resids)) + scale_colour_gradient2() + facet_wrap(~year) +   # looks pretty uncorrelated
  geom_point()

master_df$res = p$resids 
ggplot(master_df,aes(x,y,col=res)) +   # plot residuals spatially
  geom_point()    # better

# Moran's I test
moran_test <- moran.test(p$resids, weights)  
print(moran_test)   # much better...

# try plotting out residuals?

# Monte Carlo simulation for significance
# moran_mc <- moran.mc(residuals, weights, nsim = 999)
# print(moran_mc)   # tests are significant. Morans I of 0.37 indicates substantial autocorrelation

moran_plot <- moran.plot(p$resids, weights, 
                         main = "Moran's I Plot for Residuals")


## Calculate correlogram (Moran's I at different lags)
correlogram <- sp.correlogram(nb, p$resids, 
                              order = 10,          # Number of lags
                              method = "I",        # Moran's I
                              zero.policy = TRUE)

plot(correlogram, main = "Spatial Correlogram")  # ### looks better -------


# DHARMa tests 

s = simulate(fit, nsim = 500, type = "mle-mvn")
dharma_residuals(s,fit, test_dispersion = T)
res = dharma_residuals(s, fit, return_DHARMa = T)  

testResiduals(res)   # not passing all tests

# Visualizing a marginal effect ----------------------------------------

# See the visreg package or the ggeffects::ggeffect() or
# ggeffects::ggpredict() functions
# To do this manually:

nd <- data.frame(
  TREE_s = seq(min(master_df$TREE_s), max(master_df$TREE_s), length.out = 50),
  ELEV_s = 0, PFG_s = 0, SLOPE_s = 0, WATER_s = 0, year = 2010, SHRUB_s =0
)
nd$offset=0

p <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA,offset=nd$offset)
names(p)
ggplot(p, aes(TREE_s, est,
              ymin = est - 1.96 * est_se, ymax = est + 1.96 * est_se)) +
  geom_line() + geom_ribbon(alpha = 0.4)

nd <- data.frame(
  SHRUB_s = seq(min(master_df$SHRUB_s), max(master_df$SHRUB_s), length.out = 50),
  ELEV_s = 0, PFG_s = 0, TREE_s = 0, WATER_s = 0, year = 2010, SLOPE_s =0
)
nd$offset=0

p <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA,offset=nd$offset)
names(p)
ggplot(p, aes(SHRUB_s, est,
              ymin = est - 1.96 * est_se, ymax = est + 1.96 * est_se)) +
  geom_line() + geom_ribbon(alpha = 0.4)

fit

 ### use 'visreg' package

visreg(fit,"ELEV_s")    # don't run all- can take a little while..
visreg(fit,"WATER_s")
visreg(fit,"SHRUB_s")
visreg(fit,"TREE_s")
visreg(fit,"PFG_s")

##  Predict to the entire study area ----------

# ?predict.sdmTMB

   # pick a year 

covars_df3 = NULL
y=1
for(y in 1:length(allyears)){
  thisyr = allyears[y]
  thisveg = veg_byyear[[as.character(thisyr)]]
  thisvegdf = as.data.frame(thisveg,xy=T,cells=T)
  covars_df2 = left_join(covars_df,thisvegdf,by = c("x","y","cell"))
  
  t = 6
  for(t in 1:nrow(transform_df)){
    thisvar=transform_df$origvar[t]
    thisvar2=transform_df$modvar[t]
    if(is.null(covars_df2[[thisvar2]])){
      if(transform_df$transform[t]=="log"){
        covars_df2[[thisvar2]] = log(covars_df2[[thisvar]]+transform_df$offset[t])
        covars_df2[[thisvar2]] = (covars_df2[[thisvar2]] - transform_df$mean[t])/transform_df$sd[t]
      } else{
        covars_df2[[thisvar2]] = (covars_df2[[thisvar]] - transform_df$mean[t])/transform_df$sd[t]
      }
    }
    covars_df2[[thisvar2]] = pmin(pmax(covars_df2[[thisvar2]],transform_df$min2[t]),transform_df$max2[t])
    
  }
  
  covars_df2$x_s = covars_df2$x / 1000
  covars_df2$y_s = covars_df2$y / 1000
  covars_df2$offset = 0
  covars_df2$year = thisyr
  
  covars_df3 <- rbind(covars_df3,covars_df2)
  
}
nrow(covars_df3)
p2 = predict(fit,newdata=covars_df3,type="link",offset=rep(0,nrow(covars_df3)) )
ggplot(p2, aes(x, y, z = est)) + geom_contour_filled() + facet_wrap(~year)

## TODO: spatial residuals, other checks, etc.
## TODO: plot_smooth for visualizing smooths. plot random field, linear predictor,
## TODO: incorporate time-varying effect (maybe WATER_s?)
## TODO: incorporate delta/hurdle model with different formulas for each component
## TODO: incorporate spatially varying coefficient?


peren_water <- st_read("RSF_REAL/Mojave_originalfiles/Perennial_water.shp")
peren_water_proj <- st_transform(peren_water,crs_proj)
study_area <- st_read("RSF_REAL/Mojave_originalfiles/Midhills.shp")
study_area_proj <- st_transform(study_area,crs_proj)

peren_water_clip = st_intersection(peren_water_proj,study_area_proj)

p3= subset(p2,year==2015)
ggplot() + geom_contour_filled(data=p3, aes(x, y, z = est)) +
  geom_sf(data = peren_water_clip, aes(fill = NA), alpha = 0.5)


# NOTE: omega is spatial random field, epsilon is spatiotemporal random effects


# cross validation --------------

    # to run in parallel...
# library(future)
# plan(multisession, workers = 2)


## TODO: plot observed (y) vs cross validation prediction (x) to look for one to one correspondence

fit_cv <- sdmTMB_cv(
  n_points ~ poly(ELEV_s,2) + WATER_s + SLOPE_s + TREE_s + SHRUB_s + PFG_s ,
  data = master_df,
  mesh = mesh,
  offset = "offset",
  family = nbinom2(),   
  time = "year",   # include spatiotemporal and spatial random field
  k_folds = 4
)
closeAllConnections()

p = predict(fit,nsim=25,type="response",model=NA,offset=master_df$offset)

# ?sdmTMB::sdmTMB_cv

# Sum of log likelihoods of left-out data: equivalent to ELPD - useful for model comparison!
fit_cv$sum_loglik
fit_cv$fold_loglik  # log likelihood of each fold

# RMSE across entire dataset:
sqrt(mean((master_df$n_points - apply(p,1,mean) )^2))    # fitted
sqrt(mean((fit_cv$data$n_points - fit_cv$data$cv_predicted )^2))    # CV alternative




# END SCRIPT 



