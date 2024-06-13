#SETUP
setwd("C:/Users/lisah/Documents/GEDI Work/data/hotspot")
library(sf)
library(raster)
library(MODIS)
library(exactextractr)
library(viridis)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# bcr 27 boundary
#studyArea <- read_sf("gis-data2.gpkg", "PIWO_sf") %>% 
studyArea <- ebird_sf %>% #ebird_sf is generated in the GISDataLayers.R script, and seems to have trouble being saved to the gpkg, so it is easier to just run both scripts
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))

# load ebird data
ebird <- read_csv("ebird/all_ebird.csv")


#this file combines modis data with ebird locations
modis16 <- raster("data/modis/modis_mcd12q1_umd_2016.tif")
modis17<- raster("data/modis/modis_mcd12q1_umd_2017.tif")
modis18 <- raster("data/modis/modis_mcd12q1_umd_2018.tif")
modis19 <- raster("data/modis/modis_mcd12q1_umd_2019.tif")
modis20 <- raster("data/modis/modis_mcd12q1_umd_2020.tif")
#m20 <- raster("data/modis/MCD12Q1.A2020001.LC_Type2.tif")

# load the landcover data
landcover <- list.files("data/modis", "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>% 
  stack()

# label layers with year
landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)
landcover




#collapse modis data down to three cover types

rcl_collapsed <- data.frame(from = c(0:15, 255),
                            to = c(NA, 1, 1, 2, 2, 
                                   3, NA, NA, NA, NA,
                                   NA, NA, NA, NA, NA, 
                                   NA, NA))
lc_collapsed <- reclassify(landcover, rcl_collapsed)
#1 = evergreen, 2 = deciduous, 3 = mixed

max_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  max()


#rm(landcover)

#get zonal stats for ea. of the three
#multiply evergreen + 1/2 of mixed per summary distance






# Function that calculates pland
calculate_pland <- function(yr, regions, lc) {
  locs <- st_set_geometry(regions, NULL)
  exact_extract(lc[[yr]], regions, progress = FALSE) %>% 
    map(~ count(., landcover = value, wt = coverage_fraction)) %>% 
    tibble(locs, data = .) %>% 
    unnest(data)
}

buffer_diameters <- c(540, 1530, 2520)
lc_df <- data.frame()
# Add a for loop over buffer distances
for (bd in buffer_diameters) {
  
  message("Working on buffer ", bd)
  
  neighborhood_radius <- bd / 2
  ebird_buff <- ebird %>% 
    distinct(year, locality_id, latitude, longitude) %>%
    mutate(year_lc = paste0("y", year)) %>%
    # convert to spatial features
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    # transform to modis projection
    st_transform(crs = projection(lc_collapsed)) %>% 
    # buffer to create neighborhood around each point
    st_buffer(dist = neighborhood_radius) %>% 
    # nest by year
    nest(data = c(year, locality_id, geometry))
  
  
  # iterate over all years extracting landcover for all checklists in each
  lc_extract <- ebird_buff %>% 
    mutate(pland = map2(year_lc, data, calculate_pland, lc = lc_collapsed)) %>% 
    select(pland) %>% 
    unnest(cols = pland)
  
  print(Sys.time())
  
  # Add column for buffer distance and combine with main df
  lc_extract$buffer_distance <- rep(bd)
  lc_df <- rbind(lc_df, lc_extract)
  
}

# Testing to figure out how it works
#lc_test <- crop(lc_collapsed[[1]])
#test_extract <- exact_extract(lc_collapsed[[1]], ebird_buff$data[[1]])
#count(test_extract[[1]], landcover = value)
#count(test_extract[[1]], landcover = value, wt = coverage_fraction)
#test_map <- test_extract[[1:2]] %>% map(~ count(.x, landcover = value))

pland <- lc_df %>% 
  # calculate proporiton
  group_by(locality_id, year, buffer_distance) %>% #add in buffer distance
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) #%>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

# convert names to be more descriptive
lc_names <- tibble(landcover = 1:3,
                   lc_name = c("evergreen", 
                               "deciduous", 
                               "mixed"))
pland <- pland %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

#mutate evergreen and deciduous to be evergreen + 1/2 mixed (or pivot wider so that these are rows first before doing this)

# tranform to wide format, filling in implicit missing values with 0s%>% 
pland <- pland %>% 
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>%
  mutate(evergreen_combined = evergreen + 0.5 * mixed,
         deciduous_combined = deciduous + 0.5 * mixed)


# save
mname <- paste0("data/modis_pland_location-year_hotspot4.csv")
write_csv(pland, mname)

#####Prediction surface

e <- extent(studyArea)
studyAreaShp <- as(e, 'SpatialPolygons')  #create a polygon shapefile from extent of studyArea points
crs(studyAreaShp) <- paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0", "+a=6371007.181 +b=6371007.181 +units=m +no_defs")
shapefile(studyAreaShp, 'studyAreaShp.shp')
studyAreaShp_sf <- st_as_sf(studyAreaShp)

agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
r <- studyAreaShp_sf %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r, field = 1) %>% 
  # remove any empty cells at edges
  trim()
r <- writeRaster(r, filename = paste0("data/prediction-surface_hotspot.tif"), overwrite = TRUE)

# get cell centers and create neighborhoods
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% 
  transmute(id = row_number())
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)

# extract landcover values within neighborhoods, only needed most recent year
lc_extract_pred <- landcover[[paste0("y", max_lc_year)]] %>% 
  exact_extract(r_cells, progress = FALSE) %>% 
  map(~ count(., landcover = value)) %>% 
  tibble(id = r_cells$id, data = .) %>% 
  unnest(data)


# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>% 
  count(id, landcover) %>% 
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

# convert names to be more descriptive
pland_pred <- pland_pred %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s
pland_pred <- pland_pred %>% 
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>% 
  mutate(year = max_lc_year) %>% 
  select(id, year, everything())

# join in coordinates
pland_coords <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred, by = "id")

for(n in 1:nrow(pland_coords)){
  pland_coords$forest[n] <- sum(pland_coords$pland_01_evergreen_needleleaf[n], 
                                pland_coords$pland_02_evergreen_broadleaf[n], 
                                pland_coords$pland_03_deciduous_needleleaf[n], 
                                pland_coords$pland_04_deciduous_broadleaf[n], 
                                pland_coords$pland_05_mixed_forest[n])
}


forest_cover <- pland_coords %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r, field = "forest") %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

# make a map
par(mar = c(0.25, 0.25, 2, 0.25))
t <- str_glue("Proportion of Forest Cover\n",
              "{max_lc_year} MODIS Landcover")
plot(forest_cover, axes = FALSE, box = FALSE, col = viridis(10), main = t)


##ELEVATION
elev <- raster("data/elevation_1KMmd_GMTEDmd.tif")
# crop, buffer bcr by 10 km to provide a little wiggly room
elev <- studyArea %>% 
  st_buffer(dist = 10000) %>% 
  st_transform(crs = projection(elev)) %>% 
  crop(elev, .) %>% 
  projectRaster(crs = projection(landcover))

# buffer each checklist location
ebird_buff_noyear <- ebird %>% 
  distinct(locality_id, latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(elev)) %>% 
  st_buffer(dist = neighborhood_radius)

# extract elevation values and calculate median and sd
locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% 
  mutate(id = row_number())
elev_checklists <- exact_extract(elev, ebird_buff_noyear, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                   elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
  # join to lookup table to get locality_id
  bind_cols(locs, .)

# extract and calculate median and sd
elev_pred <- exact_extract(elev, r_cells, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = mean(.$value, na.rm = TRUE),
                   elevation_sd = sd(.$value, na.rm = TRUE))) %>% 
  # join to lookup table to get locality_id
  bind_cols(st_drop_geometry(r_cells), .)

# checklist covariates
pland_elev_checklist <- inner_join(pland, elev_checklists, by = "locality_id")
xname <- paste0("data/pland-elev_location-year_", sp, "_hotspot.csv")
write_csv(pland_elev_checklist, xname)

# prediction surface covariates
pland_elev_pred <- inner_join(pland_coords, elev_pred, by = "id")
psname <- paste0("data/pland-elev_prediction-surface_", sp, "_hotspot.csv")
write_csv(pland_elev_pred, psname)
glimpse(pland_elev_pred)
