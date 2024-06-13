

setwd("C:/Users/lisah/Documents/GEDI Work/data/hotspot")
library(sf)
library(exactextractr)
library(tidyverse)
library(terra)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

#add a shapefile of the study area boundary
#studyArea <- read_sf("gis-data2.gpkg", "PIWO_sf") %>% 
studyArea <- read_sf(dsn = "D:/PRISM/study_area_shp_5070", layer = "study_area_shp_5070")

er <- "MarineWCForest"
ecoRegion <- read_sf(dsn = "C:/Users/lisah/Documents/ArcGIS/Ecoregions", layer = er)

year <- c(2016, 2017, 2018, 2019, 2020)
variables <- c("fhd")
var_df <- data.frame()

for(v in variables){
  print(v)
  pattern <- paste0("^", v, "_[0-9]{4}.tif$")
  #pattern <- paste0("^fhd", "_2020", ".tif$")
  
  # load the raster data
  covar_files <- list.files("D:/spatial data_v2_w", pattern, full.names = TRUE)
  covar <- rast(covar_files)     
  
  # label layers with year
  covar <- names(covar) %>% 
    str_extract(paste0("(?<=", v,"_)[0-9]{4}")) %>% 
    paste0("y", .) %>% 
    setNames(covar, .)
  covar
  
  max_lc_year <- names(covar) %>% 
    str_extract("[0-9]{4}") %>% 
    as.integer() %>% 
    max()
}

#set buffer diameter distances (in m)
bd <- c(90)

neighborhood_radius <- bd / 2


agg_factor <- round(2 * neighborhood_radius / res(covar))
r <- rast(covar) %>% 
 aggregate(agg_factor) 
sa_poly <- ecoRegion %>% 
  st_transform(crs = crs(r)) 
v <- vect(sa_poly)

rsa <- rasterize(v, r, field = 1) %>% 
   #remove any empty cells at edges
 trim()
fn <- paste0("data/prediction-surface_", er,".tif")
r <- writeRaster(rsa, filename = fn, overwrite = TRUE)

#r <- raster("data/prediction-surface.tif")
library(raster)
r <- raster(fn)
# get cell centers and create neighborhoods
r_points <- rasterToPoints(r)
r_df <- as.data.frame(r_points)
write.csv(r_df, "data/r_df.csv")
r_sf <- st_as_sf(r_df, coords = c("x", "y"), crs = crs(r))

r_centers <- transmute(r_sf, id = row_number())
write.csv(r_centers, "data/r_centers.csv")
grp1 <- c(1:1000000)
grp2 <- c(10000001:2000000)
grp3 <- c(2000001:3000000)
grp4 <- c(3000001:4000000)
grp5 <- c(4000001:5000000)
grp6 <- c(5000001:6000000)
grp7 <- c(6000001:7000000)
grp8 <- c(7000001:8000000)
grp9 <- c(8000001:nrow(r_centers))

r_cells <- data.frame()
for(n in 1:9){
  grp <- paste0("grp", n)
  print(grp)
  grp_sub <- get(grp)
  r_sub <- subset(r_centers, r_centers$id %in% grp_sub)
  r_centers_sub <- st_buffer(r_sub, dist = neighborhood_radius)
  rcells <- rbind(r_cells, r_centers_sub)
}

r_cells <- st_buffer(r_centers, dist = neighborhood_radius)

# extract landcover values within neighborhoods, only needed most recent year
lc_extract_pred <- covar[[paste0("y", max_lc_year)]] %>% 
  exact_extract(r_cells, progress = FALSE) %>% 
  map(~ count(., covar = value)) %>% 
  tibble(id = r_cells$id, data = .) %>% 
  unnest(data)

# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>% 
  
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(covar))

# convert names to be more descriptive
pland_pred <- pland_pred %>% 
  inner_join(lc_names, by = "covar") %>% 
  arrange(covar) %>% 
  select(-covar)

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