#SETUP
setwd("C:/Users/lisah/Documents/GEDI Work/data/hotspot")
library(sf)
library(MODIS)
library(exactextractr)
library(viridis)
library(tidyverse)
library(terra)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

#add a shapefile of the study area boundary
#studyArea <- read_sf("gis-data2.gpkg", "PIWO_sf") %>% 
studyArea <- ebird_sf %>% #ebird_sf is generated in the GISDataLayers.R script, and seems to have trouble being saved to the gpkg, so it is easier to just run both scripts
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))

# load ebird data
ebird <- read_csv("ebird/all_ebird.csv")

#set buffer diameter distances (in m)
buffer_diameters <- c(90, 180, 360, 540, 1530, 2520) #(90, 

#
year <- c(2016, 2017, 2018, 2019, 2020)
variables <- c("cover", "rh98", "fhd")
var_df <- data.frame()

for(v in variables){
    print(v)
    pattern <- paste0("^", v, "_[0-9]{4}.tif$")
    
    # load the raster data
    covar <- list.files("D:/spatial data_v2_w", pattern, 
                            full.names = TRUE) %>% 
      stack()     
    
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
    
    # Function that calculates weighted mean covariate value
    calculate_mean <- function(yr, regions, covar) {
      locs <- st_set_geometry(regions, NULL)
      exact_extract(covar[[yr]], regions, progress = FALSE) %>% 
        map(~ summarise(., mean_value = weighted.mean(value, coverage_fraction, na.rm = TRUE))) %>% 
        tibble(locs, data = .) %>% 
        unnest(data)
    }
    
    
    
    # Add a for loop over buffer distances
    for (bd in buffer_diameters) {
      
      message("Working on buffer ", bd)
      
      neighborhood_radius <- bd / 2
      ebird_buff <- ebird %>% 
        distinct(year, locality_id, latitude, longitude) %>%
        mutate(year_lbl = paste0("y", year)) %>%
        # convert to spatial features
        st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
        # transform to modis projection
        st_transform(crs = projection(covar)) %>% 
        # buffer to create neighborhood around each point
        st_buffer(dist = neighborhood_radius) %>% 
        # nest by year
        nest(data = c(year, locality_id, geometry))
      
      
      # iterate over all years extracting landcover for all checklists in each
      var_extract <- ebird_buff %>% 
        mutate(covar_value = map2(year_lbl, data, calculate_mean, covar = covar)) %>% 
        select(covar_value) %>% 
        unnest(cols = covar_value)
      
      print(Sys.time())
      
      # Add column for buffer distance and combine with main df
      var_extract$buffer_distance <- rep(bd)
      var_extract$variable <- rep(v)
      var_df <- rbind(var_df, var_extract)
      
    }
      
}


# Testing to figure out how it works
#test_extract <- exact_extract(covar[[1]], ebird_buff$data[[1]])
#weighted.mean(test_extract[[1]]$value, test_extract[[1]]$coverage_fraction)
#test_map <- test_extract[[1:2]] %>% map(~ count(.x, landcover = value))

# tranform to wide format, filling in implicit missing values with 0s%>% 
wide_df <- var_df %>% 
  group_by(locality_id, year, buffer_distance) %>%
  pivot_wider(names_from = variable, 
              values_from = mean_value)


# save
mname <- paste0("data/gedi_location-year_hotspot.csv")
write_csv(wide_df, mname)



####biovars
library(dismo)
library(raster)
library(sp)
library(rgdal)

#read in PRISM data
prism_path <- "C:/Users/lisah/Documents/PRISM"
biovars_stack <- rast(file.path(prism_path, "biovars.tif"))

# Function that calculates weighted mean covariate value
calculate_mean_allyears <- function(regions, covar) {
  locs <- st_set_geometry(regions, NULL)
  extr <- exact_extract(covar, regions, progress = FALSE)
  vals <- extr %>%
    # pivot data columns to rows
    map(pivot_longer, cols = paste0("biovars_", 1:19)) %>%
    # group by name of data column (raster layer)
    map(group_by, name) %>%
    # summmarize using the weighted mean
    map(~ summarise(., value = weighted.mean(value, 
                                             coverage_fraction, 
                                             na.rm = TRUE))) %>%
    # pivot wider with one column for each data type
    map(pivot_wider) %>%
    tibble(locs, data = .) %>% 
    unnest(data)
}

# Testing
#locs <- st_set_geometry(ebird_buff$data[[1]], NULL)
#col_order <- paste0("biovars_", 1:19)
#test_extr <- exact_extract(biovars_stack, ebird_buff$data[[1]], progress = FALSE)
#test_long <- test_extr[1:10] %>%
#  map(pivot_longer, cols = paste0("biovars_", 1:19))
#test_grp <- test_long[1:10] %>%
#  map(group_by, name)
#test_sum <- test_grp %>%
#  map(~ summarise(., value = weighted.mean(value, coverage_fraction, na.rm = TRUE)))
#test_wide <- test_sum %>%
#  map(pivot_wider) %>%
#  map(relocate, all_of(col_order))
#test_ret <- test_wide %>%
#  tibble(locs[1:10,], data = .) %>% 
#  unnest(data)


#Reproject shapefile to match raster
ebird_proj <- ebird %>% 
  distinct(locality_id, latitude, longitude) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to projection of raster
  st_transform(crs = crs(biovars_stack))
saveRDS(ebird_proj, file = "data/ebird_proj.rds")



extr <- extract(biovars_stack, vect(ebird_proj), progress = TRUE)
biovar_df <- cbind(ebird_proj, extr)
# save
bname <- paste0("data/bioclim_data_800m.csv")
write_csv(biovar_df, bname)
