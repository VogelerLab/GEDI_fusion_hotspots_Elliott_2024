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
studyArea <- ebird_sf %>% #ebird_sf is generated in the GISDataLayers_hotspot.R script, and seems to have trouble being saved to the gpkg, so it is easier to just run both scripts
 st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
# load ebird data
ebird <- read_csv("ebird/all_ebird.csv")
# get list of tiles required to cover this bcr
tiles <- getTile(studyArea)
tiles@tile
#> [1]"h09v05" "h10v05" "h08v04" "h09v04" "h10v04" "h11v04""


#MODIS setup
MODIS::EarthdataLogin(usr = "lhelliott", pwd = "3DGi6Z5JybCnVZ3")
MODIS:::checkTools("GDAL")


#
# earliest year of ebird data
begin_year <- paste0(min(ebird$year), ".01.01")
# end date for ebird data
end_year <- paste0(max(ebird$year), ".12.31")
# download tiles and combine into a single raster for each year
tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01", 
                extent = studyArea %>% st_buffer(dist = 10000), 
                begin = begin_year, end = end_year, 
                outDirPath = "data", job = "modis",
                MODISserverOrder = "LAADS", config = httr::config(connecttimeout = 60)) %>% 
  pluck("MCD12Q1.006") %>% 
  unlist()

# rename tifs to have more descriptive names
new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)
file.rename(tifs, new_names)

