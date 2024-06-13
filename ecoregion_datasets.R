library(terra)
library(dplyr)
library(sf)

ecoregions <- vect("C:/Users/lisah/Documents/ArcGIS/Eco_Level1_study_area.shp")
script_dir <- ("C:/Users/lisah/Documents/GEDI Work/data/hotspot/")
data_name_ebird <- paste0(script_dir, "ebird/all_ebird.csv")
ebird <- read.csv(data_name_ebird)
  
#Reproject shapefile to match raster
ebird_proj <- ebird %>% 
  distinct(locality_id, latitude, longitude) %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to projection of raster
  st_transform(crs = crs(ecoregions))

extr <- extract(ecoregions, vect(ebird_proj), progress = TRUE)
ecoregion_locs <- cbind(ebird_proj, extr[,"NA_L1CODE"])
ecoregion_locs_list <- data.frame(locality_id = ecoregion_locs$locality_id, NA_L1CODE = as.numeric(ecoregion_locs$extr....NA_L1CODE..))
ecoregion_locs_list <- na.omit(ecoregion_locs_list)
for(r in 1:nrow(ecoregion_locs_list)){
  if(ecoregion_locs_list$NA_L1CODE[r] == 10){
    ecoregion_locs_list$ecoregion[r] = "nade" 
  } else{
    if(ecoregion_locs_list$NA_L1CODE[r] == 9){
      ecoregion_locs_list$ecoregion[r] = "grpl" 
    } else{
      if(ecoregion_locs_list$NA_L1CODE[r] == 7){
        ecoregion_locs_list$ecoregion[r] = "mwcf" 
      } else{
        if(ecoregion_locs_list$NA_L1CODE[r] == 6){
          ecoregion_locs_list$ecoregion[r] = "nwmf" 
        } 
      }
    }
  }
} 
ebird_ecoregion <- inner_join(ebird, ecoregion_locs_list, by = c("locality_id"))
name_ebird_ecoregion <- paste0(script_dir, "ebird/all_ebird_ecoregion.csv")
#save as all_ebird_ecoregion
write.csv(ebird_ecoregion, name_ebird_ecoregion)
name_ebird_ecoregion
