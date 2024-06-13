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
rasterize <- terra::rasterize


#set buffer diameter distances (in m)
buffer_diameters <- c(180)#, 360, 540, 1530, 2520) #(90, 

#
year <- c(2020)
variables <- c("fhd") #"cover", "rh98") #"fhd"
var_df <- data.frame()

for(v in variables){
    print(v)
    
    # load the raster data
    #lc_extract_pred <- readRDS("data/lc_extract_pred_bd90_fhd.rds")
    covar_30m_fn <-  paste0("D:/spatial data_v2_w/", v, "_2020.tif") 
    covar_30m_rst <- rast(covar_30m_fn)
      
    # Create 90m rast
    covar_90m_fn <- paste0("data/", v, "_90m.tif")
    if (file.exists(covar_90m_fn)) {
      covar_90m_rst <- rast(covar_90m_fn)
    } else {
      covar_90m_rst <- aggregate(covar_30m_rst, fact = 3, 
                                 fun = "mean", na.rm = TRUE,
                                 filename = covar_90m_fn)
    }
    
    # Add a for loop over buffer distances
    for (bd in buffer_diameters) {
      
      message("Working on buffer ", bd)
      fcl_rst_fn <- paste0("data/", v, "_", bd, "m.tif")
      #if (file.exists(fcl_rst_fn)) {
        #message("File already created. Skipping...")
       # next
      #}
      
      fcl_mat <- focalMat(covar_90m_rst, bd / 2, type = "circle")
      fcl_rst <- focal(covar_90m_rst, fcl_mat, 
                       filename = fcl_rst_fn, overwrite = TRUE)
      
      print(Sys.time())
      
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

covar_90m_rst <- rast("data/fhd_90m.tif")


biovars_90m <- resample(biovars_stack, covar_90m_rst, filename = "data/biovars_90.tif")

###modis
modis20 <- rast("data/modis/modis_mcd12q1_umd_2020.tif")
modis_proj <- terra::project(modis20, covar_90m_rst)
modis20_90 <- resample(modis_proj, covar_90m_rst, method = "near", filename = "data/modis20_90.tif", overwrite = TRUE)


library(raster)
modis20_90 <- raster("data/modis20_90.tif")
rcl_collapsed <- data.frame(from = c(0:15, 255),
                            to = c(0, 1, 1, 2, 2, 
                                   3, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 
                                   0, NA))
lc_collapsed <- raster::reclassify(modis20_90, rcl_collapsed)
writeRaster(lc_collapsed, "D:/hotspot/landcover.tif", overwrite = TRUE)



#1 = evergreen, 2 = deciduous, 3 = mixed

ev_rcl <- data.frame(from = c(0, 1, 2, 3),
                        to = c(0, 1, 0, 0.5))
ev_collapsed <- raster::reclassify(lc_collapsed, ev_rcl)
writeRaster(ev_collapsed, "D:/hotspot/evergreen_combined.tif", overwrite = TRUE)
de_rcl <- data.frame(from = c(0, 1, 2, 3),
                     to = c(0, 0, 1, 0.5))
de_collapsed <- raster::reclassify(lc_collapsed, de_rcl)
writeRaster(de_collapsed, "D:/hotspot/deciduous_combined.tif", overwrite = TRUE)

lc_collapsed <- rast("D:/hotspot/landcover.tif")
evergreen_combined <- rast("D:/hotspot/evergreen_combined.tif")
deciduous_combined <- rast("D:/hotspot/deciduous_combined.tif")

buffer_diameters <- c(540, 1530, 2520)#540, 
  # Add a for loop over buffer distances
  for (bd in buffer_diameters) {
    
    message("Working on buffer ", bd)
    fcl_rst_fn_ev <- paste0("D:/hotspot/evergreen_combined_", bd, ".tif")
    fcl_rst_fn_de <- paste0("D:/hotspot/deciduous_combined_", bd, ".tif")
    #if (file.exists(fcl_rst_fn)) {
    #message("File already created. Skipping...")
    # next
    #}
    
    fcl_mat_ev <- focalMat(evergreen_combined, bd / 2, type = "circle")
    fcl_rst_ev <- focal(evergreen_combined, fcl_mat_ev, fun = "mean")
    names(fcl_rst_ev) <- paste0("evergreen_combined_", bd)
    writeRaster(fcl_rst_ev, fcl_rst_fn_ev, overwrite = TRUE)
  
    fcl_mat_de <- focalMat(deciduous_combined, bd / 2, type = "circle")
    fcl_rst_de <- focal(deciduous_combined, fcl_mat_de, fun = "mean")
    names(fcl_rst_de) <- paste0("deciduous_combined_", bd)
    writeRaster( fcl_rst_de, fcl_rst_fn_de, overwrite = TRUE)
    
    print(Sys.time())
    
  }

modis20_90m <- resample(modis20, )

rst <- rast("D:/hotspot/deciduous_combined_540.tif") 
plot(rst$deciduous_combined_540)
