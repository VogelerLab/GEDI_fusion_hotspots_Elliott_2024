library(terra)
setwd("~/GEDI Work/data/hotspot/data")
buffer_diameters <- c(90, 180, 360, 540, 1530, 2520)

variables <- c("cover", "rh98", "fhd")

v_names <- paste0(rep(variables, each = length(buffer_diameters)), 
                  "_", rep(buffer_diameters, times = length(variables)))
#file_list <- paste0(getwd(), "/", v_names, "m.tif")
    
#covar_rst <- rast(file_list) 
#names(covar_rst) <- v_names

#covar_rst <- writeRaster(covar_rst, "covar_stack.tif")

#covar_df <- terra::as.data.frame(covar_rst)
 
for(f in v_names){
  print(f)
  f_name <- paste0(getwd(), "/", f, "m.tif")
  covar_rst <- rast(f_name)
  covar_df <- terra::as.data.frame(covar_rst, xy = TRUE, cells = TRUE)
  names(covar_df)[4] = f
  c_name <- paste0(getwd(), "/", f, ".rds")
  saveRDS(covar_df, c_name, compress = TRUE)
}

#combine the list of variables
file_list <- paste0(getwd(), "/", v_names, ".rds")
library(dplyr)
var_df <- readRDS(file_list[1])
for(i in 2:length(file_list)){
  print(file_list[i])
  covar_dat <- readRDS(file_list[i])
  var_df <- full_join(var_df, covar_dat, by = "cell")
}

#read in biovars
biovars_df <- 

#merge biovars with var_df
var_df <- full_join(var_df, biovars_df, by = "cell")



###############
forest_list <- c("evergreen", "deciduous")
buffs <- c(540, 1530, 2520)
modis_list <- paste0("D:/hotspot/", rep(forest_list, each = length(buffs)), 
                     "_combined_", rep(buffs, times = length(forest_list)), ".tif")
f_list <- c("covar_stack.tif", "biovars_90.tif", modis_list)
var_stack <- rast(f_list)
covar_stack_full <- writeRaster(var_stack, "D:/hotspot/data/covar_stack_full.tif", overwrite = TRUE)
##################################
#setwd("C:/Users/lisah/Documents/GEDI Work/data/hotspot/data")
library(dplyr)
library(lubridate)
library(stringr)
t_peak = 6.0

#add prediction surface filler
pred_surface <- readRDS("fhd_90.rds") 
pred_surface <- pred_surface[,c(1:3)]
# add effort covariates to prediction 
pred_surface_eff <- pred_surface %>% 
  mutate(observation_date = ymd(str_glue("2020-07-01")),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         time_observations_started = t_peak,
         duration_minutes = 60,
         number_observers = 1)
saveRDS(pred_surface_eff, "effort_pred")
rm(pred_surface)
effort_rst <- rast(pred_surface_eff)
writeRaster(pred_surface_eff, "effort.tif")


pred_stack <- c(pred_surface_eff, )
pred_surface_full <-  ("D:/hotspot/data/prediction_surface_full.tif") 
