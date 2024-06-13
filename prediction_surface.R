library(terra)

covar_stack_full <- rast("D:/hotspot/data/covar_stack_full.tif")
setwd("~/GEDI Work/data/hotspot/data")

#create prediction surface template
pred_surface <- rast("fhd_90m.tif") 
ps_proj <- project(pred_surface, covar_stack_full)


#read raster layer
ecoregions <- c("grpl", "mwcf", "nade", "nwmf")
er_name <- c("GRPL", "MarineWCForest", "NADessert", "NWMF")
vars <- c("year", "day_of_year", "time_observations_started",
          "duration_minutes", "number_observers")
vals <- c(2020, 182, 6, 60, 1)
for(er in 1:length(ecoregions)){
  print(ecoregions[er])
  ecoregion_file <- paste0("C:/Users/lisah/Documents/ArcGIS/", er_name[er] , ".shp")
  #read in MWCF layer
  ecor_vect <- vect(ecoregion_file)
  #convert mwcf to projection of covar_stack_full
  ecor_proj <- project(ecor_vect, covar_stack_full)
  ps_proj_ecor <- crop(ps_proj, ecor_proj, mask = TRUE)
  
  
  #nwmf <- rast("D:/hotspot/ps_proj_nwmf.tif")
  for(v in 1:length(vars)){
    print(vars[v])
    v_name <- paste0("D:/hotspot/", ecoregions[er], "_", vars[v], ".tif") 
    reclass <- matrix(c(0, 10, vals[v]), ncol = 3)
    var_rast <- classify(ps_proj_ecor, reclass, include.lowest = TRUE, filename = v_name, overwrite = TRUE)
    #repeat for other variables and other regions
  }
file_nm <- paste0("D:/hotspot/", rep(ecoregions[er], times = length(vars)), "_", vars, ".tif")
covar_nm <- paste0("D:/hotspot/covar_stack_", ecoregions[er], "_m.tif")
file_vect <- c(covar_nm, file_nm)
rast_nm <- paste0("D:/hotspot/", ecoregions[er], "_pred_surface.tif")
full_rst_stack <- rast(file_vect)
names(full_rst_stack)[44:48] <- vars

writeRaster(full_rst_stack, rast_nm, overwrite = TRUE)

}

rst <- rast("D:/hotspot/mwcf_pred_surface.tif")
rst
