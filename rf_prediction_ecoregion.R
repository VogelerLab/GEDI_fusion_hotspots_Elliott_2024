library(terra)
library(scam)
library(ranger)
spp <- c("nofl", "dowo",  "hawo", "rbsa", "acwo", "attw") #"piwo")#,    
#spp <- c("nofl", "whwo", "lewo", "wisa", "rnsa", "bbwo")
ecoregion <- c("mwcf")

for(sp in spp){
  print(sp)
  rf_name <- paste0("D:/hotspot/output/rf_", sp, "_", ecoregion, ".rds")
  calib_name <- paste0("D:/hotspot/output/calibration_model_", sp, "_", ecoregion, ".rds") 
  rf <- readRDS(rf_name)
  calibration_model <- readRDS(calib_name)
  ecoregions <- c("mwcf")#,"mwcf", "nade", "grpl", "nwmf")
  for(er in 1:length(ecoregions)){
    print(ecoregions[er])
    rast_nm <- paste0("D:/hotspot/", ecoregions[er], "_pred_surface.tif")
    pred_surface_eff <- rast(rast_nm)
    names(rf$variable.importance) %in% names(pred_surface_eff)
    # predict
    #pred_rf <- predict(rf, data = pred_surface_eff, type = "response")
    #pred_rf <- pred_rf$predictions[, 2]
    #pred_rf <- terra::predict(pred_surface_eff, rf, cores = 2, cpkgs = "ranger", progress = "text")
    pred_rf <- terra::predict(pred_surface_eff, rf, fun = function(model, ...) predict(model, num.threads = 1, ...)$predictions,
                              type = "response", cpkgs = "ranger", na.rm = TRUE)
    # apply calibration models
    if(names(pred_rf)[2] == "TRUE."){
      pred_rf_true <-  pred_rf$TRUE.
    } else{
      if(names(pred_rf)[1] == "lyr1"){
        pred_rf_true <- pred_rf$lyr1 
      } else {
        print("pred_rf col name?")
      }
    }
    
    
    names(pred_rf_true) <- "pred"
    writeRaster(pred_rf_true, paste0("D:/hotspot/output/ecoregion/pred_", sp, "_", ecoregions[er], ".tif"), overwrite = TRUE)
    pred_rf_cal <- terra::predict(pred_rf_true, calibration_model, 
                                  type = "response", cpkgs = "scam", na.rm = TRUE)
    pred_name <- paste0("D:/hotspot/output/ecoregion/rf-model_encounter-rate_", sp, "_", ecoregions[er], ".tif")
    writeRaster(pred_rf_cal, pred_name, overwrite = TRUE)
    r_max <- ceiling(10 * minmax(pred_rf_cal)[2,]) / 10
    print(r_max)
    plot(pred_rf_cal)
  }
} 
