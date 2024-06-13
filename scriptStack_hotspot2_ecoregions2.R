

script_dir <- "C:/Users/lisah/Documents/GEDI Work/hotspot analysis"
source(file.path(script_dir, "all_bird_data.R")) #make the individual sp_data.csv files, plus maps of ebird observations
source(file.path(script_dir, "birdData.R")) #make the file with all species combined into a single file

source(file.path(script_dir, "ecoregion_datasets.R")) #combine ebird data points with ecoregion

#source(file.path(script_dir, "MODIS_obtain.R")) #to download the modis data and save raster to computer...only should need to do this once

source(file.path(script_dir, "MODIS_combine.R"))#to combine modis data from saved raster
source(file.path(script_dir, "vars_combine2.R"))#to create gedi and bioclim csv data from saved raster
source(file.path(script_dir, "dataset_for_RF2_ecoregion.R")) #to combine all data and set up test/training data for all 



source(file.path(script_dir, "rf_prep_ecoregion.R")) #select vars for rf...this will need to be hand selected for each species
source(file.path(script_dir, "rf_analysis_ecoregion.R")) #complete rf for the final analysis

source(file.path(script_dir, "prediction_dataset.R")) #creates covar_stack_full.tif
source(file.path(script_dir, "vars_combine_pred.R"))#to create the prediction surface files including deciduous_combined_540.tif
source(file.path(script_dir, "pred_effort.R")) #to create the "D:/hotspot/covar_stack_mwcf_m.tif"
source(file.path(script_dir, "prediction_surface.R")) #makes "D:/hotspot/", ecoregions[er], "_pred_surface.tif"

source(file.path(script_dir, "rf_prediction_ecoregion.R")) #complete spatial predictions
source(file.path(script_dir, "ecoregion_quantiles.R")) #to find thresholds to get binary grids
source(file.path(script_dir, "hotspots_ecoregion.R")) #combine individual species maps to get species richness and hotspots
source(file.path(script_dir, "GAP Analysis.R"))
