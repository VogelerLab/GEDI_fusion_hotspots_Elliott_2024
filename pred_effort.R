library(terra)
library(dplyr)
library(lubridate)
library(stringr)
covar_stack_full <- rast("D:/hotspot/data/covar_stack_full.tif")
setwd("~/GEDI Work/data/hotspot/data")

#create prediction surface template
pred_surface <- rast("fhd_90m.tif") 
ps_proj <- project(pred_surface, covar_stack_full)

#read in MWCF layer
mwcf <- vect("C:/Users/lisah/Documents/ArcGIS/MarineWCForest.shp")
#convert mwcf to projection of covar_stack_full
mwcf_proj <- project(mwcf, covar_stack_full)
#use MWCF layer to crop and mask covar_stack_full
#covar_stack_mwcf <- crop(covar_stack_full, mwcf_proj, mask = TRUE, filename = "D:/hotspot/covar_stack_mwcf_m.tif", overwrite = TRUE)
covar_stack_mwcf <- crop(covar_stack_full, mwcf_proj, mask = FALSE, filename = "D:/hotspot/covar_stack_mwcf.tif", overwrite = TRUE)
covar_stack_mwcf <- mask(covar_stack_mwcf, mwcf_proj, filename = "D:/hotspot/covar_stack_mwcf_m.tif", overwrite = TRUE)
#create the effort prediction surface
ps_proj_mwcf <- crop(ps_proj, mwcf_proj, mask = TRUE)
ps_vect_mwcf <- as.data.frame(ps_proj_mwcf, xy = TRUE, cells = TRUE)
ps_vect_mwcf <- ps_vect_mwcf[,c(1:3)]
pred_surface_eff_mwcf <- ps_vect_mwcf %>% 
  mutate(observation_date = ymd(str_glue("2020-07-01")),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         time_observations_started = 6.0,
         duration_minutes = 60,
         number_observers = 1)
pred_surface_eff_mwcf <- pred_surface_eff_mwcf[,c(1:3, 5:8)]
saveRDS(pred_surface_eff_mwcf, "D:/hotspot/effort_pred_mwcf.rds")
#eff_mwcf_rst <- rast(pred_surface_eff_mwcf, type = "xyz")
eff_mwcf_rst <- rast(ps_proj_mwcf, nlyr=ncol(pred_surface_eff_mwcf))
values(eff_mwcf_rst) <- as.matrix(pred_surface_eff_mwcf)
names(eff_mwcf_rst) <- colnames(pred_surface_eff_mwcf)
eff_mwcf_rst <- writeRaster(eff_mwcf_rst, "D:/hotspot/effort_pred_mwcf.tif", overwrite = TRUE)
mwcf_rst <- rast(c("D:/hotspot/covar_stack_mwcf.tif", "D:/hotspot/effort_pred_mwcf.tif"))
pred_mwcf <- writeRaster(mwcf_rst, "D:/hotspot/pred_mwcf.tif")


nade <- vect("C:/Users/lisah/Documents/ArcGIS/NADessert.shp")
#convert mwcf to projection of covar_stack_full
nade_proj <- project(nade, covar_stack_full)
#use MWCF layer to crop and mask covar_stack_full
covar_stack_nade <- crop(covar_stack_full, nade_proj, mask = FALSE, filename = "D:/hotspot/covar_stack_nade.tif", overwrite = TRUE)
covar_stack_nade <- mask(covar_stack_nade, nade_proj, filename = "D:/hotspot/covar_stack_nade_m.tif", overwrite = TRUE)
#create the effort prediction surface
ps_proj_nade <- crop(ps_proj, nade_proj, mask = TRUE)
ps_vect_nade <- as.data.frame(ps_proj_nade, xy = TRUE, cells = TRUE)
ps_vect_nade <- ps_vect_nade[,c(1:3)]
pred_surface_eff_nade <- ps_vect_nade %>% 
  mutate(observation_date = ymd(str_glue("2020-07-01")),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         time_observations_started = 6.0,
         duration_minutes = 60,
         number_observers = 1)
pred_surface_eff_nade <- pred_surface_eff_nade[,c(1:3, 5:8)]
saveRDS(pred_surface_eff_nade, "D:/hotspot/effort_pred_nade.rds")
#eff_nade_rst <- rast(pred_surface_eff_nade, type = "xyz")
eff_nade_rst <- rast(ps_proj_nade, nlyr=ncol(pred_surface_eff_nade))
values(eff_nade_rst) <- as.matrix(pred_surface_eff_nade)
names(eff_nade_rst) <- colnames(pred_surface_eff_nade)
eff_nade_rst <- writeRaster(eff_nade_rst, "D:/hotspot/effort_pred_nade.tif", overwrite = TRUE)
nade_rst <- rast(c("D:/hotspot/covar_stack_nade.tif", "D:/hotspot/effort_pred_nade.tif"))
pred_nade <- writeRaster(nade_rst, "D:/hotspot/pred_nade.tif")

grpl <- vect("C:/Users/lisah/Documents/ArcGIS/GRPL.shp")
#convert mwcf to projection of covar_stack_full
grpl_proj <- project(grpl, covar_stack_full)
#use MWCF layer to crop and mask covar_stack_full
covar_stack_grpl <- crop(covar_stack_full, grpl_proj, mask = FALSE, filename = "D:/hotspot/covar_stack_grpl.tif", overwrite = TRUE)
covar_stack_grpl <- mask(covar_stack_grpl, grpl_proj, filename = "D:/hotspot/covar_stack_grpl_m.tif", overwrite = TRUE)
#create the effort prediction surface
ps_proj_grpl <- crop(ps_proj, grpl_proj, mask = TRUE)
ps_vect_grpl <- as.data.frame(ps_proj_grpl, xy = TRUE, cells = TRUE)
ps_vect_grpl <- ps_vect_grpl[,c(1:3)]
pred_surface_eff_grpl <- ps_vect_grpl %>% 
  mutate(observation_date = ymd(str_glue("2020-07-01")),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         time_observations_started = 6.0,
         duration_minutes = 60,
         number_observers = 1)
pred_surface_eff_grpl <- pred_surface_eff_grpl[,c(1:3, 5:8)]
saveRDS(pred_surface_eff_grpl, "D:/hotspot/effort_pred_grpl.rds")
#eff_grpl_rst <- rast(pred_surface_eff_grpl, type = "xyz")
eff_grpl_rst <- rast(ps_proj_grpl, nlyr=ncol(pred_surface_eff_grpl))
values(eff_grpl_rst) <- as.matrix(pred_surface_eff_grpl)
names(eff_grpl_rst) <- colnames(pred_surface_eff_grpl)
eff_grpl_rst <- writeRaster(eff_grpl_rst, "D:/hotspot/effort_pred_grpl.tif", overwrite = TRUE)
grpl_rst <- rast(c("D:/hotspot/covar_stack_grpl.tif", "D:/hotspot/effort_pred_grpl.tif"))
pred_grpl <- writeRaster(grpl_rst, "D:/hotspot/pred_grpl.tif")



nwmf <- vect("C:/Users/lisah/Documents/ArcGIS/NWMF.shp")
#convert mwcf to projection of covar_stack_full
nwmf_proj <- project(nwmf, covar_stack_full)
#use MWCF layer to crop and mask covar_stack_full
#covar_stack_nwmf <- crop(covar_stack_full, nwmf_proj, mask = TRUE, filename = "D:/hotspot/covar_stack_nwmf.tif", overwrite = TRUE)
covar_stack_nwmf <- crop(covar_stack_full, nwmf_proj, mask = FALSE, 
                         filename = "D:/hotspot/covar_stack_nwmf.tif", overwrite = TRUE)
#rm(covar_stack_full)
covar_stack_nwmf <- mask(covar_stack_nwmf, nwmf_proj, filename = "D:/hotspot/covar_stack_nwmf_m.tif", overwrite = TRUE)
#create the effort prediction surface
ps_proj_nwmf <- crop(ps_proj, nwmf_proj, mask = TRUE)
writeRaster(ps_proj_nwmf, "D:/hotspot/ps_proj_nwmf.tif")

ps_vect_nwmf <- as.data.frame(ps_proj_nwmf, xy = TRUE, cells = TRUE)
ps_vect_nwmf <- ps_vect_nwmf[,c(1:3)]
saveRDS(ps_vect_nwmf, "D:/hotspot/ps_vect_nwmf.rds")
pred_surface_eff_nwmf <- ps_vect_nwmf %>% 
  mutate(observation_date = ymd(str_glue("2020-07-01")),
         year = year(observation_date),
         day_of_year = yday(observation_date),
         time_observations_started = 6.0,
         duration_minutes = 60,
         number_observers = 1)
pred_surface_eff_nwmf <- pred_surface_eff_nwmf[,c(1:3, 5:8)]
saveRDS(pred_surface_eff_nwmf, "D:/hotspot/effort_pred_nwmf.rds")

library(terra)
pred_surface_eff_nwmf <-readRDS("D:/hotspot/effort_pred_nwmf.rds")
#ps_vect_nwmf <- readRDS("D:/hotspot/ps_vect_nwmf.rds")
ps_proj_nwmf <- rast("D:/hotspot/ps_proj_nwmf.tif")
#eff_nwmf_rst <- rast(pred_surface_eff_nwmf, type = "xyz")
eff_nwmf_rst <- rast(ps_proj_nwmf, nlyr=ncol(pred_surface_eff_nwmf))
#rm(pred_surface, nwmf)#, ps_proj_nwmf)
values(eff_nwmf_rst) <- as.matrix(pred_surface_eff_nwmf)
names(eff_nwmf_rst) <- colnames(pred_surface_eff_nwmf)
eff_nwmf_rst <- writeRaster(eff_nwmf_rst, "D:/hotspot/effort_pred_nwmf.tif", overwrite = TRUE)
rm(eff_nwmf_rst, ps_proj_nwmf, pred_surface_eff_nwmf, ps_vect_nwmf, ps_proj)

nwmf_rst <- rast(c("D:/hotspot/covar_stack_nwmf.tif", "D:/hotspot/effort_pred_nwmf.tif"))
pred_nwmf <- writeRaster(nwmf_rst, "D:/hotspot/pred_nwmf.tif")
