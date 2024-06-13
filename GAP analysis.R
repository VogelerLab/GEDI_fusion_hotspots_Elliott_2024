library(terra)
hotspot_rst <- rast("D:/hotspot/hotspot_mwcf.tif")
padus <- vect("D:PADUS/PADUS_mwcf.shp")
padus_rst <- rasterize(padus, hotspot_rst, field = "own_type", filename = "D:/hotspot/padus_rst.tif", overwrite =  TRUE)
padus_rst <- rast("D:/hotspot/padus_rst.tif")
hotspot_padus <- freq(padus_rst, zones = hotspot_rst)

hotspots_own <- subset(hotspot_padus, hotspot_padus$zone == 2)
notspots <- subset(hotspot_padus, hotspot_padus$zone == 1)
for(r in 1:nrow(hotspots_own)){
  hotspots_own$pct[r] <- (hotspots_own$count[r]/sum(hotspots_own$count))*100
  hotspots_own$hotspots <- "Hotspots"
}
for(r in 1:nrow(notspots)){
  notspots$pct[r] <- (notspots$count[r]/sum(notspots$count))*100
  notspots$hotspots <- "Other Land"
}


pct_own <- rbind(hotspots_own, notspots)

library(ggplot2)
library(viridis)
ggplot(data = pct_own, aes(fill=value, y=pct, x=hotspots)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  ggtitle("Ownership") +
  xlab("") +
  ylab("Percent")

hotspot_agg <- aggregate(hotspot_padus$count, list(hotspot_padus$value), FUN=sum)

for(r in 1:nrow(hotspot_agg)){
  hotspot_agg$pct[r] <- (hotspot_agg$x[r]/sum(hotspot_agg$x))*100
}

total_own <- data.frame(layer = rep(1, times = nrow(hotspot_agg)), 
                        value = hotspot_agg$Group.1, count = hotspot_agg$x,
                        zone = rep(3), pct = hotspot_agg$pct,
                        hotspots = rep("Total"))

all_own <- rbind(pct_own, total_own)
ggplot(data = all_own, aes(fill=value, y=pct, x=hotspots)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  ggtitle("Ownership") +
  xlab("") +
  ylab("Percent")

hotspots_total <- subset(all_own, all_own$hotspots == "Hotspots" | all_own$hotspots == "Total")
ggplot(data = hotspots_total, aes(fill=value, y=pct, x=hotspots)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  ggtitle("Ownership") +
  xlab("") +
  ylab("Percent")

#######
covar_stack_full <- rast("D:/hotspot/data/covar_stack_full.tif")
setwd("~/GEDI Work/data/hotspot/data")
buffer_diameters <- c(90)
gedi_variables <- c("cover", "rh98", "fhd")
bioclim_variables <- paste0("biovars_", 1:19)
comp_variables <- c("evergreen_combined_540", "deciduous_combined_540")
v_names <- c(paste0(rep(gedi_variables, each = length(buffer_diameters)), 
                  "_", rep(buffer_diameters, times = length(gedi_variables)))) #,bioclim_variables, comp_variables
file_list <- c(paste0(getwd(), "/", v_names, "m.tif"), "biovars_90.tif")

covar_rst <- rast(file_list) 
hotspot_proj <- project(hotspot_rst, covar_rst)
#use MWCF layer to crop and mask covar_stack_full
#covar_stack_mwcf <- crop(covar_stack_full, mwcf_proj, mask = TRUE, filename = "D:/hotspot/covar_stack_mwcf_m.tif", overwrite = TRUE)
covar_stack_hotspot <- crop(covar_rst, hotspot_proj, mask = FALSE, filename = "D:/hotspot/covar_stack_hotspot.tif", overwrite = TRUE)
covar_stack_hotspot <- mask(covar_stack_hotspot, hotspot_proj, filename = "D:/hotspot/covar_stack_hotspot_m.tif", overwrite = TRUE)
#create the effort prediction surface

hotspot_covars <- zonal(covar_stack_hotspot, hotspot_proj, fun = "mean", filename = "D:/hotspot/hotspot_covars.tif", overwrite = TRUE)
hotspot_covars_sd <- zonal(covar_stack_hotspot, hotspot_proj, fun = "sd", filename = "D:/hotspot/hotspot_covars.tif", overwrite = TRUE)
hotspot_covars$metric <- c("Mean")
hotspot_covars_sd$metric <- c("SD")
hotspot_summary <- rbind(hotspot_covars, hotspot_covars_sd)
#hotspot_summary <- data.frame(Location = c("Hotspot", "Not Hotspot"), 
#                              cover_90m = c(0.56, 0.51), cover_sd = c(0.08, 0.24),
#                              rh98_90m = c(26.9, 23.7), rh98_sd = c(3.8, 11.3),
#                              fhd_90m = c(2.87, 2.56), fhd_sd = c(0.14, 0.57))
write.csv(hotspot_summary, "D:/hotspot/hotspot_summary.csv")  
write.csv(hotspots_total, "D:/hotspot/hotspot_ownership.csv")


#Forest structure 
rst <- rast("D:/hotspot/covar_stack_mwcf_m.tif")
rst1 <- rst[[c("evergreen_combined_540", "deciduous_combined_540")]]
hotspot_proj <- project(hotspot_rst, rst1)
covar_stack_hotspot <- crop(rst1, hotspot_proj, mask = FALSE, filename = "D:/hotspot/covar_stack_hotspot.tif", overwrite = TRUE)
covar_stack_hotspot <- mask(covar_stack_hotspot, hotspot_proj, filename = "D:/hotspot/covar_stack_hotspot_m.tif", overwrite = TRUE)
#create the effort prediction surface
hotspot_covars <- zonal(covar_stack_hotspot, hotspot_proj, fun = "mean", filename = "D:/hotspot/hotspot_covars.tif", overwrite = TRUE)
hotspot_covars_sd <- zonal(covar_stack_hotspot, hotspot_proj, fun = "sd", filename = "D:/hotspot/hotspot_covars.tif", overwrite = TRUE)
hotspot_covars$metric <- c("Mean")
hotspot_covars_sd$metric <- c("SD")
hotspot_summary_landcover <- rbind(hotspot_covars, hotspot_covars_sd)
write.csv(hotspot_summary_landcover, "D:/hotspot/hotspot_summary_landcover.csv")  


hotspots_total$value <- rep(c("Federal Land", "Local Land", "Tribal Land", "Private Conservation Land", "Private Land", 
                          "State Land"), times = 2)
ggplot(data = hotspots_total, aes(fill=value, y=pct, x=hotspots)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  ggtitle("Ownership") +
  xlab("") +
  ylab("Percent")


#combine hotspot summary tables into a single table
t1 <- read.csv("D:/hotspot/hotspot_summary_landcover.csv")
t2 <- read.csv("D:/hotspot/hotspot_summary.csv")
t3 <- merge(t1, t2, by = c("sum", "metric"))
library(dplyr)
t4 <- t3 %>%
  select(!c(X.x, X.y, focal_sum)) %>%
  pivot_longer(cols = !c(metric, sum), names_to = "Variable", values_to = "Mean") %>%
  pivot_wider(names_from = metric, values_from = Mean) %>%
  pivot_wider(names_from = sum, values_from = c(Mean, SD)) %>%
  rename(Hotspot_mean = Mean_1, Hotspot_SD = SD_1, "Outside_mean" = Mean_0, "Outside_SD" = SD_0) %>%
  mutate(across(where(is.numeric), round, 2))
write.csv(t4, "D:/hotspot/hotspot_summary_all.csv")
