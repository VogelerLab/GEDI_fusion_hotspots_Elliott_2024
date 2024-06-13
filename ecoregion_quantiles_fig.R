library(terra)
library(dplyr)
library(sf)
library(prettymapr)
#attw
spp_list <- c("attw", "acwo", "dowo", "hawo", "nofl", "piwo", "rbsa")
ecoregion <- "mwcf"

# create a new plotting window and set the plotting area into a 1*2 array
par(mfrow = c(3, 3))

#for(sp in spp_list){
  sp <- c("attw")
  pred_rst <- rast(paste0("D:/hotspot/output/ecoregion/rf-model_encounter-rate_", sp, "_", ecoregion, ".tif"))
  quants <- global(pred_rst, quantile, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  quants[[1]]
  reclass_df <- c(quants[[1]], quants[[2]], quants[[3]], quants[[4]], quants[[5]], 
                  quants[[6]], quants[[7]], quants[[8]], quants[[9]], quants[[10]],
                  quants[[2]], quants[[3]], quants[[4]], quants[[5]], quants[[6]],
                  quants[[7]], quants[[8]], quants[[9]], quants[[10]], quants[[11]],
                  1:10)
  reclass_m <- matrix(reclass_df,
                      ncol = 3,
                      byrow = FALSE)
  quant_nm <- paste0("D:/hotspot/output/ecoregion/quants_", sp, "_", ecoregion, ".tif")
  quants_rst <- classify(pred_rst,
                         rcl = reclass_m, 
                         filename = quant_nm,
                         overwrite = TRUE)
  plot(quants_rst, bty="n", axes = FALSE, box=FALSE, legend = TRUE, main="G) ATTW")
  north(xy = "right", type =2, angle = 20)
  sbar(200000, xy = "bottomright", type = "bar", divs = 2, scaleby = 1000, below = "km")#, cex = 0.8, type = "bar", divs = 2) 
  
  
  
  
  
  
  
  
  plot_nm <- paste0("D:/hotspot/output/ecoregion/quantile_", ecoregion, "_allspp.png")
  png(plot_nm, height=nrow(quants_rst), width=ncol(quants_rst)) 
  
  dev.off()
}

par(mfrow = c(3, 3))
