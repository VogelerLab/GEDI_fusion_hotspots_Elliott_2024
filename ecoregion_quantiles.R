library(terra)
library(dplyr)
library(sf)
#attw
spp_list <- c("attw", "acwo", "dowo", "hawo", "nofl", "piwo", "rbsa")
ecoregion <- "mwcf"
for(sp in spp_list){
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
  plot(quants_rst, bty="n", axes = FALSE, box=FALSE)
  plot_nm <- paste0("D:/hotspot/output/ecoregion/quantile_", sp, "_", ecoregion, ".png")
  png(plot_nm, height=nrow(quants_rst), width=ncol(quants_rst)) 
  plot(quants_rst, bty="n", axes = FALSE, box=FALSE, legend = FALSE)
  dev.off()
  test_name <- paste0("D:/hotspot/", sp, "_", ecoregion, "_test.csv")
  bird_pts <- read.csv(test_name)
  ebird_proj <- bird_pts %>% 
    distinct(locality_id, latitude, longitude) %>%
    # convert to spatial features
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    # transform to projection of raster
    st_transform(crs = crs(quants_rst))
  
  extr <- extract(quants_rst, vect(ebird_proj), progress = TRUE)
  quants_df <- cbind(ebird_proj, extr)
  
  samps <- subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7| quants_df$pred == 6 | quants_df$pred == 5 | quants_df$pred == 4| quants_df$pred == 3 | quants_df$pred == 2 | quants_df$pred == 1) 
  pct95 <- nrow(samps)*0.75
  
  bin10 <- nrow(subset(quants_df, quants_df$pred == 10))
  bin9 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9))
  bin8 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8))
  bin7 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7))
  bin6 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7| quants_df$pred == 6))
  bin5 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7| quants_df$pred == 6 | quants_df$pred == 5))
  bin4 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7| quants_df$pred == 6 | quants_df$pred == 5| quants_df$pred == 4))
  bin3 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7| quants_df$pred == 6 | quants_df$pred == 5| quants_df$pred == 4 | quants_df$pred == 3))
  bin2 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7| quants_df$pred == 6 | quants_df$pred == 5| quants_df$pred == 4 | quants_df$pred == 3 | quants_df$pred == 2))
  bin1 <- nrow(subset(quants_df, quants_df$pred == 10 | quants_df$pred == 9 | quants_df$pred == 8 | quants_df$pred == 7| quants_df$pred == 6 | quants_df$pred == 5 | quants_df$pred == 4| quants_df$pred == 3 | quants_df$pred == 2 | quants_df$pred == 1))
  
  if(bin10 >= pct95){
    bestbin <- 10
  } else{
    if(bin9 >= pct95){
      bestbin <- 9
    } else{
      if(bin8 >= pct95){
        bestbin <- 8
      } else{
        if(bin7 >= pct95){
          bestbin <- 7
        } else{
          if(bin6 >= pct95){
            bestbin <- 6
          } else{
            if(bin5 >= pct95){
              bestbin <- 5
            } else{
              if(bin4 >= pct95){
                bestbin <- 4
              } else{
                if(bin3 >= pct95){
                  bestbin <- 3
                } else{
                  if(bin2 >= pct95){
                    bestbin <- 2
                  } else{
                    if(bin1 >= pct95){
                      bestbin <- 1
                    } else{
                      bestbin <- 0
                    }
  }}}}}}}}}
  bestbin
  reclass_m
  reclass_df2 <- c(1:10, rep(0, times = bestbin - 1), rep(1, times = 10 - bestbin + 1)) 
  reclass_m2 <- matrix(reclass_df2,
                      ncol = 2,
                      byrow = FALSE)
  binary_nm <- paste0("D:/hotspot/output/ecoregion/binary_", sp, "_", ecoregion, ".tif")
  binary_rst <- classify(quants_rst,
                         rcl = reclass_m2, 
                         filename = binary_nm,
                         overwrite = TRUE)
  plot(binary_rst)
}


