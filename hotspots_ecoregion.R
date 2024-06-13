library(terra)
#spp <- "piwo"
spp <- c("piwo", "nofl", "dowo", "hawo",  "rbsa", "acwo", "attw")
        # "whwo", "rnsa", "bbwo", "lewo", "wisa",)

ecoregions <- c("mwcf")

#Turn the binomial rasters into a raster stack
f_names <- paste0("D:/hotspot/output/ecoregion/binary_", rep(spp, each = length(ecoregions)), 
                  "_", rep(ecoregions, times = length(spp)), ".tif")
binomial_stack <- rast(f_names)
names(binomial_stack) <- spp

#extract sum to get species richness
spp_richness <- sum(binomial_stack)
plot(spp_richness, bty="n", axes = FALSE, box=FALSE)
north(xy = "right", type =2, angle = 20)
sbar(200000, xy = "bottomright", type = "bar", divs = 2, scaleby = 1000, below = "km")#, cex = 0.8, type = "bar", divs = 2) 
writeRaster(spp_richness, "D:/hotspot/output/ecoregion/spp_richness.tif", overwrite = TRUE)
plot_nm <- paste0("D:/hotspot/output/ecoregion/spp_richness_", "mwcf", ".png")
png(plot_nm, height=nrow(spp_richness), width=ncol(spp_richness)) 
plot(spp_richness, bty="n", axes = FALSE, box=FALSE, legend = TRUE)
dev.off()


spp_richness <- rast("D:/hotspot/output/ecoregion/spp_richness.tif")
#convert species richness into a raster of hotspots
reclass_df <- c(0:6, 1:7,
                0, 0, 0, 0, 0, 1, 1)
reclass_m <- matrix(reclass_df,
                    ncol = 3,
                    byrow = FALSE)
hotspot_rst <- classify(spp_richness,
                            rcl = reclass_m, 
                            filename = "D:/hotspot/output/ecoregion/hotspot_mwcf.tif",
                            overwrite = TRUE)

#rcl_collapsed <- data.frame(from = c(0:6),
#                            to = c(0, 0, 0, 0, 0, 1, 1))
#library(raster)
#hotspot_rst <- reclassify(spp_richness, rcl_collapsed)

plot(hotspot_rst)
plot_nm2 <- paste0("D:/hotspot/output/ecoregion/hotspots_", "mwcf", ".png")
png(plot_nm2, height=nrow(hotspot_rst), width=ncol(hotspot_rst)) 
plot(hotspot_rst, bty="n", axes = FALSE, box=FALSE, legend = FALSE)
north(xy = "right", type =2, angle = 20)
sbar(200000, xy = "bottomright", type = "bar", divs = 2, scaleby = 1000, below = "km")
legend(x = -1990000, y = "2700000", legend = c("Hotspots"), fill = "forestgreen", cex = 0.8)
dev.off()

fcl_mat <- focalMat(binomial_stack, 45, type = "circle")
fcl_rst <- focal(binomial_stack, fcl_mat, fun = "sum")
names(fcl_rst) <- paste0("evergreen_combined_", bd)
writeRaster(fcl_rst_ev, fcl_rst_fn_ev, overwrite = TRUE)


