#setwd("~/GEDI Work/ebird practice/data/ebird")
setwd("~/GEDI Work/data")
library(sf)
library(dplyr)
# load and project gis data
map_proj <- st_crs("ESRI:102003")
#ne_land <- read_sf("data/gis-data2.gpkg", "ne_land") %>% 
 # st_transform(crs = map_proj) %>% 
#  st_geometry()
#ne_country_lines <- read_sf("data/gis-data2.gpkg", "ne_country_lines") %>% 
 #st_transform(crs = map_proj) %>% 
  #st_geometry()
#ne_state_lines <- read_sf("data/gis-data2.gpkg", "ne_state_lines") %>% 
  #st_transform(crs = map_proj) %>% 
  #st_geometry()
#NOFL <- read_sf("data/gis-data2.gpkg", "NOFL_sf") %>% 
 # st_transform(crs = map_proj) %>% 
#  st_geometry()
ebird <- read.csv("hotspot/ebird/all_ebird_ecoregion.csv")
ebird_sf <- ebird %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj) %>% 
  dplyr::filter(ecoregion == "mwcf") %>%
  dplyr::select(dowo) 
  

SPsf <- ebird_sf

# prepare ebird data for mapping
#ebird <- read.csv("ebd_woothr_june_bcr27_zf.csv")


# map
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(st_geometry(SPsf), col = NA)
# contextual gis data
#plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
#plot(bcr, col = "#cccccc", border = NA, add = TRUE)
#plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
#plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)


ecoregion_file <- paste0("C:/Users/lisah/Documents/ArcGIS/MarineWCForest.shp")
#read in MWCF layer
ecor_vect <- vect(ecoregion_file)
#convert mwcf to projection of covar_stack_full
ecor_proj <- project(ecor_vect, SPsf)
# ebird observations
# not observed
library(ggplot2)
plot(ecor_proj, axes=FALSE, box=FALSE)
plot(st_geometry(SPsf),
     pch = 19, cex = 0.1, col = alpha("#4daf4a", 1),
     add = TRUE)
plot(foo)

# legend
legend("bottomleft", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklists", paste0(cn,"\nsightings")),
       pch = 19)
box()
par(new = TRUE, mar = c(0, 0, 3, 0))
#title("Pileated Woodpecker Stationary eBird Observations\nBreeding Season 2016-2020, eBird-Range Study Area")
mname <- paste0("maps/", sp, "_map.png")
ggsave(mname)

#####################
library(ggplot2)

plot(st_geometry(SPsf),
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)



library(raster)
library(ggplot2)
library(ggthemes)
library(ggsn)
library(terra)

mapdata <- getData("GADM", country = "usa", level = 1)
mymap <- fortify(mapdata)

mypoints <- data.frame(long = -121.6945, lat = 39.36708)

g1 <- ggplot() +
  geom_blank(data = mymap, aes(x=long, y=lat)) +
  geom_map(data = mymap, map = mymap, 
           aes(group = group, map_id = id),
           fill = "#b2b2b2", color = "black", size = 0.3) +
  geom_point(data = mypoint, aes(x = long, y = lat),
             color = "black", size = 2) +
  scale_x_continuous(limits = c(-125, -114), expand = c(0, 0)) +
  scale_y_continuous(limits = c(32.2, 42.5), expand = c(0, 0)) +
  theme_map() +
  scalebar(location = "bottomleft", dist = 200,
           dd2km = TRUE, model = 'WGS84',           
           x.min = -124.5, x.max = -114,
           y.min = 33.2, y.max = 42.5) +
  north(x.min = -115.5, x.max = -114,
        y.min = 40.5, y.max = 41.5,
        location = "toprgiht", scale = 0.1)


foo <- map_data("state") 


foo2 <- foo %>%
  filter(region == "oregon" | region == "washington") 

foo2_sf <- st_as_sf(foo2, coords = c("long", "lat"))

f2 <- vect(foo2_sf)

f2_proj <- st_transform(foo2_sf, crs = st_crs(ecor_proj))
ecor_sf <- st_as_sf(ecor_proj)
SPsf

#PLOT OF SPATIAL EXTENT OF MWCF
ggplot(ecor_sf) +
  geom_sf(fill = NA, colour = "black") +
  geom_sf(data = SPsf, colour = "forestgreen", size = 0.5) +
  geom_polygon(data = states, aes(x = x, y = y, group = group), fill = "lightgrey", color = "black") +
  coord_sf()






library(usmap)

states_of_interest <- c("OR", "WA")

#counties <- us_map(regions = "counties", include = states_of_interest)
states <- us_map(include = states_of_interest)
states_proj <- st_transform(states, crs = st_crs(SPsf))

SPsf usmap_transform(SPsf)

ggplot() +
  geom_polygon(data = states, aes(x = x, y = y, group = group), fill = "lightgrey", color = "black") +
  coord_equal()

, aes
             pch = 19, cex = 0.1, col = alpha("#4daf4a", 1),
             add = TRUE)
  ggplot() +
  geom_polygon(data = foo2,
               aes(x = long, y = lat, group = group),
               fill = "#b2b2b2", color = "black", size = 0.3) +
  st_transform(, crs = crs(ecor_proj))


g2 <- ggplotGrob(
  ggplot() +
    geom_polygon(data = foo2,
                 aes(x = long, y = lat, group = group),
                 fill = "#b2b2b2", color = "black", size = 0.3) +
    coord_map("albers",lat0=39, lat1=45) +
    
    geom_polygon(data = as.data.frame(ecor_proj), aes(x = long, y = lat),
               color = "black", size = 2) +
    coord_map("polyconic") +
    theme_map() +
    theme(panel.background = element_rect(fill = NULL))
)     

g3 <- g1 +
  annotation_custom(grob = g2, xmin = -119, xmax = -114,
                    ymin = 31.5, ymax = 36)

