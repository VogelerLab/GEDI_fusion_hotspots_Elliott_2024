setwd("~/GEDI Work/data")
library(sf)
library(rnaturalearth)
library(dplyr)

# file to save spatial data
gpkg_dir <- "data"
if (!dir.exists(gpkg_dir)) {
  dir.create(gpkg_dir)
}
f_ne <- file.path(gpkg_dir, "gis-data2.gpkg")
map_proj <- st_crs("ESRI:102003")

# political boundaries
# land border with lakes removed
#ne_land <- ne_download(scale = 50, category = "cultural",
#                       type = "admin_0_countries_lakes",
#                       returnclass = "sf") %>%
#  filter(CONTINENT == "North America") %>%
#  st_set_precision(1e6) %>%
#  st_union()
# country lines
# downloaded globally then filtered to north america with st_intersect()
#ne_country_lines <- ne_download(scale = 50, category = "cultural",
#                                type = "admin_0_boundary_lines_land",
#                                returnclass = "sf") %>% 
#  st_geometry()
#ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
#  as.logical() %>%
#  {ne_country_lines[.]}
# states, north america
#ne_state_lines <- ne_download(scale = 50, category = "cultural",
#                              type = "admin_1_states_provinces_lines",
#                              returnclass = "sf") %>%
#  filter(ADM0_A3 %in% c("USA", "CAN")) %>%
#  mutate(iso_a2 = recode(ADM0_A3, USA = "US", CAN = "CAN")) %>% 
#  select(country = ADM0_NAME, country_code = iso_a2)

#create ebird layer
ebird <- read_csv("hotspot/ebird/all_ebird.csv")
ebird_sf <- ebird %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj) %>% 
  dplyr::select(dowo)


# output
#unlink(f_ne)
#write_sf(ne_land, f_ne, "ne_land")
#write_sf(ne_country_lines, f_ne, "ne_country_lines")
#write_sf(ne_state_lines, f_ne, "ne_state_lines")
#write_sf(ebird_sf, f_ne, "PIWO_sf")

#write_sf(ebird_sf, "data", "NOFL_sf", driver = "ESRI Shapefile")


