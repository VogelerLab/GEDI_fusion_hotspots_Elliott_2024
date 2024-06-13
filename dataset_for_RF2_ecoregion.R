setwd("C:/Users/lisah/Documents/GEDI Work/data/hotspot")
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


#create ebird layer
ebird <- read.csv("D:/ebird/ebd_piwo_spr20162020_fullStudyArea_zf.csv")
ebird_sf <- ebird %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj)



library(sf)
library(raster)
library(dggridR)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(ebirdst)
library(fields)
library(gridExtra)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# set random number seed to insure fully repeatable results
set.seed(1)

# setup output directory for saved results
if (!dir.exists("output")) {
  dir.create("output")
}

# ebird data
script_dir <- ("C:/Users/lisah/Documents/GEDI Work/data/hotspot/")
data_name_ebird <- paste0(script_dir, "ebird/all_ebird_ecoregion.csv")
ebird <- read_csv(data_name_ebird)

# modis covariates
modis_dat <- read_csv("data/modis_pland_location-year_hotspot4.csv") %>% 
  mutate(year = as.integer(year))

#gedi covariates
gedi_dat <- read_csv("data/gedi_location-year_hotspot.csv")
newcols <- c("year", "locality_id", paste0("cover_", unique(gedi_dat$buffer_distance)),
             paste0("fhd_", unique(gedi_dat$buffer_distance)),
             paste0("rh98_", unique(gedi_dat$buffer_distance)))
gedi_new <- data.frame()

# tranform to wide format, filling in implicit missing values with 0s%>% 
wide_gedi <- gedi_dat %>% 
  group_by(locality_id, year) %>%
  pivot_wider(names_from = buffer_distance, 
              values_from = c(cover, fhd, rh98))

wide_modis <- modis_dat %>% 
  select(!(evergreen:mixed)) %>%
  group_by(locality_id, year) %>%
  pivot_wider(names_from = buffer_distance, 
              values_from = c(evergreen_combined, deciduous_combined))


#bioclim covariates
bioclim_dat <- read_csv("data/bioclim_data_800m.csv")


# combine ebird and habitat data
ebird_g <- inner_join(ebird, wide_gedi, by = c("locality_id", "year"))
ebird_gm <- inner_join(ebird_g, wide_modis, by = c("locality_id", "year"))
ebird_gmb <- inner_join(ebird_gm, bioclim_dat, by = c("locality_id"))



# prediction surface
pred_surface <- read_csv("data/pland-elev_prediction-surface.csv")
# latest year of landcover data
max_lc_year <- 2020
r <- raster("data/fhd_90m.tif")


#ebird_sf
gpkg_dir <- paste0("~/GEDI Work/data")
test_dir <- paste0("~/GEDI work/ebird practice/data/ebird")
library(sf)
library(rnaturalearth)
library(dplyr)
f_ne <- file.path(test_dir, "gis-data2.gpkg")
map_proj <- st_crs("ESRI:102003")
ebird_sf <- ebird_gmb %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj) %>% 
  dplyr::select(dowo)

# load gis data for making maps
map_proj <- st_crs("ESRI:102003")
setwd("D:/hotspot")
ne_land <- read_sf("data/gis-data2.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
studyArea <- ebird_sf %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data2.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- ne_download(scale = 50, category = "cultural",
                              type = "admin_1_states_provinces_lines",
                              returnclass = "sf") %>%
  filter(ADM0_A3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(ADM0_A3, USA = "US", CAN = "CAN")) %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

##Spatiotemporal subsampling:


ebird_gmb <- ebird_gmb[names(ebird_gmb)[names(ebird_gmb) != "species_observed"]]

#For each species:
# sample one checklist per grid cell per week
# sample detection/non-detection independently 
sp_list <- c("piwo", "nofl", "dowo", "whwo", "hawo", "lewo", "wisa", 
        "rbsa", "attw", "rnsa", "bbwo", "acwo")



long_bird <- ebird_gmb %>% 
  pivot_longer(cols = all_of(sp_list),
              names_to = "species",
              values_to = "species_observed")


ecoregion_list <- c("grpl", "mwcf", "nade", "nwmf")
ecoregion_list <- c("mwcf")

sp_stats <- data.frame(Species = rep(sp_list, each = length(ecoregion_list)), Ecoregion = rep(ecoregion_list, times = length(sp_list)))
for(sp in 1:length(sp_list)){
  print(sp_list[sp])
  for(er in 1:length(ecoregion_list)){
    print(ecoregion_list[er])
    obs <- long_bird %>%
      filter(species == sp_list[sp]) %>%
      filter(ecoregion == ecoregion_list[er])
    stat_row <- which(sp_stats$Species == sp_list[sp] & sp_stats$Ecoregion == ecoregion_list[er])
    sp_stats[stat_row, "n_obs"] <- nrow(obs)
    sp_stats[stat_row, "sum_obs"] <- sum(obs$species_observed) 
    if(sum(obs$species_observed) >= 30){
      #set week of year
      for(s in 1:nrow(obs)){
        obs$week[s] <- round(obs$day_of_year[s]/7, digits = 0)
      }
      # generate hexagonal grid with ~ 5 km betweeen cells
      dggs <- dgconstruct(spacing = 5)
      # get hexagonal cell id and week number for each checklist
      checklist_cell <- obs %>% 
        mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
               year, week)
      
      ebird_ss <- checklist_cell %>% 
        group_by(sp_list[sp], year, week, cell) %>% 
        sample_n(size = 1) %>% 
        ungroup()
      
      # original data
      nrow(obs)
      #67550
      #> [1] 121017
      b_samp <- count(obs, species_observed) %>% 
        mutate(percent = n / sum(n))
      sp_stats[stat_row, "b_samp"] <- b_samp$percent[2]
      #piwo      n percent
      #<lgl> <int>   <dbl>
      #  1 FALSE 65252  0.966 
      #2 TRUE   2298  0.0340
      #2 TRUE               2565  0.0212
      
      # after sampling
      nrow(ebird_ss)
      #> [1]  43920
      a_samp <- count(ebird_ss, species_observed) %>% 
        mutate(percent = n / sum(n))
      sp_stats[stat_row, "a_samp"] <- a_samp$percent[2]
      #> # A tibble: 2 x 3
      #>   piwo      n percent
      #<lgl> <int>   <dbl>
      #  1 FALSE 22240  0.929 
      #2 TRUE   1712  0.0715
      #2 TRUE              1928  0.0439
      
      
      # convert checklists to spatial features
      all_pts <- obs %>%  
        st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
        st_transform(crs = map_proj) %>% 
        select(species_observed)
      ss_pts <- ebird_ss %>%  
        st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
        st_transform(crs = map_proj) %>% 
        select(species_observed)
      both_pts <- list(before_ss = all_pts, after_ss = ss_pts)
      
      # map
      p <- par(mfrow = c(2, 1))
      for (i in seq_along(both_pts)) {
        par(mar = c(0.25, 0.25, 0.25, 0.25))
        # set up plot area
        plot(st_geometry(both_pts[[i]]), col = NA)
        # contextual gis data
        #plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
        plot(studyArea, col = "#cccccc", border = NA, add = TRUE)
        #plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
        #plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
        # ebird observations
        # not observed
        plot(st_geometry(both_pts[[i]]),
             pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
             add = TRUE)
        # observed
        plot(filter(both_pts[[i]], species_observed) %>% st_geometry(),
             pch = 19, cex = 0.3, col = alpha("#4daf4a", 0.5),
             add = TRUE)
        # legend
        legend("bottomleft", bty = "n",
               col = c("#555555", "#4daf4a"),
               legend = c("Non-detection", "Detection"),
               pch = 19)
        box()
        par(new = TRUE, mar = c(0, 0, 3, 0))
        if (names(both_pts)[i] == "before_ss") {
          title("A) Before subsampling", adj = 0)
        } else {
          title("B) After subsampling", adj = 0)
        }
      }
      par(p)
      
      
      ebird_ss$observed <- as.integer(as.logical(ebird_ss$species_observed))
      #Create training and validation datasets
      ebird_split <- ebird_ss %>% 
        # select only the columns to be used in the model
        select(state_code, locality_id, latitude, longitude,
               observation_date, species, observed, species_observed,
               year, day_of_year,
               time_observations_started, duration_minutes,
               number_observers,
               starts_with("evergreen_"),
               starts_with("deciduous_"),
               starts_with("biovars_"),
               starts_with("cover"),
               starts_with("fhd"),
               starts_with("rh98"), 
               ecoregion) %>% 
        drop_na()
      # split 80/20
      set.seed(1995)
      ebird_split <- ebird_split %>% 
        split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
      sp_stats[stat_row, "test_n"] <- map_int(ebird_split, nrow)[[1]]
      sp_stats[stat_row, "train_n"] <- map_int(ebird_split, nrow)[[2]]
      #test train 
      #3023 12112
      
      
      sp_stats[stat_row, "detection_freq"] <- mean(ebird_split$train$observed)
      # ranger requires a factor response to do classification
      ebird_split$train$observed <- factor(ebird_split$train$species_observed)
      ebird_split$test$observed <- factor(ebird_split$test$species_observed)
      
      
      #####
      #THIS IS WHERE YOU WILL SAVE THE DATASET
      #########
      setwd("D:/hotspot")
      train_name <- paste0(sp_list[sp], "_", ecoregion_list[er], "_train.csv")
      test_name <- paste0(sp_list[sp], "_", ecoregion_list[er], "_test.csv")
      write.csv(ebird_split$train, train_name, row.names = FALSE)
      write.csv(ebird_split$test, test_name, row.names = FALSE)
    } else{
        print("Too few obs, skipping")
      }  
  }  
}


write.csv(sp_stats, "sp_stats_ecoregion.csv", row.names = FALSE)
