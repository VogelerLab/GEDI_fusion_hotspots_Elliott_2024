#This script combines all of the ebird data into a single file

script_dir <- "C:/Users/lisah/Documents/GEDI Work/data/hotspot/ebird/"
#source(file.path(script_dir, "spatialDataDownload.R")) #if you need to download new versions of the GEDI-fusion data
spp <- c("piwo", "nofl", "dowo", "whwo", "hawo", "lewo", "wisa", "rbsa",
         "attw", "rnsa", "bbwo", "acwo")
cname <- c("Pileated Woodpecker", "Northern Flicker", "Downy Woodpecker",
           "White-headed Woodpecker","Hairy Woodpecker", 
           "Lewis's Woodpecker", "Williamson's Sapsucker", 
           "Red-breasted Sapsucker", "American Three-toed Woodpecker",
           "Red-naped Sapsucker", "Black-backed Woodpecker", "Acorn Woodpecker")
piwoData <- read.csv(paste0(script_dir, "piwo_data.csv"))
drop_cols <- c("scientific_name", "observation_count")
woodpecker_data <- piwoData %>%
  rename(piwo = species_observed) %>%
  select(!all_of(drop_cols)) 

for(p in 2:length(spp)){
  print(spp[p])
  sp <- spp[p]
  cn <- cname[p]
  bname <- paste0(script_dir, sp, "_data.csv")
  birdDat <- read.csv(bname)
  spDat <- data.frame(checklist_id = birdDat$checklist_id)
  spDat[, sp] <- birdDat$species_observed
  woodpecker_data <- merge(woodpecker_data, spDat, by = "checklist_id") 
} 

data_name <- paste0(script_dir, "all_ebird.csv")
write.csv(woodpecker_data, data_name, row.names = FALSE)
