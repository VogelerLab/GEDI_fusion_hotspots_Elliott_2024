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
library(car)
library(terra)
library(forcats)

#create empty dataframe
pd_all <- data.frame()
#set ecoregion
ecoregion <- c("mwcf")

spp <- c("piwo", "nofl", "dowo", "hawo",  "rbsa", "acwo", "attw")
sppCAP <- c("PIWO", "NOFL", "DOWO", "HAWO",  "RBSA", "ACWO", "ATTW")
#"rnsa", "bbwo", "whwo", "lewo", "wisa", )
vars_list <- c("rh98_90", "rh98_180",  "rh98_360", "rh98_540", "rh98_1530", "rh98_2520",
               "cover_90", "cover_180",  "cover_360", "cover_540", "cover_1530", "cover_2520",
               "fhd_90", "fhd_180",  "fhd_360", "fhd_540", "fhd_1530", "fhd_2520",
               "deciduous_combined_540", "deciduous_combined_1530", "deciduous_combined_2520",
               "evergreen_combined_540", "evergreen_combined_1530", "evergreen_combined_2520",
               "biovars_1", "biovars_2", "biovars_3", "biovars_4", "biovars_11",
               "biovars_12", "biovars_14", "biovars_15", "biovars_17", "biovars_18", "biovars_19",
               "duration_minutes", "time_observations_started", "day_of_year", "year", "number_observers"
) 


for(i in 1:length(spp)){
  sp <- spp[i]
  spCAP <- sppCAP[i]
  
  #set names
  train_name <- paste0(sp, "_", ecoregion, "_train.csv")
  test_name <- paste0(sp, "_", ecoregion, "_test.csv")
  #vars_name <- paste0(sp, "_vars_list3.csv")
  dDir <- "D:/hotspot"
  vars_name <- paste0(dDir, "/", sp, "_", ecoregion, "_vars_list3.csv")
  
  # Load  Directories
  #baseDir <- "C:/Users/lisah/Documents/GEDI Work/data/hotspot"
  baseDir <- "C:/Users/lisah/Documents/GEDI Work"
  dataDir <- file.path(baseDir, "data/hotspot/data")
  rfname <- paste0("RF_models/", sp)
  mdlDir <- file.path(baseDir, rfname)
  outputDir <- file.path(baseDir, "hotspot analysis/output2")
  
  
  #Read in data
  spDat <- read.csv(file.path(dDir, train_name))
  head(spDat)
  vars_list <- read.csv(file.path(vars_name))
  
  #subset to important data
  #spDat_subset <- spDat %>%
  # select(all_of(vars_list))
  
  
  spDat_subset <- spDat[, names(spDat) %in% vars_list$x]
  
  data <- spDat_subset
  
  
  cor(data)
  
  detection_freq <- mean(spDat_subset$species_observed)
  detection_freq #0.04339542
  
  #The actual rf analysis
  spDat_subset$species_observed <- factor(spDat_subset$species_observed)
  set.seed(625)
  rf <- ranger(formula =  species_observed ~ ., 
               data = spDat_subset,
               importance = "impurity",
               probability = TRUE,
               replace = TRUE, 
               sample.fraction = c(detection_freq, detection_freq))
  
  
  
  
  #Calibration:
  library(scam)
  library(dplyr)
  # make predictions on training data
  occ_pred <- rf$predictions[, 2]
  # convert the observered response back to a numeric value from factor
  occ_obs <- spDat$species_observed %>% 
    as.logical() %>% 
    as.integer()
  rf_pred_train <- tibble(obs = occ_obs, pred = occ_pred) #%>% 
  #drop_na()
  
  # fit calibration model
  calibration_model <- scam(obs ~ s(pred, k = 5, bs = "mpi"), 
                            gamma = 1.4,
                            data = rf_pred_train)
  
  # calculate the average observed encounter rates for different 
  # categories of estimated encounter rates 
  
  average_encounter <- rf_pred_train %>%
    mutate(pred_cat = cut(rf_pred_train$pred, breaks = seq(0, 1, by=0.02))) %>%
    group_by(pred_cat) %>%
    summarise(pred = mean(pred), obs = mean(obs), checklist_count = n()) %>%
    ungroup()
  
  # plot
  library(ggplot2)
  cal_pred <- tibble(pred = seq(0, 1, length.out = 100))
  cal_pred <- predict(calibration_model, cal_pred, type = "response") %>% 
    bind_cols(cal_pred, calibrated = .)
  cal_plot <- paste0("cal_plot_", sp, "_", ecoregion, ".png")
  ggplot(cal_pred) +
    aes(x = pred, y = calibrated) +
    geom_line() +
    theme_classic() +
    geom_point(data = average_encounter, 
               aes(x = pred, y = obs, size = sqrt(checklist_count)),
               show.legend = FALSE, shape = 1) +
    labs(x = "Estimated encounter rate",
         y = "Observed encounter rate",
         title = "Calibration model-Base Model")
  
  #ggsave(cal_plot)
  
  #Assessment:
  
  #Read in data
  spTest <- read.csv(file.path(dDir, test_name))
  spTest <- spTest[, names(spTest) %in% vars_list$x]
  
  # predict on test data using calibrated model
  p_fitted <- predict(rf, data = spTest, type = "response")
  # extract probability of detection
  p_fitted <- p_fitted$predictions[, 2]
  # calibrate
  p_calibrated <- predict(calibration_model, 
                          newdata = tibble(pred = p_fitted), 
                          type = "response")
  rf_pred_test <- data.frame(id = seq_along(p_calibrated),
                             # actual detection/non-detection
                             obs = spTest$species_observed,
                             # uncalibrated prediction
                             fit = p_fitted,
                             # calibrated prediction
                             cal = p_calibrated) %>%
    # constrain probabilities to 0-1
    mutate(cal = pmin(pmax(cal, 0), 1)) #%>% 
  #drop_na()
  unique(is.na(rf_pred_test$cal))
  # mean squared error (mse)
  mse_fit <- mean((rf_pred_test$obs - rf_pred_test$fit)^2, na.rm = TRUE)
  mse_cal <- mean((rf_pred_test$obs - rf_pred_test$cal)^2, na.rm = TRUE)
  
  # pick threshold to maximize kappa
  library(PresenceAbsence)
  select <- dplyr::select
  opt_thresh <- optimal.thresholds(rf_pred_test, opt.methods = "MaxKappa")
  
  # calculate accuracy metrics: auc, kappa, sensitivity, specificity,
  metrics_fit <- rf_pred_test %>% 
    select(id, obs, fit) %>% 
    presence.absence.accuracy(threshold = opt_thresh$fit, 
                              na.rm = TRUE, 
                              st.dev = FALSE)
  metrics_cal <- rf_pred_test %>% 
    select(id, obs, cal) %>% 
    presence.absence.accuracy(threshold = opt_thresh$cal, 
                              na.rm = TRUE, 
                              st.dev = FALSE)
  
  rf_assessment <- tibble(
    species = sp,
    model = c("RF", "Calibrated RF"),
    mse = c(mse_fit, mse_cal),
    sensitivity = c(metrics_fit$sensitivity, metrics_cal$sensitivity),
    specificity = c(metrics_fit$specificity, metrics_cal$specificity),
    auc = c(metrics_fit$AUC, metrics_cal$AUC)
  )
  knitr::kable(rf_assessment, digits = 3)
  at_name <- paste0(outputDir, "/assessment_table_", sp, "_", ecoregion, "_.csv")
  #write.csv(rf_assessment, at_name, row.names = FALSE)
  
  #Variable Importance
  library(tibble)
  library(forcats)
  pi <- enframe(rf$variable.importance, "predictor", "importance")
  # plot
  ip_name <- paste0(outputDir, "/Importance_Plot_", sp, "_", ecoregion, ".png")
  ip <- ggplot(pi) + 
    aes(x = fct_reorder(predictor, importance), y = importance) +
    geom_col() +
    geom_hline(yintercept = 0, size = 2, colour = "#555555") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    labs(x = NULL, 
         y = "Predictor Importance (Gini Index)") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.grid.major.x = element_line(colour = "#cccccc", size = 0.5)) +
    theme(plot.title.position = "plot") #NEW parameter
  #ggsave(ip_name)
  #NOTE THAT YOU WILL NEED TO REDUCE DATASET TO ONLY THOSE COLUMNS OF INTEREST!
  
  dowo_ip <- ip
  
  # top 9 predictors other than date
  top_pred <- pi %>% 
    filter(!predictor %in% c("year", "day_of_year", "state_code", "number_observers")) %>% 
    top_n(n = 9, wt = importance) %>% 
    arrange(desc(importance))
  
  #Partial Dependence
  # function to calculate partial dependence for a single predictor
  calculate_pd <- function(predictor, model, data, 
                           x_res = 25, n = 1000) {
    # create prediction grid
    rng <- range(data[[predictor]], na.rm = TRUE)
    x_grid <- seq(rng[1], rng[2], length.out = x_res)
    grid <- data.frame(covariate = predictor, x = x_grid, 
                       stringsAsFactors = FALSE)
    names(grid) <- c("covariate", predictor)
    
    # subsample training data
    n <- min(n, nrow(data))
    s <- sample(seq.int(nrow(data)), size = n, replace = FALSE)
    data <- data[s, ]
    
    # drop focal predictor from data
    data <- data[names(data) != predictor]
    grid <- merge(grid, data, all = TRUE)
    
    # predict
    p <- predict(model, data = grid)
    
    # summarize
    pd <- grid[, c("covariate", predictor)]
    names(pd) <- c("covariate", "x")
    pd$pred <- p$predictions[, 2]
    pd <- dplyr::group_by(pd, covariate, x) %>% 
      dplyr::summarise(pred = mean(pred, na.rm = TRUE)) %>% 
      dplyr::ungroup()
    
    return(pd)
  }
  
  # calculate partial dependence for each predictor
  # map is used to iteratively apply calculate_pd to each predictor
  library(tidyr)
  
  spDat2 <- spDat_subset
  
  pd <- top_pred %>%
    mutate(pd = map(predictor, calculate_pd, model = rf, 
                    data = spDat),
           pd = map(pd, ~ .[, c(2, 3)]),
           pd = map(pd, set_names, nm = c("value",  "encounter_rate"))) %>% 
    unnest(cols = pd)
  
  #pd <- top_pred
  #for (p in 1:nrow(pd)){ #!!!NEED TO LINK pd with calculate_pd
  #calculate_pd(top_pred$predictor[p], model = rf, 
  # data = noflDat)  
  #}
  
  # calibrate predictions
  pd$encounter_rate <- predict(calibration_model, 
                               newdata = tibble(pred = pd$encounter_rate), 
                               type = "response") %>% 
    as.numeric()
  
  
  #Variable Importance
  library(tibble)
  library(forcats)
  
  pd_sp <- pd %>%
    mutate("Species" = spCAP) %>%
    mutate("encounter_rate_lat" = encounter_rate / sum(encounter_rate, na.rm = TRUE))
  
  pd_all <- rbind(pd_all, pd_sp) 
}  

###############################################
# plot


pd_all_grp <- pd_all %>%
  mutate(Category = case_when(
    startsWith(predictor, "biovars") ~ "Bioclimatic",
    startsWith(predictor, "fhd") ~ "GEDI",
    startsWith(predictor, "rh") ~ "GEDI",
    startsWith(predictor, "cover") ~ "GEDI",
    startsWith(predictor, "ever") ~ "Forest Type",
    startsWith(predictor, "decid") ~ "Forest Type",
    startsWith(predictor, "year") ~ "Survey Particulars",
    startsWith(predictor, "day") ~ "Survey Particulars",
    startsWith(predictor, "time") ~ "Survey Particulars",
    startsWith(predictor, "duration") ~ "Survey Particulars",
    startsWith(predictor, "number") ~ "Survey Particulars"
    
  ))

pd_all_gedi <- pd_all_grp %>%
  filter(Category == "GEDI") %>%
  mutate(Bin = case_when(
    startsWith(predictor, "fhd") ~ "fhd (unitless)",
    startsWith(predictor, "rh") ~ "rh98 (m)",
    startsWith(predictor, "cover") ~ "cover (proportional)")) %>%
  mutate(Buffer = case_when(
    endsWith(predictor, "90") ~ "90",
    endsWith(predictor, "180") ~ "180",
    endsWith(predictor, "360") ~ "360",
    endsWith(predictor, "540") ~ "540",
    endsWith(predictor, "1530") ~ "1530",
    endsWith(predictor, "2520") ~ "2520"
    ))

pd_all_gedi$predictor <- factor(pd_all_gedi$predictor, levels= vars_list)
pd_all_gedi$Species <- factor(pd_all_gedi$Species, levels= sppCAP)
pd_all_gedi$Bin <- factor(pd_all_gedi$Bin, 
                               levels= c("rh98 (m)", "cover (proportional)", "fhd (unitless)"))
pd_all_gedi$Buffer <- factor(pd_all_gedi$Buffer, 
                          levels= c("90", "180", "360", "540", "1530", "2520"))
###
pdp_name <- paste0(outputDir, "/Partial_Dep_Plot_", sp, "_", ecoregion, ".png")

dep <- ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate, color = Species, linetype = Buffer) +
  geom_line () +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ as_factor(Bin), nrow = 3, scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  scale_colour_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                 "#44AA99")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

#v2
ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate_lat, color = Species, linetype = Buffer) +
  geom_line () +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ as_factor(Bin), nrow = 3, scales = "free") +
  labs(x = NULL, y = "Encounter Rate (Relativized)") +
  theme_minimal() +
  theme_minimal() +
  scale_colour_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))
ggsave(pdp_name)

#v3

ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate, color = Species, linetype = Buffer) +
  geom_line () +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Species ~ as_factor(Bin), scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  scale_colour_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

#v3.2
ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate, color = Species, linetype = Buffer) +
  geom_line () +
  geom_point(aes(shape = Buffer)) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Species ~ as_factor(Bin), scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  scale_colour_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

#v3.3
library(paletteer)
paletteer_c("ggthemes::Blue-Green Sequential", 7)
paletteer_c("grDevices::YlGnBu", 14)
ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate, color = Buffer, linetype = Buffer) +
  geom_line () +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Species ~ as_factor(Bin), scales = "free") +
  labs(y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  scale_colour_manual(values=c("#26185FFF", "#00579AFF", "#0076AFFF", 
                               "#00A3AEFF", "#4BC5B1FF", "#C3ECC7FF")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

#v3.3
library(paletteer)
paletteer_c("ggthemes::Blue-Green Sequential", 7)
paletteer_c("grDevices::YlGnBu", 14)
ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate, color = Buffer, linetype = Buffer) +
  geom_line () +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Species ~ as_factor(Bin), scales = "free") +
  labs(y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  scale_colour_manual(values=c("#C3ECC7FF", "#4BC5B1FF", "#00A3AEFF",  
                               "#0076AFFF", "#00579AFF","#26185FFF")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))


ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate, color = Species, linetype = Buffer) +
  geom_line () +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Species ~ as_factor(Bin), scales = "free") +#, switch = "x") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme(strip.placement = "outside") +
  theme_minimal() +
  scale_colour_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))
#v4
pd_sub <- pd_all_gedi %>%
  filter(Species == "ATTW")
my_strip_labels <- as_labeller(c(
  "rh98" = "rh98 (m)",
  "cover" = "cover (proportional)",
  "fhd" = "fhd (unitless)"
))
ggplot(pd_all_gedi) +
  aes(x = value, y = encounter_rate, color = Species, linetype = Buffer) +
  geom_line () +
  geom_point() +
  #geom_rug(data= pd_sub, sides="b") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(Species ~ as_factor(Bin), scales = "free", labeller = my_strip_labels #,
             #labeller = labeller(factor_name= as_labeller(c("rh98" = "rh98 (m)", "cover" = "cover (proportional)", "fhd" = "fhd (unitless)")))
  ) +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  scale_colour_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                               "#44AA99")) +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))


#Version 1 figure

# plot
ip_name <- paste0(outputDir, "/Importance_Plot_", sp, "_", ecoregion, ".png")

pi_lat_all_grp$predictor <- factor(pi_lat_all_grp$predictor, levels= vars_list)
pi_lat_all_grp$Category <- factor(pi_lat_all_grp$Category, 
                                  levels= c("GEDI", "Forest Type", "Bioclimatic", "Survey Particulars"))


ip <- ggplot(pi_lat_all_grp, aes(x = forcats::fct_rev(predictor), y = Importance, fill = species)) +
  geom_bar(stat='identity')+
  #geom_col() +
  geom_hline(yintercept = 0, size = 2) +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  facet_wrap(~Category, scales = "free") +
  labs(x = NULL, 
       y = "Predictor Importance") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", size = 0.5)) +
  theme(plot.title.position = "plot") #NEW parameter
#ggsave(ip_name)






#Version 2 Figure
pi_lat_all_bins <- pi_lat_all %>%
  complete(Species, predictor,
           fill = list(Gini = 0, Importance = 0)) %>%
  mutate(Bins = case_when(
    startsWith(predictor, "biovars") ~ "Non-GEDI",
    startsWith(predictor, "fhd") ~ "GEDI",
    startsWith(predictor, "rh") ~ "GEDI",
    startsWith(predictor, "cover") ~ "GEDI",
    startsWith(predictor, "ever") ~ "Non-GEDI",
    startsWith(predictor, "decid") ~ "Non-GEDI",
    startsWith(predictor, "year") ~ "Non-GEDI",
    startsWith(predictor, "day") ~ "Non-GEDI",
    startsWith(predictor, "time") ~ "Non-GEDI",
    startsWith(predictor, "duration") ~ "Non-GEDI",
    startsWith(predictor, "number") ~ "Non-GEDI"))
  
pi_lat_all_bins$predictor <- factor(pi_lat_all_bins$predictor, levels= vars_list)
pi_lat_all_bins$Species <- factor(pi_lat_all_bins$Species, levels= sppCAP)
pi_lat_all_bins$Bins <- factor(pi_lat_all_bins$Bins, 
                                  levels= c("GEDI", "Non-GEDI"))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99")#, "#999933", "#882255", "#661100", "#6699CC", "#888888")

ip <- ggplot(pi_lat_all_bins, aes(x = forcats::fct_rev(predictor), y = Importance, fill = Species)) +
  geom_bar(stat='identity', position = "dodge")+
  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99")) +
  #geom_col() +
  geom_hline(yintercept = 0, size = 2) +
  #scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  facet_wrap(~Bins, scales = "free") +
  labs(x = NULL, 
       y = "Predictor Importance") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", size = 0.5)) +
  theme(plot.title.position = "plot") #NEW parameter
#ggsave(ip_name)



# top 9 predictors other than date
top_pred <- pi %>% 
  filter(!predictor %in% c("year", "day_of_year", "state_code", "number_observers")) %>% 
  top_n(n = 9, wt = importance) %>% 
  arrange(desc(importance))

#Partial Dependence
# function to calculate partial dependence for a single predictor
calculate_pd <- function(predictor, model, data, 
                         x_res = 25, n = 1000) {
  # create prediction grid
  rng <- range(data[[predictor]], na.rm = TRUE)
  x_grid <- seq(rng[1], rng[2], length.out = x_res)
  grid <- data.frame(covariate = predictor, x = x_grid, 
                     stringsAsFactors = FALSE)
  names(grid) <- c("covariate", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  s <- sample(seq.int(nrow(data)), size = n, replace = FALSE)
  data <- data[s, ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict
  p <- predict(model, data = grid)
  
  # summarize
  pd <- grid[, c("covariate", predictor)]
  names(pd) <- c("covariate", "x")
  pd$pred <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, covariate, x) %>% 
    dplyr::summarise(pred = mean(pred, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(pd)
}

# calculate partial dependence for each predictor
# map is used to iteratively apply calculate_pd to each predictor
library(tidyr)

spDat2 <- spDat_subset

pd <- top_pred %>%
  mutate(pd = map(predictor, calculate_pd, model = rf, 
                  data = spDat),
         pd = map(pd, ~ .[, c(2, 3)]),
         pd = map(pd, set_names, nm = c("value",  "encounter_rate"))) %>% 
  unnest(cols = pd)

#pd <- top_pred
#for (p in 1:nrow(pd)){ #!!!NEED TO LINK pd with calculate_pd
  #calculate_pd(top_pred$predictor[p], model = rf, 
              # data = noflDat)  
#}

# calibrate predictions
pd$encounter_rate <- predict(calibration_model, 
                             newdata = tibble(pred = pd$encounter_rate), 
                             type = "response") %>% 
  as.numeric()

# plot
pdp_name <- paste0(outputDir, "/Partial_Dep_Plot_", sp, "_", ecoregion, ".png")

dep <- ggplot(pd) +
  aes(x = value, y = encounter_rate) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ as_factor(predictor), nrow = 3, scales = "free") +
  labs(x = NULL, y = "Encounter Rate",
       title = "A) PIWO") +
  theme_minimal() +
  theme_minimal() +
  theme(plot.title.position = "plot") +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))
ggsave(pdp_name)

piwo_dep

#rf_name <- paste0("D:/hotspot/output/rf_", sp, "_", ecoregion, ".rds")
#calib_name <- paste0("D:/hotspot/output/calibration_model_", sp, "_", ecoregion, ".rds") 
#rf <- saveRDS(rf, rf_name)
#calibration_model <- saveRDS(calibration_model, calib_name)

library(ggpubr)
ggarrange(piwo_ip, nofl_ip, dowo_ip + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)


pdp_name <- paste0(outputDir, "/Partial_Dep_Plots_", ecoregion, ".png")


,  common.legend = TRUE



########



#Calibration:
library(scam)
library(dplyr)
# make predictions on training data
occ_pred <- rf$predictions[, 2]
# convert the observered response back to a numeric value from factor
occ_obs <- spDat$species_observed %>% 
  as.logical() %>% 
  as.integer()
rf_pred_train <- tibble(obs = occ_obs, pred = occ_pred) #%>% 
#drop_na()

# fit calibration model
calibration_model <- scam(obs ~ s(pred, k = 5, bs = "mpi"), 
                          gamma = 1.4,
                          data = rf_pred_train)

# calculate the average observed encounter rates for different 
# categories of estimated encounter rates 

average_encounter <- rf_pred_train %>%
  mutate(pred_cat = cut(rf_pred_train$pred, breaks = seq(0, 1, by=0.02))) %>%
  group_by(pred_cat) %>%
  summarise(pred = mean(pred), obs = mean(obs), checklist_count = n()) %>%
  ungroup()

# plot
library(ggplot2)
cal_pred <- tibble(pred = seq(0, 1, length.out = 100))
cal_pred <- predict(calibration_model, cal_pred, type = "response") %>% 
  bind_cols(cal_pred, calibrated = .)
cal_plot <- paste0("cal_plot_", sp, "_", ecoregion, ".png")
ggplot(cal_pred) +
  aes(x = pred, y = calibrated) +
  geom_line() +
  theme_classic() +
  geom_point(data = average_encounter, 
             aes(x = pred, y = obs, size = sqrt(checklist_count)),
             show.legend = FALSE, shape = 1) +
  labs(x = "Estimated encounter rate",
       y = "Observed encounter rate",
       title = "Calibration model-Base Model")

#ggsave(cal_plot)

#Assessment:

#Read in data
spTest <- read.csv(file.path(dDir, test_name))
spTest <- spTest[, names(spTest) %in% vars_list$x]

# predict on test data using calibrated model
p_fitted <- predict(rf, data = spTest, type = "response")
# extract probability of detection
p_fitted <- p_fitted$predictions[, 2]
# calibrate
p_calibrated <- predict(calibration_model, 
                        newdata = tibble(pred = p_fitted), 
                        type = "response")
rf_pred_test <- data.frame(id = seq_along(p_calibrated),
                           # actual detection/non-detection
                           obs = spTest$species_observed,
                           # uncalibrated prediction
                           fit = p_fitted,
                           # calibrated prediction
                           cal = p_calibrated) %>%
  # constrain probabilities to 0-1
  mutate(cal = pmin(pmax(cal, 0), 1)) #%>% 
#drop_na()
unique(is.na(rf_pred_test$cal))
# mean squared error (mse)
mse_fit <- mean((rf_pred_test$obs - rf_pred_test$fit)^2, na.rm = TRUE)
mse_cal <- mean((rf_pred_test$obs - rf_pred_test$cal)^2, na.rm = TRUE)

# pick threshold to maximize kappa
library(PresenceAbsence)
select <- dplyr::select
opt_thresh <- optimal.thresholds(rf_pred_test, opt.methods = "MaxKappa")

# calculate accuracy metrics: auc, kappa, sensitivity, specificity,
metrics_fit <- rf_pred_test %>% 
  select(id, obs, fit) %>% 
  presence.absence.accuracy(threshold = opt_thresh$fit, 
                            na.rm = TRUE, 
                            st.dev = FALSE)
metrics_cal <- rf_pred_test %>% 
  select(id, obs, cal) %>% 
  presence.absence.accuracy(threshold = opt_thresh$cal, 
                            na.rm = TRUE, 
                            st.dev = FALSE)

rf_assessment <- tibble(
  species = sp,
  model = c("RF", "Calibrated RF"),
  mse = c(mse_fit, mse_cal),
  sensitivity = c(metrics_fit$sensitivity, metrics_cal$sensitivity),
  specificity = c(metrics_fit$specificity, metrics_cal$specificity),
  auc = c(metrics_fit$AUC, metrics_cal$AUC)
)
knitr::kable(rf_assessment, digits = 3)
at_name <- paste0(outputDir, "/assessment_table_", sp, "_", ecoregion, "_.csv")
#write.csv(rf_assessment, at_name, row.names = FALSE)
