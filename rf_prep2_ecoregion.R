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

#set sp
spp <- "piwo"
spp <- c("piwo", "nofl", "dowo", "whwo", "bbwo", "wisa", "attw", "acwo", "hawo", "lewo", "rbsa", "rnsa")
ecoregion_list <- c("mwcf")#, "mwcf", "grpl", "nade", "nwmf")
for(i in 10:length(spp)){
  sp <- spp[i] 
  print(sp)
  for(er in 1:length(ecoregion_list)){
    print(ecoregion_list[er])
  #set names
    train_name <- paste0(spp[i], "_", ecoregion_list[er], "_train.csv")
    test_name <- paste0(spp[i], "_", ecoregion_list[er], "_test.csv")
    
    # Load  Directories
    baseDir <- "C:/Users/lisah/Documents/GEDI Work/data/hotspot"
    dDir <- "D:/hotspot"
    dataDir <- file.path(baseDir, "data")
    rfname <- paste0("RF_models/", sp)
    mdlDir <- file.path(baseDir, rfname)
    
    if(file.exists(file.path(dDir, train_name)) == TRUE){
      #Read in data
      spDat <- read.csv(file.path(dDir, train_name))
      head(spDat)
      
      #Are any variables highly correlated?
      predList500 <- names(spDat)[14:(length(spDat)-1)] #predVars
      cormat <- cor(spDat[,predList500], use = "complete.obs")
      cor_name <- paste0(sp, "_CorMatrix_rf.csv")
      write.csv(cormat, file=paste0(dataDir, cor_name))
      corDf <-as.data.frame(cormat)
      #corList <- data.frame(row = rep(NA, nrow(corDf)*(length(corDf)-1)))
      #n=0
      #for(r in 1:nrow(corDf)){
      #  for(c in 2:ncol(corDf)){
      #    n <- n + 1
      #    corList$row[n] <- row.names(corDf[r])
      #    corList$col[n] <- colnames(corDf[c])
      #    corList$cor[n] <- corDf[r,c]
      #  }
      #}
      #corList
      
      #corDf$row <- row.names(corDf)
      long_corDf <- corDf %>%
        rownames_to_column("var1") %>%
        pivot_longer(cols = all_of(predList500),
                     names_to = c( "var2"),
                     values_to = "cor")
      
      cors80 <- subset(long_corDf, abs(long_corDf$cor) >= 0.80 & long_corDf$cor < 1)
      
      
      for(r in 1:nrow(cors80)){
        frm1 <- paste("I(log(species_observed +0.1)) ~", cors80[r,"var1"])
        frm2 <- paste("I(log(species_observed +0.1)) ~", cors80[r,"var2"])
        cors80$rAIC[r] <- AIC(lm(frm1, data=spDat))
        cors80$cAIC[r] <- AIC(lm(frm2, data=spDat))
        cors80$best[r] <- which.min(c(cors80$rAIC[r], cors80$cAIC[r]))
      }
      tab80 <- cors80[order(cors80$cor),]
      
      
      tab80_b <- tab80 %>% 
        # select only the columns to be used in the model
        filter(substr(var1,1,1)== "b" & substr(var2,1,1) == "b")%>%
        filter(best == 1)%>%
        print(n=40)
      not_list <- as.numeric(unique(substr(tab80_b$var2, 9, 10)))
      b_list <- c(1:19)
      b_list <- b_list[which(!b_list %in% not_list)]
      biovar_list <- paste0("biovars_", b_list)
      
      
      #tab80_b540 <- tab80_b %>% 
        # select only the columns to be used in the model
        #filter(substr(var1,11,13)== "540" | substr(var1,12,14)== "540") %>%
        #filter(substr(var2,11,13)== "540" | substr(var2,12,14)== "540") %>%
        #filter(best == 1) %>%
        #print(n=40)
      
      
      
      tab80_e <- tab80 %>% 
        # select only the columns to be used in the model
        filter(substr(var1,1,1)== "e")%>%
        filter(best == 1)#%>%
        #print(n=40)
      not_elist <- as.numeric(unique(substr(tab80_e$var2, 20, 23)))
      sp_elist <- as.numeric(unique(substr(tab80_e$var1, 20, 23)))
      e_list <- c(90, 180, 360, 540, 1530, 2520)
      e_list <- e_list[which(!e_list %in% not_elist)]
      e_list <- e_list [which(e_list %in% sp_elist)]
      e_list <- paste0("evergreen_combined_", e_list)
      
      tab80_d <- tab80 %>% 
        # select only the columns to be used in the model
        filter(substr(var1,1,1)== "d")%>%
        filter(best == 1)# %>%
        #print(n=40)
      not_dlist <- as.numeric(unique(substr(tab80_d$var2, 20, 23)))
      sp_dlist <- as.numeric(unique(substr(tab80_d$var1, 20, 23)))
      d_list <- c(90, 180, 360, 540, 1530, 2520)
      d_list <- d_list[which(!d_list %in% not_dlist)]
      d_list <- d_list [which(d_list %in% sp_dlist)]
      d_list <- paste0("deciduous_combined_", d_list)
      
      tab80_c <- tab80 %>% 
        # select only the columns to be used in the model
        filter(substr(var1,1,1)== "c")%>%
        filter(best == 1)#%>%
        #print(n=40)
      
      tab80_f <- tab80 %>% 
        # select only the columns to be used in the model
        filter(substr(var1,1,1)== "f") %>%
        filter(substr(var2,1,1)== "f") %>%
        filter(best == 1)
      
      
      tab80_r <- tab80 %>% 
        # select only the columns to be used in the model
        filter(substr(var1,1,1)== "r") %>%
        filter(substr(var2,1,1)== "r") %>%
        filter(best == 1)
      
        
      
      #predVars <- c(paste0("biovars_", 1:19, "_540"))
      
      
      #piwo: spDat_bioclim <- spDat[,c("species_observed", "biovars_1", "biovars_3", "biovars_4", "biovars_9", "biovars_10", "biovars_12", "biovars_15", "biovars_18")]
      spDat_bioclim <- spDat[,c("species_observed", biovar_list)]
      #spDat_bioclim <- spDat[,c("species_observed", "biovars_1_540","biovars_3_540", "biovars_5_540","biovars_9_540", "biovars_12_540", "biovars_15_540", "biovars_17_540")]
      detection_freq <- mean(spDat_bioclim$species_observed)
      detection_freq
      
      spDat_bioclim$species_observed <- factor(spDat$species_observed) # ranger requires a factor response to do classification
      
      
      # grow random forest
      library(ranger)
      set.seed(2023)
      rf <- ranger(formula =  species_observed ~ ., 
                   data = spDat_bioclim,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE, 
                   sample.fraction = c(detection_freq, detection_freq))
      
      rf$variable.importance
      #check that VIF are <10
      fun <- paste("species_observed ~ ", paste(biovar_list, collapse = " + "))
      model <- lm(fun,
                  data = spDat)
      summary(model)
      vif(model)
      print(vif(model))
      if(any(vif(model) >= 10)){
        print("vif>10")
        print(paste0("drop: ", names(which.min(rf$variable.importance))))
        biovar_list <- biovar_list[which(!biovar_list %in% names(which.min(rf$variable.importance)))]
        fun <- paste("species_observed ~ ", paste(biovar_list, collapse = " + "))
        model <- lm(fun, data = spDat)
        summary(model)
        vif(model)
        print(vif(model))
        if(any(vif(model) >= 10)){
          print("vif>10_still")
          imp <- rf$variable.importance[-c(which.min(rf$variable.importance))]
          print(paste0("drop: ", names(which.min(imp))))
          biovar_list <- biovar_list[which(!biovar_list %in% names(which.min(imp)))]
          fun <- paste("species_observed ~ ", paste(biovar_list, collapse = " + "))
          model <- lm(fun, data = spDat)
          summary(model)
          vif(model)
          print(vif(model))
          if(any(vif(model) >= 10)){
            print("vif>10_still2")
            imp2 <- imp[-c(which.min(imp))]
            print(paste0("drop: ", names(which.min(imp2))))
            biovar_list <- biovar_list[which(!biovar_list %in% names(which.min(imp2)))]
            fun <- paste("species_observed ~ ", paste(biovar_list, collapse = " + "))
            model <- lm(fun, data = spDat)
            summary(model)
            vif(model)
            print(vif(model))
            if(any(vif(model) >= 10)){
              print("vif>10_still3")
              imp3 <- imp2[-c(which.min(imp2))]
              print(paste0("drop: ", names(which.min(imp3))))
              biovar_list <- biovar_list[which(!biovar_list %in% names(which.min(imp3)))]
              fun <- paste("species_observed ~ ", paste(biovar_list, collapse = " + "))
              model <- lm(fun, data = spDat)
              summary(model)
              vif(model)
              print(vif(model))
              if(any(vif(model) >= 10)){
                print("vif>10_still4")
                imp4 <- imp3[-c(which.min(imp3))]
                print(paste0("drop: ", names(which.min(imp4))))
                biovar_list <- biovar_list[which(!biovar_list %in% names(which.min(imp4)))]
                fun <- paste("species_observed ~ ", paste(biovar_list, collapse = " + "))
                model <- lm(fun, data = spDat)
                summary(model)
                vif(model)
                print(vif(model))
                if(any(vif(model) >= 10)){
                  print("vif>10_still5")
                  imp5 <- imp4[-c(which.min(imp4))]
                  print(paste0("drop: ", names(which.min(imp5))))
                  biovar_list <- biovar_list[which(!biovar_list %in% names(which.min(imp5)))]
                  fun <- paste("species_observed ~ ", paste(biovar_list, collapse = " + "))
                  model <- lm(fun, data = spDat)
                  summary(model)
                  vif(model)
                  print(vif(model))
                }
              }
            }
          }
        }
      }
    ###########################################################################
    vars_for_all <- c("species_observed",
                      "year", "day_of_year",
                      "time_observations_started", "duration_minutes", 
                      "number_observers")
    
    ########################################################
    #write the buffer distances for each species
    if(sp == "piwo" | sp == "whwo" | sp == "bbwo" ){
      sp_buf <- c(540, 1530, 2520)
    } else{
      if(sp == "nofl" | sp == "dowo" | sp == "wisa" | sp == "attw" | sp == "acwo"){
        sp_buf <- c(180, 360, 540)
      } else{
        if(sp == "hawo" | sp== "lewo" | sp == "rbsa"){
          sp_buf <- c(90, 180, 360)
        } else {
          if(sp == "rnsa"){
            sp_buf <- c(90, 180, 360, 540)
          }
        }
      }
    }
    
    
    #######################
    #create lists of variables for cover, fhd, and rh98
    c_list <- paste0("cover_", sp_buf)
    f_list <- paste0("fhd_", sp_buf)
    r_list <- paste0("rh98_", sp_buf)
    
    ####################################################
    #create full variable list for each species
    vars_list <- c(vars_for_all, e_list, d_list, biovar_list, c_list, f_list, r_list)
    vars_name <- paste0(dDir, "/", sp, "_", ecoregion_list[er], "_vars_list3.csv")
    write.csv(vars_list, vars_name, row.names = FALSE)
  } else{
    print("No file for sp. (too few obs), skipping")
    }
  }
}
