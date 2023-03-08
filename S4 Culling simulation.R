# Culling simulation ------------

# Load required package --------
packages <- c("dplyr", "svMisc","igraph", "readxl")
lapply(packages, library, character.only = TRUE)

# load required file -----
# farm input data
dat_kernel <- readRDS("datLP.rds")
# Distance between farm matrix 
Mod_distance <- readRDS("ModLP_distance.rds")
# trade matrix
component_matrix <- readRDS("LP_component_matrix.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------
## Set the parameters for simulation -----------------------

culling_sen_cond <- function(){
  
  
  # set index cases
  index_farm <- c(409,422,84) # this is index cases from kernel_LP data
  
  # transmission parameters
  k0 =   0.00537620
  r0 =  0.17187575
  alpha = 1.49815422
  delta =  0.00061939
  
  # set rate and shape for infectious duration
  shape_inf =3.0262817
  rate_inf =0.1377417
  
  #latent period
  latent = 3
  
  # calculate the kernel transmission matrix
  Mod_kernel <- (k0/(1+((Mod_distance/r0)^alpha))) 
  diag(Mod_kernel) <- 0 # diagonal = 0
  
  
  # calculate the trade network transmission matrix  
  Mod_component_matrix_delta <- component_matrix * delta
  diag(Mod_component_matrix_delta) <- 0 # diagonal = 0
  
  
  # detection time
  culling_excute = TRUE
  detection_time_index = 7
  detectiontime_cattle = 5 # for dairy cattle
  detectiontime_beef = 7 # detection time after day since infection
  detectiontime_goat = 18
  detectiontime_pig = 8
  
  # set detection time in dat_kernel dataframe
  dat_kernel$detection_time <-NA
  dat_kernel$detection_time[dat_kernel$type3=="index"]<- detection_time_index
  dat_kernel$detection_time[dat_kernel$type3=="goat"]<-  detectiontime_goat
  dat_kernel$detection_time[dat_kernel$type3=="cattle"]<-  detectiontime_cattle
  dat_kernel$detection_time[dat_kernel$type3=="pig"]<-  detectiontime_pig
  
  list( index_farm = index_farm , 
        shape_inf = shape_inf, rate_inf = rate_inf,
        latent = latent,Mod_kernel= Mod_kernel,
        Mod_component_matrix_delta=Mod_component_matrix_delta,
        culling_excute =culling_excute,
        detection_time_index =  detection_time_index,
        detectiontime_cattle = detectiontime_cattle,
        detectiontime_beef = detectiontime_beef,
        detectiontime_goat = detectiontime_goat,
        detectiontime_pig = detectiontime_pig,
        dat_kernel=dat_kernel
        
  )
}

# save all the parameters as a list 
culling_sen_list <-culling_sen_cond()

# We will test  parameter
# 1) delay_cullling
sen_kernel <- data.frame(delay_culling= c(1, 7,14))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------
## Function to sum the outbreak simulation -------------
sum_outbreak <- function (dat) {
  
  outbreak_duration <- sapply(dat , ncol)-6
  n_susceptible <- sapply(dat , function(x) sum(x[,ncol(x)] == "S" ))
  n_latent <- sapply(dat , function(x) sum(x[,ncol(x)] == "L" ))
  n_infectious <- sapply(dat , function(x) sum(x[,ncol(x)] == "I" ))
  n_recover <- sapply(dat , function(x) sum(x[,ncol(x)] == "R" ))
  n_VI <- sapply(dat , function(x) sum(x[,ncol(x)] == "VI" ))
  n_VS <- sapply(dat , function(x) sum(x[,ncol(x)] == "VS" ))
  EMvaccinated <- sapply(dat , function(x) sum(x$day_sinceEV !=0 ))
  Quarantine <- sapply(dat , function(x) sum(x$day_sinceQ >0 ))
  n_culling <- sapply(dat , function(x) sum(x[,ncol(x)] == "C" ))
  n_UI <- sapply(dat , function(x) sum(x[,ncol(x)] == "UI" ))
  
  outbreak <- data.frame (outbreak_duration, n_susceptible, n_latent, n_infectious,
                          n_recover,  n_VI, n_VS,EMvaccinated,n_culling,Quarantine,n_UI)
  return(outbreak)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------
## Simulate outbreak for isolation----------
# before using this code don't forget to load the requirement condition files
#load distance and trade matrix

simulate_sen_culling <- function (iter, cond,sen_kernel) {
  
  # create list to keep the results
  result_sim <- list()
  # for loop parameters combination
  for (r in 1:nrow(sen_kernel)){
    progress(r, max.value = nrow(sen_kernel),progress.bar = FALSE)# to track the for loop progress
    # set condition here
    list2env(cond, envir = .GlobalEnv)
    
    # modify parameters
    delay_culling = sen_kernel$delay_culling[r] 
    
    # create list to keep results
    Mod_outbreaksim <- list()
    set.seed(1234)
    # for loop for iterated simulation
    for (n in 1:iter) {
      
      # set starter conditions
      Mod_status_df <- data.frame(ID = c(1:nrow(dat_kernel)), type = dat_kernel$type3, day_sinceinf = 0,detection = dat_kernel$detection_time,
                                  day_sinceEV = 0,day_sinceQ = 0, status_D1 = "S")
      
      #infduration 
      inf_duration <- round(rgamma(nrow(dat_kernel),shape= shape_inf, rate = rate_inf))
      inf_duration <- inf_duration+latent # plus latent period
      waning_duration = round(inf_duration + (365*2.5))
      
      
      # set index cases
      index_farm =  index_farm
      Mod_status_df$status_D1 [index_farm] <- "L" # status of index case is I
      # day since infection is equal to latent in index case as it passed latent period
      Mod_status_df$day_sinceinf [index_farm] <- 1
      
      
      # set start number of infectious farm
      n_infectious = sum(Mod_status_df$status_D1 %in% c("L","UI", "I", "Q","VI"))
      
      # simulated in one iteration until no infectious farms
      while(n_infectious > 0){
        
        
        
        #i = ncol(Mod_status_df)-5 # i = equal to column number before status_day
        
        # calculate distance dependent transmission
        status_infectious <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("I","UI"))
        status_susceptible <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("S"))
        
        
        ### normal kernel ######
        # row is ID of infectious farms, column is ID of susceptible farms
        dist_trans <- colSums(Mod_kernel*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
        
        # calculate trade network transmission
        trade_trans <- colSums(Mod_component_matrix_delta*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
        
        
        # calculate p inf
        p_inf <- 1 - exp(-1*(dist_trans+trade_trans))
        # random 0 and 1 with prob. of infection
        inf <- rbinom(length(p_inf), size = 1, prob=p_inf)
        
        # create new vector to store status next day
        Mod_status_df$status <- Mod_status_df[,ncol(Mod_status_df)]  
        
        # if farm got infected
        Mod_status_df$status[which(inf>0)] = "L"
        
        ##### update status ######
        # Culling after detection + delay_culling
        Mod_status_df$status[Mod_status_df$day_sinceinf> Mod_status_df$detection+delay_culling  & Mod_status_df$status =="I"] <- "C"
        
        # plus 1 day for farm that is not susceptible
        Mod_status_df$day_sinceinf <- if_else( Mod_status_df$status %in% c("L","UI", "I", "Q", "R","VI"), Mod_status_df$day_sinceinf+1, Mod_status_df$day_sinceinf)
        # day since EV +1 in EM vaccinated farm
        Mod_status_df$day_sinceEV <-  if_else( Mod_status_df$day_sinceEV > 0, Mod_status_df$day_sinceEV+1, Mod_status_df$day_sinceEV)
        
        # undetection time for index case, it is longer than later farms
        # if day_sinceinf > 3 & day_sinceinf <= detection time, farm become undetected infectious, and not latent anymore
        Mod_status_df$status[Mod_status_df$day_sinceinf > latent & Mod_status_df$day_sinceinf < Mod_status_df$detection  & Mod_status_df$status =="L"] <- "UI"
        
        
        # detection time for index case
        # if day_sinceinf > detection time, farm become detected infection
        Mod_status_df$status[Mod_status_df$day_sinceinf== Mod_status_df$detection & Mod_status_df$day_sinceinf <= inf_duration & Mod_status_df$status =="UI"] <- "I"
        
        
        # if day_sinceinf > inf_duration, it is recovered. 
        # All farm can be infected even before the recovered
        Mod_status_df$status [Mod_status_df$day_sinceinf > inf_duration & Mod_status_df$day_sinceinf <= waning_duration & Mod_status_df$status %in% c("I", "Q", "UI", "L","VI")] <-  "R"
        
        # if day_sinceinf > waning it become susceptible again and day since infection become 0
        Mod_status_df$status[Mod_status_df$status=="R" & Mod_status_df$day_sinceinf > waning_duration] <- "S"
        Mod_status_df$day_sinceinf[Mod_status_df$status== "S" & Mod_status_df$day_sinceinf > waning_duration] <- 0
        
        # name new column
        names(Mod_status_df)[ncol(Mod_status_df)] <- paste("status_D",ncol(Mod_status_df)-6,sep = "")
        
        # check if there are still infectious farms
        n_infectious = sum(Mod_status_df[, ncol(Mod_status_df)] %in% c("L","UI", "I", "Q", "VI"))
        
      }
      
      Mod_outbreaksim [[n]] <- Mod_status_df
    }
    #add all parameters to results
    sum_result<-sum_outbreak(dat = Mod_outbreaksim) 
    sum_result$delay_culling <-sen_kernel$delay_culling[r]
    
    
    # conclude the results and keep it in list
    result_sim[[r]] <-  sum_result
    
  }
  return(result_sim)
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~---------------------------------

## Simulation --------------
gc()
start_time <- Sys.time()
LP_sen_culling<-simulate_sen_culling(iter=500,cond=culling_sen_list,sen_kernel=sen_kernel)
saveRDS(LP_sen_culling,"LP_sen_culling.rds", compress = TRUE)
end_time <- Sys.time()
end_time - start_time