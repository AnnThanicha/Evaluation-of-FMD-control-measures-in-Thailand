# Animal movement restriction simulation ------------

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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------
## Set condition -----------------
### movement with the combination of cut off and reduce transmission ---------
Movement_comb_cond <- function(){
  
  
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
  
  
  # calculate the distance dependent transmission matrix
  Mod_kernel <- (k0/(1+((Mod_distance/r0)^alpha))) 
  diag(Mod_kernel) <- 0 # diagonal = 0
  
  
  # trade network matrix  
  Mod_component_matrix_delta <- component_matrix * delta
  diag(Mod_component_matrix_delta) <- 0 # diagonal = 0
  
  
  # parameters for IP culling
  culling_execute = FALSE # parameters to tell model whether the culling is done or not
  date_start_Culling = 7 # here we assume that culling decided first farm detection
  delay_culling_time = 1# delay culling after detection
  detectiontime_cattle = 5 # detection time after day since infection
  detectiontime_goat = 18
  detectiontime_pig = 8
  
  # parameters for emergency vaccination
  EMvaccine_execute_week1 = FALSE # set true or false to tell model if the emergency vaccination is done or not
  # here we assume that emergency vaccination decided at 6 days (detect in 5 day and 1 day do vaccine)
  date_start_EMvaccine_week1 = detectiontime_cattle+1
  
  # week2 EM vaccination 
  EMvaccine_execute_week2 = FALSE
  
  date_start_EMvaccine_week2 = date_start_EMvaccine_week1+14
  
  
  EM_outsidein = TRUE
  rate_EMvaccination =50 # we assume that emergency vaccination capacity = 30 farms per day
  EMvaccine_induce_cattle = 10 # emergency vaccine induced immunity reached protection after 4 days
  EMvaccine_induce_goat = 6
  EMvaccine_induce_pig = 13
  EMvaccine_protect_duration_pig= 210   # the duration of emergency vaccine induced immune protection
  EMvaccine_protect_duration_goat = 180
  EMvaccine_protect_duration_cattle = 180
  EMvaccine_radius = 5 #2 km radius from vaccination farm
  EM_date = 1    # use to define the list of farm that will get vaccination
  
  EM_ID_per_day1 <-c()
  EM_ID_per_day2 <-c()
  
  # parameters for movement restriction
  # When movement restriction happened, the trade is stopping and kernel is lower.
  movementR_execute = TRUE
  # date that we start movement restriction after detection the cases
  
  # we assume that alpha increase 2 time, r0 decrease 1/2, delta = 0
  r0_modify = r0
  alpha_modify = alpha
  delta_modify = 0
  k0_modify = k0
  # movementR_radius_cutoff =0.5 # after 250 m the transmission stop
  # calculate the distance dependent transmission matrix
  # Mod_kernel_movementR <- (k0_modify/(1+((Mod_distance/r0_modify )^alpha_modify)))
  
  # create modify kernel matrix, 
  # if the distance more than cutoff it will reduce by n
  # Mod_kernel_movementR[Mod_distance > movementR_radius_cutoff] <- Mod_kernel_movementR[Mod_distance > movementR_radius_cutoff]*percent_kernel
  # diag( Mod_kernel_movementR) <- 0 # diagonal = 0
  
  # trade network matrix  
  Mod_component_matrix_delta_movementR <- component_matrix * delta_modify
  diag(Mod_component_matrix_delta_movementR) <- 0 # diagonal = 0
  
  
  
  list( index_farm = index_farm , k0= k0, r0 = r0, alpha = alpha, delta = delta, 
        shape_inf = shape_inf, rate_inf = rate_inf, Mod_kernel = Mod_kernel,
        Mod_component_matrix_delta = Mod_component_matrix_delta, 
        
        EMvaccine_execute_week1 = EMvaccine_execute_week1, date_start_EMvaccine_week1 =date_start_EMvaccine_week1,
        
        EMvaccine_execute_week2 = EMvaccine_execute_week2,
        
        date_start_EMvaccine_week2 = date_start_EMvaccine_week2,
        EM_outsidein = EM_outsidein,
        rate_EMvaccination = rate_EMvaccination, 
        EMvaccine_induce_cattle = EMvaccine_induce_cattle,
        EMvaccine_induce_goat = EMvaccine_induce_goat,
        EMvaccine_induce_pig = EMvaccine_induce_pig,
        EMvaccine_protect_duration_pig= EMvaccine_protect_duration_pig,
        EMvaccine_protect_duration_goat = EMvaccine_protect_duration_goat,
        EMvaccine_protect_duration_cattle = EMvaccine_protect_duration_cattle,
        EMvaccine_radius = EMvaccine_radius,
        EM_ID_per_day1 = EM_ID_per_day1,
        EM_ID_per_day2 = EM_ID_per_day2,
        
        culling_execute = culling_execute, date_start_Culling = date_start_Culling, 
        detectiontime_cattle =  detectiontime_cattle, detectiontime_goat = detectiontime_goat,
        detectiontime_pig = detectiontime_pig, delay_culling_time =delay_culling_time,
        
        movementR_execute = movementR_execute,  
        r0_modify =  r0_modify, alpha_modify = alpha_modify ,delta_modify = delta_modify,k0_modify=k0_modify,
        Mod_component_matrix_delta_movementR=Mod_component_matrix_delta_movementR
  )
}

Movement_cond <- Movement_comb_cond()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----------------------------------

## Function to sum the outbreak simulation -------------
sum_outbreak <- function (dat) {
  # Find outbreak duration which is the number of column - 4
  outbreak_duration <- sapply(dat , ncol)-5
  n_susceptible <- sapply(dat , function(x) sum(x[,ncol(x)] == "S" ))
  n_latent <- sapply(dat , function(x) sum(x[,ncol(x)] == "L" ))
  n_infectious <- sapply(dat , function(x) sum(x[,ncol(x)] == "I" ))
  n_recover <- sapply(dat , function(x) sum(x[,ncol(x)] == "R" ))
  n_Rvaccinated <- sapply(dat , function(x) sum(x[,ncol(x)] == "RV" ))
  EMvaccinated <- sapply(dat , function(x) sum(x$day_sinceEV >0 ))
  Quarantine <- sapply(dat , function(x) sum(x$day_sinceQ >0 ))
  n_EMvaccinated <- sapply(dat , function(x) sum(x[,ncol(x)] == "EMV" ))
  n_culling <- sapply(dat , function(x) sum(x[,ncol(x)] == "C" ))
  n_UI <- sapply(dat , function(x) sum(x[,ncol(x)] == "UI" ))
  outbreak <- data.frame (outbreak_duration,EMvaccinated, n_susceptible, n_latent, n_infectious,
                          n_recover,  n_Rvaccinated, n_EMvaccinated,n_culling,Quarantine,n_UI)
  return(outbreak)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------
## Set combination parameters for movementR -------
# set radius
radius <- seq(from=0, to= 1, by = 0.1)
percent_kernel <- seq(from=0, to= 0.9, by = 0.1)
comb_par <- expand.grid(radius,percent_kernel)
# change name 
colnames(comb_par)<- c("radius","percent_kernel")
comb_par$date <- 1


## Set combination parameters for baseline movementR ----------
radius <- seq(from=0, to= 1, by = 0.1)
percent_kernel <- 1
comb_par_baseline <- expand.grid(radius,percent_kernel)
# change name 
colnames(comb_par_baseline)<- c("radius","percent_kernel")
comb_par_baseline$date <- 1
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------

## Simulate outbreak for isolation----------
# before using this code don't forget to load the requirement condition files
#load distance and trade matrix

simulate_outbreak_movementR <- function (iter,dat_kernel, cond,comb_par) {
  
  # create list to keep the results
  result_sim <- list()
  # for loop parameters combination
  for (r in 1:nrow(comb_par)){
    progress(r, max.value = nrow(comb_par),progress.bar = FALSE)# to track the for loop progress
    # set condition here
    list2env(cond, envir = .GlobalEnv)
    
    # modify parameters
    date_start_movementR = 7 + comb_par$date[r] 
    movementR_radius_cutoff <- comb_par$radius[r]
    Mod_kernel_movementR<- Mod_kernel
    Mod_kernel_movementR[Mod_distance > movementR_radius_cutoff] <- Mod_kernel_movementR[Mod_distance > movementR_radius_cutoff]*comb_par$percent_kernel[r]
    diag( Mod_kernel_movementR) <- 0 # diagonal = 0
    
    # create list to keep results
    Mod_outbreaksim <- list()
    set.seed(1234)
    # for loop for iterated simulation
    for (n in 1:iter) {
      
      # set starter conditions
      Mod_status_df <- data.frame(ID = c(1:nrow(dat_kernel)), type = dat_kernel$type2, day_sinceinf = 0,
                                  day_sinceEV = 0,day_sinceQ = 0, status_D1 = "S")
      
      
      #infduration 
      inf_duration <- round(rgamma(nrow(dat_kernel),shape= shape_inf, rate = rate_inf))
      latent = 3 # latent period 3 days
      
      inf_duration <- inf_duration+latent # plus latent period
      
      # create immunity waning from natural infection
      #1-1/1.98 # if the herd immunity below 49.5% the farm will become susceptible
      #' estimate culling rat 20% per year
      #' Then 2.5 year the farm will become susceptible
      waning_duration = round(inf_duration + (365*2.5))
      
      
      # set index cases
      index_farm =  index_farm
      Mod_status_df$status_D1 [index_farm] <- "L" # status of index case is I
      # day since infection is equal to latent in index case as it passed latent period
      Mod_status_df$day_sinceinf [index_farm] <- 1
      
      
      # set start number of infectious farm
      n_infectious = sum(Mod_status_df$status_D1 %in% c("L","UI", "I", "Q"))
      
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
        
        ###  Movement restriction  ##########
        if(movementR_execute==TRUE & date_start_movementR <= ncol(Mod_status_df)-5){ 
          # distrans for movement restriction
          dist_trans <- colSums(Mod_kernel_movementR*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
          
          # trade network transmission for movement restriction
          trade_trans <- colSums(Mod_component_matrix_delta_movementR*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
        } 
        
        # calculate p inf
        p_inf <- 1 - exp(-1*(dist_trans+trade_trans))
        # random 0 and 1 with prob. of infection
        inf <- rbinom(length(p_inf), size = 1, prob=p_inf)
        
        # create new vector to store status next day
        Mod_status_df$status <- Mod_status_df[,ncol(Mod_status_df)]  
        
        # if farm got infected
        Mod_status_df$status[which(inf>0)] = "L"
        
        ##### update status ######
        #Mod_status_df$status <- ifelse(quarantine_execute==TRUE & Mod_status_df$day_sinceinf >= date_start_quarantine , "Q", Mod_status_df$status)
        
        # plus 1 day for farm that is not susceptible
        Mod_status_df$day_sinceinf <- ifelse( Mod_status_df$status %in% c("L","UI", "I", "Q", "R"), Mod_status_df$day_sinceinf+1, Mod_status_df$day_sinceinf)
        
        
        # if day_sinceinf < 3 farm become farm is latent
        Mod_status_df$status <- ifelse(Mod_status_df$day_sinceinf == latent & Mod_status_df$day_sinceinf > 0 & Mod_status_df$status=="S", "L", Mod_status_df$status)
        
        # undetection time for index case, it is longer than later farms
        Mod_status_df$status <- ifelse(Mod_status_df$ID%in% (index_farm) & Mod_status_df$day_sinceinf >= latent & Mod_status_df$day_sinceinf < 7 & Mod_status_df$status =="L", "UI", Mod_status_df$status)
        # if day_sinceinf > 3 & day_sinceinf <= detection time, farm become undetected infectious, and not latent anymore
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))& Mod_status_df$day_sinceinf >= latent & Mod_status_df$day_sinceinf < detectiontime_cattle & Mod_status_df$type == "cattle" & Mod_status_df$status =="L", "UI", Mod_status_df$status)
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))& Mod_status_df$day_sinceinf >= latent & Mod_status_df$day_sinceinf < detectiontime_goat & Mod_status_df$type == "goat"& Mod_status_df$status =="L", "UI", Mod_status_df$status)
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))& Mod_status_df$day_sinceinf >= latent & Mod_status_df$day_sinceinf < detectiontime_pig & Mod_status_df$type == "pig"& Mod_status_df$status =="L", "UI", Mod_status_df$status)
        
        
        # detection time for index case
        Mod_status_df$status <- ifelse(Mod_status_df$ID%in% (index_farm) & Mod_status_df$day_sinceinf== 7& Mod_status_df$day_sinceinf < inf_duration & Mod_status_df$status =="UI", "I", Mod_status_df$status)
        
        # if day_sinceinf > detection time, farm become detected infection
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))&Mod_status_df$day_sinceinf== detectiontime_cattle & Mod_status_df$day_sinceinf < inf_duration & Mod_status_df$type == "cattle"& Mod_status_df$status =="UI", "I", Mod_status_df$status)
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))&Mod_status_df$day_sinceinf== detectiontime_goat & Mod_status_df$day_sinceinf < inf_duration & Mod_status_df$type == "goat"& Mod_status_df$status =="UI", "I", Mod_status_df$status)
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))&Mod_status_df$day_sinceinf== detectiontime_pig & Mod_status_df$day_sinceinf < inf_duration & Mod_status_df$type == "pig" & Mod_status_df$status =="UI", "I", Mod_status_df$status)
        
        # if day_sinceinf > inf_duration, it is recovered. 
        # All farm can be infected even before the recovered
        Mod_status_df$status <- ifelse(Mod_status_df$day_sinceinf >= inf_duration & Mod_status_df$day_sinceinf < waning_duration & Mod_status_df$status %in% c("I", "Q", "UI", "L") , "R", Mod_status_df$status)
        
        # if day_sinceinf > waning it become susceptible again and day since infection become 0
        Mod_status_df$status <- ifelse(Mod_status_df$status=="R" & Mod_status_df$day_sinceinf >= waning_duration, "S", Mod_status_df$status)
        Mod_status_df$day_sinceinf <- ifelse(Mod_status_df$status== "S" & Mod_status_df$day_sinceinf >= waning_duration, 0, Mod_status_df$day_sinceinf)
        
        # name new column
        names(Mod_status_df)[ncol(Mod_status_df)] <- paste("status_D",ncol(Mod_status_df)-5,sep = "")
        # remove parameter from this iteration before next one
        rm(dist_trans, trade_trans, p_inf, inf) 
        # check if there are still infectious farms
        n_infectious = sum(Mod_status_df[, ncol(Mod_status_df)] %in% c("L","UI", "I", "Q"))
        
      }
      
      Mod_outbreaksim [[n]] <- Mod_status_df
    }
    # conclude the results and keep it in list
    result_sim[[r]] <- sum_outbreak(dat = Mod_outbreaksim) 
    # write RDS every iterations for back up
    saveRDS(result_sim,"result_sim.rds")
    
  }
  return(result_sim)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------------------------------------------

### simulate the movement restriction ---------
gc()
set.seed(1234)
start_time <- Sys.time()
#don't forget to set date start movement in combpar
LP_sim_movement <- simulate_outbreak_movementR (iter=500,dat_kernel=dat_kernel,cond=Movement_cond,comb_par=comb_par)
# saveRDS(LP_sim_movement, "LP_sim_movement1d500iter.rds")
end_time <- Sys.time()
end_time - start_time