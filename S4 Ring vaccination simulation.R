# Ring vaccination simulation ------------

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
### EM combine -------
EM_comb_cond <- function(){
  
  
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
  EMvaccine_execute_week1 = TRUE # set true or false to tell model if the emergency vaccination is done or not
  # here we assume that emergency vaccination decided at 6 days (detect in 5 day and 1 day do vaccine)
  date_start_EMvaccine_week1 = 7+1
  
  # week2 EM vaccination 
  EMvaccine_execute_week2 = TRUE
  
  date_start_EMvaccine_week2 = date_start_EMvaccine_week1+7
  
  
  EM_outsidein = TRUE
  rate_EMvaccination =40 # we assume that emergency vaccination capacity = 30 farms per day
  EMvaccine_protect_duration_pig= 180   # the duration of emergency vaccine induced immune protection
  EMvaccine_protect_duration_goat = 180
  EMvaccine_protect_duration_cattle = 180
  EM_date = 1    # use to define the list of farm that will get vaccination
  
  EM_ID_per_day1 <-c()
  EM_ID_per_day2 <-c()
  
  # parameters for movement restriction
  # When movement restriction happened, the trade is stopping and kernel is lower.
  movementR_execute = FALSE
  # date that we start movement restriction after detection the cases
  date_start_movementR = 7 + 1
  # we assume that alpha increase 2 time, r0 decrease 1/2, delta = 0
  r0_modify = r0
  alpha_modify = alpha
  delta_modify = 0
  k0_modify = k0
  movementR_radius_cutoff =1 # after 250 m the transmission stop
  # calculate the distance dependent transmission matrix
  Mod_kernel_movementR <- (k0_modify/(1+((Mod_distance/r0_modify )^alpha_modify)))
  
  # kernel from data above cut off radius = 0
  Mod_kernel_movementR [Mod_distance > movementR_radius_cutoff] <- 0
  diag( Mod_kernel_movementR) <- 0 # diagonal = 0
  
  # trade network matrix  
  Mod_component_matrix_delta_movementR <- component_matrix * delta_modify
  diag(Mod_component_matrix_delta_movementR) <- 0 # diagonal = 0
  
  
  # parameters for quarantine
  quarantine_execute = FALSE
  # day start quarantine after detection
  
  
  list( index_farm = index_farm , k0= k0, r0 = r0, alpha = alpha, delta = delta, 
        shape_inf = shape_inf, rate_inf = rate_inf, Mod_kernel = Mod_kernel,
        Mod_component_matrix_delta = Mod_component_matrix_delta, 
        
        EMvaccine_execute_week1 = EMvaccine_execute_week1, date_start_EMvaccine_week1 =date_start_EMvaccine_week1,
        
        EMvaccine_execute_week2 = EMvaccine_execute_week2,
        
        date_start_EMvaccine_week2 = date_start_EMvaccine_week2,
        EM_outsidein = EM_outsidein,
        rate_EMvaccination = rate_EMvaccination, 
        EMvaccine_protect_duration_pig= EMvaccine_protect_duration_pig,
        EMvaccine_protect_duration_goat = EMvaccine_protect_duration_goat,
        EMvaccine_protect_duration_cattle = EMvaccine_protect_duration_cattle,
        EM_ID_per_day1 = EM_ID_per_day1,
        EM_ID_per_day2 = EM_ID_per_day2,
        
        culling_execute = culling_execute, date_start_Culling = date_start_Culling, 
        detectiontime_cattle =  detectiontime_cattle, detectiontime_goat = detectiontime_goat,
        detectiontime_pig = detectiontime_pig, delay_culling_time =delay_culling_time,
        movementR_execute = movementR_execute, date_start_movementR = date_start_movementR,  
        r0_modify =  r0_modify, alpha_modify = alpha_modify ,delta_modify = delta_modify,k0_modify=k0_modify,
        movementR_radius_cutoff = movementR_radius_cutoff,Mod_kernel_movementR=Mod_kernel_movementR,
        Mod_component_matrix_delta_movementR=Mod_component_matrix_delta_movementR,
        
        quarantine_execute = quarantine_execute
        
  )
}

EM_cond <- EM_comb_cond()

### no EM -------
EM_no_comb_cond <- function(){
  
  
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
  date_start_EMvaccine_week1 = 7+1
  
  # week2 EM vaccination 
  EMvaccine_execute_week2 = FALSE
  
  date_start_EMvaccine_week2 = date_start_EMvaccine_week1+7
  
  
  EM_outsidein = TRUE
  rate_EMvaccination =40 # we assume that emergency vaccination capacity = 30 farms per day
  EMvaccine_protect_duration_pig= 180   # the duration of emergency vaccine induced immune protection
  EMvaccine_protect_duration_goat = 180
  EMvaccine_protect_duration_cattle = 180
  EM_date = 1    # use to define the list of farm that will get vaccination
  
  EM_ID_per_day1 <-c()
  EM_ID_per_day2 <-c()
  
  # parameters for movement restriction
  # When movement restriction happened, the trade is stopping and kernel is lower.
  movementR_execute = FALSE
  # date that we start movement restriction after detection the cases
  date_start_movementR = 7 + 1
  # we assume that alpha increase 2 time, r0 decrease 1/2, delta = 0
  r0_modify = r0
  alpha_modify = alpha
  delta_modify = 0
  k0_modify = k0
  movementR_radius_cutoff =1 # after 250 m the transmission stop
  # calculate the distance dependent transmission matrix
  Mod_kernel_movementR <- (k0_modify/(1+((Mod_distance/r0_modify )^alpha_modify)))
  
  # kernel from data above cut off radius = 0
  Mod_kernel_movementR [Mod_distance > movementR_radius_cutoff] <- 0
  diag( Mod_kernel_movementR) <- 0 # diagonal = 0
  
  # trade network matrix  
  Mod_component_matrix_delta_movementR <- component_matrix * delta_modify
  diag(Mod_component_matrix_delta_movementR) <- 0 # diagonal = 0
  
  
  # parameters for quarantine
  quarantine_execute = FALSE
  # day start quarantine after detection
  
  
  list( index_farm = index_farm , k0= k0, r0 = r0, alpha = alpha, delta = delta, 
        shape_inf = shape_inf, rate_inf = rate_inf, Mod_kernel = Mod_kernel,
        Mod_component_matrix_delta = Mod_component_matrix_delta, 
        
        EMvaccine_execute_week1 = EMvaccine_execute_week1, date_start_EMvaccine_week1 =date_start_EMvaccine_week1,
        
        EMvaccine_execute_week2 = EMvaccine_execute_week2,
        
        date_start_EMvaccine_week2 = date_start_EMvaccine_week2,
        EM_outsidein = EM_outsidein,
        rate_EMvaccination = rate_EMvaccination, 
        EMvaccine_protect_duration_pig= EMvaccine_protect_duration_pig,
        EMvaccine_protect_duration_goat = EMvaccine_protect_duration_goat,
        EMvaccine_protect_duration_cattle = EMvaccine_protect_duration_cattle,
        EM_ID_per_day1 = EM_ID_per_day1,
        EM_ID_per_day2 = EM_ID_per_day2,
        
        culling_execute = culling_execute, date_start_Culling = date_start_Culling, 
        detectiontime_cattle =  detectiontime_cattle, detectiontime_goat = detectiontime_goat,
        detectiontime_pig = detectiontime_pig, delay_culling_time =delay_culling_time,
        movementR_execute = movementR_execute, date_start_movementR = date_start_movementR,  
        r0_modify =  r0_modify, alpha_modify = alpha_modify ,delta_modify = delta_modify,k0_modify=k0_modify,
        movementR_radius_cutoff = movementR_radius_cutoff,Mod_kernel_movementR=Mod_kernel_movementR,
        Mod_component_matrix_delta_movementR=Mod_component_matrix_delta_movementR,
        
        quarantine_execute = quarantine_execute
        
  )
}

EM_no_cond <- EM_no_comb_cond ()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-----------------------------------
## Function to sum the outbreak simulation -------------
sum_outbreak <- function (dat) {
  # Find outbreak duration which is the number of column - 4
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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-----------------------------------
## Set combination parameters for EM -------
# set radius
EMvaccine_induce <- seq(from = 7, to = 28, by = 7)
percent_protect <- seq(from=0.5, to= 1, by = 0.1)
radius <- 10
comb_par <- expand.grid(EMvaccine_induce ,percent_protect, radius)
# change name 
colnames(comb_par)<- c("EMvaccine_induce","percent_protect", "radius")

## Set combination parameters for EM below 50% protection -------
EMvaccine_induce <- seq(from = 7, to = 28, by = 7)
percent_protect <- seq(from=0, to= 0.4, by = 0.1)# protect 0 is equal so no need for combination
radius <- 10
comb_par_below0.5 <- expand.grid(EMvaccine_induce ,percent_protect, radius)
# change name 
colnames(comb_par_below0.5)<- c("EMvaccine_induce","percent_protect", "radius")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------
simulate_outbreak_EM <- function (iter,dat_kernel, cond,comb_par) {
  # create list to keep the results
  result_sim <- list()
  # for loop parameters combination
  for (r in 1:nrow(comb_par)){
    progress(r, max.value = nrow(comb_par),progress.bar = FALSE)# to track the for loop progress
    # set condition here
    list2env(cond, envir = .GlobalEnv)
    
    # set kernel for VI and S
    #row is ID of infectious farms, column is ID of susceptible farms
    Mod_kernel_VI_S <- Mod_kernel*(1-comb_par$percent_protect[r])
    diag( Mod_kernel_VI_S) <- 0
    Mod_component_VI_S <- component_matrix * delta *(1-comb_par$percent_protect[r])
    diag(Mod_component_VI_S) <- 0 
    
    # set kernel for I and VS
    #row is ID of infectious farms, column is ID of susceptible farms
    Mod_kernel_I_VS <- Mod_kernel*(1-comb_par$percent_protect[r])
    diag( Mod_kernel_I_VS) <- 0
    Mod_component_I_VS<- component_matrix * delta *(1-comb_par$percent_protect[r])
    diag(Mod_component_I_VS) <- 0 
    
    # set kernel for VI and VS
    #row is ID of infectious farms, column is ID of susceptible farms
    Mod_kernel_VI_VS <- Mod_kernel*(1-comb_par$percent_protect[r])*(1-comb_par$percent_protect[r])
    diag( Mod_kernel_VI_VS) <- 0
    Mod_component_VI_VS<- component_matrix * delta *(1-comb_par$percent_protect[r])*(1-comb_par$percent_protect[r])
    diag(Mod_component_VI_VS) <- 0 
    
    EMvaccine_radius <- comb_par$radius[r]
    
    EMvaccine_induce <- comb_par$EMvaccine_induce[r]
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
      n_infectious = sum(Mod_status_df$status_D1 %in% c("L","UI", "I", "Q","VI"))
      
      # simulated in one iteration until no infectious farms
      while(n_infectious > 0){
        
        
        
        #i = ncol(Mod_status_df)-5 # i = equal to column number before status_day
        
        # calculate distance dependent transmission
        status_infectious <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("I","UI"))
        status_susceptible <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("S"))
        status_VIinfectious <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("VI"))
        status_VSsusceptible <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] %in% c("VS"))
        
        ### normal kernel ######
        # row is ID of infectious farms, column is ID of susceptible farms
        dist_trans <- colSums(Mod_kernel*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
        
        # calculate trade network transmission
        trade_trans <- colSums(Mod_component_matrix_delta*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
        
        ###  transmission between VI and S  ##########
        dist_trans_VI_S <- colSums(Mod_kernel_VI_S*(as.matrix(status_VIinfectious)%*% as.matrix(t(status_susceptible)))) 
        trade_trans_VI_S <- colSums(Mod_component_VI_S*(as.matrix(status_VIinfectious)%*% as.matrix(t(status_susceptible)))) 
        
        ###  transmission between I and VS  ##########
        dist_trans_I_VS <- colSums(Mod_kernel_I_VS*(as.matrix(status_infectious)%*% as.matrix(t(status_VSsusceptible)))) 
        trade_trans_I_VS <- colSums(Mod_component_I_VS*(as.matrix(status_infectious)%*% as.matrix(t(status_VSsusceptible))))
        
        ###  transmission between VI and VS  ##########
        dist_trans_VI_VS <- colSums(Mod_kernel_VI_VS*(as.matrix(status_VIinfectious)%*% as.matrix(t(status_VSsusceptible)))) 
        trade_trans_VI_VS <- colSums(Mod_component_VI_VS*(as.matrix(status_VIinfectious)%*% as.matrix(t(status_VSsusceptible)))) 
        
        
        # calculate p inf
        p_inf <- 1 - exp(-1*(dist_trans+trade_trans+dist_trans_VI_S+trade_trans_VI_S+
                               dist_trans_I_VS+trade_trans_I_VS+dist_trans_VI_VS +trade_trans_VI_VS))
        # random 0 and 1 with prob. of infection
        inf <- rbinom(length(p_inf), size = 1, prob=p_inf)
        
        # create new vector to store status next day
        Mod_status_df$status <- Mod_status_df[,ncol(Mod_status_df)]  
        
        # if farm got infected
        Mod_status_df$status[which(inf>0)] = "L"
        
        ##### update status ######
        
        # plus 1 day for farm that is not susceptible
        Mod_status_df$day_sinceinf <- ifelse( Mod_status_df$status %in% c("L","UI", "I", "Q", "R", "VI"), Mod_status_df$day_sinceinf+1, Mod_status_df$day_sinceinf)
        Mod_status_df$day_sinceQ <- ifelse( Mod_status_df$status %in% c("Q"), Mod_status_df$day_sinceQ+1, Mod_status_df$day_sinceQ)
        # day since EV +1 in EM vaccinated farm
        Mod_status_df$day_sinceEV <-  ifelse( Mod_status_df$day_sinceEV > 0, Mod_status_df$day_sinceEV+1, Mod_status_df$day_sinceEV)
        
        # if the day since inf >= date start quarantine 
        # only the detected farm will go to quarantine
        # Mod_status_df$status <- ifelse(quarantine_execute==TRUE & Mod_status_df$status %in% c("I") & Mod_status_df$day_sinceinf >= date_start_quarantine_cattle & Mod_status_df$type == "cattle", "Q", Mod_status_df$status)
        #Mod_status_df$status <- ifelse(quarantine_execute==TRUE & Mod_status_df$status %in% c("I") & Mod_status_df$day_sinceinf >= date_start_quarantine_goat & Mod_status_df$type == "goat", "Q", Mod_status_df$status)
        # Mod_status_df$status <- ifelse(quarantine_execute==TRUE & Mod_status_df$status %in% c("I") & Mod_status_df$day_sinceinf >= date_start_quarantine_pig & Mod_status_df$type == "pig", "Q", Mod_status_df$status)
        
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
        
        
        # if day_sinceinf > detection time, and vaccinated farms become VI
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))&Mod_status_df$day_sinceinf== detectiontime_cattle & Mod_status_df$day_sinceinf < inf_duration &  Mod_status_df$day_sinceEV >= EMvaccine_induce & Mod_status_df$type == "cattle"& Mod_status_df$status =="UI", "VI", Mod_status_df$status)
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))&Mod_status_df$day_sinceinf== detectiontime_goat & Mod_status_df$day_sinceinf < inf_duration & Mod_status_df$day_sinceEV >= EMvaccine_induce &Mod_status_df$type == "goat"& Mod_status_df$status =="UI", "VI", Mod_status_df$status)
        Mod_status_df$status <- ifelse(!(Mod_status_df$ID%in% (index_farm))&Mod_status_df$day_sinceinf== detectiontime_pig & Mod_status_df$day_sinceinf < inf_duration & Mod_status_df$day_sinceEV >= EMvaccine_induce & Mod_status_df$type == "pig" & Mod_status_df$status =="UI", "VI", Mod_status_df$status)
        
        # if day_sinceinf > inf_duration, it is recovered. 
        # All farm can be infected even before the recovered
        Mod_status_df$status <- ifelse(Mod_status_df$day_sinceinf >= inf_duration & Mod_status_df$day_sinceinf < waning_duration & Mod_status_df$status %in% c("I", "Q", "UI", "L", "VI'") , "R", Mod_status_df$status)
        
        # if day_sinceinf > waning it become susceptible again and day since infection become 0
        Mod_status_df$status <- ifelse(Mod_status_df$status=="R" & Mod_status_df$day_sinceinf >= waning_duration, "S", Mod_status_df$status)
        Mod_status_df$day_sinceinf <- ifelse(Mod_status_df$status== "S" & Mod_status_df$day_sinceinf >= waning_duration, 0, Mod_status_df$day_sinceinf)
        
        
        # if the day_sinceEV >1, it means they already got EM vaccination in that case
        # the day since EM vaccine than EMvaccince_induce they change status if farm is S to VS
        # if the farm got infected or latent before emergency vaccination, it will not be protected
        # the yday that got protection depended on the species
        Mod_status_df$status <-  ifelse( Mod_status_df$day_sinceEV == EMvaccine_induce  & Mod_status_df$status %in% c("S"), "VS", Mod_status_df$status)
        
        
        # if day since EM vaccine > EMvaccine_protect_duration, 
        # the status is back to S, and day_sinceEV become -1 so it is stopped increase in next loop
        Mod_status_df$status <-  ifelse( Mod_status_df$day_sinceEV > EMvaccine_protect_duration_cattle &  Mod_status_df$status %in% c("VS")& Mod_status_df$type == "cattle", "S", Mod_status_df$status)
        Mod_status_df$day_sinceEV <-  ifelse( Mod_status_df$day_sinceEV > EMvaccine_protect_duration_cattle &  Mod_status_df$status %in% c("S")& Mod_status_df$type == "cattle", -1, Mod_status_df$day_sinceEV)
        
        Mod_status_df$status <-  ifelse( Mod_status_df$day_sinceEV > EMvaccine_protect_duration_goat &  Mod_status_df$status %in% c("VS")& Mod_status_df$type == "goat", "S", Mod_status_df$status)
        Mod_status_df$day_sinceEV <-  ifelse( Mod_status_df$day_sinceEV > EMvaccine_protect_duration_goat &  Mod_status_df$status %in% c("S")& Mod_status_df$type == "goat", -1, Mod_status_df$day_sinceEV)
        
        Mod_status_df$status <-  ifelse( Mod_status_df$day_sinceEV > EMvaccine_protect_duration_pig &  Mod_status_df$status %in% c("VS")& Mod_status_df$type == "pig", "S", Mod_status_df$status)
        Mod_status_df$day_sinceEV <-  ifelse( Mod_status_df$day_sinceEV > EMvaccine_protect_duration_pig &  Mod_status_df$status %in% c("S")& Mod_status_df$type == "pig", -1, Mod_status_df$day_sinceEV)
        
        ### List EM week1 #######
        if( EMvaccine_execute_week1 ==TRUE & ncol(Mod_status_df)-5 == date_start_EMvaccine_week1){  # set the farms for EM vaccine
          # the farm that close to the infected farms will be vaccinated first 
          # prepare radius for EM vaccination
          # get status infectious from index case
          status_infectious_EM <- as.numeric(Mod_status_df[, date_start_EMvaccine_week1+5] %in% "I")
          # because we don't know if farm being latent, undetected infected or susceptible
          # the vaccination will be done on all these three status
          status_susceptible_EM <- as.numeric(Mod_status_df[, date_start_EMvaccine_week1+5] %in% c("S","L","UI"))
          
          vaccination_radius <- (Mod_distance*(as.matrix(status_infectious_EM)%*% as.matrix(t(status_susceptible_EM))))
          # remove zero to NA
          vaccination_radius[vaccination_radius == 0 ] <- NA
          # remove farm that is not in radius of infected farm
          vaccination_radius[vaccination_radius > EMvaccine_radius ] <- NA
          
          EM_distance <- as.numeric(apply(vaccination_radius,2, function(x)min(x, na.rm = TRUE)))
          EM_ID <- data.frame(ID = 1:nrow(dat_kernel), EM_distance=EM_distance )  
          # order the farm distance by inside-out or outside-in
          if(EM_outsidein == TRUE) { EM_ID <- EM_ID[order(-EM_distance),]} else {EM_ID <- EM_ID[order(EM_distance),]}
          # select the farm that is not infinite
          EM_ID <- EM_ID$ID[which(!is.infinite(EM_ID$EM_distance))] 
          
          # The farm ID that will get vaccinated each day, it saved in a list which we can use it later
          EM_ID_per_day1 <- split(EM_ID, rep(1:ceiling(length(EM_ID)/rate_EMvaccination), each = rate_EMvaccination, length.out=length(EM_ID)))
          EM_date <-1
        }# if the EM vaccination list is not created then the EM_ID_per_day = 0
        
        
        
        ### vaccination week1 ########
        # start emergency vaccination at the set date and EM vaccine execute = TRUE 
        # and stop when run out of Farm ID to vaccinated
        
        if(EMvaccine_execute_week1 ==TRUE & ncol(Mod_status_df)-5 >=date_start_EMvaccine_week1&
           length(EM_ID_per_day1)> 0 & ncol(Mod_status_df)-5 <=date_start_EMvaccine_week2){
          # EM_ID_per_day is the list that keep ID farm that do emergency vaccination per day
          
          # farm that got vaccination 
          Mod_status_df$day_sinceEV[EM_ID_per_day1[[EM_date]]] <- 1
          EM_date <- EM_date+1 # next set of list will be get vaccinated
          # if the we run out of vaccine list or the time for week2 vaccination
          if(EM_date>length(EM_ID_per_day1)){EM_ID_per_day1 <- c()}
          
        }
        
        ### List EM week2 ############
        # This is starting when EM vaccine week2 == TRUE and  the time is more than date_start_EMvaccine_week2
        # We will update the new set of EM vaccine farms
        if(EMvaccine_execute_week2 ==TRUE & ncol(Mod_status_df)-5 == date_start_EMvaccine_week2) {
          status_infectious_EM <- as.numeric(Mod_status_df[, date_start_EMvaccine_week2+5] %in% "I")
          # because we don't know if farm being latent, undetected infected or susceptible
          # the vaccination will be done on all these three status
          #if the farms already vaccinated in the first round, it doesn't get in second round
          status_susceptible_EM <- as.numeric(Mod_status_df[, date_start_EMvaccine_week2+5] %in% c("S","L","UI") & Mod_status_df$day_sinceEV==0)
          
          
          vaccination_radius <- (Mod_distance*(as.matrix(status_infectious_EM)%*% as.matrix(t(status_susceptible_EM))))
          # remove zero to NA
          vaccination_radius[vaccination_radius == 0 ] <- NA
          
          # remove farm that is not in radius of infected farm
          vaccination_radius[vaccination_radius > EMvaccine_radius ] <- NA
          
          EM_distance <- as.numeric(apply(vaccination_radius,2, function(x)min(x, na.rm = TRUE)))
          EM_ID <- data.frame(ID = 1:nrow(dat_kernel), EM_distance=EM_distance )  
          # order the farm distance by inside-out or outside-in
          if(EM_outsidein == TRUE) { EM_ID <- EM_ID[order(-EM_distance),]} else {EM_ID <- EM_ID[order(EM_distance),]}
          
          # select the farm that is not infinite
          EM_ID <- EM_ID$ID[which(!is.infinite(EM_ID$EM_distance))] 
          
          # The farm ID that will get vaccinated each day, it saved in a list which we can use it later
          EM_ID_per_day2 <- split(EM_ID, rep(1:ceiling(length(EM_ID)/rate_EMvaccination), each = rate_EMvaccination, length.out=length(EM_ID)))
          EM_date = 1    # use to define the list of farm that will get vaccination
        }
        
        ### vaccination week2 ########
        # if all farms in radius were vaccinated in the first week this part is skipped
        if( EMvaccine_execute_week2 ==TRUE & ncol(Mod_status_df)-5>=date_start_EMvaccine_week2  
            &length(EM_ID_per_day2)>0  ){
          # EM_ID_per_day is the list that keep ID farm that do emergency vaccination per day
          
          # farm that got vaccination 
          Mod_status_df$day_sinceEV[EM_ID_per_day2[[EM_date]]] <- 1
          EM_date <- EM_date+1 # next set of list will be get vaccinated
          # if the we run out of vaccine list, the vaccination stop
          if(EM_date>length(EM_ID_per_day2)){EM_ID_per_day2 <- c()}
        }
        
        
        
        
        # name new column
        names(Mod_status_df)[ncol(Mod_status_df)] <- paste("status_D",ncol(Mod_status_df)-5,sep = "")
        # remove parameter from this iteration before next one
        rm(dist_trans, trade_trans, p_inf, inf ) 
        # check if there are still infectious farms
        n_infectious = sum(Mod_status_df[, ncol(Mod_status_df)] %in% c("L","UI", "I", "Q", "VI"))
        
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

### simulate EM -----------
gc()
set.seed(1234)
start_time <- Sys.time()
LP_simulate_EM10km <- simulate_outbreak_EM (iter=500,dat_kernel=dat_kernel,cond=EM_cond,comb_par=comb_par)
end_time <- Sys.time()
end_time - start_time
#saveRDS(LP_simulate_EM10km,"LP_simulate_EM10