
# Isolation simulation ------------

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

## Set condition -----------------

Quarantine_comb_cond <- function(){
  
  
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
  detection_time_index = 7 # detection index cases
  detectiontime_cattle = 5 # for dairy cattle
  
  
  # parameters for quarantine
  quarantine_execute = TRUE
  # day start quarantine after detection
  
  
  list( index_farm = index_farm , k0= k0, r0 = r0, alpha = alpha, delta = delta, 
        shape_inf = shape_inf, rate_inf = rate_inf, Mod_kernel = Mod_kernel,
        Mod_component_matrix_delta = Mod_component_matrix_delta, 
        
        detection_time_index = detection_time_index,
        detectiontime_cattle = detectiontime_cattle,
        quarantine_execute = quarantine_execute
        
  )
}
Quarantine_cond <- Quarantine_comb_cond()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~-----------------------------------

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

## Set combination parameters isolation-------
# set radius
radius <- seq(from=0, to= 1, by = 0.1)
percent_kernel <- seq(from=0, to= 0.9, by = 0.1)#except kernel=1 which is baseline
comb_par <- expand.grid(radius,percent_kernel)
# change name 
colnames(comb_par)<- c("radius","percent_kernel")
comb_par$date <- 1

## Set combination parameters for QR ----------
radius <- seq(from=0, to= 1, by = 0.1)
percent_kernel <- 0
comb_par_baseline <- expand.grid(radius,percent_kernel)
# change name 
colnames(comb_par_baseline)<- c("radius","percent_kernel")
comb_par_baseline$date <- 1
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------
## Simulate outbreak for quarantine----------
# before using this code don't forget to load the requirement condition files
#load distance and trade matrix

simulate_outbreak_QR <- function (iter,dat_kernel, cond,comb_par) {
  
  # create list to keep the results
  result_sim <- list()
  # for loop parameters combination
  for (r in 1:nrow(comb_par)){
    progress(r, max.value = nrow(comb_par),progress.bar = FALSE)# to track the for loop progress
    # set condition here
    list2env(cond, envir = .GlobalEnv)
    
    # modify parameters
    date_start_quarantine_cattle = detectiontime_cattle + comb_par$date[r]
    
    
    movementQ_radius_cutoff <- comb_par$radius[r]
    Mod_kernel_movementQ<- Mod_kernel
    Mod_kernel_movementQ[Mod_distance > movementQ_radius_cutoff] <- Mod_kernel_movementQ[Mod_distance > movementQ_radius_cutoff]*comb_par$percent_kernel[r]
    diag( Mod_kernel_movementQ) <- 0 # diagonal = 0
    
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
        status_susceptible <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] =="S")
        
        
        ### normal kernel ######
        # row is ID of infectious farms, column is ID of susceptible farms
        dist_trans <- colSums(Mod_kernel*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
        
        # calculate trade network transmission
        trade_trans <- colSums(Mod_component_matrix_delta*(as.matrix(status_infectious)%*% as.matrix(t(status_susceptible)))) 
        
        ### quarantine ####
        # calculate the transmission of quarantine farms 
        if(quarantine_execute==TRUE){
          status_infectious_Q <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] =="Q")
          status_susceptible_Q <- as.numeric(Mod_status_df[, ncol(Mod_status_df)] =="S")
          
          dist_trans_Q <- colSums(Mod_kernel_movementQ*(as.matrix(status_infectious_Q)%*% as.matrix(t(status_susceptible_Q)))) 
          dist_trans <- dist_trans+dist_trans_Q
          # trade transmission for quarantine is zero, so no need for add this to trade_trans
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
        
        # plus 1 day for farm that is not susceptible
        Mod_status_df$day_sinceinf <- if_else( Mod_status_df$status %in% c("L","UI", "I", "Q", "R"), Mod_status_df$day_sinceinf+1, Mod_status_df$day_sinceinf)
        Mod_status_df$day_sinceQ <- if_else( Mod_status_df$status =="Q", Mod_status_df$day_sinceQ+1, Mod_status_df$day_sinceQ)
        
        
        # if the day since inf >= date start quarantine 
        # only the detected farm will go to quarantine
        if(quarantine_execute==TRUE){
          Mod_status_df$status <- if_else( Mod_status_df$status =="I" & Mod_status_df$day_sinceinf >= date_start_quarantine_cattle , "Q", Mod_status_df$status)
        }
        
        # undetection time for index case, it is longer than later farms
        Mod_status_df$status <- if_else(Mod_status_df$ID%in% (index_farm) & Mod_status_df$day_sinceinf > latent & Mod_status_df$day_sinceinf < detection_time_index & Mod_status_df$status =="L", "UI", Mod_status_df$status)
        # if day_sinceinf > 3 & day_sinceinf <= detection time, farm become undetected infectious, and not latent anymore
        Mod_status_df$status <- if_else(!(Mod_status_df$ID%in% (index_farm)) & Mod_status_df$day_sinceinf > latent & Mod_status_df$day_sinceinf < detectiontime_cattle  & Mod_status_df$status =="L", "UI", Mod_status_df$status)
        
        
        # detection time for index case
        Mod_status_df$status <- if_else(Mod_status_df$ID%in% (index_farm) & Mod_status_df$day_sinceinf== detection_time_index & Mod_status_df$day_sinceinf <= inf_duration & Mod_status_df$status =="UI", "I", Mod_status_df$status)
        # if day_sinceinf > detection time, farm become detected infection
        Mod_status_df$status <- if_else(!(Mod_status_df$ID%in% (index_farm)) & Mod_status_df$day_sinceinf== detectiontime_cattle & Mod_status_df$day_sinceinf <= inf_duration & Mod_status_df$status =="UI", "I", Mod_status_df$status)
        
        # if day_sinceinf > inf_duration, it is recovered. 
        # All farm can be infected even before the recovered
        Mod_status_df$status <- if_else(Mod_status_df$day_sinceinf > inf_duration & Mod_status_df$day_sinceinf <= waning_duration & Mod_status_df$status %in% c("I", "Q", "UI", "L") , "R", Mod_status_df$status)
        
        # if day_sinceinf > waning it become susceptible again and day since infection become 0
        Mod_status_df$status <- if_else(Mod_status_df$status=="R" & Mod_status_df$day_sinceinf > waning_duration, "S", Mod_status_df$status)
        Mod_status_df$day_sinceinf <- if_else(Mod_status_df$status== "S" & Mod_status_df$day_sinceinf > waning_duration, 0, Mod_status_df$day_sinceinf)
        
        # name new column
        names(Mod_status_df)[ncol(Mod_status_df)] <- paste("status_D",ncol(Mod_status_df)-5,sep = "")
        # remove parameter from this iteration before next one
        rm(dist_trans, trade_trans, p_inf, inf) 
        # check if there are still infectious farms
        n_infectious = sum(Mod_status_df[, ncol(Mod_status_df)] %in% c("L","UI", "I", "Q"))
        
      }
      
      Mod_outbreaksim [[n]] <- Mod_status_df
    }
    # save RDS for raw results
    saveRDS( list(comb_par[r,],Mod_outbreaksim), file = paste('LP_QR1D',comb_par$radius[r],'km',comb_par$percent_kernel[r],'kernel','.RDS', sep=""),compress = TRUE)
    # conclude the results and keep it in list
    result_sim[[r]] <- sum_outbreak(dat = Mod_outbreaksim) 
    # clear garbage every 20 iterations
    
  }
  return(result_sim)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~------------------------------------------------------------------------------

### simulate the quarantine ------------------
gc()
#set.seed(1234)

start_time <- Sys.time()
# don't forget to set the start quarantine date in comb_par
LP_sim_quarantine <- simulate_outbreak_QR(iter=500,dat_kernel=dat_kernel,cond=Quarantine_cond,comb_par=comb_par)
#saveRDS(LP_sim_quarantine,"LP_sim_quarantine.rds")
end_time <- Sys.time()
end_time - start_time
