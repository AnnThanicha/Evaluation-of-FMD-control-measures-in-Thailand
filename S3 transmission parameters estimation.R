# ## Transmission parameters estimation ----------
##~~~~~~~~~~~~~~~~-------------------------------------------------------
# required library ---------------------------------------------------

library(bbmle)
library(optimx)
library(readr)
library(readxl)
##~~~~~~~~~~~~~~~~-------------------------------------------------------
## Import data ---------------------------------------------------
# dummy variable for susceptible status of farm for each outbreak day in area 1
status_sus_1 <- readRDS("status_infectious_1.rds") 
# dummy variable for susceptible status of farm for each outbreak day in area 2
status_sus_2 <- readRDS("status_infectious_2.rds")
# dummy variable for infected status of farm for each outbreak day in area 1
status_infected_1 <- readRDS("status_infected_1.rds")
# dummy variable for infected status of farm for each outbreak day in area 2
status_infected_2 <- readRDS("status_infected_2.rds")
# dummy variable for infectious status of farm for each outbreak day in area 1
status_infectious_1 <- readRDS("status_infectious_1.rds")
# dummy variable for infectious status of farm for each outbreak day in area 2
status_infectious_2 <- readRDS("status_infectious_2.rds")
# Matrix of distance between farms in area 1
distancematrix_1  <- readRDS("distancematrix_1.rds")
# Matrix of distance between farms in area 2
distancematrix_2  <- readRDS("distancematrix_2.rds")
# Trade network matrix for area 1
component_matrix_1  <- readRDS("component_matrix_1.rds")
# Trade network matrix for area 2
component_matrix_2 <- readRDS("component_matrix_2.rds")
##~~~~~~~~~~~~~~~~-------------------------------------------------------
## Function for kernel estimation ------------------------------------
kernel_estimate_joint <- function( k0,r0,alpha,delta){
  
  ## for the first area
  kernelmatrix_1<-(k0/(1+((distancematrix_1/r0)^alpha))) 
  diag( kernelmatrix_1) <- 0
  
  component_matrix_delta_1<-component_matrix_1*delta
  diag(component_matrix_delta_1) <- 0
  
  
  #create blank matrix to store lambda
  lamda_inf_1<-matrix(NA , ncol(status_sus_1), nrow(status_sus_1)) 
  lamda_esc_1<-matrix(NA , ncol(status_sus_1), nrow(status_sus_1))
  lamda_esc_trade_1<-matrix(NA , ncol(status_sus_1), nrow(status_sus_1))
  lamda_inf_trade_1<-matrix(NA , ncol(status_sus_1), nrow(status_sus_1))
  
  for (i in 1:ncol(status_sus_1)){
    #lambda for escape
    lamda_esc_1[i,] = rowSums(kernelmatrix_1*(as.matrix(status_sus_1[,i])%*% as.matrix(t(status_infectious_1[,i])))) 
    
    
    #lambda for infected
    lamda_inf_1[i,] = rowSums(kernelmatrix_1*(as.matrix(status_infected_1[,i])%*% as.matrix(t(status_infectious_1[,i]))))
    
    #lambda for escape trade
    
    lamda_esc_trade_1[i,] = rowSums(component_matrix_delta_1*(as.matrix(status_sus_1[,i])%*% as.matrix(t(status_infectious_1[,i])))) 
    
    
    #lambda for infected trade
    lamda_inf_trade_1[i,] = rowSums(component_matrix_delta_1*(as.matrix(status_infected_1[,i])%*% as.matrix(t(status_infectious_1[,i]))))
    
  }
  
  lamda_inf_total_1<-lamda_inf_1+lamda_inf_trade_1 # sum lambda from distance independent and distance dependent
  lamda_inf_ID_1<-colSums(lamda_inf_total_1) # sum lambda infection for each farm
  log_lamda_inf_ID_1<- ifelse(lamda_inf_ID_1>0,log(1-(exp(-1*lamda_inf_ID_1))),0) # if lambda =0 , not calculate prob inf
  logPinf_1<-sum(log_lamda_inf_ID_1)
  
  
  
  lamda_esc_ID_1<-colSums(lamda_esc_1)# sum force of infection from each ID
  logPesc_1<-sum(lamda_esc_ID_1)
  
  
  
  lamda_esc_trade_ID_1<-colSums(lamda_esc_trade_1)# sum force of infection from each ID
  logPesc_trade_1<-sum(lamda_esc_trade_ID_1)
  
  ## for the second area
  kernelmatrix_2<-(k0/(1+((distancematrix_2/r0)^alpha))) 
  diag( kernelmatrix_2) <- 0
  
  component_matrix_delta_2<-component_matrix_2*delta
  diag(component_matrix_delta_2) <- 0
  
  
  #create blank matrix to store lambda
  lamda_inf_2<-matrix(NA , ncol(status_sus_2), nrow(status_sus_2)) 
  lamda_esc_2<-matrix(NA , ncol(status_sus_2), nrow(status_sus_2))
  lamda_esc_trade_2<-matrix(NA , ncol(status_sus_2), nrow(status_sus_2))
  lamda_inf_trade_2<-matrix(NA , ncol(status_sus_2), nrow(status_sus_2))
  
  for (i in 1:ncol(status_sus_2)){
    #lambda for escape
    lamda_esc_2[i,] = rowSums(kernelmatrix_2*(as.matrix(status_sus_2[,i])%*% as.matrix(t(status_infectious_2[,i])))) 
    
    
    #lambda for infected
    lamda_inf_2[i,] = rowSums(kernelmatrix_2*(as.matrix(status_infected_2[,i])%*% as.matrix(t(status_infectious_2[,i]))))
    
    #lambda for escape trade
    
    lamda_esc_trade_2[i,] = rowSums(component_matrix_delta_2*(as.matrix(status_sus_2[,i])%*% as.matrix(t(status_infectious_2[,i])))) 
    
    
    #lambda for infected trade
    lamda_inf_trade_2[i,] = rowSums(component_matrix_delta_2*(as.matrix(status_infected_2[,i])%*% as.matrix(t(status_infectious_2[,i]))))
    
  }
  
  lamda_inf_total_2<-lamda_inf_2+lamda_inf_trade_2 # sum lambda from distance independent and distance dependent
  lamda_inf_ID_2<-colSums(lamda_inf_total_2) # sum lambda infection for each farm
  log_lamda_inf_ID_2<- ifelse(lamda_inf_ID_2>0,log(1-(exp(-1*lamda_inf_ID_2))),0) # if lambda =0 , not calculate prob inf
  logPinf_2<-sum(log_lamda_inf_ID_2)
  
  
  
  lamda_esc_ID_2<-colSums(lamda_esc_2)# sum force of infection from each ID
  logPesc_2<-sum(lamda_esc_ID_2)
  
  
  
  lamda_esc_trade_ID_2<-colSums(lamda_esc_trade_2)# sum force of infection from each ID
  logPesc_trade_2<-sum(lamda_esc_trade_ID_2)
  
  return((-1)*(logPinf_1-logPesc_1-logPesc_trade_1 +logPinf_2-logPesc_2-logPesc_trade_2))# * -1 to minimize loglikelihood
}
##~~~~~~~~~~~~~-~~~~~~~~~~~~~~~~~~~~~~~~~~---------------------------------------------
# Run maximum likelihood estimation ------
# Read bbmle package about how to set optimizer and the parameters
fit_1 <- mle2(kernel_estimate_joint, start = list(k0 = 0.003368303 ,r0 = 0.4315369,alpha= 2.803289, delta = 0.0006193897),
                         skip.hessian = FALSE,method = "nlminb",  optimizer =  "optimx",
                         lower = c(k0 = 0.0001 ,r0 = 0.0001,alpha= 0.0001, delta =0))

summary(fit_1)
