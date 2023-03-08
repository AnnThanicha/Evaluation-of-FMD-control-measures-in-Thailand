# Generate  network using the preferential attachment --------

# Required library --------------
packages <- c("dplyr", "reshape2","igraph", "readxl")
lapply(packages, library, character.only = TRUE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----------------------------------------------------
## import data -----------------
# the data contain node of
raw_vertices <- readRDS("raw_vertices.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----------------------------------------------------
## Generate trade network from preferential attachment -------------
# Function to generate network 

# raw_vertices is the dataframe containing ID and status of nodes
# prob_isolated is the probability of being isolated farms without trading

sim_trade_pref <- function(raw_vertices, prob_isolated) {
  
  raw_vertices$isolated_status<-0
  raw_vertices$isolated_status[raw_vertices$status=="farm"] <-rbinom(length(raw_vertices$isolated_status[raw_vertices$status=="farm"]), size = 1, prob= prob_isolated) 
  # trader is not isolated
  raw_vertices$isolated_status[raw_vertices$status=="trader"] <-1 
  table(raw_vertices$isolated_status[raw_vertices$status=="farm"])
  raw_vertices_trader<-raw_vertices[raw_vertices$status=="trader",]#node trader
  raw_vertices_farmer<-raw_vertices[raw_vertices$status=="farm",]#node farmer
  raw_vertices_farmeredge <-raw_vertices[raw_vertices$isolated_status==1 & raw_vertices$status=="farm",]
  # set initail edge to one for trader
  raw_vertices_trader$initialedge <- 1
  
  # add random number of edge
  raw_vertices_farmeredge$edge_number <- sample(c(1,2,3,4), nrow(raw_vertices_farmeredge), replace = TRUE, prob= c(44/63, 16/63, 2/63, 1/63))
  
  # arrange the ID so farm with 1 edges start first
  raw_vertices_farmeredge <- raw_vertices_farmeredge[order(raw_vertices_farmeredge$edge_number),]
  # create blank vector
  sim_edge <-c()
  
  for(i in 1:length(raw_vertices_farmeredge$ID)){
    # if there are still traders without the link to farmers, it will be prioritized to link with farmer first
    # to gaurantee that all trader have at least one link to the farmers
    # after all trader got at least one link it will be
    if (sum(raw_vertices_trader$initialedge==1)>0){add <- sample(raw_vertices_trader$ID[raw_vertices_trader$initialedge==1], raw_vertices_farmeredge$edge_number[i], replace = FALSE)
    #update edge for selected trader
    raw_vertices_trader$initialedge[raw_vertices_trader$ID %in% add] <-raw_vertices_trader$initialedge[raw_vertices_trader$ID %in% add]+1
    }
    else    
    { #random edge with preferential attachment
      add <- sample(raw_vertices_trader$ID, raw_vertices_farmeredge$edge_number[i], replace = FALSE, prob= (raw_vertices_trader$initialedge/sum(raw_vertices_trader$initialedge)))
      #update edge for selected trader
      raw_vertices_trader$initialedge[raw_vertices_trader$ID %in% add] <-raw_vertices_trader$initialedge[raw_vertices_trader$ID %in% add]+1
    }
    # keep add to sim_edge vecyor
    sim_edge <- c(sim_edge,add)
    
  } 
  
 # To create the trade matrix for simulation model
  sim_edge <- data.frame(from = rep(raw_vertices_farmeredge$ID, raw_vertices_farmeredge$edge_number), to = sim_edge)
  
  # group the farm with the same trader
  list_trade <- split(sim_edge, as.factor(sim_edge$to))
  list_trade  <- lapply(list_trade, "[[", 1)
  # change to numeric
  list_trade  <- lapply(list_trade, as.numeric, 1)
  #First create blank trade matrix 
  component_matrix<-matrix(0, nrow=502, ncol=502)
  
  # the farms from same component is 1 , while farms that are not in the same component is 0
  for(i in 1:length(list_trade)){
    component_matrix[list_trade[[i]],list_trade[[i]]]<-1   }
  # diagonal to 0
  diag(component_matrix)<-0 
  
  return(component_matrix)
}

##~~~~~~~~~~~~~~~~~~~~~~-~~~~~~~~~~~~~~~~~~~~~~~~~---------------------------------
 simulated_trade <-sim_trade_pref (raw_vertices = raw_vertices, prob_isolated = 0.5)
