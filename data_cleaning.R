library(igraph)

setwd("~/Desktop/BillLane_Research/metrodef-db")
dat = read.csv("clean.csv", fileEncoding="latin1")
county.pop = read.csv("county_pop.csv", fileEncoding="latin1")

county.pop$State = as.character(county.pop$State)
# Remove leading space in state names
county.pop$State = substr(county.pop$State, 2, nchar(county.pop$State))
county.pop$County = paste(as.character(county.pop$County), ", ", as.character(county.pop$State), sep="")

PadZeros = function(col, size) {
  # Convert column to string
  newCol = as.character(col)
  
  # Pad zeros
  for(i in 1:size-1){
    newCol[nchar(newCol)==i] = paste0(paste0(rep("0",size-i),collapse=""), newCol[nchar(newCol)==i])
  }
  
  return(newCol)
}

dat[,1] = PadZeros(dat[,1], 2)
dat[,2] = PadZeros(dat[,2], 3)
dat[,3] = PadZeros(dat[,3], 3)
dat[,4] = PadZeros(dat[,4], 3)

# Remove rows with Puerto Rico residence
dat_nopr = dat[dat$RES_State != "Puerto Rico",]

# Remove rows with people working overseas (e.g. military)
dat_nomil = dat_nopr[as.numeric(dat_nopr$WRK_State_FIPS_Code) <= 56,]

# Finalize data and store in properly named variable
census.data = dat_nomil
for(i in 7:ncol(census.data)){
  census.data[,i] = as.character(census.data[,i])
}
dat_nomil = NULL
census.data$RES_County = paste(census.data$RES_County, ", ", census.data$RES_State, sep="")
census.data$WRK_County = paste(census.data$WRK_County, ", ", census.data$WRK_State, sep="")
# Collect state names
state.names = unique(census.data$RES_State)

# Add percentage column to census.data
num.pct = rep(0, nrow(census.data))
moe.pct = rep(0,nrow(census.data))
population = rep(0,nrow(census.data))
for(i in 1:nrow(county.pop)){
  print(paste0("Processing Population of County ", i, " of 3143"))
  indices = intersect(which(census.data$RES_County == county.pop$County[i]), 
            which(census.data$RES_State == county.pop$State[i]))
  
  pop = county.pop$April.1.2010.Census[i]
  
  num.pct[indices] = (census.data$Number[indices] / pop)
  
  moe.pct[indices] = (census.data$MOE[indices] / pop)
  
  population[indices] = pop
}

census.data$Number_PCT = num.pct
census.data$MOE_PCT = moe.pct
census.data$RES_pop = population

# Create function that makes commuting adjacency matrix for a given state
CMatrix = function(state.name, census.data, threshold=0.15, weights=FALSE){
  
  # Gather commuting data only for particular state
  state.data = census.data[census.data$RES_State == state.name,]
  
  # Collect county names
  res.county.names = unique(state.data$RES_County)
  wrk.county.names = unique(state.data$WRK_County)
  
  # Merge work and county names
  county.names = sort(union(res.county.names, wrk.county.names))
  
  # Create adjacency matrix
  adj.mat = matrix( , nrow = length(county.names), ncol = length(county.names))
  adj.mat[,] = 0
  
  # Fill in adjacency matrix
  for(i in 1:nrow(state.data)){
    comm.pct = state.data$Number_PCT[i] - state.data$MOE_PCT[i]
    res.indx = which(county.names==state.data$RES_County[i])
    wrk.indx= which(county.names==state.data$WRK_County[i])
    if(comm.pct >= threshold && weights == FALSE) {
      adj.mat[res.indx, wrk.indx] = 1
    } else if(comm.pct >= threshold && weights == TRUE){
      adj.mat[res.indx, wrk.indx] = comm.pct
    }
  }
  
  # Strip adjacency matrix of columns/rows that are entirely zero
  connected.rows = apply(adj.mat, 2, function(x){any(x>0)})
  connected.cols = apply(adj.mat, 1, function(x){any(x>0)})
  connected = connected.rows | connected.cols
  
  # Collect subset of county.names
  labs = county.names[connected]
  
  # Plot the adjacency matrix
  plot(graph.adjacency(adj.mat[connected,connected], mode="directed", weighted=TRUE), vertex.size=10, vertex.label=labs,
       vertex.label.cex=0.6, edge.arrow.size=0.15)
  
  return(adj.mat)
}