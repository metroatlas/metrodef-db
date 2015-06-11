library(igraph)
library(ggplot2)
library(ggdendro)
library(ape)
library(rjson)
library(grid)
library(gridExtra)
require(rjson)

setwd("~/Desktop/GitHub/metrodef-db2/metrodef-db")



# =========== C L U S T E R I N G ============
# These functions, cumulativeGraph, getClusters, ClusterByState, ClusterByRegion, and Cluster, all provide functionality to transform commuting data
# into a distance matrix and then produce a hclust object. getClusters parses clustering data to get area definitions. ClusterByState and ClusterByRegion
# take a subset of data and pass it to Cluster. Cluster creates the distance matrix, runs hierarchical clustering, and then calls cumulativeGraph
# to produce an overview of the clustering behavior. The hclust object returned by Cluster can be further processed.


# FUNCTION: cumulativeGraph
# -------------------------
# Takes a dataframe and produces a graph.
# Dataframe must have labels: threshold, cumul
# Resulting graph: threshold on x axis, cumul (total number/percent of metro areas when cutting the hclust 
# object at that threshold) on y axis
cumulativeGraph <- function(cgraph.d, filen) {
  #print(cgraph.d) # for debugging only
  p <- ggplot(cgraph.d, aes(x = cgraph.d$threshold, y = cgraph.d$cumul), environment=environment()) + geom_line()
  ggsave(filename=filen, plot=p)
  return(p)
}

# FUNCTION: getClusters
# ------------------------
# This function gets a definition of clusters from a .csv file, and produces a list of clusters. 
# Output is a list that contains, for each region, a vector with all the fips codes in that region.
# EXAMPLE USAGE:
# For the large commute sheds we used in the project, we used the 3099th column (colname = 3099)
# labelled "X1.34703368193124e.06" in the .csv file named "us_avg_linkage."
# The following line can be found in the main code at the end of this file:
#       regions <- getClusters(3099, "us_avg_linkage.csv")
getClusters = function(colnum, csvname){
  data <- read.csv(csvname) 
  clustercol <- data[[colnum]] # read data from the specified column into a vector
  numclust <- max(clustercol) # total number of clusters is the maximum cluster label in that column

  # Correct fips codes with leading zeroes if applicable; not translated into numeric value. These are extracted from the first column of the
  # .csv file.
  fipscodes <- sprintf("%05d", data[[1]])
  
  clist <- list() # new list
  
  # i = current cluster number (used to identify what counties belong in which clusters); counties with the same "cluster number"
  # are in the same cluster.
  for(i in 1:numclust){ 
    # Get the indecis of the clustercol vector (the column we're extracting the cluster definitions from) that match the current cluster number.
    matches <- which(clustercol == i)
    
    # Get the corresponding fips codes from the fipscodes vector.
    allc <- fipscodes[matches]
    
    # Append this vector of matching fips codes to our list.
    clist <- c(clist, list(allc))
  }
  return(clist)
}

# FUNCTION: ClusterByState
# ------------------------
# Gather commuting data only for particular state then call Cluster with this subset of the 
# commuting data.
ClusterByState = function(state.names, census.data, type="single", county.pop){
  area.data <- census.data[census.data$RES_State %in% state.names,]
  return(Cluster(area.data, type, county.pop, state.names[1])) #TODO: concatenate state names for output filename
}

# FUNCTION: ClusterByRegion
# -------------------------
# Gather commuting data for specified residence fips codes then call Cluster with this subset
# of the commuting data.
ClusterByRegion = function(fips.codes, census.data, type="single",county.pop, filename="RegionX"){
  #print(fips.codes) # for debugging only
  
  # Gets data where the residence fips code matches one of the fips codes we're looking for.
  area.data = census.data[census.data$RES_FIPS %in% fips.codes,]
  
  # Cluster with this subset of the commuting matrix data.
  return(Cluster(area.data, type, county.pop, filename))
}

# FUNCTION: Cluster
# -----------------
# Makes commuting adjacency matrix for a given subset of the commuting data, 
# then hierarchically clusters the counties based on commuting pct.
Cluster = function(area.data, type="single", county.pop, export.file){
  # Collect county names
  res.county.names = unique(area.data$RES_County)
  wrk.county.names = unique(area.data$WRK_County)
  
  # Merge work and county names
  county.names = sort(union(res.county.names, wrk.county.names))
  
  # Create adjacency matrix
  adj.mat = matrix( , nrow = length(county.names), ncol = length(county.names))
  adj.mat[,] = 0
  
  # Fill in adjacency matrix
  for(i in 1:nrow(area.data)){
    comm.pct = area.data$Number_PCT[i] - area.data$MOE_PCT[i]
    res.indx = which(county.names==area.data$RES_County[i])
    wrk.indx= which(county.names==area.data$WRK_County[i])
    adj.mat[res.indx, wrk.indx] = comm.pct
  }
  
  # Strip adjacency matrix of columns/rows that are entirely zero (or close enough to zero)
  connected.rows = apply(adj.mat, 2, function(x){any(x>0.0001)})
  connected.cols = apply(adj.mat, 1, function(x){any(x>0.0001)})
  connected = connected.rows | connected.cols
  
  # Collect subset of county.names
  labs = county.names[connected]
  indices = 1:length(county.names)
  indices = indices[connected]
  adj.mat = adj.mat[connected,connected]
  colnames(adj.mat) = labs
  rownames(adj.mat) = labs
  
  # Size n of n by n adjacency matrix is also the number of counties
  matsize = ncol(adj.mat)
  
  # Compute symmetrical integration matrix
  adj.diag <- diag(adj.mat)
  adj.mat <- (adj.mat + t(adj.mat))

  # This line that preserves the diagonal isn't actually necessary because the matrix -> distance object
  # transformation (as.dist()) just drops the diagonal without using it anyway
  diag(adj.mat) <- adj.diag

  # Make it so that greater integration => shorter distances
  maxelement <- max(adj.mat)

  addmat <-  matrix(rep(maxelement,each=matsize*matsize),nr=matsize)
  s = addmat - adj.mat
  print(s[47])
  reversedist = as.dist(s)
  
  hci = hclust(reversedist, type)
  res <- hci$height # initialize the vector to be of length hci$height

   for(cutnum in length(hci$height):1){
    clusMember <- cutree(hci, cutnum)
    metroareas <- 0
    for(i in 1:cutnum){
      currClust <- names(clusMember[clusMember == i])
      
      pop <- sum(subset(x = county.pop, subset = County %in% currClust,select = c(April.1.2010.Census))$April.1.2010.Census)
      if(pop >= 50000){
        metroareas = metroareas + 1
      }
    }
    res[length(hci$height) - cutnum + 1] <- metroareas
  }
  print(res)
  
  
  integ <- maxelement - hci$height # reverse height to get integration

  #print("MAX HEIGHT") # for debugging only
  #print(max(hci$height)) # for debugging only
  
  # set up the dataframe to create the cumulative graph
  cgraph.d <- data.frame(integ, res)
  colnames(cgraph.d) <- c("threshold", "cumul")

  # export data for shiny app -- only use if you want to export data
  # temp.frame <- cgraph.d
  # temp.frame$Region <- export.file
  # assign("export.frame", rbind(export.frame, temp.frame), envir=.GlobalEnv)
  
  cumulativeGraph(cgraph.d, paste0(export.file, "_cumulative.pdf"))
  cgraph.d$cumul = cgraph.d$cumul / length(hci$height)
  cumulativeGraph(cgraph.d, paste0(export.file, "_cumulativeNorm.pdf"))
  
  # adjust dendrogram for visibility (while keeping the scale) and return the hclust object
  hc.mi = min(hci$height)
  submat = matrix(rep(hc.mi,each=matsize*matsize),nr=matsize)
  hc.dist = as.dist(s - submat)
  hc = hclust(hc.dist, type)
  return(hc)
}


# ======== D E N D R O G R A M S ==========
#
PlotRadial = function(hc, file.path="") {
  if(file.path != "") {
    # uses the ape library to plot a radial dendrogram
    pdf(file=file.path, width=20, height=20)
    par(mar=c(4,1,3,6))
    plot(as.phylo(hc), type="fan", cex=0.8)
    dev.off()
  } else {
    plot(as.phylo(hc), type="fan", cex=0.8)
  }
  return
}

PlotTree = function(hc, file.path="") {
  ggd <- ggdendrogram(hc, rotate=TRUE, size=1)
  if(file.path != "") {
    pdf(file=file.path, width=20, height=20)
    plot(ggd)
    par(mar=c(4,1,3,6))
    dev.off()
  } else {
    plot(ggd)
  }
  return
}

# FUNCTION: HCExport
# ------------------
# Convert output from hclust into a nested JSON file. This makes it easy to render the tree in d3.js.
HCExport<-function(hc, file_out){
  labels<-hc$labels
  merge<-data.frame(hc$merge)
  
  # Create nested node structure
  for (i in (1:nrow(merge))) {
    if (merge[i,1]<0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,1]]),list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]>0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(node", merge[i,1], ", list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]<0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,1]]), node", merge[i,2],"))")))}
    else if (merge[i,1]>0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(node",merge[i,1] , ", node" , merge[i,2]," ))")))}
  }
  eval(parse(text=paste0("JSON<-toJSON(node",nrow(merge), ")")))
  
  # Wrap nested JSON file into d3 html
  fileConn<-file(file_out)
  writeLines(paste0("data='", JSON, "'"), fileConn)
  close(fileConn)
}

# ========= M A I N ==========
# statenames = unique(census.data$RES_State) # for debugging only

# hc = ClusterByState(c('California', 'Massachusetts'), census.data, 'average', county.pop) # Cluster 2 random states - testing only
# PlotRadial(hc, paste0('radial/', 'CaliMA_radial.pdf')) # plot a radial graph with the 2 test states

regions <- getClusters(3099, "us_avg_linkage.csv") # you only need to run this once and then the regions list can be re-used 
# print(regions) # for debugging only

# export.frame <<- data.frame(height=numeric(), cumul=numeric(), Region=character()) # export for Shiny app -- global dataframe variabl
for(i in 1:length(regions)){
  if(length(regions[[i]]) > 1 ){
    hc <- ClusterByRegion(regions[[i]], census.data, 'average', county.pop, paste0("Region", i))
  }
}

#print(export.frame) # for debuggin purposes (export for Shiny app)
#colnames(export.frame) <- c("integration", "cumulative", "Region") # export for Shiny app
#write.table(export.frame, "exportcumul2.csv", sep=",") # export for Shiny app

hc = ClusterByState("California", census.data, 'average', county.pop)
hc = ClusterByState("New Jersey", census.data, 'average', county.pop)

#dir.create('linegraph_sing', showWarnings = FALSE) # data export purposes
#dir.create('linegraph_comp', showWarnings = FALSE) # data export purposes
#dir.create('linegraph_avg', showWarnings = FALSE) # data export purposes

# Cluster all of US (use the entire set of commuting data) using avg linkage. Export file => "AllUS"
hc = Cluster(census.data, 'average', county.pop, "AllUS")

for(st.name in state.names) {
  hc = ClusterByState(c(st.name), census.data, 'single', county.pop)
  hc.comp = ClusterByState(c(st.name), census.data, 'complete', county.pop)
  hc.avg = ClusterByState(c(st.name), census.data, 'average', county.pop)

  # Example usage below - not actually used for presentation
  #   PlotRadial(hc, paste0('radial_sing/', st.name,'_radial_sing.pdf'))
  #   PlotRadial(hc.comp, paste0('radial_comp/', st.name,'_radial_comp.pdf'))
  #   PlotRadial(hc.avg, paste0('radial_avg/', st.name,'_radial_avg.pdf'))
  #   PlotRadial(hc.wardd, paste0('radial_wardd/', st.name,'_radial_wardd.pdf'))
  #   PlotRadial(hc.wardd2, paste0('radial_wardd2/', st.name,'_radial_wardd2.pdf'))
  #   PlotRadial(hc.mcquitty, paste0('radial_mcquitty/', st.name,'_radial_mcquitty.pdf'))
  #   PlotRadial(hc.median, paste0('radial_median/', st.name,'_radial_median.pdf'))
  #   PlotRadial(hc.centroid, paste0('radial_cent/', st.name,'_radial_cent.pdf'))
  #
  #   HCExport(hc, paste0('hcexport/', st.name,'_hclust.json'))
  
}

