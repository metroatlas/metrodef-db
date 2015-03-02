library(igraph)
library(ggplot2)
library(ggdendro)
library(ape)
library(rjson)
require(rjson)
setwd("~/Desktop/GitHub/metrodef-db2/metrodef-db")
# Create function that makes commuting adjacency matrix for a given state, hierarchically clusters
# the counties based on commuting pct, and then displays a simple node-edge visualization of the 
# result for a desired number of clusters
Cluster = function(state.name, census.data, type="single", county.pop){
  
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

  # make it so that greater integration => shorter distances
  maxelement = max(adj.mat)
  addmsize = ncol(adj.mat)
  addmat<-  matrix(rep(maxelement,each=addmsize*addmsize),nr=addmsize)
  s = addmat - adj.mat
  reversedist = as.dist(s)
  
  #print(reversedist)
  # adjust dendrogram for visibility
  hci = hclust(reversedist, "single")
  hc.mi = min(hci$height)
  print(hc.mi)
  
  submat = matrix(rep(hc.mi,each=addmsize*addmsize),nr=addmsize)
  hc.dist = as.dist(s - submat)
  hc = hclust(hc.dist, "single")

  labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
  categories = cutree(hc,k=4)

  # uses the ape library to plot a radial dendrogram
  pdf(file=paste0('radial/', state.name,'_radial.pdf'), width=20, height=20)
  par(mar=c(4,1,3,6))
  plot(as.phylo(hc), type="fan", cex=0.8)
  dev.off()

  # uncomment this to plot a tree dendrogram 
  #ggd <- ggdendrogram(hc, rotate=TRUE, size=1)
  #plot(ggd)

  print(min(adj.mat))
  print(max(addmat-adj.mat))
  HCExport(hc, file_out="/Users/liezl/Desktop/calidata2.js")
  
  return
}

#convert output from hclust into a nested JSON file
HCExport<-function(hc, file_out){  
  labels<-hc$labels
  merge<-data.frame(hc$merge)
  
  #cophenetic distances
  cd = cophenetic(hc)
  
  #heights
  h = hc$height
  print(h)
  
  # this will get the name of the first link (Yolo County + Sacramento County)
#  i_1 = which(as.matrix(cd)==h[1], arr.ind=T)
#   print(i_1)
#   print(colnames(as.matrix(cd))[i_1[1]])
#   print(i_1[1])
#   print(min(as.matrix(cd)[42,]))

  # Create nested node structure
  for (i in (1:nrow(merge))) {
    if (merge[i,1]<0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,1]]),list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]>0 & merge[i,2]<0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(node", merge[i,1], ", list(name=labels[-merge[i,2]])))")))}
    else if (merge[i,1]<0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(list(name=labels[-merge[i,1]]), node", merge[i,2],"))")))}
    else if (merge[i,1]>0 & merge[i,2]>0) {eval(parse(text=paste0("node", i, "<-list(name=\"node", i, "\", children=list(node",merge[i,1] , ", node" , merge[i,2]," ))")))}
  }
  eval(parse(text=paste0("JSON<-toJSON(node",nrow(merge), ")")))
  
  #wrap nested JSON file into d3 html
  fileConn<-file(file_out)
  writeLines(paste0("data='", JSON, "'"), fileConn)
  close(fileConn)
}


statenames = unique(census.data$RES_State)
for(stName in statenames) {
  Cluster(stName, census.data, 'single', county.pop)
}
# Cluster('Maryland', census.data, 'single', county.pop)
# Cluster('California', census.data, 'single', county.pop)
# Cluster('Texas', census.data, 'single', county.pop)


















