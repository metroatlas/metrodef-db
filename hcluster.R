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
Cluster = function(state.names, census.data, type="single", county.pop){
  
  # Gather commuting data only for particular state
  state.data = census.data[census.data$RES_State %in% state.names,]
  
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

  # size n of nxn adjacency matrix is also the number of counties
  matsize = ncol(adj.mat)
  
  # compute symmetrical integration matrix 
  adj.diag = diag(adj.mat)
  adj.mat <- 0.5 * (adj.mat + t(adj.mat))
  diag(adj.mat) <- adj.diag 
    # preserving the diagonal isn't actually necessary because as.dist() doesn't use the diagonal anyway
  
  
  # make it so that greater integration => shorter distances
  maxelement = max(adj.mat)
  addmat<-  matrix(rep(maxelement,each=matsize*matsize),nr=matsize)
  s = addmat - adj.mat
  reversedist = as.dist(s)
  
  # adjust dendrogram for visibility
  hci = hclust(reversedist, type)
  hc.mi = min(hci$height)
  print(hc.mi)
  
  submat = matrix(rep(hc.mi,each=matsize*matsize),nr=matsize)
  hc.dist = as.dist(s - submat)
  hc = hclust(hc.dist, type)
  return(hc)
}


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
  if(filePath != "") {
    pdf(file=file.path, width=20, height=20)
    plot(ggd)    
    par(mar=c(4,1,3,6))
    dev.off()
  } else {
    plot(ggd)
  }
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
hc = Cluster(c('California', 'Massachusetts'), census.data, 'single', county.pop)
PlotRadial(hc, paste0('radial/', 'CaliMA_radial.pdf'))




for(st.name in statenames) {
  hc = Cluster(c(st.name), census.data, 'single', county.pop)
  hc.comp = Cluster(c(st.name), census.data, 'complete', county.pop)
  hc.avg = Cluster(c(st.name), census.data, 'average', county.pop)
  hc.wardd = Cluster(c(st.name), census.data, 'ward.D', county.pop)
  hc.wardd2 = Cluster(c(st.name), census.data, 'ward.D2', county.pop)
  hc.mcquitty = Cluster(c(st.name), census.data, 'single', county.pop)
  hc.median = Cluster(c(st.name), census.data, 'median', county.pop)
  hc.centroid = Cluster(c(st.name), census.data, 'centroid', county.pop)
  
  
  PlotRadial(hc, paste0('radial_sing/', st.name,'_radial_sing.pdf'))
  PlotRadial(hc.comp, paste0('radial_comp/', st.name,'_radial_comp.pdf'))
  PlotRadial(hc.avg, paste0('radial_avg/', st.name,'_radial_avg.pdf'))
  PlotRadial(hc.wardd, paste0('radial_wardd/', st.name,'_radial_wardd.pdf'))
  PlotRadial(hc.wardd2, paste0('radial_wardd2/', st.name,'_radial_wardd2.pdf'))
  PlotRadial(hc.mcquitty, paste0('radial_mcquitty/', st.name,'_radial_mcquitty.pdf'))
  PlotRadial(hc.median, paste0('radial_median/', st.name,'_radial_median.pdf'))
  PlotRadial(hc.centroid, paste0('radial_cent/', st.name,'_radial_cent.pdf'))
  
  HCExport(hc, paste0('hcexport/', st.name,'_hclust.json'))
}

# labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
# categories = cutree(hc,k=4)

















