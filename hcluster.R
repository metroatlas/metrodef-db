library(igraph)
library(ggplot2)
library(ggdendro)
library(ape)

# Create function that makes commuting adjacency matrix for a given state, hierarchically clusters
# the counties based on commuting pct, and then displays a simple node-edge visualization of the 
# result for a desired number of clusters
Cluster = function(state.name, census.data, type="single", num.clusters=10, county.pop){
  
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

  maxelement = max(adj.mat)
  addmsize = ncol(adj.mat)
  addmat<-  matrix(rep(maxelement,each=addmsize*addmsize),nr=68)
  s = addmat - adj.mat
  mi = min(s+10000*(s<=0))
  smat <-  matrix(rep(mi,each=addmsize*addmsize),nr=68)
  print('min:')
  print(mi)
  s = s - smat
  #reversedist = as.dist(addmat - adj.mat)
  reversedist = as.dist(s)
  
  #print(reversedist)
  hc = hclust(reversedist, "single")
  labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
  categories = cutree(hc,k=20)
  #plot(as.phylo(hc), type="unrooted")
  #hc <- hclust(dist(USArrests), "ave")
  x <- ggdendrogram(hc, rotate=FALSE, size=1)
  plot(x)
#   dhc <- as.dendrogram(hc)
#   # Rectangular lines
#   ddata <- dendro_data.dendrogram(dhc, type = "triangle")
#   p <- ggplot(segment(ddata)) + 
#     geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#     coord_flip() + 
#     scale_y_reverse(expand = c(0.2, 0)) +
#     theme_dendro()
#   plot(p)
  #plot(hc)
  print(min(adj.mat))
  print(max(addmat-adj.mat))
  HCExport(hc, file_out="/Users/liezl/Desktop/calidata2.js")
  
  # Hierarchically cluster according to distances
  #links = HCluster(adj.mat, type, num.clusters,labs, county.pop)
 
  # Produce plot of the clusters
  #ClusterViz(cluster.list, labs)
  
  #return(links)
  return
}

HCluster = function(adj.mat, type="single", num.clusters=21, labs, county.pop) {
  # This function takes in an adjacency matrix between counties, a list of the
  # county names, and a 'type' argument which specifies how the clustering is to be performed.
  # Distances need not be recalculated in this function, only retreived from the adjacency matrix.
  
  # Initialize clusters to be all counties
  cluster.list = list()
  for(i in 1:nrow(adj.mat)){
    cluster.list[[i]] = i
  }
  
  while(length(cluster.list)>num.clusters){
    print(paste0("Number of clusters: ", length(cluster.list)))
    dist.mat = matrix(,nrow=length(cluster.list), ncol=length(cluster.list))
    dist.mat[,] = 0
    # Calculate distance between clusters
    for(i in 2:length(cluster.list)){
      for(j in 1:(i-1)){
        # Retrieve candidate clusters
        cluster1 = cluster.list[[i]]
        cluster2 = cluster.list[[j]]
        # Find distance between all of the members of the clusters
        cluster.dist = adj.mat[unlist(cluster1),unlist(cluster2)]
        if(type == "single"){
          # In this case, our "distance" actually increases with similarity.
          # Thus, single linkage should consider the cluster-distance to be the maximum.
          dist = max(cluster.dist)
        } else if(type == "complete"){
          dist = min(cluster.dist)
        } else if(type == "average"){
          dist = mean(cluster.dist)
        }
        dist.mat[i,j] = dist
      }
    }
    
    # Find index of maximum element of the distance matrix
    max.dist = max(dist.mat)
    if(max.dist==0){return(cluster.list)}
    max.ind = which(dist.mat == max.dist, arr.ind=TRUE)
    
    # Join "nearest" clusters
    for(m in 1:nrow(max.ind)){
      cluster.list[[max.ind[m,1]]] = list(cluster.list[[max.ind[m,1]]], cluster.list[[max.ind[m,2]]])
      cluster.list[[max.ind[m,2]]] = NULL
      max.ind[max.ind >= max.ind[m,2]] = max.ind[max.ind >= max.ind[m,2]] - 1
    }
  }
  # Make links list
  links = ListLinks(cluster.list, labs, adj.mat, max.dist, county.pop)
  ClusterViz(cluster.list, labs)
  return(links)
  
}

ClusterViz = function(cluster.list, labs){
  adj.mat = matrix(,nrow=length(labs), ncol=length(labs))
  adj.mat[,] = 0
  for(i in 1:length(cluster.list)){
    curr.clust = unlist(cluster.list[[i]])
    for(j in 1:length(curr.clust)){
      for(k in 1:length(curr.clust))
        if(curr.clust[j] != curr.clust[k]){
          adj.mat[curr.clust[j],curr.clust[k]] = 1
        }
    }
  }
  plot(graph.adjacency(adj.mat, mode="undirected"), vertex.size=10, vertex.label=labs,
  vertex.label.cex=0.6, edge.arrow.size=0.15)
}

ListLinks = function(cluster.list, labs, adj.mat, dist.threshold, county.pop){
  links = list()
  pops = rep(0,length(labs))
  Source = c()
  Target = c()
  Cluster.id = c()
  Comm.pct = c()
  s.population = c()
  t.population = c()
  
  for(i in 1:length(labs)){
    pops[i] = county.pop$April.1.2010.Census[county.pop$County == labs[i]][1]
  }
  
  for(i in 1:length(cluster.list)){
    curr.clust = unlist(cluster.list[[i]])
    for(j in 1:length(curr.clust)){
      for(k in 1:length(curr.clust))
        if(curr.clust[j] != curr.clust[k]){
          pct = adj.mat[curr.clust[j],curr.clust[k]]
          if(pct >= dist.threshold){
            Source = c(Source, labs[curr.clust[j]])
            Target = c(Target, labs[curr.clust[k]])
            Cluster.id = c(Cluster.id, i)
            Comm.pct = c(Comm.pct, pct)
            s.population = c(s.population, pops[curr.clust[j]])
            t.population = c(t.population, pops[curr.clust[k]])
          }
      }
    }
  }
  
  links$"source" = Source
  links$"target" = Target
  links$"cluster" = Cluster.id
  links$"pct" = Comm.pct
  links$"srcpop" = s.population
  links$"tarpop" = t.population
  
  return(links)
}

library(rjson)
require(rjson)
#convert output from hclust into a nested JSON file

HCExport<-function(hc, file_out){
  
  labels<-hc$labels
  merge<-data.frame(hc$merge)
  cd = cophenetic(hc)
  #print(cd)
  h = hc$height
  print(h)
  #FAILED - plot(cut(as.dendrogram(hc),h=0.41))
  i_1 = which(as.matrix(cd)==h[1], arr.ind=T)
#   print(i_1)
#   print(colnames(as.matrix(cd))[i_1[1]])
#   print(i_1[1])
#   print(min(as.matrix(cd)[42,]))
  
   

  
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






Cluster('California', census.data, 'single', 42, county.pop)



















