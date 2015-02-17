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
  # Hierarchically cluster according to distances
  links = HCluster(adj.mat, type, num.clusters,labs, county.pop)
 
  # Produce plot of the clusters
  # ClusterViz(cluster.list, labs)
  
  return(links)
}

HCluster = function(adj.mat, type="single", num.clusters=10, labs, county.pop) {
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
        cluster1 = unlist(cluster.list[[i]])
        cluster2 = unlist(cluster.list[[j]])
        # Find distance between all of the members of the clusters
        cluster.dist = matrix(,nrow=length(cluster1), ncol=length(cluster2))
        cluster.dist[,] = 0
        for(k in length(cluster1)){
          for(t in length(cluster2)){
            cluster.dist[k,t] = adj.mat[cluster1[k],cluster2[t]]
          }
        }
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
    print(max.dist)
    print(length(cluster.list))
    if(max.dist==0){
      links = ListLinks(cluster.list, labs, adj.mat, max.dist, county.pop)
      return(links)
    }
    max.ind = which(dist.mat == max.dist, arr.ind=TRUE)
    print(max.ind)
    print(max.dist)
    # Join "nearest" clusters
    for(m in 1:nrow(max.ind)){
      cluster.list[[max.ind[m,1]]] = list(cluster.list[[max.ind[m,1]]], cluster.list[[max.ind[m,2]]])
      cluster.list[[max.ind[m,2]]] = NULL
      max.ind[max.ind >= max.ind[m,2]] = max.ind[max.ind >= max.ind[m,2]] - 1
    }
  }
  
  # Make links list
  links = ListLinks(cluster.list, labs, adj.mat, max.dist, county.pop)
  return(links)
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
        if(TRUE){
          pct = adj.mat[curr.clust[j],curr.clust[k]]
          if(TRUE){
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