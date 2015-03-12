makeMaps <- function() {
  
  data <- read.csv("data/cluster_table.csv")
  colnames(data) <- c("geoid", paste("c", 1:68, sep=""))
  data$geoid <- sprintf("%05d", data$geoid)
  data <- data[grep("^06", data$geoid),]
  
  for(i in 3:70){
    data[,i] <- as.factor
  }
  
  vMakeMapCat(2:69)
}

makeMapCol <- function(col) {
  library(ggplot2)
  library(maptools)
  library(RColorBrewer)
  library(plyr)
  
  county <- readShapeSpatial("geo/ca_county/ca_county.shp")
  gpclibPermit()
  county <- fortify(county, region="GEOID")
  
  distcenters <- ddply(county, .(id), summarize, clat = mean(lat), clong = mean(long))
  
  distcenters$cluster <- data[,"c20"]
  
  plot <- ggplot() + geom_map(data = data, aes(map_id = geoid, fill = as.factor(c20 %% 12)), 
                      map = county) + expand_limits(x = county$long, y = county$lat) +
    coord_map(projection="mercator") +
    scale_fill_brewer(palette = "Paired") +
    geom_text(data = distcenters, aes(x = clong, y = clat, label = cluster, size = 0.2))
  
  filename  <- paste("exports/", col, ".pdf", sep="", collapse="")
  ggsave(filename = filename,
         plot = plot,
         width = 21,
         height = 29.7,
         unit = 'cm')
}

vMakeMapCol <- Vectorize(makeMapCol)

makeMapCol("c20")