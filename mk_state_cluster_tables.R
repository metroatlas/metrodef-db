setwd("~/Desktop/BillLane_Research/metrodef-db")
source("data_cleaning.R")
source("hcluster.R")
US = c("Alabama", "Arizona", "Arkansas", "California",
       "Colorado", "Connecticut", "Delaware", "Florida", "Georgia", "Idaho",
       "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana",
       "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", 
       "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada",
       "New Hampshire", "New Jersey", "New Mexico", "New York", 
       "North Carolina" ,"North Dakota", "Ohio", "Oklahoma", "Oregon",
       "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee",
       "Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia",
       "Wisconsin", "Wyoming")

for(i in 1:length(US)){
  single = Cluster(US[i], census.data, type="single", num.clusters=1, county.pop, output="cluster.col")
  average = Cluster(US[i], census.data, type="average", num.clusters=1, county.pop, output="cluster.col")
  write.csv(single,file=paste0("clusterMap/data/",US[i],"_single_linkage.csv"))
  write.csv(average,file=paste0("clusterMap/data/",US[i],"_average_linkage.csv"))
}
