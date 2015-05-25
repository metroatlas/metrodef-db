setwd("~/Desktop/BillLane_Research/metrodef-db")
source("data_cleaning.R")
source("link_list_avg.R")
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
  average = Cluster(US[i], census.data, type="average", num.clusters=1, county.pop, output="links")
  write.csv(average,file=paste0("networkData/",US[i],"_average_linkage.csv"))
}
