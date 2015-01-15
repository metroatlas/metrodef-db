C_CommuterFlow_2010 <- function(d=TRUE)

  #Download data spreadsheet
  if(d) {
    fileUrl <- "http://www.census.gov/population/metro/files/commuting/Table1.xlsx"
    download.file(fileUrl, destfile="data/C_CommuterFlow_2010.xlsx", method="curl")
    dateDownloaded <- date()
    write(dateDownloaded,file="data/C_CommuterFlow_2010.xlsx.date.txt")
  }

  #Load data
  library(xlsx)
  cLinks <- read.xlsx("data/C_CommuterFlow_2010.xlsx",
    sheetIndex=1,
    colIndex=1:10,
    rowIndex=4:136799,
    colClasses="character",
    header=TRUE)

  #Transform names
  names(cLinks) <- gsub("\\.", "", names(cLinks))
  names(cLinks)[2] <- "ResCountyFIPS"
  names(cLinks)[4] <- "WorkCountyFIPS"
  cLinks[c(1:10)] <- lapply(cLinks[c(1:10)], as.character)

  #Remove Puerto Rico/ Foreign Countries (Iraq, At Sea, Africa, etc.): where FIPS Code > 56
  cLinks$ResFIPSNum <- lapply(cLinks$ResCountyFIPS, as.integer)
  cLinks$WorkFIPSNum <- lapply(cLinks$WorkCuntyFIPS, as.integer)
  cLinks <- cLinks[(cLinks$ResFIPSNum <= 56 | cLinks$WorkFIPSNum <= 56),]

  #Drop extra columns used for processing
  cLinks$ResFIPSNum <- NULL
  cLinks$WorkFIPSNum <- NULL

