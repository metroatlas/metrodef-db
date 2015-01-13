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
	cLinks[c(1:10)] <- lapply(cLinks[c(1:10)], as.character)

	#Remove Puerto Rico/ Foreign Countries (Iraq, At Sea, Africa, etc.)
	# TODO