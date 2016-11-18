##################
# Load Libraries
##################

library(zipcode)
library(ggplot2)
library(utils)
library(geosphere)
library(randomForest)
library(plyr)

##################
# Load data
##################

data(zipcode)
dma.data <- read.csv("//Users//psullivan01//Documents//dmazip.csv", 
                     header=TRUE, sep=",")
zip.pop <- read.csv("//Users//psullivan01//Documents//Popbyzip.csv", 
                    header=TRUE, sep=",")	#zip code and population (33,092 rows)
dma.data2 <- read.csv("//Users//psullivan01//Documents//ziptodma.csv", 
                      header=TRUE, sep=",")

##################
# DMA Cleanup
##################

# deduplicate
dma.data2 <- dma.data2[!duplicated(dma.data2), ]
dma.data2$numeraire <- 1

# remove remaining zips with 2 entries - split across DMAs
zipcount <- aggregate(numeraire ~ nzip + LONGITUDE + LATITUDE, data=dma.data2, FUN=sum)
zipcount <- subset(zipcount, numeraire>1)
dma.data2 <- merge(dma.data2, zipcount, by="nzip", all.x=TRUE)
dma.data2 <- dma.data2[rowSums(is.na(dma.data2)) > 0,]
dma.data2<- dma.data2[1:5]
colnames(dma.data2) <- c("zip", "latitude", "longitude", "dma", "dmaname")

# build model off of zips with no overlap
dma.data2$dmafactor <- as.factor(dma.data2$dma)
zip.RF <- randomForest(dmafactor ~ longitude + latitude, data = dma.data2)

# apply model to zips with multiple DMAs
colnames(zipcount) <- c("zip", "longitude", "latitude", "numeraire")
zipcount$dma <- predict(zip.RF, newdata=zipcount)
zipcount <- subset(zipcount, select=c('zip', 'longitude', 'latitude', 'dma'))

# create dmaname lookup table and merge with zipcount
dmalookup <- subset(dma.data2, select=c("dma", "dmaname"))
dmalookup <- dmalookup[!duplicated(dmalookup), ]
zipcount <- merge(zipcount, dmalookup, by="dma")

# Fold predictions back into original dataset
dma.data2 <- dma.data2[1:5]
dma.data2 <- rbind(dma.data2, zipcount)

##################
# Format data
##################

dma.rank <- subset(dma.data, select=c("rank", "dma"))
dma.rank <- dma.rank[!duplicated(dma.rank), ]

# filter for top 25 DMAs				
colnames(zip.pop) <- c("zip", "pop")
dma.data2 <- merge(dma.data2, dma.rank, by="dma")

# add population data	
zip.pop <- zip.pop[!duplicated(zip.pop), ]
zip.pop <- aggregate(pop ~ zip, data=zip.pop, FUN=sum)		
zip.data <- merge(dma.data2, zip.pop, by="zip")

# Add leading zeros to zips shorter than 5 characters
zip.data$zip <- as.character(zip.data$zip)
zip.data$zip <- ifelse(nchar(zip.data$zip, type="chars")==4, 
                       paste0("0", zip.data$zip), zip.data$zip)
zip.data$zip <- ifelse(nchar(zip.data$zip, type="chars")==3, 
                       paste0("00", zip.data$zip), zip.data$zip)

# Add city and state
zip.data <- merge(zip.data, zipcode, by="zip")
zip.data <- zip.data[1:9]
colnames(zip.data)[3:4] <- c("latitude", "longitude")

# Calculate DMA populations 
dmapop <- aggregate(pop ~ dmaname, data=zip.data, sum)
colnames(dmapop) <- c("dmaname", "dmapop")
zip.data <- merge(zip.data, dmapop, by="dmaname")
zip.data$optimalsplit <- as.integer(zip.data$dmapop/1500000)+1
zip.data$optimalpercent <- 1/zip.data$optimalsplit

#Update Wheatfield, IN with correct coordinates
zip.data[1755,4:5] <- c(41.02044, -87.1024)

########################
# Create objects - North
########################

zip.data.n <- zip.data				
remainder.data.n <- 0
zip.data.n$zone <- 0
zip.data.n$latrank <- 0

dma.data.n <- list()
dmadata.n <- list()
city.data.n <- list()

#####################
# Begin loop - North
#####################

for (n in 1:25) #by dma rank
{
  my.data.n <- subset(zip.data, rank==n) #subsets by DMA rank
  my.data.n$ziprank <- rank(-my.data.n$pop, ties.method = c("first"), 0)
  citydata.n <- my.data.n
  
  for (c in 1:citydata.n$optimalsplit) #by city rank
  {
    citydata.n <- citydata.n[order(-citydata.n$latitude),]
    citydata.n$latrank <- rank(-citydata.n$latitude, ties.method = c("first"), 0)
    northmost <- subset(citydata.n, latrank==1, select=c("dmaname", "zip", 
                                                         "longitude", "latitude"))
    colnames(northmost) <- c("dmaname", "nzip", "nlongitude", "nlatitude")
    allzips.n <- merge(northmost, citydata.n, by="dmaname")
    distance.n <- 0
    result.data.n <- list()
    
    for (i in 1:nrow(allzips.n)) #calculate distances
    {
      record.n <- subset(allzips.n, as.numeric(rownames(allzips.n)) == i)
      distance.n[i] <- distm(record.n[,c('nlongitude', 'nlatitude')],
                             record.n[,c('longitude', 'latitude')],
                             fun = distHaversine)
      record.n$distance <- distance.n[i]
      result.data.n[[i]] <- record.n
    }
    allzips.n <- do.call(rbind, result.data.n)
    allzips.n <- allzips.n[order(allzips.n$distance),]
    
    allzips.n$runpop <- cumsum(allzips.n$pop)
    allzips.n$runpercent <- allzips.n$runpop/allzips.n$dmapop
    allzips.n$zone <-	ifelse(allzips.n$runpercent < allzips.n$optimalpercent, c, 0)
    
    city.data.n[[c]] <- subset(allzips.n, zone > 0)
    citydata.n <- subset(citydata.n, !(zip %in% city.data.n[[c]]$zip))
  }
  mydmadata.n <- do.call(rbind, city.data.n)
  dmadata.n[[n]] <- mydmadata.n
  my.data.n <- subset(my.data.n, !(zip %in% dmadata.n[[n]]$zip))
}
finaldata.n <- do.call(rbind, dmadata.n)
remainder.data.n <- subset(zip.data, !(zip %in% finaldata.n$zip))
remainder.data.n$distance <- 0
remainder.data.n$zone <- remainder.data.n$optimalsplit

########################
# Format results - North
########################

finaldata.n <- subset(finaldata.n, select=c('dmaname', 'zip', 'rank', 'dma' , 
                                            'city', 'pop', 'state', 'latitude', 'longitude', 
                                            'dmapop', 'optimalsplit', 'optimalpercent', 'zone', 'distance'))
remainder.data.n$distance <- 0
finaldata.n <- rbind(finaldata.n, remainder.data.n)
finaldata.n <- finaldata.n[!duplicated(finaldata.n), ]
finaldata.n$zone <- ifelse(nchar(finaldata.n$zone, type="chars")==1, paste0("0", 
                                                                            finaldata.n$zone), finaldata.n$zone)
finaldata.n$zoneid <- paste(finaldata.n$dma, finaldata.n$zone, sep="-")

#Calculate Zone Population
zonepop.n <- aggregate(pop ~ zoneid, data=finaldata.n, sum)
colnames(zonepop.n) <- c("zoneid", "zonepop")
finaldata.n <- merge(finaldata.n, zonepop.n, by="zoneid")

finaldata.n <- finaldata.n[ , !(names(finaldata.n) %in% c('city'))]
finaldata.n <- finaldata.n[!duplicated(finaldata.n), ]
finaldata.n$method <- "Auto"
finaldata.n$direction <- "North"

##################
# Create objects - East
##################	

zip.data.e <- zip.data				
remainder.data.e <- 0
zip.data.e$zone <- 0
zip.data.e$latrank <- 0

dma.data.e <- list()
dmadata.e <- list()
city.data.e <- list()

#####################
# Begin loop - East
#####################

for (n in 1:25) #by dma rank
{
  my.data.e <- subset(zip.data, rank==n) #subsets by DMA rank
  my.data.e$ziprank <- rank(-my.data.e$pop, ties.method = c("first"), 0)
  citydata.e <- my.data.e
  
  for (c in 1:citydata.e$optimalsplit) #by city rank
  {
    citydata.e <- citydata.e[order(-citydata.e$longitude),]
    citydata.e$latrank <- rank(-citydata.e$longitude, ties.method = c("first"), 0)
    northmost <- subset(citydata.e, latrank==1, select=c("dmaname", "zip", 
                                                         "longitude", "latitude"))
    colnames(northmost) <- c("dmaname", "nzip", "nlongitude", "nlatitude")
    allzips.e <- merge(northmost, citydata.e, by="dmaname")
    distance.e <- 0
    result.data.e <- list()
    
    for (i in 1:nrow(allzips.e)) #calculate distances
    {
      record.e <- subset(allzips.e, as.numeric(rownames(allzips.e)) == i)
      distance.e[i] <- distm(record.e[,c('nlongitude', 'nlatitude')],
                             record.e[,c('longitude', 'latitude')],
                             fun = distHaversine)
      record.e$distance <- distance.e[i]
      result.data.e[[i]] <- record.e
    }
    allzips.e <- do.call(rbind, result.data.e)
    allzips.e <- allzips.e[order(allzips.e$distance),]
    
    allzips.e$runpop <- cumsum(allzips.e$pop)
    allzips.e$runpercent <- allzips.e$runpop/allzips.e$dmapop
    allzips.e$zone <-	ifelse(allzips.e$runpercent < allzips.e$optimalpercent, c, 0)
    
    city.data.e[[c]] <- subset(allzips.e, zone > 0)
    citydata.e <- subset(citydata.e, !(zip %in% city.data.e[[c]]$zip))
  }
  mydmadata.e <- do.call(rbind, city.data.e)
  dmadata.e[[n]] <- mydmadata.e
  my.data.e <- subset(my.data.e, !(zip %in% dmadata.e[[n]]$zip))
}
finaldata.e <- do.call(rbind, dmadata.e)
remainder.data.e <- subset(zip.data, !(zip %in% finaldata.e$zip))
remainder.data.n$distance <- 0
remainder.data.e$zone <- remainder.data.e$optimalsplit

########################
# Format results - East
########################

finaldata.e <- subset(finaldata.e, select=c('dmaname', 'zip', 'rank', 'dma' , 
                                            'city', 'pop', 'state', 'latitude', 'longitude', 
                                            'dmapop', 'optimalsplit', 'optimalpercent', 'zone', 'distance'))
remainder.data.e$distance <- 0
finaldata.e <- rbind(finaldata.e, remainder.data.e)
finaldata.e <- finaldata.e[!duplicated(finaldata.e), ]
finaldata.e$zone <- ifelse(nchar(finaldata.e$zone, type="chars")==1, paste0("0", 
                                                                            finaldata.e$zone), finaldata.e$zone)
finaldata.e$zoneid <- paste(finaldata.e$dma, finaldata.e$zone, sep="-")

#Calculate Zone Population
zonepop.e <- aggregate(pop ~ zoneid, data=finaldata.e, sum)
colnames(zonepop.e) <- c("zoneid", "zonepop")
finaldata.e <- merge(finaldata.e, zonepop.e, by="zoneid")

finaldata.e <- finaldata.e[ , !(names(finaldata.e) %in% c('city'))]
finaldata.e <- finaldata.e[!duplicated(finaldata.e), ]

finaldata.e$method <- "Auto"
finaldata.e$direction <- "East"

########################
# Create objects - South
########################

zip.data.s <- zip.data				
remainder.data.s <- 0
zip.data.s$zone <- 0
zip.data.s$latrank <- 0

dma.data.s <- list()
dmadata.s <- list()
city.data.s <- list()

#####################
# Begin loop - South
#####################

for (n in 1:25) #by dma rank
{
  my.data.s <- subset(zip.data, rank==n) #subsets by DMA rank
  my.data.s$ziprank <- rank(-my.data.s$pop, ties.method = c("first"), 0)
  citydata.s <- my.data.s
  
  for (c in 1:citydata.s$optimalsplit) #by city rank
  {
    citydata.s <- citydata.s[order(citydata.s$latitude),]
    citydata.s$latrank <- rank(citydata.s$latitude, ties.method = c("first"), 0)
    northmost <- subset(citydata.s, latrank==1, select=c("dmaname", "zip", 
                                                         "longitude", "latitude"))
    colnames(northmost) <- c("dmaname", "nzip", "nlongitude", "nlatitude")
    allzips.s <- merge(northmost, citydata.s, by="dmaname")
    distance.s <- 0
    result.data.s <- list()
    
    for (i in 1:nrow(allzips.s)) #calculate distances
    {
      record.s <- subset(allzips.s, as.numeric(rownames(allzips.s)) == i)
      distance.s[i] <- distm(record.s[,c('nlongitude', 'nlatitude')],
                             record.s[,c('longitude', 'latitude')],
                             fun = distHaversine)
      record.s$distance <- distance.s[i]
      result.data.s[[i]] <- record.s
    }
    allzips.s <- do.call(rbind, result.data.s)
    allzips.s <- allzips.s[order(allzips.s$distance),]
    
    allzips.s$runpop <- cumsum(allzips.s$pop)
    allzips.s$runpercent <- allzips.s$runpop/allzips.s$dmapop
    allzips.s$zone <-	ifelse(allzips.s$runpercent < allzips.s$optimalpercent, c, 0)
    
    city.data.s[[c]] <- subset(allzips.s, zone > 0)
    citydata.s <- subset(citydata.s, !(zip %in% city.data.s[[c]]$zip))
  }
  mydmadata.s <- do.call(rbind, city.data.s)
  dmadata.s[[n]] <- mydmadata.s
  my.data.s <- subset(my.data.s, !(zip %in% dmadata.s[[n]]$zip))
}
finaldata.s <- do.call(rbind, dmadata.s)
remainder.data.s <- subset(zip.data, !(zip %in% finaldata.s$zip))
remainder.data.n$distance <- 0
remainder.data.s$zone <- remainder.data.s$optimalsplit

########################
# Format results - South
########################

finaldata.s <- subset(finaldata.s, select=c('dmaname', 'zip', 'rank', 'dma' , 
                                            'city', 'pop', 'state', 'latitude', 'longitude', 
                                            'dmapop', 'optimalsplit', 'optimalpercent', 'zone', 'distance'))
remainder.data.s$distance <- 0
finaldata.s <- rbind(finaldata.s, remainder.data.s)
finaldata.s <- finaldata.s[!duplicated(finaldata.s), ]
finaldata.s$zone <- ifelse(nchar(finaldata.s$zone, type="chars")==1, paste0("0", 
                                                                            finaldata.s$zone), finaldata.s$zone)
finaldata.s$zoneid <- paste(finaldata.s$dma, finaldata.s$zone, sep="-")

#Calculate Zone Population
zonepop.s <- aggregate(pop ~ zoneid, data=finaldata.s, sum)
colnames(zonepop.s) <- c("zoneid", "zonepop")
finaldata.s <- merge(finaldata.s, zonepop.s, by="zoneid")

finaldata.s <- finaldata.s[ , !(names(finaldata.s) %in% c('city'))]
finaldata.s <- finaldata.s[!duplicated(finaldata.s), ]

finaldata.s$method <- "Auto"
finaldata.s$direction <- "South"

#######################
# Create objects - West
#######################

zip.data.w <- zip.data					
remainder.data.w <- 0
zip.data.w$zone <- 0
zip.data.w$latrank <- 0

dma.data.w <- list()
dmadata.w <- list()
city.data.w <- list()

#####################
# Begin loop - West
#####################

for (n in 1:25) #by dma rank
{
  my.data.w <- subset(zip.data, rank==n) #subsets by DMA rank
  my.data.w$ziprank <- rank(-my.data.w$pop, ties.method = c("first"), 0)
  citydata.w <- my.data.w
  
  for (c in 1:citydata.w$optimalsplit) #by city rank
  {
    citydata.w <- citydata.w[order(citydata.w$longitude),]
    citydata.w$latrank <- rank(citydata.w$longitude, ties.method = c("first"), 0)
    northmost <- subset(citydata.w, latrank==1, select=c("dmaname", "zip", 
                                                         "longitude", "latitude"))
    colnames(northmost) <- c("dmaname", "nzip", "nlongitude", "nlatitude")
    allzips.w <- merge(northmost, citydata.w, by="dmaname")
    distance.w <- 0
    result.data.w <- list()
    
    for (i in 1:nrow(allzips.w)) #calculate distances
    {
      record.w <- subset(allzips.w, as.numeric(rownames(allzips.w)) == i)
      distance.w[i] <- distm(record.w[,c('nlongitude', 'nlatitude')],
                             record.w[,c('longitude', 'latitude')],
                             fun = distHaversine)
      record.w$distance <- distance.w[i]
      result.data.w[[i]] <- record.w
    }
    allzips.w <- do.call(rbind, result.data.w)
    allzips.w <- allzips.w[order(allzips.w$distance),]
    
    allzips.w$runpop <- cumsum(allzips.w$pop)
    allzips.w$runpercent <- allzips.w$runpop/allzips.w$dmapop
    allzips.w$zone <-	ifelse(allzips.w$runpercent < allzips.w$optimalpercent, c, 0)
    
    city.data.w[[c]] <- subset(allzips.w, zone > 0)
    citydata.w <- subset(citydata.w, !(zip %in% city.data.w[[c]]$zip))
  }
  mydmadata.w <- do.call(rbind, city.data.w)
  dmadata.w[[n]] <- mydmadata.w
  my.data.w <- subset(my.data.w, !(zip %in% dmadata.w[[n]]$zip))
}
finaldata.w <- do.call(rbind, dmadata.w)
remainder.data.w <- subset(zip.data, !(zip %in% finaldata.w$zip))
remainder.data.n$distance <- 0
remainder.data.w$zone <- remainder.data.w$optimalsplit

########################
# Format results - West
########################

finaldata.w <- subset(finaldata.w, select=c('dmaname', 'zip', 'rank', 'dma' , 
                                            'city', 'pop', 'state', 'latitude', 'longitude', 
                                            'dmapop', 'optimalsplit', 'optimalpercent', 'zone', 'distance'))
remainder.data.w$distance <- 0
finaldata.w <- rbind(finaldata.w, remainder.data.w)
finaldata.w <- finaldata.w[!duplicated(finaldata.w), ]
finaldata.w$zone <- ifelse(nchar(finaldata.w$zone, type="chars")==1, paste0("0", 
                                                                            finaldata.w$zone), finaldata.w$zone)
finaldata.w$zoneid <- paste(finaldata.w$dma, finaldata.w$zone, sep="-")

#Calculate Zone Population
zonepop.w <- aggregate(pop ~ zoneid, data=finaldata.w, sum)
colnames(zonepop.w) <- c("zoneid", "zonepop")
finaldata.w <- merge(finaldata.w, zonepop.w, by="zoneid")

finaldata.w <- finaldata.w[ , !(names(finaldata.w) %in% c('city'))]
finaldata.w <- finaldata.w[!duplicated(finaldata.w), ]

finaldata.w$method <- "Auto"
finaldata.w$direction <- "West"

#############################
# Manually Assign
#############################

finaldata.m <- rbind(subset(finaldata.w, rank==1), subset(finaldata.s, rank==2),
                     subset(finaldata.e, rank==3), subset(finaldata.w, rank==4),
                     subset(finaldata.n, rank==5), subset(finaldata.e, rank==6),
                     subset(finaldata.n, rank==7), subset(finaldata.n, rank==8),
                     subset(finaldata.s, rank==9), subset(finaldata.n, rank==10),
                     subset(finaldata.n, rank==11), subset(finaldata.s, rank==12),
                     subset(finaldata.s, rank==13), subset(finaldata.s, rank==14),
                     subset(finaldata.s, rank==15), subset(finaldata.s, rank==16),
                     subset(finaldata.s, rank==17), subset(finaldata.e, rank==18),
                     subset(finaldata.n, rank==19), subset(finaldata.s, rank==20),
                     subset(finaldata.n, rank==21), subset(finaldata.w, rank==22),
                     subset(finaldata.s, rank==23), subset(finaldata.n, rank==24),
                     subset(finaldata.n, rank==25))
finaldata.m$method <- "Manual"

# NY 1
zips <- c("06907", "06906", "06905", "06903", "06902", "06901", "06897", "06890", 
          "06883", "06880", "06856", "06855", "06854", "06853", "06851", "06850", 
          "06840", "06825", "06824", "06820", "06615", "06614", "06612", "06611", 
          "06610", "06608", "06607", "06606", "06605", "06604", "06484", "06482", 
          "06468")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "501-03"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "03"

# NY 2
zips <- c("12483")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "501-03"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "03"

# CHI 1
zips <- c("47964", "47963", "47951", "47922", "46349", "47943")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "602-01"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "01"

# CHI 2
zips <- c("60403")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "602-02"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "02"

# CHI 3
zips <- c("60431")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "602-06"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "06"

# PHI 1
zips <- c("08327")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "504-06"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "06"

# TSP 1
zips <- c("33857", "33876")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "539-01"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "01"

# TSP 1
zips <- c("33178")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "528-02"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "02"

# AZ 1
zips <- c("85392")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "753-02"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "02"

# CO 1
zips <- c("81325")
zips <- data.frame(zips)
finaldata.m[finaldata.m$zip %in% zips$zips, 1] <- "751-01"
finaldata.m[finaldata.m$zip %in% zips$zips, 13] <- "01"

finaldata <- rbind(finaldata.n, finaldata.e, finaldata.s, finaldata.w, finaldata.m)

write.table(finaldata, file = "//Users//psullivan01//Documents//DMAZones9.csv",
            row.names=FALSE, sep=",")
