set.seed(1111)
setwd("/home/edgarandres/modelos/especies_2/")

library(rgdal)
library(raster)
library(maptools) 
library(rgeos)
library(data.table)
library(parallel)
library(filesstrings)
resultado <- setNames(data.table(matrix(data=character(),nrow = 10000, ncol = 2)), c("scientificName","10km"))#, "0.20","0.41","0.62","0.83" "1km","2.5km","5km",

archivos <- data.table(list.files(path="/home/edgarandres/modelos/especies_1/",pattern =".csv",full.name=F))
a083 <- data.table(list.files(path="/home/edgarandres/modelos/especies_1/",pattern =".csv",full.name=T))

for (j in 1:nrow(archivos)) {
  a <- read.csv(paste0(a083[j]),header = T)
  f <- nrow(a)-300
  for (k in 1:(f)) {  
    e <- data.table(matrix(data=double(),nrow = nrow(a)-1, ncol = nrow(a)))
    z <- data.table(matrix(data=double(),nrow = nrow(a), ncol = 2))
      for (i in 1:nrow(a)) {
      b <- cbind(a[i,1],a[i,2])
      c <- cbind(a[,1],a[,2])
      c <- c[-i,]
      e[[i]] <- pointDistance(b,c,lonlat = TRUE)
    } 
    for (l in 1:nrow(a)) {
      z[l,1] <- min(e[[l]])
      z[l,2] <- which.min(e[[l]])
    }
    if(sum(e[[as.numeric(z[which.min(z[[1]]),2])]]) > sum(e[[which.min(z[[1]])]])){
      a <- a[-which.min(z[[1]]),]  
    } else {
      a <- a[-as.numeric(z[which.min(z[[1]]),2]),]
    }
  }
  write.csv(a, paste0(archivos[j]),row.names = F)
}
