setwd("/home/")

library(raster)
library(data.table)

species <- data.table(list.files(path="/home/especies/",pattern =".csv",full.name=F))
species_path <- data.table(list.files(path="/home/especies/",pattern =".csv",full.name=T))

for (j in 1:nrow(species)) {
  a <- read.csv(paste0(species_path[j]),header = T)
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
  write.csv(a, paste0(species[j]),row.names = F)
}
