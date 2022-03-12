set.seed(1111)
setwd("/Users/edgarsanchez/Mamiferos/poligonos_a_puntos/puntos/")

library(rgdal)
library(raster)
library(maptools) 
library(rgeos)
library(data.table)
library(parallel)
library(filesstrings)
resultado <- setNames(data.table(matrix(data=character(),nrow = 10000, ncol = 6)), c("scientificName","0.083", "0.20","0.41","0.62","0.83"))

archivos <- data.table(list.files(path="/Users/edgarsanchez/Mamiferos/poligonos_a_puntos/shapes/end/",pattern =".csv",full.name=F))
a083 <- data.table(list.files(path="/Users/edgarsanchez/Mamiferos/poligonos_a_puntos/shapes/end/",pattern =".csv",full.name=T))

for (j in 1:nrow(archivos)) {
  a1 <- data.table(read.csv(paste0(a083[j]),header = T))
  a2 <- a1[sample(nrow(a1),300),]#round(300*(1-exp(-nrow(a1)/1873.32)))
  write.csv(a2,paste0(archivos[j]),row.names = F)
}
 
