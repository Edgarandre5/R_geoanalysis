set.seed(1111)
setwd("/home/")

library(data.table)


species <- data.table(list.files(path="/home/species_presences/",pattern =".csv",full.name=F))
species_path <- data.table(list.files(path="/home/species_presences/",pattern =".csv",full.name=T))

for (j in 1:nrow(species)) {
  a1 <- data.table(read.csv(paste0(species_path[j]),header = T))
  write.csv(a1[sample(nrow(a1),300),] ,paste0(species[j]),row.names = F)
}
 
