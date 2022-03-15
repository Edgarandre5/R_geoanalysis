setwd("/home/")

library(rgdal)
library(raster)
library(rgeos)
library(dplyr)
library(data.table)


anp_s <- shapefile("/home/anp_mexico.shp")
hurricanes <- shapefile("/home/hurricanes.shp")
fires <- shapefile("/home/fires.shp")
vegetationp <- shapefile("/home/vegetation.shp")
u30 <- shapefile("/home/landuse_without_30.shp")
u50 <- shapefile("/home/landuse_without_50.shp")
u70 <- shapefile("/home/landuse_without_70.shp")
uc30 <- shapefile("/home/landuse_with_30.shp")
uc50 <- shapefile("/home/landuse_with_50.shp")
uc70 <- shapefile("/home/landuse_with_70.shp")
b1_30 <- raster("/home/b1_585_30.tif")
b1_50 <- raster("/home/b1_585_50.tif")
b1_70 <- raster("/home/b1_585_70.tif")
b12_30 <- raster("/home/b12_585_30.tif")
b12_50 <- raster("/home/b12_585_50.tif")
b12_70 <- raster("/home/b12_585_70.tif")


species <- data.table(list.files(path = "/home/present", full.names = TRUE, recursive = FALSE))
species_names <- data.table(list.files(path = "/home/present", full.names = FALSE, recursive = FALSE))
species1 <- data.table(list.files(path = "/home/ssp245_BCC_2030_mx", full.names = TRUE, recursive = FALSE))
species_names1 <- data.table(list.files(path = "/home/ssp245_BCC_2030_mx", full.names = FALSE, recursive = FALSE))
species2 <- data.table(list.files(path = "/home/ssp245_BCC_2050_mx", full.names = TRUE, recursive = FALSE))
species_names2 <- data.table(list.files(path = "/home/ssp245_BCC_2050_mx", full.names = FALSE, recursive = FALSE))
species3 <- data.table(list.files(path = "/home/ssp245_BCC_2070_mx", full.names = TRUE, recursive = FALSE))
species_names3 <- data.table(list.files(path = "/home/ssp245_BCC_2070_mx", full.names = FALSE, recursive = FALSE))

result <- setNames(data.table(matrix(data=character(),nrow = nrow(especies), ncol = 45)), c("scientificName","presente","ssp245_2030","ssp245_2050","ssp245_2070","Change_ssp245_2030","Change_ssp245_2050","Change_ssp245_2070","Anp_overlap_present","Anp_overlap_ssp245_30","Anp_overlap_ssp245_50","Anp_overlap_ssp245_70","Anp_overlap_change_ssp245_2030","Anp_overlap_change_ssp245_2050","Anp_overlap_change_ssp245_2070","fires_present","fires_ssp245_2030","fires_ssp245_2050","fires_ssp245_2070","hurricanes_present","hurricanes_ssp245_2030","hurricanes_ssp245_2050","hurricanes_ssp245_2070","U30_overlap_ssp245","U50_overlap_ssp245","U70_overlap_ssp245","UC30_overlap_ssp245","UC50_overlap_ssp245","UC70_overlap_ssp245","T_VP_BOSQUE_DE_CONIFERAS_Y_ENCINOS","T_VP_BOSQUE_ESPINOSO","T_VP_BOSQUE_MESOFILO_DE_MONTANA","T_VP_BOSQUE_TROPICAL_CADUCIFOLIO","T_VP_BOSQUE_TROPICAL_PERENNIFOLIO","T_VP_BOSQUE_TROPICAL_SUBCADUCIFOLIO","T_VP_CUERPOS_DE_AGUA","T_VP_MATORRAL_XEROFILO","T_VP_PASTIZAL","T_VP_vegetation_ACUATICA_Y_SUBACUATICA","Average_increase_Bio1_585_30","Average_increase_Bio1_585_50","Average_increase_Bio1_585_70","Average_increase_Bio12_585_30","Average_increase_Bio12_585_50","Average_increase_Bio12_585_70"))

for (i in 1:nrow(species)) {
  
  result[i,1] <- gsub(".tif","",gsub("presente_","",paste0(species_names[i]))) 
  t1n <- raster(paste0(especies[i]))
  ones <- which(extract(t1n == 1,fires) == 1)
  result[i,16] <- nrow(fires[ones,]) 
  t1n[is.na(t1n)] <- 0
  rm(f)
  f <- freq(t1n)
  result[i,2] <- f[2,2]
  t2n <- mask(t1n,anp_s)
  b <- 0
  for (l in 1:139){
    t3n <- mask(t1n,subset(hurricanes, Intensidad==paste0(hurricanes$Intensidad[l])))
    if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
      b <- b + 0
    } else {
      rm(h)
      h <- freq(t3n)
      b <- b + (h[which.max(h[,1]),2]*as.numeric(hurricanes$Intensidad[l])/as.numeric(result[i,2]))
    }
  }
  result[i,20] <- b
  if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
    result[i,9] <- 0
  } else {
    rm(g)
    g <- freq(t2n)
    result[i,9] <- g[which.max(g[,1]),2]
  }
  #############vegetation_POTENCIAL#####################
  for (l in 1:10){
    v <- mask(t1n,subset(vegetationp, TIPOS==paste0(vegetationp$TIPOS[l])))
    if (v@data@min == v@data@max | maxValue(v)==0 | is.na(maxValue(v))){
      result[i,29+l] <- 0 
    } else {
      rm(fv)
      fv <- freq(v)
      result[i,29+l] <- fv[which.max(fv[,1]),2]
    }
  }
  
  
  ######################################2030######################################
  t1n <- raster(paste0(species1[i]))
  ones <- which(extract(t1n == 1,fires) == 1)
  result[i,17] <- nrow(fires[ones,]) 
  t1n[is.na(t1n)] <- 0
  if (t1n@data@min == t1n@data@max | maxValue(t1n)==0 | is.na(maxValue(t1n))){
    result[i,3] <- 0
    result[i,10] <- 0
    result[i,21] <- 0
    result[i,24] <- 0
    result[i,27] <- 0
    result[i,40] <- 0
    result[i,43] <- 0
  } else {
    t4n <- mask(t1n,u30)###cambio_uso_suelo_s
    if (t4n@data@min == t4n@data@max | maxValue(t4n)==0 | is.na(maxValue(t4n))){
      result[i,24] <- 0
    } else {
      rm(fu)
      fu <- freq(t4n)
      result[i,24] <- fu[which.max(fu[,1]),2]
    }
    t5n <- mask(t1n,uc30)###cambio_uso_suelo_c
    if (t5n@data@min == t5n@data@max | maxValue(t5n)==0 | is.na(maxValue(t5n))){
      result[i,27] <- 0
    } else {
      rm(gu)
      gu <- freq(t5n)
      result[i,27] <- gu[which.max(gu[,1]),2]
    }
    rm(f)
    f <- freq(t1n)
    result[i,3] <- f[2,2]
    ############Promedio_resta_bio1
    at <- crop(b1_30,t1n)
    ax <- overlay(t1n, at, fun=function(x,y){x*y})
    ax[ax==0] <- NA
    ay <- cellStats(ax,stat="mean")
    result[i,40] <- ay
    ############Promedio_resta_bio12
    bt <- crop(b12_30,t1n)
    bx <- overlay(t1n, bt, fun=function(x,y){x*y})
    bx[bx==0] <- NA
    by <- cellStats(bx,stat="mean")
    result[i,43] <- by
    ####anp
    t2n <- mask(t1n,anp_s)
    if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
      result[i,10] <- 0
    } else {
      rm(g)
      g <- freq(t2n)
      result[i,10] <- g[2,2]
    }
    b <- 0
    for (l in 1:139){
      t3n <- mask(t1n,subset(hurricanes, Intensidad==paste0(hurricanes$Intensidad[l])))
      if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
        b <- b + 0
      } else {
        rm(h)
        h <- freq(t3n)
        b <- b + (h[which.max(h[,1]),2]*as.numeric(hurricanes$Intensidad[l])/as.numeric(result[i,3]))
      }
    }
    result[i,21] <- b
  }
  ######################################2050######################################
  t1n <- raster(paste0(species2[i]))
  ones <- which(extract(t1n == 1,fires) == 1)
  result[i,18] <- nrow(fires[ones,]) 
  t1n[is.na(t1n)] <- 0
  if (t1n@data@min == t1n@data@max | maxValue(t1n)==0 | is.na(maxValue(t1n))){
    result[i,4] <- 0
    result[i,11] <- 0
    result[i,22] <- 0
    result[i,25] <- 0
    result[i,41] <- 0
    result[i,44] <- 0
  } else {
    t4n <- mask(t1n,u50)###cambio_uso_suelo_s
    if (t4n@data@min == t4n@data@max | maxValue(t4n)==0 | is.na(maxValue(t4n))){
      result[i,25] <- 0
    } else {
      rm(fu)
      fu <- freq(t4n)
      result[i,25] <- fu[which.max(fu[,1]),2]
    }
    t5n <- mask(t1n,uc50)###cambio_uso_suelo_c
    if (t5n@data@min == t5n@data@max | maxValue(t5n)==0 | is.na(maxValue(t5n))){
      result[i,28] <- 0
    } else {
      rm(gu)
      gu <- freq(t5n)
      result[i,28] <- gu[which.max(gu[,1]),2]
    }
    rm(f)
    f <- freq(t1n)
    result[i,4] <- f[2,2]
    ############Promedio_resta_bio1
    at <- crop(b1_50,t1n)
    ax <- overlay(t1n, at, fun=function(x,y){x*y})
    ax[ax==0] <- NA
    ay <- cellStats(ax,stat="mean")
    result[i,41] <- ay
    ############Promedio_resta_bio12
    bt <- crop(b12_50,t1n)
    bx <- overlay(t1n, bt, fun=function(x,y){x*y})
    bx[bx==0] <- NA
    by <- cellStats(bx,stat="mean")
    result[i,44] <- by
    ####anp
    t2n <- mask(t1n,anp_s)
    if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
      result[i,11] <- 0
    } else {
      rm(g)
      g <- freq(t2n)
      result[i,11] <- g[2,2]
    }            
    b <- 0
    for (l in 1:139){
      t3n <- mask(t1n,subset(hurricanes, Intensidad==paste0(hurricanes$Intensidad[l])))
      if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
        b <- b + 0
      } else {
        rm(h)
        h <- freq(t3n)
        b <- b + (h[which.max(h[,1]),2]*as.numeric(hurricanes$Intensidad[l])/as.numeric(result[i,4]))
      }
    }
    result[i,22] <- b
  }
  ######################################2070######################################
  t1n <- raster(paste0(species3[i]))
  ones <- which(extract(t1n == 1,fires) == 1)
  result[i,19] <- nrow(fires[ones,]) 
  t1n[is.na(t1n)] <- 0
  if (t1n@data@min == t1n@data@max | maxValue(t1n)==0 | is.na(maxValue(t1n))){
    result[i,5] <- 0
    result[i,12] <- 0
    result[i,23] <- 0
    result[i,26] <- 0
    result[i,42] <- 0
    result[i,45] <- 0
  } else {
    t4n <- mask(t1n,u70)###cambio_uso_suelo
    if (t4n@data@min == t4n@data@max | maxValue(t4n)==0 | is.na(maxValue(t4n))){
      result[i,26] <- 0
    } else {
      rm(fu)
      fu <- freq(t4n)
      result[i,26] <- fu[which.max(fu[,1]),2]
    }
    t5n <- mask(t1n,uc70)###cambio_uso_suelo_c
    if (t5n@data@min == t5n@data@max | maxValue(t5n)==0 | is.na(maxValue(t5n))){
      result[i,29] <- 0
    } else {
      rm(gu)
      gu <- freq(t5n)
      result[i,29] <- gu[which.max(gu[,1]),2]
    }
    rm(f)
    f <- freq(t1n)
    result[i,5] <- f[2,2]
    ############Promedio_resta_bio1
    at <- crop(b1_70,t1n)
    ax <- overlay(t1n, at, fun=function(x,y){x*y})
    ax[ax==0] <- NA
    ay <- cellStats(ax,stat="mean")
    result[i,42] <- ay
    ############Promedio_resta_bio12
    bt <- crop(b12_70,t1n)
    bx <- overlay(t1n, bt, fun=function(x,y){x*y})
    bx[bx==0] <- NA
    by <- cellStats(bx,stat="mean")
    result[i,45] <- by
    ####anp
    t2n <- mask(t1n,anp_s)
    if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
      result[i,12] <- 0
    } else {
      rm(g)
      g <- freq(t2n)
      result[i,12] <- g[2,2]
    }   
    b <- 0
    for (l in 1:139){
      t3n <- mask(t1n,subset(hurricanes, Intensidad==paste0(hurricanes$Intensidad[l])))
      if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
        b <- b + 0
      } else {
        rm(h)
        h <- freq(t3n)
        b <- b + (h[which.max(h[,1]),2]*as.numeric(hurricanes$Intensidad[l])/as.numeric(result[i,5]))
      }
    }
    result[i,23] <- b
  }
  
}


write.csv(result,"/home/traslapes_245BCC.csv", row.names = TRUE,col.names = TRUE)
