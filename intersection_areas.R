setwd("/home/edgarandres/modelos/")

library(rgdal)
library(raster)
library(maptools) 
library(rgeos)
library(dplyr)
library(data.table)
library(parallel)
library(filesstrings)
#library(fasterize)

anp_s <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/anp_mexico.shp")
huracanes <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/huracanes_21marzo.shp")
incendios <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/incendios_final.shp")
vegetacionp <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/vegetacion_potencial.shp")
u30 <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/uso_sin_30.shp")
u50 <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/uso_sin_50.shp")
u70 <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/uso_sin_70.shp")
uc30 <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/uso_con_30.shp")
uc50 <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/uso_con_50.shp")
uc70 <- shapefile("/home/edgarandres/modelos/variables/capas_traslapes/uso_con_70.shp")
b1_30 <- raster("/home/edgarandres/modelos/variables/capas_traslapes/b1_585_30.tif")
b1_50 <- raster("/home/edgarandres/modelos/variables/capas_traslapes/b1_585_50.tif")
b1_70 <- raster("/home/edgarandres/modelos/variables/capas_traslapes/b1_585_70.tif")
b12_30 <- raster("/home/edgarandres/modelos/variables/capas_traslapes/b12_585_30.tif")
b12_50 <- raster("/home/edgarandres/modelos/variables/capas_traslapes/b12_585_50.tif")
b12_70 <- raster("/home/edgarandres/modelos/variables/capas_traslapes/b12_585_70.tif")


especies <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/presente_mx", full.names = TRUE, recursive = FALSE))
nombres <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/presente_mx", full.names = FALSE, recursive = FALSE))
especies1 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2030_mx", full.names = TRUE, recursive = FALSE))
nombres1 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2030_mx", full.names = FALSE, recursive = FALSE))
especies2 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2050_mx", full.names = TRUE, recursive = FALSE))
nombres2 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2050_mx", full.names = FALSE, recursive = FALSE))
especies3 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2070_mx", full.names = TRUE, recursive = FALSE))
nombres3 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2070_mx", full.names = FALSE, recursive = FALSE))

# especies1 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2030", full.names = TRUE, recursive = FALSE))
# nombres1 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2030", full.names = FALSE, recursive = FALSE))
# especies2 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2050", full.names = TRUE, recursive = FALSE))
# nombres2 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2050", full.names = FALSE, recursive = FALSE))
# especies3 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2070", full.names = TRUE, recursive = FALSE))
# nombres3 <- data.table(list.files(path = "/home/edgarandres/modelos/resultados/ssp245_BCC_2070", full.names = FALSE, recursive = FALSE))

resultado <- setNames(data.table(matrix(data=character(),nrow = nrow(especies), ncol = 45)), c("scientificName","presente","ssp585_2030","ssp585_2050","ssp585_2070","Cambio_porcentual_ssp585_2030","Cambio_porcentual_ssp585_2050","Cambio_porcentual_ssp585_2070","Traslape_anp_presente","Traslape_anp_presente_ssp585_30","Traslape_anp_presente_ssp585_50","Traslape_anp_presente_ssp585_70","Cambio_porcentual_traslape_anp_ssp585_2030","Cambio_porcentual_traslape_anp_ssp585_2050","Cambio_porcentual_traslape_anp_ssp585_2070","Incendios_presente","Incendios_ssp585_2030","Incendios_ssp585_2050","Incendios_ssp585_2070","Huracanes_presente","Huracanes_ssp585_2030","Huracanes_ssp585_2050","Huracanes_ssp585_2070","Traslape_cus_spp_ssp585_2030","Traslape_cus_spp_ssp585_2050","Traslape_cus_spp_ssp585_2070","Traslape_cus_cpp_ssp585_2030","Traslape_cus_cpp_ssp585_2050","Traslape_cus_cpp_ssp585_2070","T_VP_BOSQUE_DE_CONIFERAS_Y_ENCINOS","T_VP_BOSQUE_ESPINOSO","T_VP_BOSQUE_MESOFILO_DE_MONTANA","T_VP_BOSQUE_TROPICAL_CADUCIFOLIO","T_VP_BOSQUE_TROPICAL_PERENNIFOLIO","T_VP_BOSQUE_TROPICAL_SUBCADUCIFOLIO","T_VP_CUERPOS_DE_AGUA","T_VP_MATORRAL_XEROFILO","T_VP_PASTIZAL","T_VP_VEGETACION_ACUATICA_Y_SUBACUATICA","Aumento_promedio_Bio1_585_30","Aumento_promedio_Bio1_585_50","Aumento_promedio_Bio1_585_70","Aumento_promedio_Bio12_585_30","Aumento_promedio_Bio12_585_50","Aumento_promedio_Bio12_585_70"))#, "0.20","0.41","0.62","0.83" "1km","2.5km","5km",
resultado_t <- setNames(data.table(matrix(data=character(),nrow = nrow(especies), ncol = 45)), c("scientificName","presente","ssp585_2030","ssp585_2050","ssp585_2070","Cambio_porcentual_ssp585_2030","Cambio_porcentual_ssp585_2050","Cambio_porcentual_ssp585_2070","Traslape_anp_presente","Traslape_anp_presente_ssp585_30","Traslape_anp_presente_ssp585_50","Traslape_anp_presente_ssp585_70","Cambio_porcentual_traslape_anp_ssp585_2030","Cambio_porcentual_traslape_anp_ssp585_2050","Cambio_porcentual_traslape_anp_ssp585_2070","Incendios_presente","Incendios_ssp585_2030","Incendios_ssp585_2050","Incendios_ssp585_2070","Huracanes_presente","Huracanes_ssp585_2030","Huracanes_ssp585_2050","Huracanes_ssp585_2070","Traslape_cus_spp_ssp585_2030","Traslape_cus_spp_ssp585_2050","Traslape_cus_spp_ssp585_2070","Traslape_cus_cpp_ssp585_2030","Traslape_cus_cpp_ssp585_2050","Traslape_cus_cpp_ssp585_2070","T_VP_BOSQUE_DE_CONIFERAS_Y_ENCINOS","T_VP_BOSQUE_ESPINOSO","T_VP_BOSQUE_MESOFILO_DE_MONTANA","T_VP_BOSQUE_TROPICAL_CADUCIFOLIO","T_VP_BOSQUE_TROPICAL_PERENNIFOLIO","T_VP_BOSQUE_TROPICAL_SUBCADUCIFOLIO","T_VP_CUERPOS_DE_AGUA","T_VP_MATORRAL_XEROFILO","T_VP_PASTIZAL","T_VP_VEGETACION_ACUATICA_Y_SUBACUATICA","Aumento_promedio_Bio1_585_30","Aumento_promedio_Bio1_585_50","Aumento_promedio_Bio1_585_70","Aumento_promedio_Bio12_585_30","Aumento_promedio_Bio12_585_50","Aumento_promedio_Bio12_585_70"))#, "0.20","0.41","0.62","0.83" "1km","2.5km","5km",

for (i in 1:100) {
  
  resultado[i,1] <- gsub(".tif","",gsub("presente_","",paste0(nombres[i]))) 
  t1n <- raster(paste0(especies[i]))
  ones <- which(extract(t1n == 1,incendios) == 1)
  resultado[i,16] <- nrow(incendios[ones,]) 
  t1n[is.na(t1n)] <- 0
  rm(f)
  f <- freq(t1n)
  resultado[i,2] <- f[2,2]
  t2n <- mask(t1n,anp_s)
  b <- 0
  for (l in 1:139){
    t3n <- mask(t1n,subset(huracanes, Intensidad==paste0(huracanes$Intensidad[l])))
    if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
      b <- b + 0
    } else {
      rm(h)
      h <- freq(t3n)
      b <- b + (h[which.max(h[,1]),2]*as.numeric(huracanes$Intensidad[l])/as.numeric(resultado[i,2]))
    }
  }
  resultado[i,20] <- b
  if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
    resultado[i,9] <- 0
  } else {
    rm(g)
    g <- freq(t2n)
    resultado[i,9] <- g[which.max(g[,1]),2]
  }
  #############VEGETACION_POTENCIAL#####################
  for (l in 1:10){
    v <- mask(t1n,subset(vegetacionp, TIPOS==paste0(vegetacionp$TIPOS[l])))
    if (v@data@min == v@data@max | maxValue(v)==0 | is.na(maxValue(v))){
      resultado[i,29+l] <- 0 
    } else {
      rm(fv)
      fv <- freq(v)
      resultado[i,29+l] <- fv[which.max(fv[,1]),2]
    }
  }
  

######################################2030######################################
t1n <- raster(paste0(especies1[i]))
ones <- which(extract(t1n == 1,incendios) == 1)
resultado[i,17] <- nrow(incendios[ones,]) 
t1n[is.na(t1n)] <- 0
if (t1n@data@min == t1n@data@max | maxValue(t1n)==0 | is.na(maxValue(t1n))){
  resultado[i,3] <- 0
  resultado[i,10] <- 0
  resultado[i,21] <- 0
  resultado[i,24] <- 0
  resultado[i,27] <- 0
  resultado[i,40] <- 0
  resultado[i,43] <- 0
} else {
  t4n <- mask(t1n,u30)###cambio_uso_suelo_s
  if (t4n@data@min == t4n@data@max | maxValue(t4n)==0 | is.na(maxValue(t4n))){
    resultado[i,24] <- 0
  } else {
    rm(fu)
    fu <- freq(t4n)
    resultado[i,24] <- fu[which.max(fu[,1]),2]
  }
  t5n <- mask(t1n,uc30)###cambio_uso_suelo_c
  if (t5n@data@min == t5n@data@max | maxValue(t5n)==0 | is.na(maxValue(t5n))){
    resultado[i,27] <- 0
  } else {
    rm(gu)
    gu <- freq(t5n)
    resultado[i,27] <- gu[which.max(gu[,1]),2]
  }
  rm(f)
  f <- freq(t1n)
  resultado[i,3] <- f[2,2]
  ############Promedio_resta_bio1
  at <- crop(b1_30,t1n)
  ax <- overlay(t1n, at, fun=function(x,y){x*y})
  ax[ax==0] <- NA
  ay <- cellStats(ax,stat="mean")
  resultado[i,40] <- ay
  ############Promedio_resta_bio12
  bt <- crop(b12_30,t1n)
  bx <- overlay(t1n, bt, fun=function(x,y){x*y})
  bx[bx==0] <- NA
  by <- cellStats(bx,stat="mean")
  resultado[i,43] <- by
  ####anp
  t2n <- mask(t1n,anp_s)
  if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
    resultado[i,10] <- 0
  } else {
    rm(g)
    g <- freq(t2n)
    resultado[i,10] <- g[2,2]
  }
  b <- 0
  for (l in 1:139){
    t3n <- mask(t1n,subset(huracanes, Intensidad==paste0(huracanes$Intensidad[l])))
    if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
      b <- b + 0
    } else {
      rm(h)
      h <- freq(t3n)
      b <- b + (h[which.max(h[,1]),2]*as.numeric(huracanes$Intensidad[l])/as.numeric(resultado[i,3]))
    }
  }
  resultado[i,21] <- b
}
######################################2050######################################
t1n <- raster(paste0(especies2[i]))
ones <- which(extract(t1n == 1,incendios) == 1)
resultado[i,18] <- nrow(incendios[ones,]) 
t1n[is.na(t1n)] <- 0
if (t1n@data@min == t1n@data@max | maxValue(t1n)==0 | is.na(maxValue(t1n))){
  resultado[i,4] <- 0
  resultado[i,11] <- 0
  resultado[i,22] <- 0
  resultado[i,25] <- 0
  resultado[i,41] <- 0
  resultado[i,44] <- 0
} else {
  t4n <- mask(t1n,u50)###cambio_uso_suelo_s
  if (t4n@data@min == t4n@data@max | maxValue(t4n)==0 | is.na(maxValue(t4n))){
    resultado[i,25] <- 0
  } else {
    rm(fu)
    fu <- freq(t4n)
    resultado[i,25] <- fu[which.max(fu[,1]),2]
  }
  t5n <- mask(t1n,uc50)###cambio_uso_suelo_c
  if (t5n@data@min == t5n@data@max | maxValue(t5n)==0 | is.na(maxValue(t5n))){
    resultado[i,28] <- 0
  } else {
    rm(gu)
    gu <- freq(t5n)
    resultado[i,28] <- gu[which.max(gu[,1]),2]
  }
  rm(f)
  f <- freq(t1n)
  resultado[i,4] <- f[2,2]
  ############Promedio_resta_bio1
  at <- crop(b1_50,t1n)
  ax <- overlay(t1n, at, fun=function(x,y){x*y})
  ax[ax==0] <- NA
  ay <- cellStats(ax,stat="mean")
  resultado[i,41] <- ay
  ############Promedio_resta_bio12
  bt <- crop(b12_50,t1n)
  bx <- overlay(t1n, bt, fun=function(x,y){x*y})
  bx[bx==0] <- NA
  by <- cellStats(bx,stat="mean")
  resultado[i,44] <- by
  ####anp
  t2n <- mask(t1n,anp_s)
  if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
    resultado[i,11] <- 0
  } else {
    rm(g)
    g <- freq(t2n)
    resultado[i,11] <- g[2,2]
  }            
  b <- 0
  for (l in 1:139){
    t3n <- mask(t1n,subset(huracanes, Intensidad==paste0(huracanes$Intensidad[l])))
    if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
      b <- b + 0
    } else {
      rm(h)
      h <- freq(t3n)
      b <- b + (h[which.max(h[,1]),2]*as.numeric(huracanes$Intensidad[l])/as.numeric(resultado[i,4]))
    }
  }
  resultado[i,22] <- b
}
######################################2070######################################
t1n <- raster(paste0(especies3[i]))
ones <- which(extract(t1n == 1,incendios) == 1)
resultado[i,19] <- nrow(incendios[ones,]) 
t1n[is.na(t1n)] <- 0
if (t1n@data@min == t1n@data@max | maxValue(t1n)==0 | is.na(maxValue(t1n))){
  resultado[i,5] <- 0
  resultado[i,12] <- 0
  resultado[i,23] <- 0
  resultado[i,26] <- 0
  resultado[i,42] <- 0
  resultado[i,45] <- 0
} else {
  t4n <- mask(t1n,u70)###cambio_uso_suelo
  if (t4n@data@min == t4n@data@max | maxValue(t4n)==0 | is.na(maxValue(t4n))){
    resultado[i,26] <- 0
  } else {
    rm(fu)
    fu <- freq(t4n)
    resultado[i,26] <- fu[which.max(fu[,1]),2]
  }
  t5n <- mask(t1n,uc70)###cambio_uso_suelo_c
  if (t5n@data@min == t5n@data@max | maxValue(t5n)==0 | is.na(maxValue(t5n))){
    resultado[i,29] <- 0
  } else {
    rm(gu)
    gu <- freq(t5n)
    resultado[i,29] <- gu[which.max(gu[,1]),2]
  }
  rm(f)
  f <- freq(t1n)
  resultado[i,5] <- f[2,2]
  ############Promedio_resta_bio1
  at <- crop(b1_70,t1n)
  ax <- overlay(t1n, at, fun=function(x,y){x*y})
  ax[ax==0] <- NA
  ay <- cellStats(ax,stat="mean")
  resultado[i,42] <- ay
  ############Promedio_resta_bio12
  bt <- crop(b12_70,t1n)
  bx <- overlay(t1n, bt, fun=function(x,y){x*y})
  bx[bx==0] <- NA
  by <- cellStats(bx,stat="mean")
  resultado[i,45] <- by
  ####anp
  t2n <- mask(t1n,anp_s)
  if (t2n@data@min == t2n@data@max | maxValue(t2n)==0 | is.na(maxValue(t2n))){
    resultado[i,12] <- 0
  } else {
    rm(g)
    g <- freq(t2n)
    resultado[i,12] <- g[2,2]
  }   
  b <- 0
  for (l in 1:139){
    t3n <- mask(t1n,subset(huracanes, Intensidad==paste0(huracanes$Intensidad[l])))
    if (t3n@data@min == t3n@data@max | maxValue(t3n)==0 | is.na(maxValue(t3n))){
      b <- b + 0
    } else {
      rm(h)
      h <- freq(t3n)
      b <- b + (h[which.max(h[,1]),2]*as.numeric(huracanes$Intensidad[l])/as.numeric(resultado[i,5]))
    }
  }
  resultado[i,23] <- b
}

}


write.csv(resultado,"/home/edgarandres/modelos/traslapes_2B_mx_1.csv", row.names = TRUE,col.names = TRUE)


