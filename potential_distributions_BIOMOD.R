set.seed(1111)
setwd("/home/")
library("sp")
library("raster")
library("dplyr")
library("stringr")
library("dismo")
library("readr")
library("rgdal")
library("fuzzySim")
library("rgeos")
library("magrittr")
library("tools")
library("tidyverse")
library("biomod2")


Sys.setenv(JAVA_HOME='/opt/jdk1.8.0_202/jre/bin/')
shapePath <- '/home/polygon/'
shapeLayer <- "wwf_terr_ecos"
rgnzn <- rgdal::readOGR(shapePath, shapeLayer)

crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

bioclima<-stack(list.files(path="/home/present2.5/", pattern = "*.tif$", full.names=TRUE)) 

#--------------------------------- 2030
covarDataFolder_ssp245_BCC_2030<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2021-2040.tif') 
covarDataFolder_ssp245_CAN_2030<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2021-2040.tif')   
covarDataFolder_ssp585_BCC_2030<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2021-2040.tif') 
covarDataFolder_ssp585_CAN_2030<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2021-2040.tif') 
#--------------------------------- 2050
covarDataFolder_ssp245_BCC_2050<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2041-2060.tif') 
covarDataFolder_ssp245_CAN_2050<-stack('/home/wc2.1_2.5m_bioc_CanESM5_ssp245_2041-2060.tif')   
covarDataFolder_ssp585_BCC_2050<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp585_2041-2060.tif') 
covarDataFolder_ssp585_CAN_2050<-stack('/home/wc2.1_2.5m_bioc_CanESM5_ssp585_2041-2060.tif') 
#--------------------------------- 2070
covarDataFolder_ssp245_BCC_2070<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp245_2061-2080.tif') 
covarDataFolder_ssp245_CAN_2070<-stack('/home/wc2.1_2.5m_bioc_CanESM5_ssp245_2061-2080.tif')   
covarDataFolder_ssp585_BCC_2070<-stack('/home/wc2.1_2.5m_bioc_BCC-CSM2-MR_ssp585_2061-2080.tif') 
covarDataFolder_ssp585_CAN_2070<-stack('/home/wc2.1_2.5m_bioc_CanESM5_ssp585_2061-2080.tif') 

for (j in 1:19){
  names(covarDataFolder_ssp245_BCC_2030)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp245_CAN_2030)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp585_BCC_2030)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp585_CAN_2030)[j]=paste0("Bio",j)
}

for (j in 1:19){
  names(covarDataFolder_ssp245_BCC_2050)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp245_CAN_2050)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp585_BCC_2050)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp585_CAN_2050)[j]=paste0("Bio",j)
}

for (j in 1:19){
  names(covarDataFolder_ssp245_BCC_2070)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp245_CAN_2070)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp585_BCC_2070)[j]=paste0("Bio",j)
  names(covarDataFolder_ssp585_CAN_2070)[j]=paste0("Bio",j)
}

args <- list.files("/home/species", pattern = ".csv",full.names = TRUE)

for (i in 1:length(args)) {
  
inputDataFile <- args[i]

outputFolder <- inputDataFile %>%
  basename %>%
  file_path_sans_ext
 
if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}
 
crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84")
occsData<-readr::read_csv(inputDataFile)
 
sp::coordinates(occsData) <- c("x", "y")
sp::proj4string(occsData) <- crs.wgs84

covarData <- raster::extract(bioclima, occsData)
covarData <- cbind(occsData, covarData)

completeDataCases <- covarData@data %>%
  dplyr::select(.dots=names(bioclima)) %>%
  complete.cases 
covarData <- covarData[completeDataCases, ]

speciesCol <- raster::match(outputFolder, names(occsData))
speciesCol <- speciesCol + 1
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData@data,
  sp.cols =  speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)
select_var <- correlacion$selected.vars

enviromentalVariables <- bioclima[[select_var]]

selectedVariables <- enviromentalVariables[[select_var]]

#--------------------------------- 2030
m_covarDataFolder_ssp245_BCC_2030<-covarDataFolder_ssp245_BCC_2030[[select_var]]
m_covarDataFolder_ssp245_CAN_2030<-covarDataFolder_ssp245_CAN_2030[[select_var]]
m_covarDataFolder_ssp585_BCC_2030<-covarDataFolder_ssp585_BCC_2030[[select_var]]
m_covarDataFolder_ssp585_CAN_2030<-covarDataFolder_ssp585_CAN_2030[[select_var]]

#--------------------------------- 2050
m_covarDataFolder_ssp245_BCC_2050<-covarDataFolder_ssp245_BCC_2050[[select_var]]
m_covarDataFolder_ssp245_CAN_2050<-covarDataFolder_ssp245_CAN_2050[[select_var]]
m_covarDataFolder_ssp585_BCC_2050<-covarDataFolder_ssp585_BCC_2050[[select_var]]
m_covarDataFolder_ssp585_CAN_2050<-covarDataFolder_ssp585_CAN_2050[[select_var]]

#--------------------------------- 2070
m_covarDataFolder_ssp245_BCC_2070<-covarDataFolder_ssp245_BCC_2070[[select_var]]
m_covarDataFolder_ssp245_CAN_2070<-covarDataFolder_ssp245_CAN_2070[[select_var]]
m_covarDataFolder_ssp585_BCC_2070<-covarDataFolder_ssp585_BCC_2070[[select_var]]
m_covarDataFolder_ssp585_CAN_2070<-covarDataFolder_ssp585_CAN_2070[[select_var]]

rgnzn <- spTransform(rgnzn,crs.wgs84)
ecoregionsOfInterest <- sp::over(occsData, rgnzn) %>%
  filter(!is.na(ECO_ID))

idsEcoRegions <- unique(ecoregionsOfInterest$ECO_ID)
polygonsOfInterest <- rgnzn[rgnzn$ECO_ID %in% idsEcoRegions, ]
pts_b <- gBuffer(occsData, width=3)
pts_b <- as(pts_b, 'SpatialPolygonsDataFrame')

polygonsOfInterest<-gIntersection(pts_b, polygonsOfInterest, drop_lower_td = T)
polygonsOfInterest<-gBuffer(polygonsOfInterest, width=2)

polyTransferencia<-polygonsOfInterest 

selectedVariablesCrop <- raster::crop(selectedVariables, polygonsOfInterest)
myExpl <- raster::mask(selectedVariablesCrop,polygonsOfInterest) 
myExpl<-stack(myExpl)

#ssp245_BCC_2030
env_ssp245_BCC_2030a <- raster::crop(m_covarDataFolder_ssp245_BCC_2030, polyTransferencia)
env_ssp245_BCC_2030 <- raster::mask(env_ssp245_BCC_2030a ,  polyTransferencia)
env_ssp245_BCC_2030<-stack(env_ssp245_BCC_2030)
rm(m_covarDataFolder_ssp245_BCC_2030, env_ssp245_BCC_2030a)
#ssp245_CAN_2030
env_ssp245_CAN_2030a <- raster::crop(m_covarDataFolder_ssp245_CAN_2030, polyTransferencia)
env_ssp245_CAN_2030 <- raster::mask(env_ssp245_CAN_2030a ,  polyTransferencia)
env_ssp245_CAN_2030<-stack(env_ssp245_CAN_2030)
rm(m_covarDataFolder_ssp245_CAN_2030, env_ssp245_CAN_2030a)
#ssp585_BCC_2030
env_ssp585_BCC_2030a <- raster::crop(m_covarDataFolder_ssp585_BCC_2030, polyTransferencia)
env_ssp585_BCC_2030 <- raster::mask(env_ssp585_BCC_2030a ,  polyTransferencia)
env_ssp585_BCC_2030<-stack(env_ssp585_BCC_2030)
rm(m_covarDataFolder_ssp585_BCC_2030, env_ssp585_BCC_2030a)
#ssp585_CAN_2030
env_ssp585_CAN_2030a <- raster::crop(m_covarDataFolder_ssp585_CAN_2030, polyTransferencia)
env_ssp585_CAN_2030 <- raster::mask(env_ssp585_CAN_2030a ,  polyTransferencia)
env_ssp585_CAN_2030<-stack(env_ssp585_CAN_2030)
rm(m_covarDataFolder_ssp585_CAN_2030, env_ssp585_CAN_2030a)

#ssp245_BCC_2050
env_ssp245_BCC_2050a <- raster::crop(m_covarDataFolder_ssp245_BCC_2050, polyTransferencia)
env_ssp245_BCC_2050 <- raster::mask(env_ssp245_BCC_2050a ,  polyTransferencia)
env_ssp245_BCC_2050<-stack(env_ssp245_BCC_2050)
rm(m_covarDataFolder_ssp245_BCC_2050, env_ssp245_BCC_2050a)
#ssp245_CAN_2050
env_ssp245_CAN_2050a <- raster::crop(m_covarDataFolder_ssp245_CAN_2050, polyTransferencia)
env_ssp245_CAN_2050 <- raster::mask(env_ssp245_CAN_2050a ,  polyTransferencia)
env_ssp245_CAN_2050<-stack(env_ssp245_CAN_2050)
rm(m_covarDataFolder_ssp245_CAN_2050, env_ssp245_CAN_2050a)
#ssp585_BCC_2050
env_ssp585_BCC_2050a <- raster::crop(m_covarDataFolder_ssp585_BCC_2050, polyTransferencia)
env_ssp585_BCC_2050 <- raster::mask(env_ssp585_BCC_2050a ,  polyTransferencia)
env_ssp585_BCC_2050<-stack(env_ssp585_BCC_2050)
rm(m_covarDataFolder_ssp585_BCC_2050, env_ssp585_BCC_2050a)
#ssp585_CAN_2050
env_ssp585_CAN_2050a <- raster::crop(m_covarDataFolder_ssp585_CAN_2050, polyTransferencia)
env_ssp585_CAN_2050 <- raster::mask(env_ssp585_CAN_2050a ,  polyTransferencia)
env_ssp585_CAN_2050<-stack(env_ssp585_CAN_2050)
rm(m_covarDataFolder_ssp585_CAN_2050, env_ssp585_CAN_2050a)

#ssp245_BCC_2070
env_ssp245_BCC_2070a <- raster::crop(m_covarDataFolder_ssp245_BCC_2070, polyTransferencia)
env_ssp245_BCC_2070 <- raster::mask(env_ssp245_BCC_2070a ,  polyTransferencia)
env_ssp245_BCC_2070<-stack(env_ssp245_BCC_2070)
rm(m_covarDataFolder_ssp245_BCC_2070, env_ssp245_BCC_2070a)
#ssp245_CAN_2070
env_ssp245_CAN_2070a <- raster::crop(m_covarDataFolder_ssp245_CAN_2070, polyTransferencia)
env_ssp245_CAN_2070 <- raster::mask(env_ssp245_CAN_2070a ,  polyTransferencia)
env_ssp245_CAN_2070<-stack(env_ssp245_CAN_2070)
rm(m_covarDataFolder_ssp245_CAN_2070, env_ssp245_CAN_2070a)
#ssp585_BCC_2070
env_ssp585_BCC_2070a <- raster::crop(m_covarDataFolder_ssp585_BCC_2070, polyTransferencia)
env_ssp585_BCC_2070 <- raster::mask(env_ssp585_BCC_2070a ,  polyTransferencia)
env_ssp585_BCC_2070<-stack(env_ssp585_BCC_2070)
rm(m_covarDataFolder_ssp585_BCC_2070, env_ssp585_BCC_2070a)
#ssp585_CAN_2070
env_ssp585_CAN_2070a <- raster::crop(m_covarDataFolder_ssp585_CAN_2070, polyTransferencia)
env_ssp585_CAN_2070 <- raster::mask(env_ssp585_CAN_2070a ,  polyTransferencia)
env_ssp585_CAN_2070<-stack(env_ssp585_CAN_2070)
rm(m_covarDataFolder_ssp585_CAN_2070, env_ssp585_CAN_2070a)

presencias<-data.frame(occsData)
names(presencias)[names(presencias)=="outputFolder"] <- outputFolder
presencias<-dplyr::select(presencias, c("x", "y",outputFolder))
names(presencias)

DataSpecies<-presencias

myRespName <- outputFolder
myResp <- as.numeric(DataSpecies[,myRespName])
myRespCoord = DataSpecies[c("x", "y")]
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 1000,
                                     PA.strategy = 'random')
myBiomodData
  
myBiomodOption <-BIOMOD_ModelingOptions()

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM', 'GAM', 'CTA','RF', "MAXENT.Phillips"),
  models.options = myBiomodOption,
  NbRunEval=10,
  DataSplit=70,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(outputFolder))

myBiomodModelEval <- get_evaluations(myBiomodModelOut)
write.csv(myBiomodModelEval, file = file.path(outputFolder, "myBiomodModelEval.csv"),
          row.names = FALSE)

myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = outputFolder,
  selected.models = 'all',
  binary.meth = "TSS",
  compress = 'xz',
  build.clamping.mask = TRUE,
  output.format = '.grd')

myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  prob.cv = T,
  prob.ci = T, 
  prob.ci.alpha = 0.05,
  eval.metric.quality.threshold = c(0.7),
  prob.mean.weight = T, 
  VarImport = 1)

myVarImportEM<-data.frame(get_variables_importance(myBiomodEM))
myVarImportEM<-myVarImportEM[5]
write.csv(myVarImportEM, file = file.path(outputFolder, "myVarImportEM.csv"),
          row.names = T)

myBiomodEMEval<-get_evaluations(myBiomodEM)
write.csv(myBiomodEMEval, file = file.path(outputFolder, "myBiomodEMEval.csv"),
          row.names = FALSE)

myBiomodEM_proj <-BIOMOD_EnsembleForecasting(EM.output  = myBiomodEM,
                                             projection.output = myBiomodProj,
                                             selected.models = 'all',
                                             proj.name = outputFolder,
                                             binary.meth = "TSS")

outpath<-file.path(outputFolder,paste0("proj_",outputFolder))
currentPred <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder,"_ensemble_TSSbin.grd",sep="")))
writeRaster(currentPred,
            file.path(outpath,"TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

currentPred_c <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder,"_ensemble.grd",sep="")))
writeRaster(currentPred_c,
            file.path(outpath,"c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

ClampMask <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_ClampingMask.grd",sep="")))
writeRaster(ClampMask,file.path(outpath,paste0(outputFolder,"_ClampingMask.tif")),overwrite= TRUE)

myBiomodEM_ssp245_BCC_2030<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp245_BCC_2030,
                                                        selected.models = 'all',
                                                        proj.name = "ssp245_BCC_2030",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp245_BCC_2030/","proj_ssp245_BCC_2030_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp245_BCC_2030/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp245_BCC_2030_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp245_BCC_2030/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp245_BCC_2030/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp245_BCC_2030)

myBiomodEM_ssp585_BCC_2030<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp585_BCC_2030,
                                                        selected.models = 'all',
                                                        proj.name = "ssp585_BCC_2030",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp585_BCC_2030/","proj_ssp585_BCC_2030_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp585_BCC_2030/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp585_BCC_2030_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp585_BCC_2030/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp585_BCC_2030/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp585_BCC_2030)

myBiomodEM_ssp245_CAN_2030<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp245_CAN_2030,
                                                        selected.models = 'all',
                                                        proj.name = "ssp245_CAN_2030",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp245_CAN_2030/","proj_ssp245_CAN_2030_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp245_CAN_2030/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp245_CAN_2030_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp245_CAN_2030/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp245_CAN_2030/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp245_CAN_2030)

myBiomodEM_ssp585_CAN_2030<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp585_CAN_2030,
                                                        selected.models = 'all',
                                                        proj.name = "ssp585_CAN_2030",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp585_CAN_2030/","proj_ssp585_CAN_2030_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp585_CAN_2030/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp585_CAN_2030_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp585_CAN_2030/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp585_CAN_2030/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp585_CAN_2030)

myBiomodEM_ssp245_BCC_2050<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp245_BCC_2050,
                                                        selected.models = 'all',
                                                        proj.name = "ssp245_BCC_2050",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp245_BCC_2050/","proj_ssp245_BCC_2050_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp245_BCC_2050/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp245_BCC_2050_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp245_BCC_2050/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp245_BCC_2050/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp245_BCC_2030)

myBiomodEM_ssp585_BCC_2050<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp585_BCC_2050,
                                                        selected.models = 'all',
                                                        proj.name = "ssp585_BCC_2050",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp585_BCC_2050/","proj_ssp585_BCC_2050_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp585_BCC_2050/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp585_BCC_2050_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp585_BCC_2050/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp585_BCC_2050/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp585_BCC_2050)

myBiomodEM_ssp245_CAN_2050<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp245_CAN_2050,
                                                        selected.models = 'all',
                                                        proj.name = "ssp245_CAN_2050",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp245_CAN_2050/","proj_ssp245_CAN_2050_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp245_CAN_2050/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp245_CAN_2050_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp245_CAN_2050/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp245_CAN_2050/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp245_CAN_2030)

myBiomodEM_ssp585_CAN_2050<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp585_CAN_2050,
                                                        selected.models = 'all',
                                                        proj.name = "ssp585_CAN_2050",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp585_CAN_2050/","proj_ssp585_CAN_2050_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp585_CAN_2050/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), #outputFolder debe ir entre comillas
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp585_CAN_2050_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp585_CAN_2050/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp585_CAN_2050/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp585_CAN_2050)

myBiomodEM_ssp245_BCC_2070<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp245_BCC_2070,
                                                        selected.models = 'all',
                                                        proj.name = "ssp245_BCC_2070",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp245_BCC_2070/","proj_ssp245_BCC_2070_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp245_BCC_2070/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), #outputFolder debe ir entre comillas
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp245_BCC_2070_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp245_BCC_2070/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp245_BCC_2070/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp245_BCC_2070)

myBiomodEM_ssp585_BCC_2070<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp585_BCC_2070,
                                                        selected.models = 'all',
                                                        proj.name = "ssp585_BCC_2070",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp585_BCC_2070/","proj_ssp585_BCC_2070_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp585_BCC_2070/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp585_BCC_2070_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp585_BCC_2070/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp585_BCC_2070/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp585_BCC_2070)

myBiomodEM_ssp245_CAN_2070<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp245_CAN_2070,
                                                        selected.models = 'all',
                                                        proj.name = "ssp245_CAN_2070",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp245_CAN_2070/","proj_ssp245_CAN_2070_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp245_CAN_2070/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), 
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp245_CAN_2070_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp245_CAN_2070/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp245_CAN_2070/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp245_CAN_2070)

myBiomodEM_ssp585_CAN_2070<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_ssp585_CAN_2070,
                                                        selected.models = 'all',
                                                        proj.name = "ssp585_CAN_2070",
                                                        binary.meth = "TSS")

x <- paste0(outputFolder,"/proj_ssp585_CAN_2070/","proj_ssp585_CAN_2070_", outputFolder,"_ensemble_TSSbin.grd")

FutureProj <- stack(file.path(x))

y <- paste0(outputFolder,"/proj_ssp585_CAN_2070/TSSbin.tif")

writeRaster(FutureProj, 
            file.path(y), #outputFolder debe ir entre comillas
            suffix='names',
            bylayer=TRUE, 
            overwrite= TRUE)

myBiomodRangeSize<-BIOMOD_RangeSize(currentPred, FutureProj,  
                                    SpChange.Save="ssp585_CAN_2070_rs")

write.csv(myBiomodRangeSize$Compt.By.Models,  
          file = file.path(paste0(outputFolder,"/proj_ssp585_CAN_2070/rangesize.csv")), 
          row.names = TRUE)
writeRaster(myBiomodRangeSize$Diff.By.Pixel, 
            file.path(paste0(outputFolder, "/proj_ssp585_CAN_2070/RS.tif")), 
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)
rm(FutureProj, myBiomodRangeSize, myBiomodEM_ssp585_CAN_2070)

unlink(file.path(outputFolder, paste("models/")), recursive = T)
unlink(list.files(outputFolder, pattern = "*.grd$", full.names=TRUE, recursive = T))
unlink(list.files(outputFolder, pattern = "*.gri$", full.names=TRUE, recursive = T))

 
}
