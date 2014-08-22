# Copyright notice: ----
# This script is provided with a Creative Commons - Attribution license, as defined on:
# http://creativecommons.org/licenses/by/3.0/us/
#
#
# Author Contact:
# Daniel Godwin
# danielg7@gmail.com
# Savanna Ecology Lab
# Division of Biological Sciences
# University of Missouri - Columbia
#
# Script Intent: ---
# This script downloads and assembles a series of useful raster datasets for Kruger National Park
#
# Completeness: Incomplete
#
# Inputs: ----
# 
#
#
# Outputs: ----
# 
# 
#
# TODO:  ----
#        
#        
#
#        
# 
# Load required packages ----
library(raster)
library(weatherData)
# Define working directories ----
DataFolder <- "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data"
ScratchData <- "/Users/danielgodwin/Documents/"
# Functions ----
projectAndCrop <- function(image,crs,maskOutline){
  projection(image) <- crs
  
  scratch <- crop(image,extent(maskOutline))
  scratch <- mask(scratch,maskOutline)
  
  return(scratch)
}
# Read in data ----
krugerPrec1 <- getData(name="worldclim",var="prec",res=.5,lat=-24,lon=28)
krugerPrec2 <- getData(name="worldclim",var="prec",res=.5,lat=-24,lon=31)
krugertmin1 <- getData(name="worldclim",var="tmin",res=.5,lat=-24,lon=28)
krugertmin2 <- getData(name="worldclim",var="tmin",res=.5,lat=-24,lon=31)
# Mosaic data ----
krugerPrec <- mosaic(krugerPrec1,krugerPrec2,fun=max)
krugertmin <- mosaic(krugertmin1,krugertmin2,fun=max)

crs.k <- CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
crs.m <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

krugerOutline <- readOGR(dsn="/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/",layer="boundary_kruger")

landtypes <- readOGR(dsn="/Users/danielg7/Documents/KNP-GIS",layer="LANDC_LANDTYPES")
LTYPE <- c("Malelane","Skukuza","Satara","Vutome","Sabiepoort","Phalaborwa","Letaba","Bulweni","Klipkoppies","Nwambiya","Pafuri")
MFRI <- c(5.3,4.2,5.0,5.8,4.1,7.1,4.4,3.7,6.1,3.6,2.7)
landtypesMFRI <- data.frame(LTYPE = LTYPE, MFRI = MFRI)
landtypes <- merge(landtypes,landtypesMFRI, by="LTYPE")
landtypes <- spTransform(x=landtypes, CRSobj=crs.m)

ext <- extent(30.89167, 32.03333 ,-25.53333 , -22.325 )
xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))
n <- 10
landtypes_raster <- raster(ext, ncol=xy[1]*5, nrow=xy[2]*5)
landtypes_raster <- rasterize(landtypes,landtypes_raster,MFRI)

krugerOutline_UTM <- spTransform(x=krugerOutline, CRSobj=crs.k)

krugerPrec_sub <- projectRaster(krugerPrec_sub,crs=crs.k)

krugerPrec_sub <- projectAndCrop(krugerPrec_sub,crs.k,krugerOutline_UTM)
krugerMAP <- sum(krugerPrec_sub$layer.1,krugerPrec_sub$layer.2,krugerPrec_sub$layer.3,krugerPrec_sub$layer.4,krugerPrec_sub$layer.5,krugerPrec_sub$layer.6,krugerPrec_sub$layer.7,krugerPrec_sub$layer.8,krugerPrec_sub$layer.9,krugerPrec_sub$layer.10,krugerPrec_sub$layer.11,krugerPrec_sub$layer.12)



setwd(DataFolder)

km <- writeRaster(krugerMAP,filename="krugerMAP_UTM")

#km <- writeRaster(krugertmin_sub,filename="krugerTmin_UTM")
#krugertmin <- projectRaster(krugertmin,crs=crs.k)
krugertmin_sub <- projectAndCrop(krugertmin,crs.m,krugerOutline_proj)
krugerMeanTmin <- mean(krugertmin_sub)/10
krugerMeanTmin$MeanMinTemp <- krugerMeanTmin$layer
krugerMeanTmin <- projectRaster(krugerMeanTmin,crs=crs.k)

km <- writeRaster(krugerMeanTmin$MeanMinTemp,filename="krugerAvgTmin_UTM",overwrite=TRUE)

#krugerMeanRH <- predict(krugerMeanTmin,lm_rhbymintemp)/100

# Clean up ----
rm(krugerPrec_sub)
rm(krugerPrec)
rm(krugerPrec1)
rm(krugerPrec2)
rm(krugertmin)
rm(krugertmin1)
rm(krugertmin2)



# 
# 
# MAP <- c(450,550,650,750,850)
# MFRI <- c(5.0,5.3,5.2,3.5,3.5)
# test <- data.frame(MAP,MFRI)
# lm_MAP <- lm(MFRI~MAP,test)
# 
# krugerMAP$MAP <- krugerMAP$layer
# 
# krugerMFRIfromMAP <- predict(krugerMAP,lm_MAP)
# 
# landtypes_raster2 <- resample(x=landtypes_raster,y=krugerMAP)
# krugerBrick <- brick(krugerMAP$MAP,landtypes_raster2,krugerMFRIfromMAP,krugerMeanRH)
# names(krugerBrick) <- c("MAP","TrollopeMAP","MFRIfromMAP","RH")
# 
# krugerFuelLoad <- 382.9 + 3.3 * krugerBrick$MAP + 979.4 * krugerBrick$TrollopeMAP - 0.001 * krugerBrick$MAP^2 + 0.37*krugerBrick$MAP*krugerBrick$TrollopeMAP - 161.8*krugerBrick$TrollopeMAP^2

#RelativeHumidity <- .25 # Range values taken from Trollope 2002
#FuelMoisture <- .06 # Fuels assumed to take the value of RH immediately
#WindSpeed <- 5

#krugerFirelineIntensity <- 2729 + 0.8684*krugerFuelLoad - 530*sqrt(FuelMoisture) - 0.1907*krugerBrick$RH^2 - 5961/WindSpeed

#krugerIntensityLog <- log(krugerFirelineIntensity)

#krugerTopKillProb <- 4.14e-05 * krugerIntensity + 4.43e-01

#krugerTopKillProb <- exp(4.3 - 5.003*log(2) + 0.004408*sqrt(krugerFirelineIntensity)) / (1 + exp(4.3 - 5.003*log(2) + 0.004408*sqrt(krugerFirelineIntensity)))

