## @knitr loadEverything

library(gstat)
library(ggplot2)
library(raster)
library(lattice)
library(rgdal)
library(lubridate)

crs.k <- CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
crs.m <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

FIRMS <- readOGR(dsn="/Users/danielg7/Documents/FIRMS/",layer="firms138621406915198_MCD14ML")
krugerOutline <- readOGR(dsn="/Users/danielg7/Documents/LandsatImagery/Boundary/",layer="boundary_kruger")
krugerAvgTmin_UTM <- raster(x="/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/krugerAvgTmin_UTM")
krugerMAP_UTM <- raster(x="/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/krugerMAP_UTM")
krugerFirelineIntensity_UTM <- raster(x="/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/krugerFireLineIntensity_UTM")
krugerGly <- readOGR(dsn="Data/",layer="KNP_GraniticAndBasaltic")

krugerWoodyCover <- raster(x="Data/WoodyCover/wcp_map_fin.tif")

krugerMAP_UTM <- resample(krugerMAP_UTM,krugerAvgTmin_UTM)
krugerFirelineIntensity_UTM <- resample(krugerFirelineIntensity_UTM,krugerAvgTmin_UTM)
krugerWoodyCover_UTM <- projectRaster(krugerWoodyCover,krugerAvgTmin_UTM)

rm(krugerWoodyCover)

names(krugerMAP_UTM) <- "MAP"
names(krugerFirelineIntensity_UTM) <- "Fireline_Intensity"
names(krugerWoodyCover_UTM) <- "WoodyCover"

FIRMS <- spTransform(FIRMS,crs.k)
krugerOutline <- spTransform(krugerOutline,crs.k)
krugerGly_UTM <- spTransform(krugerGly,crs.k)

## @knitr SubsetAndExtract
FIRMS_highConfidence <- subset(FIRMS,CONFIDENCE >= 95)

FIRMS_highConfidence$ACQ_DATE <- ymd(as.character(FIRMS_highConfidence$ACQ_DATE))
FIRMS_highConfidence$Month <- month(FIRMS_highConfidence$ACQ_DATE)
FIRMS_Kruger_DrySeason <- subset(FIRMS_highConfidence,Month >= 7 & Month < 10)
FIRMS_Kruger_WetSeason <- subset(FIRMS_highConfidence,Month < 7 | Month > 10)

krugerBrick <- brick(krugerMAP_UTM,krugerAvgTmin_UTM,krugerFirelineIntensity_UTM,krugerWoodyCover_UTM)

DrySeasonMAP <- extract(krugerBrick,FIRMS_Kruger_DrySeason,method='bilinear',df=TRUE,sp=TRUE)
WetSeasonMAP <- extract(krugerBrick,FIRMS_Kruger_WetSeason,method='bilinear',df=TRUE,sp=TRUE)

DrySeasonMAP <- as.data.frame(DrySeasonMAP)
WetSeasonMAP <- as.data.frame(WetSeasonMAP)

DrySeasonMAP$Season <- "Dry"
WetSeasonMAP$Season <- "Wet"

FRP_Variables <- rbind(DrySeasonMAP,WetSeasonMAP)

FRP_Variables_subsetWC <- subset(FRP_Variables,WoodyCover >= 50)
