library(rgdal)
library(maptools)
library(lubridate)

source("Scripts/WeatherProcessor.R")

Aggregated_FireScars <- readOGR(dsn="Data/BurnScars/", layer="agg2")
crs.k <- CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

Aggregated_FireScars$STARTDATE <- ymd(as.character(Aggregated_FireScars$STARTDATE))
Aggregated_FireScars$DATE_START <- ymd(as.character(Aggregated_FireScars$DATE_START))
Aggregated_FireScars$DATE_END <- ymd(as.character(Aggregated_FireScars$DATE_END))
Aggregated_FireScars <- spTransform(Aggregated_FireScars,crs.k)

FireScars_WoodyImpact <- subset(Aggregated_FireScars,!is.na(WOODYIMPAC))

krugerMAP_UTM <- raster(x="/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/krugerMAP_UTM")
krugerWoodyCover <- raster(x="Data/WoodyCover/wcp_map_fin.tif")
krugerWoodyCover_UTM <- projectRaster(krugerWoodyCover,krugerMAP_UTM)
names(krugerWoodyCover_UTM) <- "WoodyCover"


krugerFireScarAnalysis_brick <- brick(krugerMAP_UTM,krugerWoodyCover_UTM)


extractedFireScars <- extract(x = krugerFireScarAnalysis_brick,
                              y = FireScars_WoodyImpact,
                              fun = mean,
                              sp = TRUE)

extractedFireScars_df <- as.data.frame(extractedFireScars)

extractedFireScarsAgg_df <- ddply(extractedFireScars_df,.(FIREID),transform,
                                  IGNITIONSE = IGNITIONSE,
                                  STARTDATE = STARTDATE,
                                  HERBACEOUS = HERBACEOUS,
                                  GENERALOBS = GENERALOBS,
                                  AREA_HA = sum(AREA_HA),
                                  MAP = mean(layer,na.rm = TRUE),
                                  WoodyCover = mean(WoodyCover,na.rm = TRUE))
extractedFireScarsAgg_df <-subset(extractedFireScarsAgg_df, !is.na(layer))
extractedFireScarsAgg_df$layer <- NULL
extractedFireScarsAgg_df <- extractedFireScarsAgg_df[!duplicated(extractedFireScarsAgg_df),]

extractedFireScars_wx_df <- merge(extractedFireScarsAgg_df,Kruger_Wx_Combined,by.x = "IGNITIONSE", by.y = "Station_Long")

extractedFireScars_wx_df$Year <- year(as.Date(as.character(extractedFireScars_wx_df$Year),format="%Y"))


testPlied <- ddply(extractedFireScars_wx_df,
                   .(FIREID),
                   mutate,
                   PreviousYears = year(STARTDATE) - 2)

testPlied <- ddply(testPlied,
                   .(FIREID),
                   subset,
                   Year > PreviousYears & Year <= year(STARTDATE)
                  )
aggPlied <- ddply(testPlied,"FIREID",summarise,
                   PreviousMAP = sum(AnnualPrecip))

firescar_final <- merge(testPlied,aggPlied,by="FIREID")

firescar_final$WOODYIMPAC <- factor(firescar_final$WOODYIMPAC,levels(firescar_final$WOODYIMPAC)[c(2,4,1,3,5)])
levels(firescar_final$HERBACEOUS) <- c("Clean","Clean","Moderately clean","Patchy","Very patchy")
firescar_final$HERBACEOUS <- factor(firescar_final$HERBACEOUS,levels(firescar_final$HERBACEOUS)[c(2,4,1,3,5)])

firescar_final_dry <- subset(firescar_final,month(STARTDATE) >= 7 & month(STARTDATE) < 10)
firescar_final_wet <- subset(firescar_final,month(STARTDATE) >= 7 | month(STARTDATE) > 10)

firescar_final_dry$Season <- "Dry"
firescar_final_wet$Season <- "Wet"

firescar_final <- rbind(firescar_final_dry,firescar_final_wet)
firescar_final_wet$Season <- as.factor(firescar_final_wet$Season)

rm(firescar_final_dry)
rm(firescar_final_wet)
rm(aggPlied)
rm(testPlied)
rm(extractedFireScars)
rm(extractedFireScars_df)
rm(extractedFireScars_wx_df)
rm(FireScars_WoodyImpact)
rm(FireScars_WoodyImpact_df)
rm(krugerFireScarAnalysis_brick)