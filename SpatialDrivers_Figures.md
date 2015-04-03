# Spatial Drivers of Fire Intensity





We don't know much about what drives fire intensity across different scales. To address this, fire radiative power data were obtained from [FIRMS](https://earthdata.nasa.gov/data/near-real-time-data/firms) for the period 1 Jan 2004 to 1 Jan 2014.

Data were subset to only include fire detections of > 95% confidence.

Fire detections were then associated with:
* Mean annual precipitation
* Geologic parent material
* Woody cover

**Methods**

Run the analyses related to the fire radiative power detections.


```r
source("Scripts/FIRMS_Map.R",echo=TRUE,verbose = FALSE)
```

```
## 
## > library(ggplot2)
## 
## > library(raster)
## 
## > library(rgdal)
## 
## > library(lubridate)
## 
## > library(plyr)
## 
## > crs.k <- CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
## 
## > crs.m <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
## 
## > FIRMS <- readOGR(dsn = "/Users/danielg7/Documents/FIRMS/", 
## +     layer = "firms138621406915198_MCD14ML")
## OGR data source with driver: ESRI Shapefile 
## Source: "/Users/danielg7/Documents/FIRMS/", layer: "firms138621406915198_MCD14ML"
## with 10370 features and 12 fields
## Feature type: wkbPoint with 2 dimensions
## 
## > krugerOutline <- readOGR(dsn = "/Users/danielg7/Documents/LandsatImagery/Boundary/", 
## +     layer = "boundary_kruger")
## OGR data source with driver: ESRI Shapefile 
## Source: "/Users/danielg7/Documents/LandsatImagery/Boundary/", layer: "boundary_kruger"
## with 1 features and 7 fields
## Feature type: wkbPolygon with 2 dimensions
## 
## > krugerAvgTmin_UTM <- raster(x = "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Dat ..." ... [TRUNCATED] 
## 
## > krugerMAP_UTM <- raster(x = "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/kr ..." ... [TRUNCATED] 
## 
## > krugerFirelineIntensity_UTM <- raster(x = "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofT ..." ... [TRUNCATED] 
## 
## > krugerGly <- readOGR(dsn = "Data/", layer = "KNP_GraniticAndBasaltic")
## OGR data source with driver: ESRI Shapefile 
## Source: "Data/", layer: "KNP_GraniticAndBasaltic"
## with 2 features and 1 fields
## Feature type: wkbPolygon with 2 dimensions
## 
## > krugerWoodyCover <- raster(x = "Data/WoodyCover/wcp_map_fin.tif")
## 
## > krugerMAP_UTM <- resample(krugerMAP_UTM, krugerAvgTmin_UTM)
## 
## > krugerGly_UTM <- spTransform(krugerGly, crs.k)
## 
## > krugerFirelineIntensity_UTM <- resample(krugerFirelineIntensity_UTM, 
## +     krugerAvgTmin_UTM)
## 
## > krugerWoodyCover_UTM <- projectRaster(krugerWoodyCover, 
## +     krugerAvgTmin_UTM)
## 
## > krugerOverlayRaster <- raster(krugerWoodyCover_UTM)
## 
## > krugerGly_UTM <- spTransform(krugerGly, crs.k)
## 
## > krugerGlyRaster <- rasterize(krugerGly_UTM, krugerOverlayRaster)
## Found 2 region(s) and 14 polygon(s)
## 
## > names(krugerMAP_UTM) <- "MAP"
## 
## > names(krugerFirelineIntensity_UTM) <- "Fireline_Intensity"
## 
## > names(krugerWoodyCover_UTM) <- "WoodyCover"
## 
## > FIRMS <- spTransform(FIRMS, crs.k)
## 
## > krugerOutline <- spTransform(krugerOutline, crs.k)
## 
## > FIRMS_highConfidence <- subset(FIRMS, CONFIDENCE >= 
## +     95)
## 
## > FIRMS_highConfidence$ACQ_DATE <- ymd(as.character(FIRMS_highConfidence$ACQ_DATE))
## 
## > FIRMS_highConfidence$Month <- month(FIRMS_highConfidence$ACQ_DATE)
## 
## > FIRMS_Kruger_DrySeason <- subset(FIRMS_highConfidence, 
## +     Month >= 7 & Month < 10)
## 
## > FIRMS_Kruger_WetSeason <- subset(FIRMS_highConfidence, 
## +     Month < 7 | Month > 10)
## 
## > krugerBrick <- brick(krugerMAP_UTM, krugerAvgTmin_UTM, 
## +     krugerFirelineIntensity_UTM, krugerWoodyCover_UTM, krugerGlyRaster$layer)
## 
## > names(krugerBrick)[5] <- "Geology"
## 
## > DrySeasonMAP <- raster::extract(krugerBrick, FIRMS_Kruger_DrySeason, 
## +     method = "bilinear", df = TRUE, sp = TRUE)
## 
## > WetSeasonMAP <- raster::extract(krugerBrick, FIRMS_Kruger_WetSeason, 
## +     method = "bilinear", df = TRUE, sp = TRUE)
## 
## > DrySeasonMAP <- as.data.frame(DrySeasonMAP)
## 
## > WetSeasonMAP <- as.data.frame(WetSeasonMAP)
## 
## > DrySeasonMAP$Season <- "Dry"
## 
## > WetSeasonMAP$Season <- "Wet"
## 
## > FRP_Variables <- rbind(DrySeasonMAP, WetSeasonMAP)
## 
## > FRP_Variables$Geology <- as.factor(FRP_Variables$Geology)
## 
## > levels(FRP_Variables$Geology)[9] <- "Basaltic"
## 
## > levels(FRP_Variables$Geology)[1] <- "Granitic"
## 
## > levels(FRP_Variables$Geology)[c(2, 3, 4, 5, 6, 7, 
## +     8)] <- "Other"
## 
## > FRP_Variables$MAP_cut <- cut(FRP_Variables$MAP, breaks = seq(400, 
## +     950, 50))
## 
## > levels(FRP_Variables$MAP_cut) <- seq(400, 950, 50)
## 
## > FRP_Variables$MAP_cut <- as.character(FRP_Variables$MAP_cut)
## 
## > FRP_Variables$MAP_cut <- as.numeric(FRP_Variables$MAP_cut)
## 
## > FRP_Variables_agg <- ddply(FRP_Variables, .(MAP_cut), 
## +     summarize, WC = mean(WoodyCover, na.rm = TRUE), WC_SE = sd(WoodyCover)/sqrt(length(Wood .... [TRUNCATED]
```
Run the model as it relates to predicting intensity.

```r
source("Scripts/IntensityMap.R",echo=TRUE,verbose = FALSE)
```

```
## 
## > library(ggplot2)
## 
## > library(raster)
## 
## > library(rgdal)
## 
## > library(lubridate)
## 
## > crs.k <- CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
## 
## > krugerAvgTmin_UTM <- raster(x = "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Dat ..." ... [TRUNCATED] 
## 
## > krugerMAP_UTM <- raster(x = "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/kr ..." ... [TRUNCATED] 
## 
## > minRH_KNP_UTM <- raster(x = "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/mi ..." ... [TRUNCATED] 
## 
## > krugerFRI <- readOGR(dsn = "Data/FRI/", layer = "fire return interval_1941to2006")
## OGR data source with driver: ESRI Shapefile 
## Source: "Data/FRI/", layer: "fire return interval_1941to2006"
## with 34856 features and 4 fields
## Feature type: wkbPolygon with 2 dimensions
## 
## > krugerFRI_UTM <- spTransform(krugerFRI, crs.k)
## 
## > krugerFRI_UTM <- subset(krugerFRI_UTM, MFRI <= 6)
## 
## > krugerFRIRaster <- raster(krugerMAP_UTM)
## 
## > source("Scripts/GMED.R")
## OGR data source with driver: ESRI Shapefile 
## Source: "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/", layer: "boundary_kruger"
## with 1 features and 7 fields
## Feature type: wkbPolygon with 2 dimensions
## 
## > krugerMAP_UTM <- resample(x = gmap_ZA_KNP, y = krugerFRIRaster)
## 
## > minRH_KNP_UTM <- resample(x = minRH_KNP_UTM, y = krugerFRIRaster)
## 
## > krugerFRIRaster <- rasterize(krugerFRI_UTM, krugerFRIRaster, 
## +     field = krugerFRI_UTM$MFRI)
## Found 29349 region(s) and 29487 polygon(s)
## 
## > krugerBrick <- brick(krugerMAP_UTM, krugerFRIRaster, 
## +     minRH_KNP_UTM)
## 
## > names(krugerBrick) <- c("MAP", "MFRI", "RH")
## 
## > krugerFuelLoad <- 382.9 + 3.3 * krugerBrick$MAP + 
## +     979.4 * krugerBrick$MFRI - 0.001 * krugerBrick$MFRI^2 + 0.37 * 
## +     krugerBrick$MAP * kru .... [TRUNCATED] 
## 
## > RelativeHumidity <- c(4.2, 36.6, 82)
## 
## > FuelMoisture <- c(7.5, 32.1, 68.8)
## 
## > WindSpeed <- c(6.7, 2.6, 0.3)
## 
## > krugerFirelineIntensity <- 2729 + 0.8684 * krugerFuelLoad - 
## +     530 * sqrt(FuelMoisture[1]) - 0.1907 * krugerBrick$RH^2 - 
## +     5961/WindSpeed[1 .... [TRUNCATED] 
## 
## > krugerFirelineIntensity$Average <- 2729 + 0.8684 * 
## +     krugerFuelLoad - 530 * sqrt(FuelMoisture[2]) - 0.1907 * RelativeHumidity[2]^2 - 
## +     596 .... [TRUNCATED] 
## 
## > names(krugerFirelineIntensity)[1] <- "High"
## 
## > krugerWoodyCover <- raster(x = "Data/WoodyCover/wcp_map_fin.tif")
## 
## > krugerWoodyCover_UTM <- projectRaster(krugerWoodyCover, 
## +     krugerAvgTmin_UTM)
## 
## > names(krugerWoodyCover_UTM) <- "WoodyCover"
## 
## > krugerWoodyCover_newExtent <- resample(krugerWoodyCover_UTM, 
## +     krugerFirelineIntensity)
## 
## > krugerGlyRaster_newExtent <- resample(krugerGlyRaster, 
## +     krugerFirelineIntensity)
## 
## > krugerIntensityInvestigation <- brick(krugerWoodyCover_newExtent, 
## +     krugerFirelineIntensity, krugerGlyRaster_newExtent, krugerMAP_UTM)
## 
## > krugerIntensityInvestigationDF <- na.omit(as.data.frame(krugerIntensityInvestigation))
## 
## > names(krugerIntensityInvestigationDF) <- c("WoodyCover", 
## +     "FirelineIntensity_High", "FirelineIntensity_Average", "Geology", 
## +     "MAP")
## 
## > krugerIntensityInvestigationDF$Geology <- as.factor(krugerIntensityInvestigationDF$Geology)
## 
## > levels(krugerIntensityInvestigationDF$Geology) <- c("Granitic", 
## +     "Basaltic")
```
Run the Fire Scar Analysis

```r
source("Scripts/FireScar_Analysis.R",echo=TRUE,verbose = FALSE)
```

```
## 
## > library(rgdal)
## 
## > library(lubridate)
## 
## > source("Scripts/WeatherProcessor.R")
## 
## > Aggregated_FireScars <- readOGR(dsn = "Data/BurnScars/", 
## +     layer = "agg2")
## OGR data source with driver: ESRI Shapefile 
## Source: "Data/BurnScars/", layer: "agg2"
## with 3533 features and 12 fields
## Feature type: wkbPolygon with 2 dimensions
## 
## > crs.k <- CRS("+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
## 
## > Aggregated_FireScars$STARTDATE <- ymd(as.character(Aggregated_FireScars$STARTDATE))
## 
## > Aggregated_FireScars$DATE_START <- ymd(as.character(Aggregated_FireScars$DATE_START))
## 
## > Aggregated_FireScars$DATE_END <- ymd(as.character(Aggregated_FireScars$DATE_END))
## 
## > Aggregated_FireScars <- spTransform(Aggregated_FireScars, 
## +     crs.k)
## 
## > FireScars_WoodyImpact <- subset(Aggregated_FireScars, 
## +     !is.na(WOODYIMPAC))
## 
## > krugerMAP_UTM <- raster(x = "/Users/danielgodwin/Dropbox/Graduate School/Dissertation/Chapter 2 - Ignition Variation/SpatialDriversofTopKill/Data/kr ..." ... [TRUNCATED] 
## 
## > krugerWoodyCover <- raster(x = "Data/WoodyCover/wcp_map_fin.tif")
## 
## > krugerWoodyCover_UTM <- projectRaster(krugerWoodyCover, 
## +     krugerMAP_UTM)
## 
## > names(krugerWoodyCover_UTM) <- "WoodyCover"
## 
## > krugerFireScarAnalysis_brick <- brick(krugerMAP_UTM, 
## +     krugerWoodyCover_UTM)
## 
## > extractedFireScars <- raster::extract(x = krugerFireScarAnalysis_brick, 
## +     y = FireScars_WoodyImpact, fun = mean, sp = TRUE)
## 
## > extractedFireScars_df <- as.data.frame(extractedFireScars)
## 
## > extractedFireScars_df <- subset(extractedFireScars_df, 
## +     WoodyCover >= 0)
## 
## > extractedFireScarsAgg_df <- ddply(extractedFireScars_df, 
## +     .(FIREID), transform, IGNITIONSE = IGNITIONSE, STARTDATE = STARTDATE, 
## +     HERBACE .... [TRUNCATED] 
## 
## > extractedFireScarsAgg_df <- subset(extractedFireScarsAgg_df, 
## +     !is.na(layer))
## 
## > extractedFireScarsAgg_df$layer <- NULL
## 
## > extractedFireScarsAgg_df <- extractedFireScarsAgg_df[!duplicated(extractedFireScarsAgg_df), 
## +     ]
## 
## > extractedFireScars_wx_df <- merge(extractedFireScarsAgg_df, 
## +     Kruger_Wx_Combined, by.x = "IGNITIONSE", by.y = "Station_Long")
## 
## > extractedFireScars_wx_df$Year <- year(as.Date(as.character(extractedFireScars_wx_df$Year), 
## +     format = "%Y"))
## 
## > testPlied <- ddply(extractedFireScars_wx_df, .(FIREID), 
## +     mutate, PreviousYears = year(STARTDATE) - 2)
## 
## > testPlied <- ddply(testPlied, .(FIREID), subset, Year > 
## +     PreviousYears & Year <= year(STARTDATE))
## 
## > aggPlied <- ddply(testPlied, "FIREID", summarise, 
## +     PreviousMAP = mean(AnnualPrecip))
## 
## > firescar_final <- merge(extractedFireScarsAgg_df, 
## +     aggPlied, by = "FIREID")
## 
## > firescar_final$Year <- year(as.Date(as.character(firescar_final$STARTDATE), 
## +     format = "%Y"))
## 
## > firescar_final$WOODYIMPAC <- factor(firescar_final$WOODYIMPAC, 
## +     levels(firescar_final$WOODYIMPAC)[c(2, 4, 1, 3, 5)])
## 
## > levels(firescar_final$HERBACEOUS) <- c("Clean", "Clean", 
## +     "Moderately clean", "Patchy", "Very patchy")
## 
## > firescar_final_dry <- subset(firescar_final, month(STARTDATE) <= 
## +     7 & month(STARTDATE) < 10)
## 
## > firescar_final_wet <- subset(firescar_final, month(STARTDATE) >= 
## +     7 | month(STARTDATE) > 10)
## 
## > firescar_final_dry$Season <- "Dry"
## 
## > firescar_final_wet$Season <- "Wet"
## 
## > firescar_final <- rbind(firescar_final_dry, firescar_final_wet)
## 
## > firescar_final$Season <- as.factor(firescar_final$Season)
## 
## > firescar_final$WoodyImpactLevel <- firescar_final$WOODYIMPAC
## 
## > levels(firescar_final$WoodyImpactLevel) <- c(0, 1, 
## +     2, 3, 4)
## 
## > firescar_final$WoodyImpactLevel <- as.numeric(firescar_final$WoodyImpactLevel)
## 
## > Impact_null <- glm(WoodyImpactLevel ~ 1, data = firescar_final)
## 
## > Impact_season <- glm(WoodyImpactLevel ~ Season, data = firescar_final)
## 
## > Impact_MAP <- glm(WoodyImpactLevel ~ MAP, data = firescar_final)
## 
## > Impact_CAUSE <- glm(WoodyImpactLevel ~ CAUSE, data = firescar_final)
## 
## > Impact_WC <- glm(WoodyImpactLevel ~ WoodyCover, data = firescar_final)
## 
## > rm(firescar_final_dry)
## 
## > rm(firescar_final_wet)
## 
## > (aggPlied)
##          FIREID PreviousMAP
## 1  20020408MOO1       424.5
## 2  20020417KRO1       523.5
## 3  20020424MAH1       442.7
## 4  20020502LET1       399.2
## 5  20020502PUN1       641.9
## 6  20020503OLI1       462.1
## 7  20020506PUN1       641.9
## 8  20020507PUN1       641.9
## 9  20020509SHI1       334.8
## 10 20020511PAF1       424.6
## 11 20020513PRE1       534.8
## 12 20020513VLA1       414.8
## 13 20020519VLA1       414.8
## 14 20020522OSA1       645.0
## 15 20020523MAL1       463.8
## 16 20020523TSH1       462.1
## 17 20020524LET1       399.2
## 18 20020524TSH1       462.1
## 19 20020530MAH1       442.7
## 20 20020601PAF1       424.6
## 21 20020703OSA1       645.0
## 22 20020705OSA1       645.0
## 23 20020707OLI1       462.1
## 24 20020708VLA1       414.8
## 25 20020711KRO1       523.5
## 26 20020716PRE1       534.8
## 27 20020716VLA1       414.8
## 28 20020718KRO1       523.5
## 29 20020723SKZ1       519.2
## 30 20020727OLI1       462.1
## 31 20020729PRE1       534.8
## 32 20020731OLI1       462.1
## 33 20020802PUN1       641.9
## 34 20020802TSH1       462.1
## 35 20020812PRE1       534.8
## 36 20020814VLA1       414.8
## 37 20020816PAF1       424.6
## 38 20020817VLA1       414.8
## 39 20020821TSH1       462.1
## 40 20020822OSA1       645.0
## 41 20020822SHA1       511.2
## 42 20020825KRO1       523.5
## 43 20020826PAF1       424.6
## 44 20020827KRO1       523.5
## 45 20020830OLI1       462.1
## 46 20020902PAF1       424.6
## 47 20020904PAF1       424.6
## 48 20020911SAT1       353.9
## 49 20020913PAF1       424.6
## 50 20020913PRE1       534.8
## 51 20020918SAT1       353.9
## 52 20020922MOO1       424.5
## 53 20020924PRE1       534.8
## 54 20020925PAF1       424.6
## 55 20020925PHA1       437.4
## 56 20020926VLA1       414.8
## 57 20021006MAL1       463.8
## 58 20021018PHA1       437.4
## 59 20021024VLA1       414.8
## 60 20021123TSH1       462.1
## 61 20021209OLI1       462.1
## 62 20021222SHI1       334.8
## 63 20021224VLA1       414.8
## 
## > (testPlied)
##           IGNITIONSE       FIREID  STARTDATE      CAUSE     AGENT
## 1          Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 2          Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 3   Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 4   Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 5         Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 6         Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 7             Letaba 20020502LET1 2002-05-02 Management    Ranger
## 8             Letaba 20020502LET1 2002-05-02 Management    Ranger
## 9        Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 10       Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 11          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 12          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 13       Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 14       Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 15       Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 16       Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 17        Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 18        Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 19            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 20            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 21      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 22      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 23       Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 24       Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 25       Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 26       Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 27       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 28       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 29          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 30          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 31         Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 32         Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 33            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 34            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 35         Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 36         Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 37        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 38        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 39            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 40            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 41       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 42       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 43       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 44       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 45          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 46          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 47       Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 48       Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 49  Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 50  Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 51      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 52      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 53       Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 54       Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 55  Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 56  Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 57           Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 58           Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 59          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 60          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 61      Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 62      Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 63          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 64          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 65       Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 66       Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 67         Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 68         Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 69      Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 70      Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 71       Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 72       Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 73            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 74            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 75       Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 76       Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 77         Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 78         Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 79       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 80       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 81          Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 82          Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 83  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 84  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 85            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 86            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 87  Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 88  Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 89          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 90          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 91            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 92            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 93            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 94            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 95            Satara 20020911SAT1 2002-09-11      Other     Other
## 96            Satara 20020911SAT1 2002-09-11      Other     Other
## 97            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 98            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 99      Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 100     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 101           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 102           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 103        Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 104        Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 105     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 106     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 107           Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 108           Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 109       Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 110       Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 111      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 112      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 113         Malelane 20021006MAL1 2002-10-06      Other     Other
## 114         Malelane 20021006MAL1 2002-10-06      Other     Other
## 115       Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 116       Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 117      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 118      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 119        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 120        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 121         Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 122         Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 123       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 124       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 125      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 126      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
##           HERBACEOUS  WOODYIMPAC
## 1   Moderately clean    Moderate
## 2   Moderately clean    Moderate
## 3             Patchy      Slight
## 4             Patchy      Slight
## 5              Clean      Severe
## 6              Clean      Severe
## 7              Clean      Severe
## 8              Clean      Severe
## 9   Moderately clean    Moderate
## 10  Moderately clean    Moderate
## 11  Moderately clean    Moderate
## 12  Moderately clean    Moderate
## 13  Moderately clean    Moderate
## 14  Moderately clean    Moderate
## 15  Moderately clean    Moderate
## 16  Moderately clean    Moderate
## 17  Moderately clean      Slight
## 18  Moderately clean      Slight
## 19       Very patchy        None
## 20       Very patchy        None
## 21             Clean      Severe
## 22             Clean      Severe
## 23             Clean      Severe
## 24             Clean      Severe
## 25             Clean      Severe
## 26             Clean      Severe
## 27  Moderately clean    Moderate
## 28  Moderately clean    Moderate
## 29             Clean      Slight
## 30             Clean      Slight
## 31             Clean    Moderate
## 32             Clean    Moderate
## 33             Clean      Severe
## 34             Clean      Severe
## 35            Patchy      Slight
## 36            Patchy      Slight
## 37  Moderately clean    Moderate
## 38  Moderately clean    Moderate
## 39  Moderately clean    Moderate
## 40  Moderately clean    Moderate
## 41  Moderately clean      Severe
## 42  Moderately clean      Severe
## 43  Moderately clean    Moderate
## 44  Moderately clean    Moderate
## 45  Moderately clean    Moderate
## 46  Moderately clean    Moderate
## 47             Clean Very severe
## 48             Clean Very severe
## 49  Moderately clean      Slight
## 50  Moderately clean      Slight
## 51            Patchy      Slight
## 52            Patchy      Slight
## 53  Moderately clean      Severe
## 54  Moderately clean      Severe
## 55  Moderately clean      Severe
## 56  Moderately clean      Severe
## 57       Very patchy    Moderate
## 58       Very patchy    Moderate
## 59  Moderately clean    Moderate
## 60  Moderately clean    Moderate
## 61             Clean      Slight
## 62             Clean      Slight
## 63       Very patchy    Moderate
## 64       Very patchy    Moderate
## 65  Moderately clean        None
## 66  Moderately clean        None
## 67       Very patchy        None
## 68       Very patchy        None
## 69             Clean      Slight
## 70             Clean      Slight
## 71  Moderately clean      Severe
## 72  Moderately clean      Severe
## 73       Very patchy        None
## 74       Very patchy        None
## 75  Moderately clean      Severe
## 76  Moderately clean      Severe
## 77  Moderately clean    Moderate
## 78  Moderately clean    Moderate
## 79             Clean      Severe
## 80             Clean      Severe
## 81  Moderately clean    Moderate
## 82  Moderately clean    Moderate
## 83            Patchy    Moderate
## 84            Patchy    Moderate
## 85  Moderately clean    Moderate
## 86  Moderately clean    Moderate
## 87             Clean      Severe
## 88             Clean      Severe
## 89  Moderately clean    Moderate
## 90  Moderately clean    Moderate
## 91  Moderately clean    Moderate
## 92  Moderately clean    Moderate
## 93  Moderately clean    Moderate
## 94  Moderately clean    Moderate
## 95             Clean Very severe
## 96             Clean Very severe
## 97  Moderately clean    Moderate
## 98  Moderately clean    Moderate
## 99             Clean      Slight
## 100            Clean      Slight
## 101            Clean Very severe
## 102            Clean Very severe
## 103 Moderately clean    Moderate
## 104 Moderately clean    Moderate
## 105       Clean burn    Moderate
## 106       Clean burn    Moderate
## 107 Moderately clean    Moderate
## 108 Moderately clean    Moderate
## 109           Patchy    Moderate
## 110           Patchy    Moderate
## 111 Moderately clean    Moderate
## 112 Moderately clean    Moderate
## 113            Clean Very severe
## 114            Clean Very severe
## 115            Clean      Severe
## 116            Clean      Severe
## 117 Moderately clean    Moderate
## 118 Moderately clean    Moderate
## 119 Moderately clean      Slight
## 120 Moderately clean      Slight
## 121           Patchy      Slight
## 122           Patchy      Slight
## 123 Moderately clean    Moderate
## 124 Moderately clean    Moderate
## 125           Patchy    Moderate
## 126           Patchy    Moderate
##                                                                                                                                                                                                                                                             GENERALOBS
## 1                                                                                                                                                                                                       First patch burn since the introduction of the new fire policy
## 2                                                                                                                                                                                                       First patch burn since the introduction of the new fire policy
## 3                                                                                                                                                                                                                                                                 <NA>
## 4                                                                                                                                                                                                                                                                 <NA>
## 5                                                                                                                                                                                                                                                                 <NA>
## 6                                                                                                                                                                                                                                                                 <NA>
## 7                                                                                                                                                                                                                                                       Second attempt
## 8                                                                                                                                                                                                                                                       Second attempt
## 9                                                                                                                                                                                                                                                                 <NA>
## 10                                                                                                                                                                                                                                                                <NA>
## 11                                                                                                                                                                                                                                                                <NA>
## 12                                                                                                                                                                                                                                                                <NA>
## 13                                                                                                                                                                                                                   Please see original fire report for further notes
## 14                                                                                                                                                                                                                   Please see original fire report for further notes
## 15                                                                                                                                                                                                                                                                <NA>
## 16                                                                                                                                                                                                                                                                <NA>
## 17                                                                                                                                                                                                                                                                <NA>
## 18                                                                                                                                                                                                                                                                <NA>
## 19                                                                                                                                                                                                                                                                <NA>
## 20                                                                                                                                                                                                                                                                <NA>
## 21                                                                                                                                                                                                                                                     Nhlangwine camp
## 22                                                                                                                                                                                                                                                     Nhlangwine camp
## 23                                                                                                                                                                                                                                                           Roan Camp
## 24                                                                                                                                                                                                                                                           Roan Camp
## 25                                                                                                                                                                                                                                                                <NA>
## 26                                                                                                                                                                                                                                                                <NA>
## 27                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 28                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 29                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 30                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 31                                                                                                                                                                                                                                                                <NA>
## 32                                                                                                                                                                                                                                                                <NA>
## 33                                                                                                                                                                                                                                                          See report
## 34                                                                                                                                                                                                                                                          See report
## 35                                                                                                                                                                                                                                                                <NA>
## 36                                                                                                                                                                                                                                                                <NA>
## 37                                                                                                                                                                                                                                                                <NA>
## 38                                                                                                                                                                                                                                                                <NA>
## 39                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 40                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 41                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 42                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 43                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 44                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 45                                                                                                                                                                                                                                                                <NA>
## 46                                                                                                                                                                                                                                                                <NA>
## 47                                                                                                                                                                                                                                                       Very hot burn
## 48                                                                                                                                                                                                                                                       Very hot burn
## 49                                                                                                                                                                                                                                                                <NA>
## 50                                                                                                                                                                                                                                                                <NA>
## 51                                                                                                                                                                                                                                                                <NA>
## 52                                                                                                                                                                                                                                                                <NA>
## 53                                                                                                                                                                                                                                                                <NA>
## 54                                                                                                                                                                                                                                                                <NA>
## 55                                                                                                                                                                                                                                                     Concession burn
## 56                                                                                                                                                                                                                                                     Concession burn
## 57                                                                                                           The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 58                                                                                                           The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 59                                                                                                                                                                                                                                                                <NA>
## 60                                                                                                                                                                                                                                                                <NA>
## 61                                                                                                                                                                                     The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 62                                                                                                                                                                                     The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 63                                                                                                                                                                                                                                                                <NA>
## 64                                                                                                                                                                                                                                                                <NA>
## 65                                                                                                                                                                                                                   Please see original fire report for further notes
## 66                                                                                                                                                                                                                   Please see original fire report for further notes
## 67                                                                                                                                                                                                                   Please see original fire report for further notes
## 68                                                                                                                                                                                                                   Please see original fire report for further notes
## 69                                                                                                                                                                           A number of plots catched fire the fire that started from block S042 with westerly winds.
## 70                                                                                                                                                                           A number of plots catched fire the fire that started from block S042 with westerly winds.
## 71                                                                                                                                                                                                         Very high litter load\r\nBlock did not burn for a few years
## 72                                                                                                                                                                                                         Very high litter load\r\nBlock did not burn for a few years
## 73                                                                                                                                                                                                                                                    Good Patchy burn
## 74                                                                                                                                                                                                                                                    Good Patchy burn
## 75                                                                                                                                                                                                                                                                <NA>
## 76                                                                                                                                                                                                                                                                <NA>
## 77                                                                                                                                             Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 78                                                                                                                                             Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 79                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 80                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 81                                                                                                                                                                                                                                                                <NA>
## 82                                                                                                                                                                                                                                                                <NA>
## 83                                                                                                                                                                                                                                                                <NA>
## 84                                                                                                                                                                                                                                                                <NA>
## 85                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 86                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 87                                                                                                                                                                                                                                                                <NA>
## 88                                                                                                                                                                                                                                                                <NA>
## 89                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 90                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 91                                                                                                                                                                                                                                                                <NA>
## 92                                                                                                                                                                                                                                                                <NA>
## 93                                                                                                                                                                                                                                                                <NA>
## 94                                                                                                                                                                                                                                                                <NA>
## 95                                                                                                                                                Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 96                                                                                                                                                Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 97                                                                                                                                                                                                                                                                <NA>
## 98                                                                                                                                                                                                                                                                <NA>
## 99                                                                                                                                                                                                                               The fire burnt along the railway line
## 100                                                                                                                                                                                                                              The fire burnt along the railway line
## 101                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 102                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 103 Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 104 Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 105                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 106                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 107                                                                                                                                                                                                                                                               <NA>
## 108                                                                                                                                                                                                                                                               <NA>
## 109                                                                                                                                       The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 110                                                                                                                                       The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 111                                                                                                                                                                                                                                                               <NA>
## 112                                                                                                                                                                                                                                                               <NA>
## 113                                                                                                                                                                                                                  Please see original fire report for further notes
## 114                                                                                                                                                                                                                  Please see original fire report for further notes
## 115                                                                                                                                                                                                                                                               <NA>
## 116                                                                                                                                                                                                                                                               <NA>
## 117                                                                                                                                                                                                                                                               <NA>
## 118                                                                                                                                                                                                                                                               <NA>
## 119                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 120                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 121                                                                                                                                                                                                                                                               <NA>
## 122                                                                                                                                                                                                                                                               <NA>
## 123                                                                                                                                                                                                                                                               <NA>
## 124                                                                                                                                                                                                                                                               <NA>
## 125                                                                                                                                                                                                                                                               <NA>
## 126                                                                                                                                                                                                                                                               <NA>
##      AREA_HA DATE_START DATE_END FIRE_YEAR WoodyCover   MAP Station Year
## 1   10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2002
## 2   10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2001
## 3    1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2001
## 4    1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2002
## 5     149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2002
## 6     149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2001
## 7    7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2001
## 8    7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2002
## 9      43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2001
## 10     43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2002
## 11   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2001
## 12   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2002
## 13     16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2001
## 14     16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2002
## 15   5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2001
## 16   5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2002
## 17   1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2002
## 18   1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2001
## 19    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2002
## 20    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2001
## 21     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2002
## 22     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2001
## 23     62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2001
## 24     62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2002
## 25   6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2001
## 26   6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2002
## 27  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2002
## 28  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2001
## 29    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2001
## 30    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2002
## 31   2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2002
## 32   2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2001
## 33    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2001
## 34    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2002
## 35    143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2002
## 36    143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2001
## 37    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2002
## 38    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2001
## 39    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2002
## 40    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2001
## 41   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2002
## 42   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2001
## 43   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2002
## 44   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2001
## 45    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2001
## 46    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2002
## 47   7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2001
## 48   7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2002
## 49    184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2001
## 50    184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2002
## 51   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2002
## 52   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2001
## 53  17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2001
## 54  17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2002
## 55   2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2001
## 56   2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2002
## 57  53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2002
## 58  53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2001
## 59   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2001
## 60   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2002
## 61   2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2002
## 62   2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2001
## 63   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2001
## 64   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2002
## 65   6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2001
## 66   6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2002
## 67   2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2002
## 68   2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2001
## 69   4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2002
## 70   4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2001
## 71   8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2001
## 72   8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2002
## 73    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2002
## 74    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2001
## 75   3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2001
## 76   3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2002
## 77   2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2002
## 78   2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2001
## 79  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2002
## 80  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2001
## 81   1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2001
## 82   1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2002
## 83   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2001
## 84   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2002
## 85  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2002
## 86  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2001
## 87  17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2001
## 88  17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2002
## 89   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2001
## 90   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2002
## 91    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2002
## 92    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2001
## 93   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2002
## 94   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2001
## 95    228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2001
## 96    228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2002
## 97    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2002
## 98    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2001
## 99    157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2002
## 100   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2001
## 101  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2001
## 102  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2002
## 103 21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2002
## 104 21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2001
## 105 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2002
## 106 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2001
## 107  2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2002
## 108  2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2001
## 109  3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2001
## 110  3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2002
## 111  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2001
## 112  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2002
## 113 17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2001
## 114 17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2002
## 115  7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2001
## 116  7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2002
## 117  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2001
## 118  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2002
## 119  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2002
## 120  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2001
## 121   120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2001
## 122   120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2002
## 123   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2002
## 124   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2001
## 125 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2001
## 126 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2002
##     AnnualPrecip PreviousYears
## 1          265.4          2000
## 2          583.6          2000
## 3          738.3          2000
## 4          308.6          2000
## 5          258.4          2000
## 6          627.0          2000
## 7          663.1          2000
## 8          135.3          2000
## 9          886.5          2000
## 10         397.2          2000
## 11         632.9          2000
## 12         291.4          2000
## 13         886.5          2000
## 14         397.2          2000
## 15         886.5          2000
## 16         397.2          2000
## 17         164.1          2000
## 18         505.5          2000
## 19         217.5          2000
## 20         631.6          2000
## 21         460.9          2000
## 22         608.8          2000
## 23         610.8          2000
## 24         218.9          2000
## 25         610.8          2000
## 26         218.9          2000
## 27         436.3          2000
## 28         853.8          2000
## 29         534.3          2000
## 30         393.4          2000
## 31         203.7          2000
## 32         720.5          2000
## 33         663.1          2000
## 34         135.3          2000
## 35         203.7          2000
## 36         720.5          2000
## 37         258.4          2000
## 38         627.0          2000
## 39         217.5          2000
## 40         631.6          2000
## 41         436.3          2000
## 42         853.8          2000
## 43         436.3          2000
## 44         853.8          2000
## 45         632.9          2000
## 46         291.4          2000
## 47         610.8          2000
## 48         218.9          2000
## 49         738.3          2000
## 50         308.6          2000
## 51         460.9          2000
## 52         608.8          2000
## 53         610.8          2000
## 54         218.9          2000
## 55         738.3          2000
## 56         308.6          2000
## 57         318.9          2000
## 58         719.6          2000
## 59         632.9          2000
## 60         291.4          2000
## 61         460.9          2000
## 62         608.8          2000
## 63         632.9          2000
## 64         291.4          2000
## 65         886.5          2000
## 66         397.2          2000
## 67         203.7          2000
## 68         720.5          2000
## 69         460.9          2000
## 70         608.8          2000
## 71         610.8          2000
## 72         218.9          2000
## 73         217.5          2000
## 74         631.6          2000
## 75         610.8          2000
## 76         218.9          2000
## 77         203.7          2000
## 78         720.5          2000
## 79         436.3          2000
## 80         853.8          2000
## 81         726.8          2000
## 82         295.6          2000
## 83         738.3          2000
## 84         308.6          2000
## 85         217.5          2000
## 86         631.6          2000
## 87         738.3          2000
## 88         308.6          2000
## 89         632.9          2000
## 90         291.4          2000
## 91         217.5          2000
## 92         631.6          2000
## 93         217.5          2000
## 94         631.6          2000
## 95         591.2          2000
## 96         116.5          2000
## 97         217.5          2000
## 98         631.6          2000
## 99         460.9          2000
## 100        608.8          2000
## 101        591.2          2000
## 102        116.5          2000
## 103        265.4          2000
## 104        583.6          2000
## 105        460.9          2000
## 106        608.8          2000
## 107        217.5          2000
## 108        631.6          2000
## 109        711.5          2000
## 110        163.3          2000
## 111        610.8          2000
## 112        218.9          2000
## 113        534.3          2000
## 114        393.4          2000
## 115        711.5          2000
## 116        163.3          2000
## 117        610.8          2000
## 118        218.9          2000
## 119        203.7          2000
## 120        720.5          2000
## 121        632.9          2000
## 122        291.4          2000
## 123        164.1          2000
## 124        505.5          2000
## 125        610.8          2000
## 126        218.9          2000
## 
## > rm(extractedFireScars)
## 
## > rm(extractedFireScars_df)
## 
## > (extractedFireScars_wx_df)
##            IGNITIONSE       FIREID  STARTDATE      CAUSE     AGENT
## 1    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 2    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 3    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 4    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 5    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 6    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 7    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 8    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 9    Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 10   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 11   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 12   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 13   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 14   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 15   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 16   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 17   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 18   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 19   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 20   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 21   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 22   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 23   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 24   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 25   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 26   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 27   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 28   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 29   Crocodile Bridge 20020711KRO1 2002-07-11      Arson Neighbour
## 30   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 31   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 32   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 33   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 34   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 35   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 36   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 37   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 38   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 39   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 40   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 41   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 42   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 43   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 44   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 45   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 46   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 47   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 48   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 49   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 50   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 51   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 52   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 53   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 54   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 55   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 56   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 57   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 58   Crocodile Bridge 20020827KRO1 2002-08-27      Arson Immigrant
## 59   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 60   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 61   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 62   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 63   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 64   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 65   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 66   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 67   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 68   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 69   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 70   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 71   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 72   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 73   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 74   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 75   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 76   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 77   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 78   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 79   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 80   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 81   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 82   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 83   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 84   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 85   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 86   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 87   Crocodile Bridge 20020718KRO1 2002-07-18 Management    Ranger
## 88   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 89   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 90   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 91   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 92   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 93   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 94   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 95   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 96   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 97   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 98   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 99   Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 100  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 101  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 102  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 103  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 104  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 105  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 106  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 107  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 108  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 109  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 110  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 111  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 112  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 113  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 114  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 115  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 116  Crocodile Bridge 20020825KRO1 2002-08-25      Arson Immigrant
## 117  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 118  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 119  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 120  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 121  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 122  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 123  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 124  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 125  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 126  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 127  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 128  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 129  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 130  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 131  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 132  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 133  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 134  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 135  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 136  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 137  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 138  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 139  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 140  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 141  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 142  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 143  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 144  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 145  Crocodile Bridge 20020417KRO1 2002-04-17 Management    Ranger
## 146            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 147            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 148            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 149            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 150            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 151            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 152            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 153            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 154            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 155            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 156            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 157            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 158            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 159            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 160            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 161            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 162            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 163            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 164            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 165            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 166            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 167            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 168            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 169            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 170            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 171            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 172            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 173            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 174            Letaba 20020524LET1 2002-05-24      Arson Immigrant
## 175            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 176            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 177            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 178            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 179            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 180            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 181            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 182            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 183            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 184            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 185            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 186            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 187            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 188            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 189            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 190            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 191            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 192            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 193            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 194            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 195            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 196            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 197            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 198            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 199            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 200            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 201            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 202            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 203            Letaba 20020502LET1 2002-05-02 Management    Ranger
## 204       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 205       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 206       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 207       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 208       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 209       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 210       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 211       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 212       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 213       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 214       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 215       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 216       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 217       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 218       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 219       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 220       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 221       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 222       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 223       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 224       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 225       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 226       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 227       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 228       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 229       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 230       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 231       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 232       Lower Sabie 20020703OSA1 2002-07-03 Management    Ranger
## 233       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 234       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 235       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 236       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 237       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 238       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 239       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 240       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 241       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 242       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 243       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 244       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 245       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 246       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 247       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 248       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 249       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 250       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 251       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 252       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 253       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 254       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 255       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 256       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 257       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 258       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 259       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 260       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 261       Lower Sabie 20020705OSA1 2002-07-05 Management    Ranger
## 262       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 263       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 264       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 265       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 266       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 267       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 268       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 269       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 270       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 271       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 272       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 273       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 274       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 275       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 276       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 277       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 278       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 279       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 280       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 281       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 282       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 283       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 284       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 285       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 286       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 287       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 288       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 289       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 290       Lower Sabie 20020522OSA1 2002-05-22 Management    Ranger
## 291       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 292       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 293       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 294       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 295       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 296       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 297       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 298       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 299       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 300       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 301       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 302       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 303       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 304       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 305       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 306       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 307       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 308       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 309       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 310       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 311       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 312       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 313       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 314       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 315       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 316       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 317       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 318       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 319       Lower Sabie 20020822OSA1 2002-08-22      Arson   Tourist
## 320        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 321        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 322        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 323        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 324        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 325        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 326        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 327        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 328        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 329        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 330        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 331        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 332        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 333        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 334        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 335        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 336        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 337        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 338        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 339        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 340        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 341        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 342        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 343        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 344        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 345        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 346        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 347        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 348        Mahlangeni 20020530MAH1 2002-05-30      Arson Immigrant
## 349        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 350        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 351        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 352        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 353        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 354        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 355        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 356        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 357        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 358        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 359        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 360        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 361        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 362        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 363        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 364        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 365        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 366        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 367        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 368        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 369        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 370        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 371        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 372        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 373        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 374        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 375        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 376        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 377        Mahlangeni 20020424MAH1 2002-04-24 Management    Ranger
## 378          Malelane 20021006MAL1 2002-10-06      Other     Other
## 379          Malelane 20021006MAL1 2002-10-06      Other     Other
## 380          Malelane 20021006MAL1 2002-10-06      Other     Other
## 381          Malelane 20021006MAL1 2002-10-06      Other     Other
## 382          Malelane 20021006MAL1 2002-10-06      Other     Other
## 383          Malelane 20021006MAL1 2002-10-06      Other     Other
## 384          Malelane 20021006MAL1 2002-10-06      Other     Other
## 385          Malelane 20021006MAL1 2002-10-06      Other     Other
## 386          Malelane 20021006MAL1 2002-10-06      Other     Other
## 387          Malelane 20021006MAL1 2002-10-06      Other     Other
## 388          Malelane 20021006MAL1 2002-10-06      Other     Other
## 389          Malelane 20021006MAL1 2002-10-06      Other     Other
## 390          Malelane 20021006MAL1 2002-10-06      Other     Other
## 391          Malelane 20021006MAL1 2002-10-06      Other     Other
## 392          Malelane 20021006MAL1 2002-10-06      Other     Other
## 393          Malelane 20021006MAL1 2002-10-06      Other     Other
## 394          Malelane 20021006MAL1 2002-10-06      Other     Other
## 395          Malelane 20021006MAL1 2002-10-06      Other     Other
## 396          Malelane 20021006MAL1 2002-10-06      Other     Other
## 397          Malelane 20021006MAL1 2002-10-06      Other     Other
## 398          Malelane 20021006MAL1 2002-10-06      Other     Other
## 399          Malelane 20021006MAL1 2002-10-06      Other     Other
## 400          Malelane 20021006MAL1 2002-10-06      Other     Other
## 401          Malelane 20021006MAL1 2002-10-06      Other     Other
## 402          Malelane 20021006MAL1 2002-10-06      Other     Other
## 403          Malelane 20021006MAL1 2002-10-06      Other     Other
## 404          Malelane 20021006MAL1 2002-10-06      Other     Other
## 405          Malelane 20021006MAL1 2002-10-06      Other     Other
## 406          Malelane 20021006MAL1 2002-10-06      Other     Other
## 407          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 408          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 409          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 410          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 411          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 412          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 413          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 414          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 415          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 416          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 417          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 418          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 419          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 420          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 421          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 422          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 423          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 424          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 425          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 426          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 427          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 428          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 429          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 430          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 431          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 432          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 433          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 434          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 435          Malelane 20020523MAL1 2002-05-23      Arson     Other
## 436         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 437         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 438         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 439         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 440         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 441         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 442         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 443         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 444         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 445         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 446         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 447         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 448         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 449         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 450         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 451         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 452         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 453         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 454         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 455         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 456         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 457         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 458         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 459         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 460         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 461         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 462         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 463         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 464         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 465         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 466         Mooiplaas 20020408MOO1 2002-04-08 Management    Ranger
## 467         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 468         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 469         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 470         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 471         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 472         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 473         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 474         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 475         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 476         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 477         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 478         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 479         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 480         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 481         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 482         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 483         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 484         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 485         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 486         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 487         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 488         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 489         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 490         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 491         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 492         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 493         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 494         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 495         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 496         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 497         Mooiplaas 20020922MOO1 2002-09-22      Arson Immigrant
## 498          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 499          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 500          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 501          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 502          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 503          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 504          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 505          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 506          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 507          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 508          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 509          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 510          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 511          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 512          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 513          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 514          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 515          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 516          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 517          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 518          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 519          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 520          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 521          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 522          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 523          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 524          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 525          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 526          Olifants 20020503OLI1 2002-05-03 Management    Ranger
## 527          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 528          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 529          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 530          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 531          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 532          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 533          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 534          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 535          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 536          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 537          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 538          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 539          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 540          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 541          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 542          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 543          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 544          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 545          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 546          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 547          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 548          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 549          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 550          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 551          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 552          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 553          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 554          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 555          Olifants 20020830OLI1 2002-08-30      Arson Immigrant
## 556          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 557          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 558          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 559          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 560          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 561          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 562          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 563          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 564          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 565          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 566          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 567          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 568          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 569          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 570          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 571          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 572          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 573          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 574          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 575          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 576          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 577          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 578          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 579          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 580          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 581          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 582          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 583          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 584          Olifants 20020707OLI1 2002-07-07      Arson Immigrant
## 585          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 586          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 587          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 588          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 589          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 590          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 591          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 592          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 593          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 594          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 595          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 596          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 597          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 598          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 599          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 600          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 601          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 602          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 603          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 604          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 605          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 606          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 607          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 608          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 609          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 610          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 611          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 612          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 613          Olifants 20021209OLI1 2002-12-09      Arson Immigrant
## 614          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 615          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 616          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 617          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 618          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 619          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 620          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 621          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 622          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 623          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 624          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 625          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 626          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 627          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 628          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 629          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 630          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 631          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 632          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 633          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 634          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 635          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 636          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 637          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 638          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 639          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 640          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 641          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 642          Olifants 20020727OLI1 2002-07-27      Arson Immigrant
## 643          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 644          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 645          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 646          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 647          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 648          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 649          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 650          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 651          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 652          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 653          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 654          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 655          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 656          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 657          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 658          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 659          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 660          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 661          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 662          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 663          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 664          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 665          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 666          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 667          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 668          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 669          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 670          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 671          Olifants 20020731OLI1 2002-07-31      Arson Immigrant
## 672            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 673            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 674            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 675            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 676            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 677            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 678            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 679            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 680            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 681            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 682            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 683            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 684            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 685            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 686            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 687            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 688            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 689            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 690            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 691            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 692            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 693            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 694            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 695            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 696            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 697            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 698            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 699            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 700            Pafuri 20020816PAF1 2002-08-16      Arson Immigrant
## 701            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 702            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 703            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 704            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 705            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 706            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 707            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 708            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 709            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 710            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 711            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 712            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 713            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 714            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 715            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 716            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 717            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 718            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 719            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 720            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 721            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 722            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 723            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 724            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 725            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 726            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 727            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 728            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 729            Pafuri 20020902PAF1 2002-09-02      Arson Immigrant
## 730            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 731            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 732            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 733            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 734            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 735            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 736            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 737            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 738            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 739            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 740            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 741            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 742            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 743            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 744            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 745            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 746            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 747            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 748            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 749            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 750            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 751            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 752            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 753            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 754            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 755            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 756            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 757            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 758            Pafuri 20020826PAF1 2002-08-26      Arson Immigrant
## 759            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 760            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 761            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 762            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 763            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 764            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 765            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 766            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 767            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 768            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 769            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 770            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 771            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 772            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 773            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 774            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 775            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 776            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 777            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 778            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 779            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 780            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 781            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 782            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 783            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 784            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 785            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 786            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 787            Pafuri 20020925PAF1 2002-09-25      Arson Neighbour
## 788            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 789            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 790            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 791            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 792            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 793            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 794            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 795            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 796            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 797            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 798            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 799            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 800            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 801            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 802            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 803            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 804            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 805            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 806            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 807            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 808            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 809            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 810            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 811            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 812            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 813            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 814            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 815            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 816            Pafuri 20020511PAF1 2002-05-11      Arson Immigrant
## 817            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 818            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 819            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 820            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 821            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 822            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 823            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 824            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 825            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 826            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 827            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 828            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 829            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 830            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 831            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 832            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 833            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 834            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 835            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 836            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 837            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 838            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 839            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 840            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 841            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 842            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 843            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 844            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 845            Pafuri 20020904PAF1 2002-09-04      Arson Immigrant
## 846            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 847            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 848            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 849            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 850            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 851            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 852            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 853            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 854            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 855            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 856            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 857            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 858            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 859            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 860            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 861            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 862            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 863            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 864            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 865            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 866            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 867            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 868            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 869            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 870            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 871            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 872            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 873            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 874            Pafuri 20020913PAF1 2002-09-13      Arson Neighbour
## 875            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 876            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 877            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 878            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 879            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 880            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 881            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 882            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 883            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 884            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 885            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 886            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 887            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 888            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 889            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 890            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 891            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 892            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 893            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 894            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 895            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 896            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 897            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 898            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 899            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 900            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 901            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 902            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 903            Pafuri 20020601PAF1 2002-06-01      Arson Neighbour
## 904        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 905        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 906        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 907        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 908        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 909        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 910        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 911        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 912        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 913        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 914        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 915        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 916        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 917        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 918        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 919        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 920        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 921        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 922        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 923        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 924        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 925        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 926        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 927        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 928        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 929        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 930        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 931        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 932        Phalaborwa 20021018PHA1 2002-10-18      Arson Immigrant
## 933        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 934        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 935        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 936        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 937        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 938        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 939        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 940        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 941        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 942        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 943        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 944        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 945        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 946        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 947        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 948        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 949        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 950        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 951        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 952        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 953        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 954        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 955        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 956        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 957        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 958        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 959        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 960        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 961        Phalaborwa 20020925PHA1 2002-09-25      Arson Immigrant
## 962      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 963      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 964      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 965      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 966      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 967      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 968      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 969      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 970      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 971      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 972      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 973      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 974      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 975      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 976      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 977      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 978      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 979      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 980      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 981      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 982      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 983      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 984      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 985      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 986      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 987      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 988      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 989      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 990      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 991      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 992      Pretoriuskop 20020513PRE1 2002-05-13 Management    Ranger
## 993      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 994      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 995      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 996      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 997      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 998      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 999      Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1000     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1001     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1002     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1003     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1004     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1005     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1006     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1007     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1008     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1009     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1010     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1011     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1012     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1013     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1014     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1015     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1016     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1017     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1018     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1019     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1020     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1021     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1022     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1023     Pretoriuskop 20020716PRE1 2002-07-16      Arson Neighbour
## 1024     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1025     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1026     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1027     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1028     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1029     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1030     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1031     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1032     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1033     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1034     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1035     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1036     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1037     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1038     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1039     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1040     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1041     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1042     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1043     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1044     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1045     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1046     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1047     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1048     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1049     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1050     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1051     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1052     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1053     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1054     Pretoriuskop 20020812PRE1 2002-08-12      Arson   Poacher
## 1055     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1056     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1057     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1058     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1059     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1060     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1061     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1062     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1063     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1064     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1065     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1066     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1067     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1068     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1069     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1070     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1071     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1072     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1073     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1074     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1075     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1076     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1077     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1078     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1079     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1080     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1081     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1082     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1083     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1084     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1085     Pretoriuskop 20020924PRE1 2002-09-24      Arson     Other
## 1086     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1087     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1088     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1089     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1090     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1091     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1092     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1093     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1094     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1095     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1096     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1097     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1098     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1099     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1100     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1101     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1102     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1103     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1104     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1105     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1106     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1107     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1108     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1109     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1110     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1111     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1112     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1113     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1114     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1115     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1116     Pretoriuskop 20020913PRE1 2002-09-13      Arson Neighbour
## 1117     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1118     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1119     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1120     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1121     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1122     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1123     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1124     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1125     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1126     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1127     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1128     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1129     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1130     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1131     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1132     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1133     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1134     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1135     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1136     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1137     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1138     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1139     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1140     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1141     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1142     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1143     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1144     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1145     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1146     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1147     Pretoriuskop 20020729PRE1 2002-07-29      Other    Nature
## 1148      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1149      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1150      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1151      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1152      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1153      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1154      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1155      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1156      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1157      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1158      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1159      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1160      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1161      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1162      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1163      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1164      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1165      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1166      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1167      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1168      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1169      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1170      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1171      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1172      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1173      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1174      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1175      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1176      Punda Maria 20020502PUN1 2002-05-02 Management    Ranger
## 1177      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1178      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1179      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1180      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1181      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1182      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1183      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1184      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1185      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1186      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1187      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1188      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1189      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1190      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1191      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1192      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1193      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1194      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1195      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1196      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1197      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1198      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1199      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1200      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1201      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1202      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1203      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1204      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1205      Punda Maria 20020802PUN1 2002-08-02      Arson Immigrant
## 1206      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1207      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1208      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1209      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1210      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1211      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1212      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1213      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1214      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1215      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1216      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1217      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1218      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1219      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1220      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1221      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1222      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1223      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1224      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1225      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1226      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1227      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1228      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1229      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1230      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1231      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1232      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1233      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1234      Punda Maria 20020507PUN1 2002-05-07 Management    Ranger
## 1235      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1236      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1237      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1238      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1239      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1240      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1241      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1242      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1243      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1244      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1245      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1246      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1247      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1248      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1249      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1250      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1251      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1252      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1253      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1254      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1255      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1256      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1257      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1258      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1259      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1260      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1261      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1262      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1263      Punda Maria 20020506PUN1 2002-05-06 Management    Ranger
## 1264           Satara 20020911SAT1 2002-09-11      Other     Other
## 1265           Satara 20020911SAT1 2002-09-11      Other     Other
## 1266           Satara 20020911SAT1 2002-09-11      Other     Other
## 1267           Satara 20020911SAT1 2002-09-11      Other     Other
## 1268           Satara 20020911SAT1 2002-09-11      Other     Other
## 1269           Satara 20020911SAT1 2002-09-11      Other     Other
## 1270           Satara 20020911SAT1 2002-09-11      Other     Other
## 1271           Satara 20020911SAT1 2002-09-11      Other     Other
## 1272           Satara 20020911SAT1 2002-09-11      Other     Other
## 1273           Satara 20020911SAT1 2002-09-11      Other     Other
## 1274           Satara 20020911SAT1 2002-09-11      Other     Other
## 1275           Satara 20020911SAT1 2002-09-11      Other     Other
## 1276           Satara 20020911SAT1 2002-09-11      Other     Other
## 1277           Satara 20020911SAT1 2002-09-11      Other     Other
## 1278           Satara 20020911SAT1 2002-09-11      Other     Other
## 1279           Satara 20020911SAT1 2002-09-11      Other     Other
## 1280           Satara 20020911SAT1 2002-09-11      Other     Other
## 1281           Satara 20020911SAT1 2002-09-11      Other     Other
## 1282           Satara 20020911SAT1 2002-09-11      Other     Other
## 1283           Satara 20020911SAT1 2002-09-11      Other     Other
## 1284           Satara 20020911SAT1 2002-09-11      Other     Other
## 1285           Satara 20020911SAT1 2002-09-11      Other     Other
## 1286           Satara 20020911SAT1 2002-09-11      Other     Other
## 1287           Satara 20020911SAT1 2002-09-11      Other     Other
## 1288           Satara 20020911SAT1 2002-09-11      Other     Other
## 1289           Satara 20020911SAT1 2002-09-11      Other     Other
## 1290           Satara 20020911SAT1 2002-09-11      Other     Other
## 1291           Satara 20020911SAT1 2002-09-11      Other     Other
## 1292           Satara 20020911SAT1 2002-09-11      Other     Other
## 1293           Satara 20020911SAT1 2002-09-11      Other     Other
## 1294           Satara 20020911SAT1 2002-09-11      Other     Other
## 1295           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1296           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1297           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1298           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1299           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1300           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1301           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1302           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1303           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1304           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1305           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1306           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1307           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1308           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1309           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1310           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1311           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1312           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1313           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1314           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1315           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1316           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1317           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1318           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1319           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1320           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1321           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1322           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1323           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1324           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1325           Satara 20020918SAT1 2002-09-18      Arson Immigrant
## 1326         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1327         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1328         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1329         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1330         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1331         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1332         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1333         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1334         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1335         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1336         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1337         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1338         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1339         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1340         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1341         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1342         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1343         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1344         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1345         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1346         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1347         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1348         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1349         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1350         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1351         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1352         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1353         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1354         Shangoni 20020822SHA1 2002-08-22 Management    Ranger
## 1355       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1356       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1357       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1358       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1359       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1360       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1361       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1362       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1363       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1364       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1365       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1366       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1367       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1368       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1369       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1370       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1371       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1372       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1373       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1374       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1375       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1376       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1377       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1378       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1379       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1380       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1381       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1382       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1383       Shingwedzi 20020509SHI1 2002-05-09 Management    Ranger
## 1384       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1385       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1386       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1387       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1388       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1389       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1390       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1391       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1392       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1393       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1394       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1395       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1396       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1397       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1398       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1399       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1400       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1401       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1402       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1403       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1404       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1405       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1406       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1407       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1408       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1409       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1410       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1411       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1412       Shingwedzi 20021222SHI1 2002-12-22      Arson Immigrant
## 1413          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1414          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1415          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1416          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1417          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1418          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1419          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1420          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1421          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1422          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1423          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1424          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1425          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1426          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1427          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1428          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1429          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1430          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1431          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1432          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1433          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1434          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1435          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1436          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1437          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1438          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1439          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1440          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1441          Skukuza 20020723SKZ1 2002-07-23 Management    Ranger
## 1442        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1443        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1444        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1445        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1446        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1447        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1448        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1449        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1450        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1451        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1452        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1453        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1454        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1455        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1456        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1457        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1458        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1459        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1460        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1461        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1462        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1463        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1464        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1465        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1466        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1467        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1468        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1469        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1470        Tshokwane 20020523TSH1 2002-05-23 Management    Ranger
## 1471        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1472        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1473        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1474        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1475        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1476        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1477        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1478        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1479        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1480        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1481        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1482        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1483        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1484        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1485        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1486        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1487        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1488        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1489        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1490        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1491        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1492        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1493        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1494        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1495        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1496        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1497        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1498        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1499        Tshokwane 20020524TSH1 2002-05-24      Arson Neighbour
## 1500        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1501        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1502        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1503        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1504        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1505        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1506        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1507        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1508        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1509        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1510        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1511        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1512        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1513        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1514        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1515        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1516        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1517        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1518        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1519        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1520        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1521        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1522        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1523        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1524        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1525        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1526        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1527        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1528        Tshokwane 20021123TSH1 2002-11-23  Lightning    Nature
## 1529        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1530        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1531        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1532        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1533        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1534        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1535        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1536        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1537        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1538        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1539        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1540        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1541        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1542        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1543        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1544        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1545        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1546        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1547        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1548        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1549        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1550        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1551        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1552        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1553        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1554        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1555        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1556        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1557        Tshokwane 20020802TSH1 2002-08-02 Management    Ranger
## 1558        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1559        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1560        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1561        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1562        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1563        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1564        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1565        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1566        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1567        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1568        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1569        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1570        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1571        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1572        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1573        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1574        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1575        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1576        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1577        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1578        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1579        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1580        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1581        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1582        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1583        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1584        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1585        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1586        Tshokwane 20020821TSH1 2002-08-21      Arson   Tourist
## 1587      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1588      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1589      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1590      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1591      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1592      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1593      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1594      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1595      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1596      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1597      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1598      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1599      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1600      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1601      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1602      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1603      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1604      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1605      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1606      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1607      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1608      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1609      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1610      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1611      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1612      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1613      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1614      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1615      Vlakteplaas 20020817VLA1 2002-08-17      Arson   Tourist
## 1616      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1617      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1618      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1619      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1620      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1621      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1622      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1623      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1624      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1625      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1626      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1627      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1628      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1629      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1630      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1631      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1632      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1633      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1634      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1635      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1636      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1637      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1638      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1639      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1640      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1641      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1642      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1643      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1644      Vlakteplaas 20020716VLA1 2002-07-16      Arson Immigrant
## 1645      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1646      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1647      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1648      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1649      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1650      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1651      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1652      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1653      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1654      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1655      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1656      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1657      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1658      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1659      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1660      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1661      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1662      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1663      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1664      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1665      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1666      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1667      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1668      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1669      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1670      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1671      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1672      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1673      Vlakteplaas 20020513VLA1 2002-05-13 Management    Ranger
## 1674      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1675      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1676      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1677      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1678      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1679      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1680      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1681      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1682      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1683      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1684      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1685      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1686      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1687      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1688      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1689      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1690      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1691      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1692      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1693      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1694      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1695      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1696      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1697      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1698      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1699      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1700      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1701      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1702      Vlakteplaas 20021224VLA1 2002-12-24  Lightning    Nature
## 1703      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1704      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1705      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1706      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1707      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1708      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1709      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1710      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1711      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1712      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1713      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1714      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1715      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1716      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1717      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1718      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1719      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1720      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1721      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1722      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1723      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1724      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1725      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1726      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1727      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1728      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1729      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1730      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1731      Vlakteplaas 20020708VLA1 2002-07-08      Arson Immigrant
## 1732      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1733      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1734      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1735      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1736      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1737      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1738      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1739      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1740      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1741      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1742      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1743      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1744      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1745      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1746      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1747      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1748      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1749      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1750      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1751      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1752      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1753      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1754      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1755      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1756      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1757      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1758      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1759      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1760      Vlakteplaas 20020926VLA1 2002-09-26      Arson Immigrant
## 1761      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1762      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1763      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1764      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1765      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1766      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1767      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1768      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1769      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1770      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1771      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1772      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1773      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1774      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1775      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1776      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1777      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1778      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1779      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1780      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1781      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1782      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1783      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1784      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1785      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1786      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1787      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1788      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1789      Vlakteplaas 20020814VLA1 2002-08-14      Arson Immigrant
## 1790      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1791      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1792      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1793      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1794      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1795      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1796      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1797      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1798      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1799      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1800      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1801      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1802      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1803      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1804      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1805      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1806      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1807      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1808      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1809      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1810      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1811      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1812      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1813      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1814      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1815      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1816      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1817      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1818      Vlakteplaas 20020519VLA1 2002-05-19      Arson Immigrant
## 1819      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1820      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1821      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1822      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1823      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1824      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1825      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1826      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1827      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1828      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1829      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1830      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1831      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1832      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1833      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1834      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1835      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1836      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1837      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1838      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1839      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1840      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1841      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1842      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1843      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1844      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1845      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1846      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
## 1847      Vlakteplaas 20021024VLA1 2002-10-24  Lightning    Nature
##            HERBACEOUS  WOODYIMPAC
## 1    Moderately clean      Slight
## 2    Moderately clean      Slight
## 3    Moderately clean      Slight
## 4    Moderately clean      Slight
## 5    Moderately clean      Slight
## 6    Moderately clean      Slight
## 7    Moderately clean      Slight
## 8    Moderately clean      Slight
## 9    Moderately clean      Slight
## 10   Moderately clean      Slight
## 11   Moderately clean      Slight
## 12   Moderately clean      Slight
## 13   Moderately clean      Slight
## 14   Moderately clean      Slight
## 15   Moderately clean      Slight
## 16   Moderately clean      Slight
## 17   Moderately clean      Slight
## 18   Moderately clean      Slight
## 19   Moderately clean      Slight
## 20   Moderately clean      Slight
## 21   Moderately clean      Slight
## 22   Moderately clean      Slight
## 23   Moderately clean      Slight
## 24   Moderately clean      Slight
## 25   Moderately clean      Slight
## 26   Moderately clean      Slight
## 27   Moderately clean      Slight
## 28   Moderately clean      Slight
## 29   Moderately clean      Slight
## 30              Clean      Severe
## 31              Clean      Severe
## 32              Clean      Severe
## 33              Clean      Severe
## 34              Clean      Severe
## 35              Clean      Severe
## 36              Clean      Severe
## 37              Clean      Severe
## 38              Clean      Severe
## 39              Clean      Severe
## 40              Clean      Severe
## 41              Clean      Severe
## 42              Clean      Severe
## 43              Clean      Severe
## 44              Clean      Severe
## 45              Clean      Severe
## 46              Clean      Severe
## 47              Clean      Severe
## 48              Clean      Severe
## 49              Clean      Severe
## 50              Clean      Severe
## 51              Clean      Severe
## 52              Clean      Severe
## 53              Clean      Severe
## 54              Clean      Severe
## 55              Clean      Severe
## 56              Clean      Severe
## 57              Clean      Severe
## 58              Clean      Severe
## 59   Moderately clean      Severe
## 60   Moderately clean      Severe
## 61   Moderately clean      Severe
## 62   Moderately clean      Severe
## 63   Moderately clean      Severe
## 64   Moderately clean      Severe
## 65   Moderately clean      Severe
## 66   Moderately clean      Severe
## 67   Moderately clean      Severe
## 68   Moderately clean      Severe
## 69   Moderately clean      Severe
## 70   Moderately clean      Severe
## 71   Moderately clean      Severe
## 72   Moderately clean      Severe
## 73   Moderately clean      Severe
## 74   Moderately clean      Severe
## 75   Moderately clean      Severe
## 76   Moderately clean      Severe
## 77   Moderately clean      Severe
## 78   Moderately clean      Severe
## 79   Moderately clean      Severe
## 80   Moderately clean      Severe
## 81   Moderately clean      Severe
## 82   Moderately clean      Severe
## 83   Moderately clean      Severe
## 84   Moderately clean      Severe
## 85   Moderately clean      Severe
## 86   Moderately clean      Severe
## 87   Moderately clean      Severe
## 88             Patchy    Moderate
## 89             Patchy    Moderate
## 90             Patchy    Moderate
## 91             Patchy    Moderate
## 92             Patchy    Moderate
## 93             Patchy    Moderate
## 94             Patchy    Moderate
## 95             Patchy    Moderate
## 96             Patchy    Moderate
## 97             Patchy    Moderate
## 98             Patchy    Moderate
## 99             Patchy    Moderate
## 100            Patchy    Moderate
## 101            Patchy    Moderate
## 102            Patchy    Moderate
## 103            Patchy    Moderate
## 104            Patchy    Moderate
## 105            Patchy    Moderate
## 106            Patchy    Moderate
## 107            Patchy    Moderate
## 108            Patchy    Moderate
## 109            Patchy    Moderate
## 110            Patchy    Moderate
## 111            Patchy    Moderate
## 112            Patchy    Moderate
## 113            Patchy    Moderate
## 114            Patchy    Moderate
## 115            Patchy    Moderate
## 116            Patchy    Moderate
## 117            Patchy      Slight
## 118            Patchy      Slight
## 119            Patchy      Slight
## 120            Patchy      Slight
## 121            Patchy      Slight
## 122            Patchy      Slight
## 123            Patchy      Slight
## 124            Patchy      Slight
## 125            Patchy      Slight
## 126            Patchy      Slight
## 127            Patchy      Slight
## 128            Patchy      Slight
## 129            Patchy      Slight
## 130            Patchy      Slight
## 131            Patchy      Slight
## 132            Patchy      Slight
## 133            Patchy      Slight
## 134            Patchy      Slight
## 135            Patchy      Slight
## 136            Patchy      Slight
## 137            Patchy      Slight
## 138            Patchy      Slight
## 139            Patchy      Slight
## 140            Patchy      Slight
## 141            Patchy      Slight
## 142            Patchy      Slight
## 143            Patchy      Slight
## 144            Patchy      Slight
## 145            Patchy      Slight
## 146             Clean      Severe
## 147             Clean      Severe
## 148             Clean      Severe
## 149             Clean      Severe
## 150             Clean      Severe
## 151             Clean      Severe
## 152             Clean      Severe
## 153             Clean      Severe
## 154             Clean      Severe
## 155             Clean      Severe
## 156             Clean      Severe
## 157             Clean      Severe
## 158             Clean      Severe
## 159             Clean      Severe
## 160             Clean      Severe
## 161             Clean      Severe
## 162             Clean      Severe
## 163             Clean      Severe
## 164             Clean      Severe
## 165             Clean      Severe
## 166             Clean      Severe
## 167             Clean      Severe
## 168             Clean      Severe
## 169             Clean      Severe
## 170             Clean      Severe
## 171             Clean      Severe
## 172             Clean      Severe
## 173             Clean      Severe
## 174             Clean      Severe
## 175             Clean      Severe
## 176             Clean      Severe
## 177             Clean      Severe
## 178             Clean      Severe
## 179             Clean      Severe
## 180             Clean      Severe
## 181             Clean      Severe
## 182             Clean      Severe
## 183             Clean      Severe
## 184             Clean      Severe
## 185             Clean      Severe
## 186             Clean      Severe
## 187             Clean      Severe
## 188             Clean      Severe
## 189             Clean      Severe
## 190             Clean      Severe
## 191             Clean      Severe
## 192             Clean      Severe
## 193             Clean      Severe
## 194             Clean      Severe
## 195             Clean      Severe
## 196             Clean      Severe
## 197             Clean      Severe
## 198             Clean      Severe
## 199             Clean      Severe
## 200             Clean      Severe
## 201             Clean      Severe
## 202             Clean      Severe
## 203             Clean      Severe
## 204  Moderately clean      Severe
## 205  Moderately clean      Severe
## 206  Moderately clean      Severe
## 207  Moderately clean      Severe
## 208  Moderately clean      Severe
## 209  Moderately clean      Severe
## 210  Moderately clean      Severe
## 211  Moderately clean      Severe
## 212  Moderately clean      Severe
## 213  Moderately clean      Severe
## 214  Moderately clean      Severe
## 215  Moderately clean      Severe
## 216  Moderately clean      Severe
## 217  Moderately clean      Severe
## 218  Moderately clean      Severe
## 219  Moderately clean      Severe
## 220  Moderately clean      Severe
## 221  Moderately clean      Severe
## 222  Moderately clean      Severe
## 223  Moderately clean      Severe
## 224  Moderately clean      Severe
## 225  Moderately clean      Severe
## 226  Moderately clean      Severe
## 227  Moderately clean      Severe
## 228  Moderately clean      Severe
## 229  Moderately clean      Severe
## 230  Moderately clean      Severe
## 231  Moderately clean      Severe
## 232  Moderately clean      Severe
## 233  Moderately clean    Moderate
## 234  Moderately clean    Moderate
## 235  Moderately clean    Moderate
## 236  Moderately clean    Moderate
## 237  Moderately clean    Moderate
## 238  Moderately clean    Moderate
## 239  Moderately clean    Moderate
## 240  Moderately clean    Moderate
## 241  Moderately clean    Moderate
## 242  Moderately clean    Moderate
## 243  Moderately clean    Moderate
## 244  Moderately clean    Moderate
## 245  Moderately clean    Moderate
## 246  Moderately clean    Moderate
## 247  Moderately clean    Moderate
## 248  Moderately clean    Moderate
## 249  Moderately clean    Moderate
## 250  Moderately clean    Moderate
## 251  Moderately clean    Moderate
## 252  Moderately clean    Moderate
## 253  Moderately clean    Moderate
## 254  Moderately clean    Moderate
## 255  Moderately clean    Moderate
## 256  Moderately clean    Moderate
## 257  Moderately clean    Moderate
## 258  Moderately clean    Moderate
## 259  Moderately clean    Moderate
## 260  Moderately clean    Moderate
## 261  Moderately clean    Moderate
## 262  Moderately clean    Moderate
## 263  Moderately clean    Moderate
## 264  Moderately clean    Moderate
## 265  Moderately clean    Moderate
## 266  Moderately clean    Moderate
## 267  Moderately clean    Moderate
## 268  Moderately clean    Moderate
## 269  Moderately clean    Moderate
## 270  Moderately clean    Moderate
## 271  Moderately clean    Moderate
## 272  Moderately clean    Moderate
## 273  Moderately clean    Moderate
## 274  Moderately clean    Moderate
## 275  Moderately clean    Moderate
## 276  Moderately clean    Moderate
## 277  Moderately clean    Moderate
## 278  Moderately clean    Moderate
## 279  Moderately clean    Moderate
## 280  Moderately clean    Moderate
## 281  Moderately clean    Moderate
## 282  Moderately clean    Moderate
## 283  Moderately clean    Moderate
## 284  Moderately clean    Moderate
## 285  Moderately clean    Moderate
## 286  Moderately clean    Moderate
## 287  Moderately clean    Moderate
## 288  Moderately clean    Moderate
## 289  Moderately clean    Moderate
## 290  Moderately clean    Moderate
## 291             Clean      Severe
## 292             Clean      Severe
## 293             Clean      Severe
## 294             Clean      Severe
## 295             Clean      Severe
## 296             Clean      Severe
## 297             Clean      Severe
## 298             Clean      Severe
## 299             Clean      Severe
## 300             Clean      Severe
## 301             Clean      Severe
## 302             Clean      Severe
## 303             Clean      Severe
## 304             Clean      Severe
## 305             Clean      Severe
## 306             Clean      Severe
## 307             Clean      Severe
## 308             Clean      Severe
## 309             Clean      Severe
## 310             Clean      Severe
## 311             Clean      Severe
## 312             Clean      Severe
## 313             Clean      Severe
## 314             Clean      Severe
## 315             Clean      Severe
## 316             Clean      Severe
## 317             Clean      Severe
## 318             Clean      Severe
## 319             Clean      Severe
## 320  Moderately clean    Moderate
## 321  Moderately clean    Moderate
## 322  Moderately clean    Moderate
## 323  Moderately clean    Moderate
## 324  Moderately clean    Moderate
## 325  Moderately clean    Moderate
## 326  Moderately clean    Moderate
## 327  Moderately clean    Moderate
## 328  Moderately clean    Moderate
## 329  Moderately clean    Moderate
## 330  Moderately clean    Moderate
## 331  Moderately clean    Moderate
## 332  Moderately clean    Moderate
## 333  Moderately clean    Moderate
## 334  Moderately clean    Moderate
## 335  Moderately clean    Moderate
## 336  Moderately clean    Moderate
## 337  Moderately clean    Moderate
## 338  Moderately clean    Moderate
## 339  Moderately clean    Moderate
## 340  Moderately clean    Moderate
## 341  Moderately clean    Moderate
## 342  Moderately clean    Moderate
## 343  Moderately clean    Moderate
## 344  Moderately clean    Moderate
## 345  Moderately clean    Moderate
## 346  Moderately clean    Moderate
## 347  Moderately clean    Moderate
## 348  Moderately clean    Moderate
## 349             Clean      Severe
## 350             Clean      Severe
## 351             Clean      Severe
## 352             Clean      Severe
## 353             Clean      Severe
## 354             Clean      Severe
## 355             Clean      Severe
## 356             Clean      Severe
## 357             Clean      Severe
## 358             Clean      Severe
## 359             Clean      Severe
## 360             Clean      Severe
## 361             Clean      Severe
## 362             Clean      Severe
## 363             Clean      Severe
## 364             Clean      Severe
## 365             Clean      Severe
## 366             Clean      Severe
## 367             Clean      Severe
## 368             Clean      Severe
## 369             Clean      Severe
## 370             Clean      Severe
## 371             Clean      Severe
## 372             Clean      Severe
## 373             Clean      Severe
## 374             Clean      Severe
## 375             Clean      Severe
## 376             Clean      Severe
## 377             Clean      Severe
## 378             Clean Very severe
## 379             Clean Very severe
## 380             Clean Very severe
## 381             Clean Very severe
## 382             Clean Very severe
## 383             Clean Very severe
## 384             Clean Very severe
## 385             Clean Very severe
## 386             Clean Very severe
## 387             Clean Very severe
## 388             Clean Very severe
## 389             Clean Very severe
## 390             Clean Very severe
## 391             Clean Very severe
## 392             Clean Very severe
## 393             Clean Very severe
## 394             Clean Very severe
## 395             Clean Very severe
## 396             Clean Very severe
## 397             Clean Very severe
## 398             Clean Very severe
## 399             Clean Very severe
## 400             Clean Very severe
## 401             Clean Very severe
## 402             Clean Very severe
## 403             Clean Very severe
## 404             Clean Very severe
## 405             Clean Very severe
## 406             Clean Very severe
## 407             Clean      Slight
## 408             Clean      Slight
## 409             Clean      Slight
## 410             Clean      Slight
## 411             Clean      Slight
## 412             Clean      Slight
## 413             Clean      Slight
## 414             Clean      Slight
## 415             Clean      Slight
## 416             Clean      Slight
## 417             Clean      Slight
## 418             Clean      Slight
## 419             Clean      Slight
## 420             Clean      Slight
## 421             Clean      Slight
## 422             Clean      Slight
## 423             Clean      Slight
## 424             Clean      Slight
## 425             Clean      Slight
## 426             Clean      Slight
## 427             Clean      Slight
## 428             Clean      Slight
## 429             Clean      Slight
## 430             Clean      Slight
## 431             Clean      Slight
## 432             Clean      Slight
## 433             Clean      Slight
## 434             Clean      Slight
## 435             Clean      Slight
## 436  Moderately clean    Moderate
## 437  Moderately clean    Moderate
## 438  Moderately clean    Moderate
## 439  Moderately clean    Moderate
## 440  Moderately clean    Moderate
## 441  Moderately clean    Moderate
## 442  Moderately clean    Moderate
## 443  Moderately clean    Moderate
## 444  Moderately clean    Moderate
## 445  Moderately clean    Moderate
## 446  Moderately clean    Moderate
## 447  Moderately clean    Moderate
## 448  Moderately clean    Moderate
## 449  Moderately clean    Moderate
## 450  Moderately clean    Moderate
## 451  Moderately clean    Moderate
## 452  Moderately clean    Moderate
## 453  Moderately clean    Moderate
## 454  Moderately clean    Moderate
## 455  Moderately clean    Moderate
## 456  Moderately clean    Moderate
## 457  Moderately clean    Moderate
## 458  Moderately clean    Moderate
## 459  Moderately clean    Moderate
## 460  Moderately clean    Moderate
## 461  Moderately clean    Moderate
## 462  Moderately clean    Moderate
## 463  Moderately clean    Moderate
## 464  Moderately clean    Moderate
## 465  Moderately clean    Moderate
## 466  Moderately clean    Moderate
## 467  Moderately clean    Moderate
## 468  Moderately clean    Moderate
## 469  Moderately clean    Moderate
## 470  Moderately clean    Moderate
## 471  Moderately clean    Moderate
## 472  Moderately clean    Moderate
## 473  Moderately clean    Moderate
## 474  Moderately clean    Moderate
## 475  Moderately clean    Moderate
## 476  Moderately clean    Moderate
## 477  Moderately clean    Moderate
## 478  Moderately clean    Moderate
## 479  Moderately clean    Moderate
## 480  Moderately clean    Moderate
## 481  Moderately clean    Moderate
## 482  Moderately clean    Moderate
## 483  Moderately clean    Moderate
## 484  Moderately clean    Moderate
## 485  Moderately clean    Moderate
## 486  Moderately clean    Moderate
## 487  Moderately clean    Moderate
## 488  Moderately clean    Moderate
## 489  Moderately clean    Moderate
## 490  Moderately clean    Moderate
## 491  Moderately clean    Moderate
## 492  Moderately clean    Moderate
## 493  Moderately clean    Moderate
## 494  Moderately clean    Moderate
## 495  Moderately clean    Moderate
## 496  Moderately clean    Moderate
## 497  Moderately clean    Moderate
## 498  Moderately clean    Moderate
## 499  Moderately clean    Moderate
## 500  Moderately clean    Moderate
## 501  Moderately clean    Moderate
## 502  Moderately clean    Moderate
## 503  Moderately clean    Moderate
## 504  Moderately clean    Moderate
## 505  Moderately clean    Moderate
## 506  Moderately clean    Moderate
## 507  Moderately clean    Moderate
## 508  Moderately clean    Moderate
## 509  Moderately clean    Moderate
## 510  Moderately clean    Moderate
## 511  Moderately clean    Moderate
## 512  Moderately clean    Moderate
## 513  Moderately clean    Moderate
## 514  Moderately clean    Moderate
## 515  Moderately clean    Moderate
## 516  Moderately clean    Moderate
## 517  Moderately clean    Moderate
## 518  Moderately clean    Moderate
## 519  Moderately clean    Moderate
## 520  Moderately clean    Moderate
## 521  Moderately clean    Moderate
## 522  Moderately clean    Moderate
## 523  Moderately clean    Moderate
## 524  Moderately clean    Moderate
## 525  Moderately clean    Moderate
## 526  Moderately clean    Moderate
## 527  Moderately clean    Moderate
## 528  Moderately clean    Moderate
## 529  Moderately clean    Moderate
## 530  Moderately clean    Moderate
## 531  Moderately clean    Moderate
## 532  Moderately clean    Moderate
## 533  Moderately clean    Moderate
## 534  Moderately clean    Moderate
## 535  Moderately clean    Moderate
## 536  Moderately clean    Moderate
## 537  Moderately clean    Moderate
## 538  Moderately clean    Moderate
## 539  Moderately clean    Moderate
## 540  Moderately clean    Moderate
## 541  Moderately clean    Moderate
## 542  Moderately clean    Moderate
## 543  Moderately clean    Moderate
## 544  Moderately clean    Moderate
## 545  Moderately clean    Moderate
## 546  Moderately clean    Moderate
## 547  Moderately clean    Moderate
## 548  Moderately clean    Moderate
## 549  Moderately clean    Moderate
## 550  Moderately clean    Moderate
## 551  Moderately clean    Moderate
## 552  Moderately clean    Moderate
## 553  Moderately clean    Moderate
## 554  Moderately clean    Moderate
## 555  Moderately clean    Moderate
## 556  Moderately clean    Moderate
## 557  Moderately clean    Moderate
## 558  Moderately clean    Moderate
## 559  Moderately clean    Moderate
## 560  Moderately clean    Moderate
## 561  Moderately clean    Moderate
## 562  Moderately clean    Moderate
## 563  Moderately clean    Moderate
## 564  Moderately clean    Moderate
## 565  Moderately clean    Moderate
## 566  Moderately clean    Moderate
## 567  Moderately clean    Moderate
## 568  Moderately clean    Moderate
## 569  Moderately clean    Moderate
## 570  Moderately clean    Moderate
## 571  Moderately clean    Moderate
## 572  Moderately clean    Moderate
## 573  Moderately clean    Moderate
## 574  Moderately clean    Moderate
## 575  Moderately clean    Moderate
## 576  Moderately clean    Moderate
## 577  Moderately clean    Moderate
## 578  Moderately clean    Moderate
## 579  Moderately clean    Moderate
## 580  Moderately clean    Moderate
## 581  Moderately clean    Moderate
## 582  Moderately clean    Moderate
## 583  Moderately clean    Moderate
## 584  Moderately clean    Moderate
## 585            Patchy      Slight
## 586            Patchy      Slight
## 587            Patchy      Slight
## 588            Patchy      Slight
## 589            Patchy      Slight
## 590            Patchy      Slight
## 591            Patchy      Slight
## 592            Patchy      Slight
## 593            Patchy      Slight
## 594            Patchy      Slight
## 595            Patchy      Slight
## 596            Patchy      Slight
## 597            Patchy      Slight
## 598            Patchy      Slight
## 599            Patchy      Slight
## 600            Patchy      Slight
## 601            Patchy      Slight
## 602            Patchy      Slight
## 603            Patchy      Slight
## 604            Patchy      Slight
## 605            Patchy      Slight
## 606            Patchy      Slight
## 607            Patchy      Slight
## 608            Patchy      Slight
## 609            Patchy      Slight
## 610            Patchy      Slight
## 611            Patchy      Slight
## 612            Patchy      Slight
## 613            Patchy      Slight
## 614  Moderately clean    Moderate
## 615  Moderately clean    Moderate
## 616  Moderately clean    Moderate
## 617  Moderately clean    Moderate
## 618  Moderately clean    Moderate
## 619  Moderately clean    Moderate
## 620  Moderately clean    Moderate
## 621  Moderately clean    Moderate
## 622  Moderately clean    Moderate
## 623  Moderately clean    Moderate
## 624  Moderately clean    Moderate
## 625  Moderately clean    Moderate
## 626  Moderately clean    Moderate
## 627  Moderately clean    Moderate
## 628  Moderately clean    Moderate
## 629  Moderately clean    Moderate
## 630  Moderately clean    Moderate
## 631  Moderately clean    Moderate
## 632  Moderately clean    Moderate
## 633  Moderately clean    Moderate
## 634  Moderately clean    Moderate
## 635  Moderately clean    Moderate
## 636  Moderately clean    Moderate
## 637  Moderately clean    Moderate
## 638  Moderately clean    Moderate
## 639  Moderately clean    Moderate
## 640  Moderately clean    Moderate
## 641  Moderately clean    Moderate
## 642  Moderately clean    Moderate
## 643       Very patchy    Moderate
## 644       Very patchy    Moderate
## 645       Very patchy    Moderate
## 646       Very patchy    Moderate
## 647       Very patchy    Moderate
## 648       Very patchy    Moderate
## 649       Very patchy    Moderate
## 650       Very patchy    Moderate
## 651       Very patchy    Moderate
## 652       Very patchy    Moderate
## 653       Very patchy    Moderate
## 654       Very patchy    Moderate
## 655       Very patchy    Moderate
## 656       Very patchy    Moderate
## 657       Very patchy    Moderate
## 658       Very patchy    Moderate
## 659       Very patchy    Moderate
## 660       Very patchy    Moderate
## 661       Very patchy    Moderate
## 662       Very patchy    Moderate
## 663       Very patchy    Moderate
## 664       Very patchy    Moderate
## 665       Very patchy    Moderate
## 666       Very patchy    Moderate
## 667       Very patchy    Moderate
## 668       Very patchy    Moderate
## 669       Very patchy    Moderate
## 670       Very patchy    Moderate
## 671       Very patchy    Moderate
## 672       Very patchy        None
## 673       Very patchy        None
## 674       Very patchy        None
## 675       Very patchy        None
## 676       Very patchy        None
## 677       Very patchy        None
## 678       Very patchy        None
## 679       Very patchy        None
## 680       Very patchy        None
## 681       Very patchy        None
## 682       Very patchy        None
## 683       Very patchy        None
## 684       Very patchy        None
## 685       Very patchy        None
## 686       Very patchy        None
## 687       Very patchy        None
## 688       Very patchy        None
## 689       Very patchy        None
## 690       Very patchy        None
## 691       Very patchy        None
## 692       Very patchy        None
## 693       Very patchy        None
## 694       Very patchy        None
## 695       Very patchy        None
## 696       Very patchy        None
## 697       Very patchy        None
## 698       Very patchy        None
## 699       Very patchy        None
## 700       Very patchy        None
## 701  Moderately clean    Moderate
## 702  Moderately clean    Moderate
## 703  Moderately clean    Moderate
## 704  Moderately clean    Moderate
## 705  Moderately clean    Moderate
## 706  Moderately clean    Moderate
## 707  Moderately clean    Moderate
## 708  Moderately clean    Moderate
## 709  Moderately clean    Moderate
## 710  Moderately clean    Moderate
## 711  Moderately clean    Moderate
## 712  Moderately clean    Moderate
## 713  Moderately clean    Moderate
## 714  Moderately clean    Moderate
## 715  Moderately clean    Moderate
## 716  Moderately clean    Moderate
## 717  Moderately clean    Moderate
## 718  Moderately clean    Moderate
## 719  Moderately clean    Moderate
## 720  Moderately clean    Moderate
## 721  Moderately clean    Moderate
## 722  Moderately clean    Moderate
## 723  Moderately clean    Moderate
## 724  Moderately clean    Moderate
## 725  Moderately clean    Moderate
## 726  Moderately clean    Moderate
## 727  Moderately clean    Moderate
## 728  Moderately clean    Moderate
## 729  Moderately clean    Moderate
## 730  Moderately clean    Moderate
## 731  Moderately clean    Moderate
## 732  Moderately clean    Moderate
## 733  Moderately clean    Moderate
## 734  Moderately clean    Moderate
## 735  Moderately clean    Moderate
## 736  Moderately clean    Moderate
## 737  Moderately clean    Moderate
## 738  Moderately clean    Moderate
## 739  Moderately clean    Moderate
## 740  Moderately clean    Moderate
## 741  Moderately clean    Moderate
## 742  Moderately clean    Moderate
## 743  Moderately clean    Moderate
## 744  Moderately clean    Moderate
## 745  Moderately clean    Moderate
## 746  Moderately clean    Moderate
## 747  Moderately clean    Moderate
## 748  Moderately clean    Moderate
## 749  Moderately clean    Moderate
## 750  Moderately clean    Moderate
## 751  Moderately clean    Moderate
## 752  Moderately clean    Moderate
## 753  Moderately clean    Moderate
## 754  Moderately clean    Moderate
## 755  Moderately clean    Moderate
## 756  Moderately clean    Moderate
## 757  Moderately clean    Moderate
## 758  Moderately clean    Moderate
## 759  Moderately clean    Moderate
## 760  Moderately clean    Moderate
## 761  Moderately clean    Moderate
## 762  Moderately clean    Moderate
## 763  Moderately clean    Moderate
## 764  Moderately clean    Moderate
## 765  Moderately clean    Moderate
## 766  Moderately clean    Moderate
## 767  Moderately clean    Moderate
## 768  Moderately clean    Moderate
## 769  Moderately clean    Moderate
## 770  Moderately clean    Moderate
## 771  Moderately clean    Moderate
## 772  Moderately clean    Moderate
## 773  Moderately clean    Moderate
## 774  Moderately clean    Moderate
## 775  Moderately clean    Moderate
## 776  Moderately clean    Moderate
## 777  Moderately clean    Moderate
## 778  Moderately clean    Moderate
## 779  Moderately clean    Moderate
## 780  Moderately clean    Moderate
## 781  Moderately clean    Moderate
## 782  Moderately clean    Moderate
## 783  Moderately clean    Moderate
## 784  Moderately clean    Moderate
## 785  Moderately clean    Moderate
## 786  Moderately clean    Moderate
## 787  Moderately clean    Moderate
## 788       Very patchy        None
## 789       Very patchy        None
## 790       Very patchy        None
## 791       Very patchy        None
## 792       Very patchy        None
## 793       Very patchy        None
## 794       Very patchy        None
## 795       Very patchy        None
## 796       Very patchy        None
## 797       Very patchy        None
## 798       Very patchy        None
## 799       Very patchy        None
## 800       Very patchy        None
## 801       Very patchy        None
## 802       Very patchy        None
## 803       Very patchy        None
## 804       Very patchy        None
## 805       Very patchy        None
## 806       Very patchy        None
## 807       Very patchy        None
## 808       Very patchy        None
## 809       Very patchy        None
## 810       Very patchy        None
## 811       Very patchy        None
## 812       Very patchy        None
## 813       Very patchy        None
## 814       Very patchy        None
## 815       Very patchy        None
## 816       Very patchy        None
## 817  Moderately clean    Moderate
## 818  Moderately clean    Moderate
## 819  Moderately clean    Moderate
## 820  Moderately clean    Moderate
## 821  Moderately clean    Moderate
## 822  Moderately clean    Moderate
## 823  Moderately clean    Moderate
## 824  Moderately clean    Moderate
## 825  Moderately clean    Moderate
## 826  Moderately clean    Moderate
## 827  Moderately clean    Moderate
## 828  Moderately clean    Moderate
## 829  Moderately clean    Moderate
## 830  Moderately clean    Moderate
## 831  Moderately clean    Moderate
## 832  Moderately clean    Moderate
## 833  Moderately clean    Moderate
## 834  Moderately clean    Moderate
## 835  Moderately clean    Moderate
## 836  Moderately clean    Moderate
## 837  Moderately clean    Moderate
## 838  Moderately clean    Moderate
## 839  Moderately clean    Moderate
## 840  Moderately clean    Moderate
## 841  Moderately clean    Moderate
## 842  Moderately clean    Moderate
## 843  Moderately clean    Moderate
## 844  Moderately clean    Moderate
## 845  Moderately clean    Moderate
## 846  Moderately clean    Moderate
## 847  Moderately clean    Moderate
## 848  Moderately clean    Moderate
## 849  Moderately clean    Moderate
## 850  Moderately clean    Moderate
## 851  Moderately clean    Moderate
## 852  Moderately clean    Moderate
## 853  Moderately clean    Moderate
## 854  Moderately clean    Moderate
## 855  Moderately clean    Moderate
## 856  Moderately clean    Moderate
## 857  Moderately clean    Moderate
## 858  Moderately clean    Moderate
## 859  Moderately clean    Moderate
## 860  Moderately clean    Moderate
## 861  Moderately clean    Moderate
## 862  Moderately clean    Moderate
## 863  Moderately clean    Moderate
## 864  Moderately clean    Moderate
## 865  Moderately clean    Moderate
## 866  Moderately clean    Moderate
## 867  Moderately clean    Moderate
## 868  Moderately clean    Moderate
## 869  Moderately clean    Moderate
## 870  Moderately clean    Moderate
## 871  Moderately clean    Moderate
## 872  Moderately clean    Moderate
## 873  Moderately clean    Moderate
## 874  Moderately clean    Moderate
## 875  Moderately clean    Moderate
## 876  Moderately clean    Moderate
## 877  Moderately clean    Moderate
## 878  Moderately clean    Moderate
## 879  Moderately clean    Moderate
## 880  Moderately clean    Moderate
## 881  Moderately clean    Moderate
## 882  Moderately clean    Moderate
## 883  Moderately clean    Moderate
## 884  Moderately clean    Moderate
## 885  Moderately clean    Moderate
## 886  Moderately clean    Moderate
## 887  Moderately clean    Moderate
## 888  Moderately clean    Moderate
## 889  Moderately clean    Moderate
## 890  Moderately clean    Moderate
## 891  Moderately clean    Moderate
## 892  Moderately clean    Moderate
## 893  Moderately clean    Moderate
## 894  Moderately clean    Moderate
## 895  Moderately clean    Moderate
## 896  Moderately clean    Moderate
## 897  Moderately clean    Moderate
## 898  Moderately clean    Moderate
## 899  Moderately clean    Moderate
## 900  Moderately clean    Moderate
## 901  Moderately clean    Moderate
## 902  Moderately clean    Moderate
## 903  Moderately clean    Moderate
## 904             Clean      Severe
## 905             Clean      Severe
## 906             Clean      Severe
## 907             Clean      Severe
## 908             Clean      Severe
## 909             Clean      Severe
## 910             Clean      Severe
## 911             Clean      Severe
## 912             Clean      Severe
## 913             Clean      Severe
## 914             Clean      Severe
## 915             Clean      Severe
## 916             Clean      Severe
## 917             Clean      Severe
## 918             Clean      Severe
## 919             Clean      Severe
## 920             Clean      Severe
## 921             Clean      Severe
## 922             Clean      Severe
## 923             Clean      Severe
## 924             Clean      Severe
## 925             Clean      Severe
## 926             Clean      Severe
## 927             Clean      Severe
## 928             Clean      Severe
## 929             Clean      Severe
## 930             Clean      Severe
## 931             Clean      Severe
## 932             Clean      Severe
## 933            Patchy    Moderate
## 934            Patchy    Moderate
## 935            Patchy    Moderate
## 936            Patchy    Moderate
## 937            Patchy    Moderate
## 938            Patchy    Moderate
## 939            Patchy    Moderate
## 940            Patchy    Moderate
## 941            Patchy    Moderate
## 942            Patchy    Moderate
## 943            Patchy    Moderate
## 944            Patchy    Moderate
## 945            Patchy    Moderate
## 946            Patchy    Moderate
## 947            Patchy    Moderate
## 948            Patchy    Moderate
## 949            Patchy    Moderate
## 950            Patchy    Moderate
## 951            Patchy    Moderate
## 952            Patchy    Moderate
## 953            Patchy    Moderate
## 954            Patchy    Moderate
## 955            Patchy    Moderate
## 956            Patchy    Moderate
## 957            Patchy    Moderate
## 958            Patchy    Moderate
## 959            Patchy    Moderate
## 960            Patchy    Moderate
## 961            Patchy    Moderate
## 962             Clean      Severe
## 963             Clean      Severe
## 964             Clean      Severe
## 965             Clean      Severe
## 966             Clean      Severe
## 967             Clean      Severe
## 968             Clean      Severe
## 969             Clean      Severe
## 970             Clean      Severe
## 971             Clean      Severe
## 972             Clean      Severe
## 973             Clean      Severe
## 974             Clean      Severe
## 975             Clean      Severe
## 976             Clean      Severe
## 977             Clean      Severe
## 978             Clean      Severe
## 979             Clean      Severe
## 980             Clean      Severe
## 981             Clean      Severe
## 982             Clean      Severe
## 983             Clean      Severe
## 984             Clean      Severe
## 985             Clean      Severe
## 986             Clean      Severe
## 987             Clean      Severe
## 988             Clean      Severe
## 989             Clean      Severe
## 990             Clean      Severe
## 991             Clean      Severe
## 992             Clean      Severe
## 993            Patchy      Slight
## 994            Patchy      Slight
## 995            Patchy      Slight
## 996            Patchy      Slight
## 997            Patchy      Slight
## 998            Patchy      Slight
## 999            Patchy      Slight
## 1000           Patchy      Slight
## 1001           Patchy      Slight
## 1002           Patchy      Slight
## 1003           Patchy      Slight
## 1004           Patchy      Slight
## 1005           Patchy      Slight
## 1006           Patchy      Slight
## 1007           Patchy      Slight
## 1008           Patchy      Slight
## 1009           Patchy      Slight
## 1010           Patchy      Slight
## 1011           Patchy      Slight
## 1012           Patchy      Slight
## 1013           Patchy      Slight
## 1014           Patchy      Slight
## 1015           Patchy      Slight
## 1016           Patchy      Slight
## 1017           Patchy      Slight
## 1018           Patchy      Slight
## 1019           Patchy      Slight
## 1020           Patchy      Slight
## 1021           Patchy      Slight
## 1022           Patchy      Slight
## 1023           Patchy      Slight
## 1024            Clean      Slight
## 1025            Clean      Slight
## 1026            Clean      Slight
## 1027            Clean      Slight
## 1028            Clean      Slight
## 1029            Clean      Slight
## 1030            Clean      Slight
## 1031            Clean      Slight
## 1032            Clean      Slight
## 1033            Clean      Slight
## 1034            Clean      Slight
## 1035            Clean      Slight
## 1036            Clean      Slight
## 1037            Clean      Slight
## 1038            Clean      Slight
## 1039            Clean      Slight
## 1040            Clean      Slight
## 1041            Clean      Slight
## 1042            Clean      Slight
## 1043            Clean      Slight
## 1044            Clean      Slight
## 1045            Clean      Slight
## 1046            Clean      Slight
## 1047            Clean      Slight
## 1048            Clean      Slight
## 1049            Clean      Slight
## 1050            Clean      Slight
## 1051            Clean      Slight
## 1052            Clean      Slight
## 1053            Clean      Slight
## 1054            Clean      Slight
## 1055       Clean burn    Moderate
## 1056       Clean burn    Moderate
## 1057       Clean burn    Moderate
## 1058       Clean burn    Moderate
## 1059       Clean burn    Moderate
## 1060       Clean burn    Moderate
## 1061       Clean burn    Moderate
## 1062       Clean burn    Moderate
## 1063       Clean burn    Moderate
## 1064       Clean burn    Moderate
## 1065       Clean burn    Moderate
## 1066       Clean burn    Moderate
## 1067       Clean burn    Moderate
## 1068       Clean burn    Moderate
## 1069       Clean burn    Moderate
## 1070       Clean burn    Moderate
## 1071       Clean burn    Moderate
## 1072       Clean burn    Moderate
## 1073       Clean burn    Moderate
## 1074       Clean burn    Moderate
## 1075       Clean burn    Moderate
## 1076       Clean burn    Moderate
## 1077       Clean burn    Moderate
## 1078       Clean burn    Moderate
## 1079       Clean burn    Moderate
## 1080       Clean burn    Moderate
## 1081       Clean burn    Moderate
## 1082       Clean burn    Moderate
## 1083       Clean burn    Moderate
## 1084       Clean burn    Moderate
## 1085       Clean burn    Moderate
## 1086            Clean      Slight
## 1087            Clean      Slight
## 1088            Clean      Slight
## 1089            Clean      Slight
## 1090            Clean      Slight
## 1091            Clean      Slight
## 1092            Clean      Slight
## 1093            Clean      Slight
## 1094            Clean      Slight
## 1095            Clean      Slight
## 1096            Clean      Slight
## 1097            Clean      Slight
## 1098            Clean      Slight
## 1099            Clean      Slight
## 1100            Clean      Slight
## 1101            Clean      Slight
## 1102            Clean      Slight
## 1103            Clean      Slight
## 1104            Clean      Slight
## 1105            Clean      Slight
## 1106            Clean      Slight
## 1107            Clean      Slight
## 1108            Clean      Slight
## 1109            Clean      Slight
## 1110            Clean      Slight
## 1111            Clean      Slight
## 1112            Clean      Slight
## 1113            Clean      Slight
## 1114            Clean      Slight
## 1115            Clean      Slight
## 1116            Clean      Slight
## 1117            Clean      Slight
## 1118            Clean      Slight
## 1119            Clean      Slight
## 1120            Clean      Slight
## 1121            Clean      Slight
## 1122            Clean      Slight
## 1123            Clean      Slight
## 1124            Clean      Slight
## 1125            Clean      Slight
## 1126            Clean      Slight
## 1127            Clean      Slight
## 1128            Clean      Slight
## 1129            Clean      Slight
## 1130            Clean      Slight
## 1131            Clean      Slight
## 1132            Clean      Slight
## 1133            Clean      Slight
## 1134            Clean      Slight
## 1135            Clean      Slight
## 1136            Clean      Slight
## 1137            Clean      Slight
## 1138            Clean      Slight
## 1139            Clean      Slight
## 1140            Clean      Slight
## 1141            Clean      Slight
## 1142            Clean      Slight
## 1143            Clean      Slight
## 1144            Clean      Slight
## 1145            Clean      Slight
## 1146            Clean      Slight
## 1147            Clean      Slight
## 1148 Moderately clean    Moderate
## 1149 Moderately clean    Moderate
## 1150 Moderately clean    Moderate
## 1151 Moderately clean    Moderate
## 1152 Moderately clean    Moderate
## 1153 Moderately clean    Moderate
## 1154 Moderately clean    Moderate
## 1155 Moderately clean    Moderate
## 1156 Moderately clean    Moderate
## 1157 Moderately clean    Moderate
## 1158 Moderately clean    Moderate
## 1159 Moderately clean    Moderate
## 1160 Moderately clean    Moderate
## 1161 Moderately clean    Moderate
## 1162 Moderately clean    Moderate
## 1163 Moderately clean    Moderate
## 1164 Moderately clean    Moderate
## 1165 Moderately clean    Moderate
## 1166 Moderately clean    Moderate
## 1167 Moderately clean    Moderate
## 1168 Moderately clean    Moderate
## 1169 Moderately clean    Moderate
## 1170 Moderately clean    Moderate
## 1171 Moderately clean    Moderate
## 1172 Moderately clean    Moderate
## 1173 Moderately clean    Moderate
## 1174 Moderately clean    Moderate
## 1175 Moderately clean    Moderate
## 1176 Moderately clean    Moderate
## 1177 Moderately clean        None
## 1178 Moderately clean        None
## 1179 Moderately clean        None
## 1180 Moderately clean        None
## 1181 Moderately clean        None
## 1182 Moderately clean        None
## 1183 Moderately clean        None
## 1184 Moderately clean        None
## 1185 Moderately clean        None
## 1186 Moderately clean        None
## 1187 Moderately clean        None
## 1188 Moderately clean        None
## 1189 Moderately clean        None
## 1190 Moderately clean        None
## 1191 Moderately clean        None
## 1192 Moderately clean        None
## 1193 Moderately clean        None
## 1194 Moderately clean        None
## 1195 Moderately clean        None
## 1196 Moderately clean        None
## 1197 Moderately clean        None
## 1198 Moderately clean        None
## 1199 Moderately clean        None
## 1200 Moderately clean        None
## 1201 Moderately clean        None
## 1202 Moderately clean        None
## 1203 Moderately clean        None
## 1204 Moderately clean        None
## 1205 Moderately clean        None
## 1206 Moderately clean    Moderate
## 1207 Moderately clean    Moderate
## 1208 Moderately clean    Moderate
## 1209 Moderately clean    Moderate
## 1210 Moderately clean    Moderate
## 1211 Moderately clean    Moderate
## 1212 Moderately clean    Moderate
## 1213 Moderately clean    Moderate
## 1214 Moderately clean    Moderate
## 1215 Moderately clean    Moderate
## 1216 Moderately clean    Moderate
## 1217 Moderately clean    Moderate
## 1218 Moderately clean    Moderate
## 1219 Moderately clean    Moderate
## 1220 Moderately clean    Moderate
## 1221 Moderately clean    Moderate
## 1222 Moderately clean    Moderate
## 1223 Moderately clean    Moderate
## 1224 Moderately clean    Moderate
## 1225 Moderately clean    Moderate
## 1226 Moderately clean    Moderate
## 1227 Moderately clean    Moderate
## 1228 Moderately clean    Moderate
## 1229 Moderately clean    Moderate
## 1230 Moderately clean    Moderate
## 1231 Moderately clean    Moderate
## 1232 Moderately clean    Moderate
## 1233 Moderately clean    Moderate
## 1234 Moderately clean    Moderate
## 1235 Moderately clean    Moderate
## 1236 Moderately clean    Moderate
## 1237 Moderately clean    Moderate
## 1238 Moderately clean    Moderate
## 1239 Moderately clean    Moderate
## 1240 Moderately clean    Moderate
## 1241 Moderately clean    Moderate
## 1242 Moderately clean    Moderate
## 1243 Moderately clean    Moderate
## 1244 Moderately clean    Moderate
## 1245 Moderately clean    Moderate
## 1246 Moderately clean    Moderate
## 1247 Moderately clean    Moderate
## 1248 Moderately clean    Moderate
## 1249 Moderately clean    Moderate
## 1250 Moderately clean    Moderate
## 1251 Moderately clean    Moderate
## 1252 Moderately clean    Moderate
## 1253 Moderately clean    Moderate
## 1254 Moderately clean    Moderate
## 1255 Moderately clean    Moderate
## 1256 Moderately clean    Moderate
## 1257 Moderately clean    Moderate
## 1258 Moderately clean    Moderate
## 1259 Moderately clean    Moderate
## 1260 Moderately clean    Moderate
## 1261 Moderately clean    Moderate
## 1262 Moderately clean    Moderate
## 1263 Moderately clean    Moderate
## 1264            Clean Very severe
## 1265            Clean Very severe
## 1266            Clean Very severe
## 1267            Clean Very severe
## 1268            Clean Very severe
## 1269            Clean Very severe
## 1270            Clean Very severe
## 1271            Clean Very severe
## 1272            Clean Very severe
## 1273            Clean Very severe
## 1274            Clean Very severe
## 1275            Clean Very severe
## 1276            Clean Very severe
## 1277            Clean Very severe
## 1278            Clean Very severe
## 1279            Clean Very severe
## 1280            Clean Very severe
## 1281            Clean Very severe
## 1282            Clean Very severe
## 1283            Clean Very severe
## 1284            Clean Very severe
## 1285            Clean Very severe
## 1286            Clean Very severe
## 1287            Clean Very severe
## 1288            Clean Very severe
## 1289            Clean Very severe
## 1290            Clean Very severe
## 1291            Clean Very severe
## 1292            Clean Very severe
## 1293            Clean Very severe
## 1294            Clean Very severe
## 1295            Clean Very severe
## 1296            Clean Very severe
## 1297            Clean Very severe
## 1298            Clean Very severe
## 1299            Clean Very severe
## 1300            Clean Very severe
## 1301            Clean Very severe
## 1302            Clean Very severe
## 1303            Clean Very severe
## 1304            Clean Very severe
## 1305            Clean Very severe
## 1306            Clean Very severe
## 1307            Clean Very severe
## 1308            Clean Very severe
## 1309            Clean Very severe
## 1310            Clean Very severe
## 1311            Clean Very severe
## 1312            Clean Very severe
## 1313            Clean Very severe
## 1314            Clean Very severe
## 1315            Clean Very severe
## 1316            Clean Very severe
## 1317            Clean Very severe
## 1318            Clean Very severe
## 1319            Clean Very severe
## 1320            Clean Very severe
## 1321            Clean Very severe
## 1322            Clean Very severe
## 1323            Clean Very severe
## 1324            Clean Very severe
## 1325            Clean Very severe
## 1326 Moderately clean    Moderate
## 1327 Moderately clean    Moderate
## 1328 Moderately clean    Moderate
## 1329 Moderately clean    Moderate
## 1330 Moderately clean    Moderate
## 1331 Moderately clean    Moderate
## 1332 Moderately clean    Moderate
## 1333 Moderately clean    Moderate
## 1334 Moderately clean    Moderate
## 1335 Moderately clean    Moderate
## 1336 Moderately clean    Moderate
## 1337 Moderately clean    Moderate
## 1338 Moderately clean    Moderate
## 1339 Moderately clean    Moderate
## 1340 Moderately clean    Moderate
## 1341 Moderately clean    Moderate
## 1342 Moderately clean    Moderate
## 1343 Moderately clean    Moderate
## 1344 Moderately clean    Moderate
## 1345 Moderately clean    Moderate
## 1346 Moderately clean    Moderate
## 1347 Moderately clean    Moderate
## 1348 Moderately clean    Moderate
## 1349 Moderately clean    Moderate
## 1350 Moderately clean    Moderate
## 1351 Moderately clean    Moderate
## 1352 Moderately clean    Moderate
## 1353 Moderately clean    Moderate
## 1354 Moderately clean    Moderate
## 1355 Moderately clean      Slight
## 1356 Moderately clean      Slight
## 1357 Moderately clean      Slight
## 1358 Moderately clean      Slight
## 1359 Moderately clean      Slight
## 1360 Moderately clean      Slight
## 1361 Moderately clean      Slight
## 1362 Moderately clean      Slight
## 1363 Moderately clean      Slight
## 1364 Moderately clean      Slight
## 1365 Moderately clean      Slight
## 1366 Moderately clean      Slight
## 1367 Moderately clean      Slight
## 1368 Moderately clean      Slight
## 1369 Moderately clean      Slight
## 1370 Moderately clean      Slight
## 1371 Moderately clean      Slight
## 1372 Moderately clean      Slight
## 1373 Moderately clean      Slight
## 1374 Moderately clean      Slight
## 1375 Moderately clean      Slight
## 1376 Moderately clean      Slight
## 1377 Moderately clean      Slight
## 1378 Moderately clean      Slight
## 1379 Moderately clean      Slight
## 1380 Moderately clean      Slight
## 1381 Moderately clean      Slight
## 1382 Moderately clean      Slight
## 1383 Moderately clean      Slight
## 1384 Moderately clean    Moderate
## 1385 Moderately clean    Moderate
## 1386 Moderately clean    Moderate
## 1387 Moderately clean    Moderate
## 1388 Moderately clean    Moderate
## 1389 Moderately clean    Moderate
## 1390 Moderately clean    Moderate
## 1391 Moderately clean    Moderate
## 1392 Moderately clean    Moderate
## 1393 Moderately clean    Moderate
## 1394 Moderately clean    Moderate
## 1395 Moderately clean    Moderate
## 1396 Moderately clean    Moderate
## 1397 Moderately clean    Moderate
## 1398 Moderately clean    Moderate
## 1399 Moderately clean    Moderate
## 1400 Moderately clean    Moderate
## 1401 Moderately clean    Moderate
## 1402 Moderately clean    Moderate
## 1403 Moderately clean    Moderate
## 1404 Moderately clean    Moderate
## 1405 Moderately clean    Moderate
## 1406 Moderately clean    Moderate
## 1407 Moderately clean    Moderate
## 1408 Moderately clean    Moderate
## 1409 Moderately clean    Moderate
## 1410 Moderately clean    Moderate
## 1411 Moderately clean    Moderate
## 1412 Moderately clean    Moderate
## 1413      Very patchy    Moderate
## 1414      Very patchy    Moderate
## 1415      Very patchy    Moderate
## 1416      Very patchy    Moderate
## 1417      Very patchy    Moderate
## 1418      Very patchy    Moderate
## 1419      Very patchy    Moderate
## 1420      Very patchy    Moderate
## 1421      Very patchy    Moderate
## 1422      Very patchy    Moderate
## 1423      Very patchy    Moderate
## 1424      Very patchy    Moderate
## 1425      Very patchy    Moderate
## 1426      Very patchy    Moderate
## 1427      Very patchy    Moderate
## 1428      Very patchy    Moderate
## 1429      Very patchy    Moderate
## 1430      Very patchy    Moderate
## 1431      Very patchy    Moderate
## 1432      Very patchy    Moderate
## 1433      Very patchy    Moderate
## 1434      Very patchy    Moderate
## 1435      Very patchy    Moderate
## 1436      Very patchy    Moderate
## 1437      Very patchy    Moderate
## 1438      Very patchy    Moderate
## 1439      Very patchy    Moderate
## 1440      Very patchy    Moderate
## 1441      Very patchy    Moderate
## 1442            Clean    Moderate
## 1443            Clean    Moderate
## 1444            Clean    Moderate
## 1445            Clean    Moderate
## 1446            Clean    Moderate
## 1447            Clean    Moderate
## 1448            Clean    Moderate
## 1449            Clean    Moderate
## 1450            Clean    Moderate
## 1451            Clean    Moderate
## 1452            Clean    Moderate
## 1453            Clean    Moderate
## 1454            Clean    Moderate
## 1455            Clean    Moderate
## 1456            Clean    Moderate
## 1457            Clean    Moderate
## 1458            Clean    Moderate
## 1459            Clean    Moderate
## 1460            Clean    Moderate
## 1461            Clean    Moderate
## 1462            Clean    Moderate
## 1463            Clean    Moderate
## 1464            Clean    Moderate
## 1465            Clean    Moderate
## 1466            Clean    Moderate
## 1467            Clean    Moderate
## 1468            Clean    Moderate
## 1469            Clean    Moderate
## 1470            Clean    Moderate
## 1471           Patchy      Slight
## 1472           Patchy      Slight
## 1473           Patchy      Slight
## 1474           Patchy      Slight
## 1475           Patchy      Slight
## 1476           Patchy      Slight
## 1477           Patchy      Slight
## 1478           Patchy      Slight
## 1479           Patchy      Slight
## 1480           Patchy      Slight
## 1481           Patchy      Slight
## 1482           Patchy      Slight
## 1483           Patchy      Slight
## 1484           Patchy      Slight
## 1485           Patchy      Slight
## 1486           Patchy      Slight
## 1487           Patchy      Slight
## 1488           Patchy      Slight
## 1489           Patchy      Slight
## 1490           Patchy      Slight
## 1491           Patchy      Slight
## 1492           Patchy      Slight
## 1493           Patchy      Slight
## 1494           Patchy      Slight
## 1495           Patchy      Slight
## 1496           Patchy      Slight
## 1497           Patchy      Slight
## 1498           Patchy      Slight
## 1499           Patchy      Slight
## 1500 Moderately clean      Slight
## 1501 Moderately clean      Slight
## 1502 Moderately clean      Slight
## 1503 Moderately clean      Slight
## 1504 Moderately clean      Slight
## 1505 Moderately clean      Slight
## 1506 Moderately clean      Slight
## 1507 Moderately clean      Slight
## 1508 Moderately clean      Slight
## 1509 Moderately clean      Slight
## 1510 Moderately clean      Slight
## 1511 Moderately clean      Slight
## 1512 Moderately clean      Slight
## 1513 Moderately clean      Slight
## 1514 Moderately clean      Slight
## 1515 Moderately clean      Slight
## 1516 Moderately clean      Slight
## 1517 Moderately clean      Slight
## 1518 Moderately clean      Slight
## 1519 Moderately clean      Slight
## 1520 Moderately clean      Slight
## 1521 Moderately clean      Slight
## 1522 Moderately clean      Slight
## 1523 Moderately clean      Slight
## 1524 Moderately clean      Slight
## 1525 Moderately clean      Slight
## 1526 Moderately clean      Slight
## 1527 Moderately clean      Slight
## 1528 Moderately clean      Slight
## 1529      Very patchy        None
## 1530      Very patchy        None
## 1531      Very patchy        None
## 1532      Very patchy        None
## 1533      Very patchy        None
## 1534      Very patchy        None
## 1535      Very patchy        None
## 1536      Very patchy        None
## 1537      Very patchy        None
## 1538      Very patchy        None
## 1539      Very patchy        None
## 1540      Very patchy        None
## 1541      Very patchy        None
## 1542      Very patchy        None
## 1543      Very patchy        None
## 1544      Very patchy        None
## 1545      Very patchy        None
## 1546      Very patchy        None
## 1547      Very patchy        None
## 1548      Very patchy        None
## 1549      Very patchy        None
## 1550      Very patchy        None
## 1551      Very patchy        None
## 1552      Very patchy        None
## 1553      Very patchy        None
## 1554      Very patchy        None
## 1555      Very patchy        None
## 1556      Very patchy        None
## 1557      Very patchy        None
## 1558 Moderately clean    Moderate
## 1559 Moderately clean    Moderate
## 1560 Moderately clean    Moderate
## 1561 Moderately clean    Moderate
## 1562 Moderately clean    Moderate
## 1563 Moderately clean    Moderate
## 1564 Moderately clean    Moderate
## 1565 Moderately clean    Moderate
## 1566 Moderately clean    Moderate
## 1567 Moderately clean    Moderate
## 1568 Moderately clean    Moderate
## 1569 Moderately clean    Moderate
## 1570 Moderately clean    Moderate
## 1571 Moderately clean    Moderate
## 1572 Moderately clean    Moderate
## 1573 Moderately clean    Moderate
## 1574 Moderately clean    Moderate
## 1575 Moderately clean    Moderate
## 1576 Moderately clean    Moderate
## 1577 Moderately clean    Moderate
## 1578 Moderately clean    Moderate
## 1579 Moderately clean    Moderate
## 1580 Moderately clean    Moderate
## 1581 Moderately clean    Moderate
## 1582 Moderately clean    Moderate
## 1583 Moderately clean    Moderate
## 1584 Moderately clean    Moderate
## 1585 Moderately clean    Moderate
## 1586 Moderately clean    Moderate
## 1587 Moderately clean      Severe
## 1588 Moderately clean      Severe
## 1589 Moderately clean      Severe
## 1590 Moderately clean      Severe
## 1591 Moderately clean      Severe
## 1592 Moderately clean      Severe
## 1593 Moderately clean      Severe
## 1594 Moderately clean      Severe
## 1595 Moderately clean      Severe
## 1596 Moderately clean      Severe
## 1597 Moderately clean      Severe
## 1598 Moderately clean      Severe
## 1599 Moderately clean      Severe
## 1600 Moderately clean      Severe
## 1601 Moderately clean      Severe
## 1602 Moderately clean      Severe
## 1603 Moderately clean      Severe
## 1604 Moderately clean      Severe
## 1605 Moderately clean      Severe
## 1606 Moderately clean      Severe
## 1607 Moderately clean      Severe
## 1608 Moderately clean      Severe
## 1609 Moderately clean      Severe
## 1610 Moderately clean      Severe
## 1611 Moderately clean      Severe
## 1612 Moderately clean      Severe
## 1613 Moderately clean      Severe
## 1614 Moderately clean      Severe
## 1615 Moderately clean      Severe
## 1616 Moderately clean      Severe
## 1617 Moderately clean      Severe
## 1618 Moderately clean      Severe
## 1619 Moderately clean      Severe
## 1620 Moderately clean      Severe
## 1621 Moderately clean      Severe
## 1622 Moderately clean      Severe
## 1623 Moderately clean      Severe
## 1624 Moderately clean      Severe
## 1625 Moderately clean      Severe
## 1626 Moderately clean      Severe
## 1627 Moderately clean      Severe
## 1628 Moderately clean      Severe
## 1629 Moderately clean      Severe
## 1630 Moderately clean      Severe
## 1631 Moderately clean      Severe
## 1632 Moderately clean      Severe
## 1633 Moderately clean      Severe
## 1634 Moderately clean      Severe
## 1635 Moderately clean      Severe
## 1636 Moderately clean      Severe
## 1637 Moderately clean      Severe
## 1638 Moderately clean      Severe
## 1639 Moderately clean      Severe
## 1640 Moderately clean      Severe
## 1641 Moderately clean      Severe
## 1642 Moderately clean      Severe
## 1643 Moderately clean      Severe
## 1644 Moderately clean      Severe
## 1645            Clean      Severe
## 1646            Clean      Severe
## 1647            Clean      Severe
## 1648            Clean      Severe
## 1649            Clean      Severe
## 1650            Clean      Severe
## 1651            Clean      Severe
## 1652            Clean      Severe
## 1653            Clean      Severe
## 1654            Clean      Severe
## 1655            Clean      Severe
## 1656            Clean      Severe
## 1657            Clean      Severe
## 1658            Clean      Severe
## 1659            Clean      Severe
## 1660            Clean      Severe
## 1661            Clean      Severe
## 1662            Clean      Severe
## 1663            Clean      Severe
## 1664            Clean      Severe
## 1665            Clean      Severe
## 1666            Clean      Severe
## 1667            Clean      Severe
## 1668            Clean      Severe
## 1669            Clean      Severe
## 1670            Clean      Severe
## 1671            Clean      Severe
## 1672            Clean      Severe
## 1673            Clean      Severe
## 1674           Patchy    Moderate
## 1675           Patchy    Moderate
## 1676           Patchy    Moderate
## 1677           Patchy    Moderate
## 1678           Patchy    Moderate
## 1679           Patchy    Moderate
## 1680           Patchy    Moderate
## 1681           Patchy    Moderate
## 1682           Patchy    Moderate
## 1683           Patchy    Moderate
## 1684           Patchy    Moderate
## 1685           Patchy    Moderate
## 1686           Patchy    Moderate
## 1687           Patchy    Moderate
## 1688           Patchy    Moderate
## 1689           Patchy    Moderate
## 1690           Patchy    Moderate
## 1691           Patchy    Moderate
## 1692           Patchy    Moderate
## 1693           Patchy    Moderate
## 1694           Patchy    Moderate
## 1695           Patchy    Moderate
## 1696           Patchy    Moderate
## 1697           Patchy    Moderate
## 1698           Patchy    Moderate
## 1699           Patchy    Moderate
## 1700           Patchy    Moderate
## 1701           Patchy    Moderate
## 1702           Patchy    Moderate
## 1703            Clean Very severe
## 1704            Clean Very severe
## 1705            Clean Very severe
## 1706            Clean Very severe
## 1707            Clean Very severe
## 1708            Clean Very severe
## 1709            Clean Very severe
## 1710            Clean Very severe
## 1711            Clean Very severe
## 1712            Clean Very severe
## 1713            Clean Very severe
## 1714            Clean Very severe
## 1715            Clean Very severe
## 1716            Clean Very severe
## 1717            Clean Very severe
## 1718            Clean Very severe
## 1719            Clean Very severe
## 1720            Clean Very severe
## 1721            Clean Very severe
## 1722            Clean Very severe
## 1723            Clean Very severe
## 1724            Clean Very severe
## 1725            Clean Very severe
## 1726            Clean Very severe
## 1727            Clean Very severe
## 1728            Clean Very severe
## 1729            Clean Very severe
## 1730            Clean Very severe
## 1731            Clean Very severe
## 1732 Moderately clean    Moderate
## 1733 Moderately clean    Moderate
## 1734 Moderately clean    Moderate
## 1735 Moderately clean    Moderate
## 1736 Moderately clean    Moderate
## 1737 Moderately clean    Moderate
## 1738 Moderately clean    Moderate
## 1739 Moderately clean    Moderate
## 1740 Moderately clean    Moderate
## 1741 Moderately clean    Moderate
## 1742 Moderately clean    Moderate
## 1743 Moderately clean    Moderate
## 1744 Moderately clean    Moderate
## 1745 Moderately clean    Moderate
## 1746 Moderately clean    Moderate
## 1747 Moderately clean    Moderate
## 1748 Moderately clean    Moderate
## 1749 Moderately clean    Moderate
## 1750 Moderately clean    Moderate
## 1751 Moderately clean    Moderate
## 1752 Moderately clean    Moderate
## 1753 Moderately clean    Moderate
## 1754 Moderately clean    Moderate
## 1755 Moderately clean    Moderate
## 1756 Moderately clean    Moderate
## 1757 Moderately clean    Moderate
## 1758 Moderately clean    Moderate
## 1759 Moderately clean    Moderate
## 1760 Moderately clean    Moderate
## 1761 Moderately clean      Severe
## 1762 Moderately clean      Severe
## 1763 Moderately clean      Severe
## 1764 Moderately clean      Severe
## 1765 Moderately clean      Severe
## 1766 Moderately clean      Severe
## 1767 Moderately clean      Severe
## 1768 Moderately clean      Severe
## 1769 Moderately clean      Severe
## 1770 Moderately clean      Severe
## 1771 Moderately clean      Severe
## 1772 Moderately clean      Severe
## 1773 Moderately clean      Severe
## 1774 Moderately clean      Severe
## 1775 Moderately clean      Severe
## 1776 Moderately clean      Severe
## 1777 Moderately clean      Severe
## 1778 Moderately clean      Severe
## 1779 Moderately clean      Severe
## 1780 Moderately clean      Severe
## 1781 Moderately clean      Severe
## 1782 Moderately clean      Severe
## 1783 Moderately clean      Severe
## 1784 Moderately clean      Severe
## 1785 Moderately clean      Severe
## 1786 Moderately clean      Severe
## 1787 Moderately clean      Severe
## 1788 Moderately clean      Severe
## 1789 Moderately clean      Severe
## 1790            Clean      Severe
## 1791            Clean      Severe
## 1792            Clean      Severe
## 1793            Clean      Severe
## 1794            Clean      Severe
## 1795            Clean      Severe
## 1796            Clean      Severe
## 1797            Clean      Severe
## 1798            Clean      Severe
## 1799            Clean      Severe
## 1800            Clean      Severe
## 1801            Clean      Severe
## 1802            Clean      Severe
## 1803            Clean      Severe
## 1804            Clean      Severe
## 1805            Clean      Severe
## 1806            Clean      Severe
## 1807            Clean      Severe
## 1808            Clean      Severe
## 1809            Clean      Severe
## 1810            Clean      Severe
## 1811            Clean      Severe
## 1812            Clean      Severe
## 1813            Clean      Severe
## 1814            Clean      Severe
## 1815            Clean      Severe
## 1816            Clean      Severe
## 1817            Clean      Severe
## 1818            Clean      Severe
## 1819 Moderately clean    Moderate
## 1820 Moderately clean    Moderate
## 1821 Moderately clean    Moderate
## 1822 Moderately clean    Moderate
## 1823 Moderately clean    Moderate
## 1824 Moderately clean    Moderate
## 1825 Moderately clean    Moderate
## 1826 Moderately clean    Moderate
## 1827 Moderately clean    Moderate
## 1828 Moderately clean    Moderate
## 1829 Moderately clean    Moderate
## 1830 Moderately clean    Moderate
## 1831 Moderately clean    Moderate
## 1832 Moderately clean    Moderate
## 1833 Moderately clean    Moderate
## 1834 Moderately clean    Moderate
## 1835 Moderately clean    Moderate
## 1836 Moderately clean    Moderate
## 1837 Moderately clean    Moderate
## 1838 Moderately clean    Moderate
## 1839 Moderately clean    Moderate
## 1840 Moderately clean    Moderate
## 1841 Moderately clean    Moderate
## 1842 Moderately clean    Moderate
## 1843 Moderately clean    Moderate
## 1844 Moderately clean    Moderate
## 1845 Moderately clean    Moderate
## 1846 Moderately clean    Moderate
## 1847 Moderately clean    Moderate
##                                                                                                                                                                                                                                                              GENERALOBS
## 1                                                                                                                                                                                                                                                                  <NA>
## 2                                                                                                                                                                                                                                                                  <NA>
## 3                                                                                                                                                                                                                                                                  <NA>
## 4                                                                                                                                                                                                                                                                  <NA>
## 5                                                                                                                                                                                                                                                                  <NA>
## 6                                                                                                                                                                                                                                                                  <NA>
## 7                                                                                                                                                                                                                                                                  <NA>
## 8                                                                                                                                                                                                                                                                  <NA>
## 9                                                                                                                                                                                                                                                                  <NA>
## 10                                                                                                                                                                                                                                                                 <NA>
## 11                                                                                                                                                                                                                                                                 <NA>
## 12                                                                                                                                                                                                                                                                 <NA>
## 13                                                                                                                                                                                                                                                                 <NA>
## 14                                                                                                                                                                                                                                                                 <NA>
## 15                                                                                                                                                                                                                                                                 <NA>
## 16                                                                                                                                                                                                                                                                 <NA>
## 17                                                                                                                                                                                                                                                                 <NA>
## 18                                                                                                                                                                                                                                                                 <NA>
## 19                                                                                                                                                                                                                                                                 <NA>
## 20                                                                                                                                                                                                                                                                 <NA>
## 21                                                                                                                                                                                                                                                                 <NA>
## 22                                                                                                                                                                                                                                                                 <NA>
## 23                                                                                                                                                                                                                                                                 <NA>
## 24                                                                                                                                                                                                                                                                 <NA>
## 25                                                                                                                                                                                                                                                                 <NA>
## 26                                                                                                                                                                                                                                                                 <NA>
## 27                                                                                                                                                                                                                                                                 <NA>
## 28                                                                                                                                                                                                                                                                 <NA>
## 29                                                                                                                                                                                                                                                                 <NA>
## 30                                                                                                                                                                                                                                                                 <NA>
## 31                                                                                                                                                                                                                                                                 <NA>
## 32                                                                                                                                                                                                                                                                 <NA>
## 33                                                                                                                                                                                                                                                                 <NA>
## 34                                                                                                                                                                                                                                                                 <NA>
## 35                                                                                                                                                                                                                                                                 <NA>
## 36                                                                                                                                                                                                                                                                 <NA>
## 37                                                                                                                                                                                                                                                                 <NA>
## 38                                                                                                                                                                                                                                                                 <NA>
## 39                                                                                                                                                                                                                                                                 <NA>
## 40                                                                                                                                                                                                                                                                 <NA>
## 41                                                                                                                                                                                                                                                                 <NA>
## 42                                                                                                                                                                                                                                                                 <NA>
## 43                                                                                                                                                                                                                                                                 <NA>
## 44                                                                                                                                                                                                                                                                 <NA>
## 45                                                                                                                                                                                                                                                                 <NA>
## 46                                                                                                                                                                                                                                                                 <NA>
## 47                                                                                                                                                                                                                                                                 <NA>
## 48                                                                                                                                                                                                                                                                 <NA>
## 49                                                                                                                                                                                                                                                                 <NA>
## 50                                                                                                                                                                                                                                                                 <NA>
## 51                                                                                                                                                                                                                                                                 <NA>
## 52                                                                                                                                                                                                                                                                 <NA>
## 53                                                                                                                                                                                                                                                                 <NA>
## 54                                                                                                                                                                                                                                                                 <NA>
## 55                                                                                                                                                                                                                                                                 <NA>
## 56                                                                                                                                                                                                                                                                 <NA>
## 57                                                                                                                                                                                                                                                                 <NA>
## 58                                                                                                                                                                                                                                                                 <NA>
## 59                                                                                                                                                                                                                                                      Concession burn
## 60                                                                                                                                                                                                                                                      Concession burn
## 61                                                                                                                                                                                                                                                      Concession burn
## 62                                                                                                                                                                                                                                                      Concession burn
## 63                                                                                                                                                                                                                                                      Concession burn
## 64                                                                                                                                                                                                                                                      Concession burn
## 65                                                                                                                                                                                                                                                      Concession burn
## 66                                                                                                                                                                                                                                                      Concession burn
## 67                                                                                                                                                                                                                                                      Concession burn
## 68                                                                                                                                                                                                                                                      Concession burn
## 69                                                                                                                                                                                                                                                      Concession burn
## 70                                                                                                                                                                                                                                                      Concession burn
## 71                                                                                                                                                                                                                                                      Concession burn
## 72                                                                                                                                                                                                                                                      Concession burn
## 73                                                                                                                                                                                                                                                      Concession burn
## 74                                                                                                                                                                                                                                                      Concession burn
## 75                                                                                                                                                                                                                                                      Concession burn
## 76                                                                                                                                                                                                                                                      Concession burn
## 77                                                                                                                                                                                                                                                      Concession burn
## 78                                                                                                                                                                                                                                                      Concession burn
## 79                                                                                                                                                                                                                                                      Concession burn
## 80                                                                                                                                                                                                                                                      Concession burn
## 81                                                                                                                                                                                                                                                      Concession burn
## 82                                                                                                                                                                                                                                                      Concession burn
## 83                                                                                                                                                                                                                                                      Concession burn
## 84                                                                                                                                                                                                                                                      Concession burn
## 85                                                                                                                                                                                                                                                      Concession burn
## 86                                                                                                                                                                                                                                                      Concession burn
## 87                                                                                                                                                                                                                                                      Concession burn
## 88                                                                                                                                                                                                                                                                 <NA>
## 89                                                                                                                                                                                                                                                                 <NA>
## 90                                                                                                                                                                                                                                                                 <NA>
## 91                                                                                                                                                                                                                                                                 <NA>
## 92                                                                                                                                                                                                                                                                 <NA>
## 93                                                                                                                                                                                                                                                                 <NA>
## 94                                                                                                                                                                                                                                                                 <NA>
## 95                                                                                                                                                                                                                                                                 <NA>
## 96                                                                                                                                                                                                                                                                 <NA>
## 97                                                                                                                                                                                                                                                                 <NA>
## 98                                                                                                                                                                                                                                                                 <NA>
## 99                                                                                                                                                                                                                                                                 <NA>
## 100                                                                                                                                                                                                                                                                <NA>
## 101                                                                                                                                                                                                                                                                <NA>
## 102                                                                                                                                                                                                                                                                <NA>
## 103                                                                                                                                                                                                                                                                <NA>
## 104                                                                                                                                                                                                                                                                <NA>
## 105                                                                                                                                                                                                                                                                <NA>
## 106                                                                                                                                                                                                                                                                <NA>
## 107                                                                                                                                                                                                                                                                <NA>
## 108                                                                                                                                                                                                                                                                <NA>
## 109                                                                                                                                                                                                                                                                <NA>
## 110                                                                                                                                                                                                                                                                <NA>
## 111                                                                                                                                                                                                                                                                <NA>
## 112                                                                                                                                                                                                                                                                <NA>
## 113                                                                                                                                                                                                                                                                <NA>
## 114                                                                                                                                                                                                                                                                <NA>
## 115                                                                                                                                                                                                                                                                <NA>
## 116                                                                                                                                                                                                                                                                <NA>
## 117                                                                                                                                                                                                                                                                <NA>
## 118                                                                                                                                                                                                                                                                <NA>
## 119                                                                                                                                                                                                                                                                <NA>
## 120                                                                                                                                                                                                                                                                <NA>
## 121                                                                                                                                                                                                                                                                <NA>
## 122                                                                                                                                                                                                                                                                <NA>
## 123                                                                                                                                                                                                                                                                <NA>
## 124                                                                                                                                                                                                                                                                <NA>
## 125                                                                                                                                                                                                                                                                <NA>
## 126                                                                                                                                                                                                                                                                <NA>
## 127                                                                                                                                                                                                                                                                <NA>
## 128                                                                                                                                                                                                                                                                <NA>
## 129                                                                                                                                                                                                                                                                <NA>
## 130                                                                                                                                                                                                                                                                <NA>
## 131                                                                                                                                                                                                                                                                <NA>
## 132                                                                                                                                                                                                                                                                <NA>
## 133                                                                                                                                                                                                                                                                <NA>
## 134                                                                                                                                                                                                                                                                <NA>
## 135                                                                                                                                                                                                                                                                <NA>
## 136                                                                                                                                                                                                                                                                <NA>
## 137                                                                                                                                                                                                                                                                <NA>
## 138                                                                                                                                                                                                                                                                <NA>
## 139                                                                                                                                                                                                                                                                <NA>
## 140                                                                                                                                                                                                                                                                <NA>
## 141                                                                                                                                                                                                                                                                <NA>
## 142                                                                                                                                                                                                                                                                <NA>
## 143                                                                                                                                                                                                                                                                <NA>
## 144                                                                                                                                                                                                                                                                <NA>
## 145                                                                                                                                                                                                                                                                <NA>
## 146                                                                                                                                                                                                                                                          See report
## 147                                                                                                                                                                                                                                                          See report
## 148                                                                                                                                                                                                                                                          See report
## 149                                                                                                                                                                                                                                                          See report
## 150                                                                                                                                                                                                                                                          See report
## 151                                                                                                                                                                                                                                                          See report
## 152                                                                                                                                                                                                                                                          See report
## 153                                                                                                                                                                                                                                                          See report
## 154                                                                                                                                                                                                                                                          See report
## 155                                                                                                                                                                                                                                                          See report
## 156                                                                                                                                                                                                                                                          See report
## 157                                                                                                                                                                                                                                                          See report
## 158                                                                                                                                                                                                                                                          See report
## 159                                                                                                                                                                                                                                                          See report
## 160                                                                                                                                                                                                                                                          See report
## 161                                                                                                                                                                                                                                                          See report
## 162                                                                                                                                                                                                                                                          See report
## 163                                                                                                                                                                                                                                                          See report
## 164                                                                                                                                                                                                                                                          See report
## 165                                                                                                                                                                                                                                                          See report
## 166                                                                                                                                                                                                                                                          See report
## 167                                                                                                                                                                                                                                                          See report
## 168                                                                                                                                                                                                                                                          See report
## 169                                                                                                                                                                                                                                                          See report
## 170                                                                                                                                                                                                                                                          See report
## 171                                                                                                                                                                                                                                                          See report
## 172                                                                                                                                                                                                                                                          See report
## 173                                                                                                                                                                                                                                                          See report
## 174                                                                                                                                                                                                                                                          See report
## 175                                                                                                                                                                                                                                                      Second attempt
## 176                                                                                                                                                                                                                                                      Second attempt
## 177                                                                                                                                                                                                                                                      Second attempt
## 178                                                                                                                                                                                                                                                      Second attempt
## 179                                                                                                                                                                                                                                                      Second attempt
## 180                                                                                                                                                                                                                                                      Second attempt
## 181                                                                                                                                                                                                                                                      Second attempt
## 182                                                                                                                                                                                                                                                      Second attempt
## 183                                                                                                                                                                                                                                                      Second attempt
## 184                                                                                                                                                                                                                                                      Second attempt
## 185                                                                                                                                                                                                                                                      Second attempt
## 186                                                                                                                                                                                                                                                      Second attempt
## 187                                                                                                                                                                                                                                                      Second attempt
## 188                                                                                                                                                                                                                                                      Second attempt
## 189                                                                                                                                                                                                                                                      Second attempt
## 190                                                                                                                                                                                                                                                      Second attempt
## 191                                                                                                                                                                                                                                                      Second attempt
## 192                                                                                                                                                                                                                                                      Second attempt
## 193                                                                                                                                                                                                                                                      Second attempt
## 194                                                                                                                                                                                                                                                      Second attempt
## 195                                                                                                                                                                                                                                                      Second attempt
## 196                                                                                                                                                                                                                                                      Second attempt
## 197                                                                                                                                                                                                                                                      Second attempt
## 198                                                                                                                                                                                                                                                      Second attempt
## 199                                                                                                                                                                                                                                                      Second attempt
## 200                                                                                                                                                                                                                                                      Second attempt
## 201                                                                                                                                                                                                                                                      Second attempt
## 202                                                                                                                                                                                                                                                      Second attempt
## 203                                                                                                                                                                                                                                                      Second attempt
## 204                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 205                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 206                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 207                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 208                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 209                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 210                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 211                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 212                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 213                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 214                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 215                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 216                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 217                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 218                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 219                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 220                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 221                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 222                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 223                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 224                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 225                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 226                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 227                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 228                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 229                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 230                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 231                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 232                                                                                                                                                                                       This was the first block that a point ignition fire was set in. Moderate burn
## 233                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 234                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 235                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 236                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 237                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 238                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 239                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 240                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 241                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 242                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 243                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 244                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 245                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 246                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 247                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 248                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 249                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 250                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 251                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 252                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 253                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 254                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 255                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 256                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 257                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 258                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 259                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 260                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 261                                                                                                                                                                          Grass cover was moderate even if there was overgrazing and trampling. Area burnt was clean
## 262                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 263                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 264                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 265                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 266                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 267                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 268                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 269                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 270                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 271                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 272                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 273                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 274                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 275                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 276                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 277                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 278                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 279                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 280                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 281                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 282                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 283                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 284                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 285                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 286                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 287                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 288                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 289                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 290                                                                                                                           The fire burnt nice and clean on the open plains and patchly within the valleys and marshes. Dry material burnt 100%, especially grasses.
## 291                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 292                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 293                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 294                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 295                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 296                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 297                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 298                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 299                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 300                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 301                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 302                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 303                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 304                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 305                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 306                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 307                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 308                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 309                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 310                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 311                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 312                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 313                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 314                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 315                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 316                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 317                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 318                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 319                                                                      This fire took the whole week to burn block 95 with the help of wind changing to all directions on a daily basis.  Block 103 only burnt a piece and extinguished inbetween two adjacent rivers
## 320                                                                                                                                                                                                                                                                <NA>
## 321                                                                                                                                                                                                                                                                <NA>
## 322                                                                                                                                                                                                                                                                <NA>
## 323                                                                                                                                                                                                                                                                <NA>
## 324                                                                                                                                                                                                                                                                <NA>
## 325                                                                                                                                                                                                                                                                <NA>
## 326                                                                                                                                                                                                                                                                <NA>
## 327                                                                                                                                                                                                                                                                <NA>
## 328                                                                                                                                                                                                                                                                <NA>
## 329                                                                                                                                                                                                                                                                <NA>
## 330                                                                                                                                                                                                                                                                <NA>
## 331                                                                                                                                                                                                                                                                <NA>
## 332                                                                                                                                                                                                                                                                <NA>
## 333                                                                                                                                                                                                                                                                <NA>
## 334                                                                                                                                                                                                                                                                <NA>
## 335                                                                                                                                                                                                                                                                <NA>
## 336                                                                                                                                                                                                                                                                <NA>
## 337                                                                                                                                                                                                                                                                <NA>
## 338                                                                                                                                                                                                                                                                <NA>
## 339                                                                                                                                                                                                                                                                <NA>
## 340                                                                                                                                                                                                                                                                <NA>
## 341                                                                                                                                                                                                                                                                <NA>
## 342                                                                                                                                                                                                                                                                <NA>
## 343                                                                                                                                                                                                                                                                <NA>
## 344                                                                                                                                                                                                                                                                <NA>
## 345                                                                                                                                                                                                                                                                <NA>
## 346                                                                                                                                                                                                                                                                <NA>
## 347                                                                                                                                                                                                                                                                <NA>
## 348                                                                                                                                                                                                                                                                <NA>
## 349                                                                                                                                                                                                                                                                <NA>
## 350                                                                                                                                                                                                                                                                <NA>
## 351                                                                                                                                                                                                                                                                <NA>
## 352                                                                                                                                                                                                                                                                <NA>
## 353                                                                                                                                                                                                                                                                <NA>
## 354                                                                                                                                                                                                                                                                <NA>
## 355                                                                                                                                                                                                                                                                <NA>
## 356                                                                                                                                                                                                                                                                <NA>
## 357                                                                                                                                                                                                                                                                <NA>
## 358                                                                                                                                                                                                                                                                <NA>
## 359                                                                                                                                                                                                                                                                <NA>
## 360                                                                                                                                                                                                                                                                <NA>
## 361                                                                                                                                                                                                                                                                <NA>
## 362                                                                                                                                                                                                                                                                <NA>
## 363                                                                                                                                                                                                                                                                <NA>
## 364                                                                                                                                                                                                                                                                <NA>
## 365                                                                                                                                                                                                                                                                <NA>
## 366                                                                                                                                                                                                                                                                <NA>
## 367                                                                                                                                                                                                                                                                <NA>
## 368                                                                                                                                                                                                                                                                <NA>
## 369                                                                                                                                                                                                                                                                <NA>
## 370                                                                                                                                                                                                                                                                <NA>
## 371                                                                                                                                                                                                                                                                <NA>
## 372                                                                                                                                                                                                                                                                <NA>
## 373                                                                                                                                                                                                                                                                <NA>
## 374                                                                                                                                                                                                                                                                <NA>
## 375                                                                                                                                                                                                                                                                <NA>
## 376                                                                                                                                                                                                                                                                <NA>
## 377                                                                                                                                                                                                                                                                <NA>
## 378                                                                                                                                                                                                                   Please see original fire report for further notes
## 379                                                                                                                                                                                                                   Please see original fire report for further notes
## 380                                                                                                                                                                                                                   Please see original fire report for further notes
## 381                                                                                                                                                                                                                   Please see original fire report for further notes
## 382                                                                                                                                                                                                                   Please see original fire report for further notes
## 383                                                                                                                                                                                                                   Please see original fire report for further notes
## 384                                                                                                                                                                                                                   Please see original fire report for further notes
## 385                                                                                                                                                                                                                   Please see original fire report for further notes
## 386                                                                                                                                                                                                                   Please see original fire report for further notes
## 387                                                                                                                                                                                                                   Please see original fire report for further notes
## 388                                                                                                                                                                                                                   Please see original fire report for further notes
## 389                                                                                                                                                                                                                   Please see original fire report for further notes
## 390                                                                                                                                                                                                                   Please see original fire report for further notes
## 391                                                                                                                                                                                                                   Please see original fire report for further notes
## 392                                                                                                                                                                                                                   Please see original fire report for further notes
## 393                                                                                                                                                                                                                   Please see original fire report for further notes
## 394                                                                                                                                                                                                                   Please see original fire report for further notes
## 395                                                                                                                                                                                                                   Please see original fire report for further notes
## 396                                                                                                                                                                                                                   Please see original fire report for further notes
## 397                                                                                                                                                                                                                   Please see original fire report for further notes
## 398                                                                                                                                                                                                                   Please see original fire report for further notes
## 399                                                                                                                                                                                                                   Please see original fire report for further notes
## 400                                                                                                                                                                                                                   Please see original fire report for further notes
## 401                                                                                                                                                                                                                   Please see original fire report for further notes
## 402                                                                                                                                                                                                                   Please see original fire report for further notes
## 403                                                                                                                                                                                                                   Please see original fire report for further notes
## 404                                                                                                                                                                                                                   Please see original fire report for further notes
## 405                                                                                                                                                                                                                   Please see original fire report for further notes
## 406                                                                                                                                                                                                                   Please see original fire report for further notes
## 407                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 408                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 409                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 410                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 411                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 412                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 413                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 414                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 415                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 416                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 417                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 418                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 419                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 420                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 421                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 422                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 423                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 424                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 425                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 426                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 427                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 428                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 429                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 430                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 431                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 432                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 433                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 434                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 435                                                                                                                                                                                                                            Area did burn cleanly te previous season
## 436                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 437                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 438                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 439                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 440                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 441                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 442                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 443                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 444                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 445                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 446                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 447                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 448                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 449                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 450                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 451                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 452                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 453                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 454                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 455                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 456                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 457                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 458                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 459                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 460                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 461                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 462                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 463                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 464                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 465                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 466                                                                                                                                                                                                      First patch burn since the introduction of the new fire policy
## 467  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 468  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 469  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 470  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 471  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 472  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 473  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 474  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 475  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 476  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 477  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 478  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 479  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 480  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 481  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 482  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 483  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 484  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 485  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 486  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 487  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 488  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 489  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 490  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 491  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 492  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 493  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 494  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 495  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 496  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 497  Original burn started on the night of the 22/09/2002 and say of 23/09/2002. Very moderate burn. Hot north - west winds on the 22/09/2002 caused hot burns over the rest of the block. \r\n\r\nThe fire in Block N138 was a management (firebreak) that was put in.
## 498                                                                                                                                                                                                                                                                <NA>
## 499                                                                                                                                                                                                                                                                <NA>
## 500                                                                                                                                                                                                                                                                <NA>
## 501                                                                                                                                                                                                                                                                <NA>
## 502                                                                                                                                                                                                                                                                <NA>
## 503                                                                                                                                                                                                                                                                <NA>
## 504                                                                                                                                                                                                                                                                <NA>
## 505                                                                                                                                                                                                                                                                <NA>
## 506                                                                                                                                                                                                                                                                <NA>
## 507                                                                                                                                                                                                                                                                <NA>
## 508                                                                                                                                                                                                                                                                <NA>
## 509                                                                                                                                                                                                                                                                <NA>
## 510                                                                                                                                                                                                                                                                <NA>
## 511                                                                                                                                                                                                                                                                <NA>
## 512                                                                                                                                                                                                                                                                <NA>
## 513                                                                                                                                                                                                                                                                <NA>
## 514                                                                                                                                                                                                                                                                <NA>
## 515                                                                                                                                                                                                                                                                <NA>
## 516                                                                                                                                                                                                                                                                <NA>
## 517                                                                                                                                                                                                                                                                <NA>
## 518                                                                                                                                                                                                                                                                <NA>
## 519                                                                                                                                                                                                                                                                <NA>
## 520                                                                                                                                                                                                                                                                <NA>
## 521                                                                                                                                                                                                                                                                <NA>
## 522                                                                                                                                                                                                                                                                <NA>
## 523                                                                                                                                                                                                                                                                <NA>
## 524                                                                                                                                                                                                                                                                <NA>
## 525                                                                                                                                                                                                                                                                <NA>
## 526                                                                                                                                                                                                                                                                <NA>
## 527                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 528                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 529                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 530                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 531                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 532                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 533                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 534                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 535                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 536                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 537                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 538                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 539                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 540                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 541                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 542                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 543                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 544                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 545                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 546                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 547                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 548                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 549                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 550                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 551                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 552                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 553                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 554                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 555                                                                                                                                                                                                          The blocks that were involved were not burnt for 3 seasons
## 556                                                                                                                                                                                                                                                                <NA>
## 557                                                                                                                                                                                                                                                                <NA>
## 558                                                                                                                                                                                                                                                                <NA>
## 559                                                                                                                                                                                                                                                                <NA>
## 560                                                                                                                                                                                                                                                                <NA>
## 561                                                                                                                                                                                                                                                                <NA>
## 562                                                                                                                                                                                                                                                                <NA>
## 563                                                                                                                                                                                                                                                                <NA>
## 564                                                                                                                                                                                                                                                                <NA>
## 565                                                                                                                                                                                                                                                                <NA>
## 566                                                                                                                                                                                                                                                                <NA>
## 567                                                                                                                                                                                                                                                                <NA>
## 568                                                                                                                                                                                                                                                                <NA>
## 569                                                                                                                                                                                                                                                                <NA>
## 570                                                                                                                                                                                                                                                                <NA>
## 571                                                                                                                                                                                                                                                                <NA>
## 572                                                                                                                                                                                                                                                                <NA>
## 573                                                                                                                                                                                                                                                                <NA>
## 574                                                                                                                                                                                                                                                                <NA>
## 575                                                                                                                                                                                                                                                                <NA>
## 576                                                                                                                                                                                                                                                                <NA>
## 577                                                                                                                                                                                                                                                                <NA>
## 578                                                                                                                                                                                                                                                                <NA>
## 579                                                                                                                                                                                                                                                                <NA>
## 580                                                                                                                                                                                                                                                                <NA>
## 581                                                                                                                                                                                                                                                                <NA>
## 582                                                                                                                                                                                                                                                                <NA>
## 583                                                                                                                                                                                                                                                                <NA>
## 584                                                                                                                                                                                                                                                                <NA>
## 585                                                                                                                                                                                                                                                                <NA>
## 586                                                                                                                                                                                                                                                                <NA>
## 587                                                                                                                                                                                                                                                                <NA>
## 588                                                                                                                                                                                                                                                                <NA>
## 589                                                                                                                                                                                                                                                                <NA>
## 590                                                                                                                                                                                                                                                                <NA>
## 591                                                                                                                                                                                                                                                                <NA>
## 592                                                                                                                                                                                                                                                                <NA>
## 593                                                                                                                                                                                                                                                                <NA>
## 594                                                                                                                                                                                                                                                                <NA>
## 595                                                                                                                                                                                                                                                                <NA>
## 596                                                                                                                                                                                                                                                                <NA>
## 597                                                                                                                                                                                                                                                                <NA>
## 598                                                                                                                                                                                                                                                                <NA>
## 599                                                                                                                                                                                                                                                                <NA>
## 600                                                                                                                                                                                                                                                                <NA>
## 601                                                                                                                                                                                                                                                                <NA>
## 602                                                                                                                                                                                                                                                                <NA>
## 603                                                                                                                                                                                                                                                                <NA>
## 604                                                                                                                                                                                                                                                                <NA>
## 605                                                                                                                                                                                                                                                                <NA>
## 606                                                                                                                                                                                                                                                                <NA>
## 607                                                                                                                                                                                                                                                                <NA>
## 608                                                                                                                                                                                                                                                                <NA>
## 609                                                                                                                                                                                                                                                                <NA>
## 610                                                                                                                                                                                                                                                                <NA>
## 611                                                                                                                                                                                                                                                                <NA>
## 612                                                                                                                                                                                                                                                                <NA>
## 613                                                                                                                                                                                                                                                                <NA>
## 614                                                                                                                                                                                                                                                                <NA>
## 615                                                                                                                                                                                                                                                                <NA>
## 616                                                                                                                                                                                                                                                                <NA>
## 617                                                                                                                                                                                                                                                                <NA>
## 618                                                                                                                                                                                                                                                                <NA>
## 619                                                                                                                                                                                                                                                                <NA>
## 620                                                                                                                                                                                                                                                                <NA>
## 621                                                                                                                                                                                                                                                                <NA>
## 622                                                                                                                                                                                                                                                                <NA>
## 623                                                                                                                                                                                                                                                                <NA>
## 624                                                                                                                                                                                                                                                                <NA>
## 625                                                                                                                                                                                                                                                                <NA>
## 626                                                                                                                                                                                                                                                                <NA>
## 627                                                                                                                                                                                                                                                                <NA>
## 628                                                                                                                                                                                                                                                                <NA>
## 629                                                                                                                                                                                                                                                                <NA>
## 630                                                                                                                                                                                                                                                                <NA>
## 631                                                                                                                                                                                                                                                                <NA>
## 632                                                                                                                                                                                                                                                                <NA>
## 633                                                                                                                                                                                                                                                                <NA>
## 634                                                                                                                                                                                                                                                                <NA>
## 635                                                                                                                                                                                                                                                                <NA>
## 636                                                                                                                                                                                                                                                                <NA>
## 637                                                                                                                                                                                                                                                                <NA>
## 638                                                                                                                                                                                                                                                                <NA>
## 639                                                                                                                                                                                                                                                                <NA>
## 640                                                                                                                                                                                                                                                                <NA>
## 641                                                                                                                                                                                                                                                                <NA>
## 642                                                                                                                                                                                                                                                                <NA>
## 643                                                                                                                                                                                                                                                                <NA>
## 644                                                                                                                                                                                                                                                                <NA>
## 645                                                                                                                                                                                                                                                                <NA>
## 646                                                                                                                                                                                                                                                                <NA>
## 647                                                                                                                                                                                                                                                                <NA>
## 648                                                                                                                                                                                                                                                                <NA>
## 649                                                                                                                                                                                                                                                                <NA>
## 650                                                                                                                                                                                                                                                                <NA>
## 651                                                                                                                                                                                                                                                                <NA>
## 652                                                                                                                                                                                                                                                                <NA>
## 653                                                                                                                                                                                                                                                                <NA>
## 654                                                                                                                                                                                                                                                                <NA>
## 655                                                                                                                                                                                                                                                                <NA>
## 656                                                                                                                                                                                                                                                                <NA>
## 657                                                                                                                                                                                                                                                                <NA>
## 658                                                                                                                                                                                                                                                                <NA>
## 659                                                                                                                                                                                                                                                                <NA>
## 660                                                                                                                                                                                                                                                                <NA>
## 661                                                                                                                                                                                                                                                                <NA>
## 662                                                                                                                                                                                                                                                                <NA>
## 663                                                                                                                                                                                                                                                                <NA>
## 664                                                                                                                                                                                                                                                                <NA>
## 665                                                                                                                                                                                                                                                                <NA>
## 666                                                                                                                                                                                                                                                                <NA>
## 667                                                                                                                                                                                                                                                                <NA>
## 668                                                                                                                                                                                                                                                                <NA>
## 669                                                                                                                                                                                                                                                                <NA>
## 670                                                                                                                                                                                                                                                                <NA>
## 671                                                                                                                                                                                                                                                                <NA>
## 672                                                                                                                                                                                                                                                    Good Patchy burn
## 673                                                                                                                                                                                                                                                    Good Patchy burn
## 674                                                                                                                                                                                                                                                    Good Patchy burn
## 675                                                                                                                                                                                                                                                    Good Patchy burn
## 676                                                                                                                                                                                                                                                    Good Patchy burn
## 677                                                                                                                                                                                                                                                    Good Patchy burn
## 678                                                                                                                                                                                                                                                    Good Patchy burn
## 679                                                                                                                                                                                                                                                    Good Patchy burn
## 680                                                                                                                                                                                                                                                    Good Patchy burn
## 681                                                                                                                                                                                                                                                    Good Patchy burn
## 682                                                                                                                                                                                                                                                    Good Patchy burn
## 683                                                                                                                                                                                                                                                    Good Patchy burn
## 684                                                                                                                                                                                                                                                    Good Patchy burn
## 685                                                                                                                                                                                                                                                    Good Patchy burn
## 686                                                                                                                                                                                                                                                    Good Patchy burn
## 687                                                                                                                                                                                                                                                    Good Patchy burn
## 688                                                                                                                                                                                                                                                    Good Patchy burn
## 689                                                                                                                                                                                                                                                    Good Patchy burn
## 690                                                                                                                                                                                                                                                    Good Patchy burn
## 691                                                                                                                                                                                                                                                    Good Patchy burn
## 692                                                                                                                                                                                                                                                    Good Patchy burn
## 693                                                                                                                                                                                                                                                    Good Patchy burn
## 694                                                                                                                                                                                                                                                    Good Patchy burn
## 695                                                                                                                                                                                                                                                    Good Patchy burn
## 696                                                                                                                                                                                                                                                    Good Patchy burn
## 697                                                                                                                                                                                                                                                    Good Patchy burn
## 698                                                                                                                                                                                                                                                    Good Patchy burn
## 699                                                                                                                                                                                                                                                    Good Patchy burn
## 700                                                                                                                                                                                                                                                    Good Patchy burn
## 701                                                                                                                                                                                                                                                                <NA>
## 702                                                                                                                                                                                                                                                                <NA>
## 703                                                                                                                                                                                                                                                                <NA>
## 704                                                                                                                                                                                                                                                                <NA>
## 705                                                                                                                                                                                                                                                                <NA>
## 706                                                                                                                                                                                                                                                                <NA>
## 707                                                                                                                                                                                                                                                                <NA>
## 708                                                                                                                                                                                                                                                                <NA>
## 709                                                                                                                                                                                                                                                                <NA>
## 710                                                                                                                                                                                                                                                                <NA>
## 711                                                                                                                                                                                                                                                                <NA>
## 712                                                                                                                                                                                                                                                                <NA>
## 713                                                                                                                                                                                                                                                                <NA>
## 714                                                                                                                                                                                                                                                                <NA>
## 715                                                                                                                                                                                                                                                                <NA>
## 716                                                                                                                                                                                                                                                                <NA>
## 717                                                                                                                                                                                                                                                                <NA>
## 718                                                                                                                                                                                                                                                                <NA>
## 719                                                                                                                                                                                                                                                                <NA>
## 720                                                                                                                                                                                                                                                                <NA>
## 721                                                                                                                                                                                                                                                                <NA>
## 722                                                                                                                                                                                                                                                                <NA>
## 723                                                                                                                                                                                                                                                                <NA>
## 724                                                                                                                                                                                                                                                                <NA>
## 725                                                                                                                                                                                                                                                                <NA>
## 726                                                                                                                                                                                                                                                                <NA>
## 727                                                                                                                                                                                                                                                                <NA>
## 728                                                                                                                                                                                                                                                                <NA>
## 729                                                                                                                                                                                                                                                                <NA>
## 730                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 731                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 732                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 733                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 734                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 735                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 736                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 737                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 738                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 739                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 740                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 741                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 742                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 743                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 744                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 745                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 746                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 747                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 748                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 749                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 750                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 751                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 752                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 753                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 754                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 755                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 756                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 757                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 758                                                                                                                                           No top kill. Fire lasted for a few days with differenct intensities\r\nFire burnt from the 29/08/2002 till the 02/09/2002
## 759                                                                                                                                                                                                                                                                <NA>
## 760                                                                                                                                                                                                                                                                <NA>
## 761                                                                                                                                                                                                                                                                <NA>
## 762                                                                                                                                                                                                                                                                <NA>
## 763                                                                                                                                                                                                                                                                <NA>
## 764                                                                                                                                                                                                                                                                <NA>
## 765                                                                                                                                                                                                                                                                <NA>
## 766                                                                                                                                                                                                                                                                <NA>
## 767                                                                                                                                                                                                                                                                <NA>
## 768                                                                                                                                                                                                                                                                <NA>
## 769                                                                                                                                                                                                                                                                <NA>
## 770                                                                                                                                                                                                                                                                <NA>
## 771                                                                                                                                                                                                                                                                <NA>
## 772                                                                                                                                                                                                                                                                <NA>
## 773                                                                                                                                                                                                                                                                <NA>
## 774                                                                                                                                                                                                                                                                <NA>
## 775                                                                                                                                                                                                                                                                <NA>
## 776                                                                                                                                                                                                                                                                <NA>
## 777                                                                                                                                                                                                                                                                <NA>
## 778                                                                                                                                                                                                                                                                <NA>
## 779                                                                                                                                                                                                                                                                <NA>
## 780                                                                                                                                                                                                                                                                <NA>
## 781                                                                                                                                                                                                                                                                <NA>
## 782                                                                                                                                                                                                                                                                <NA>
## 783                                                                                                                                                                                                                                                                <NA>
## 784                                                                                                                                                                                                                                                                <NA>
## 785                                                                                                                                                                                                                                                                <NA>
## 786                                                                                                                                                                                                                                                                <NA>
## 787                                                                                                                                                                                                                                                                <NA>
## 788                                                                                                                                                                                                                                                                <NA>
## 789                                                                                                                                                                                                                                                                <NA>
## 790                                                                                                                                                                                                                                                                <NA>
## 791                                                                                                                                                                                                                                                                <NA>
## 792                                                                                                                                                                                                                                                                <NA>
## 793                                                                                                                                                                                                                                                                <NA>
## 794                                                                                                                                                                                                                                                                <NA>
## 795                                                                                                                                                                                                                                                                <NA>
## 796                                                                                                                                                                                                                                                                <NA>
## 797                                                                                                                                                                                                                                                                <NA>
## 798                                                                                                                                                                                                                                                                <NA>
## 799                                                                                                                                                                                                                                                                <NA>
## 800                                                                                                                                                                                                                                                                <NA>
## 801                                                                                                                                                                                                                                                                <NA>
## 802                                                                                                                                                                                                                                                                <NA>
## 803                                                                                                                                                                                                                                                                <NA>
## 804                                                                                                                                                                                                                                                                <NA>
## 805                                                                                                                                                                                                                                                                <NA>
## 806                                                                                                                                                                                                                                                                <NA>
## 807                                                                                                                                                                                                                                                                <NA>
## 808                                                                                                                                                                                                                                                                <NA>
## 809                                                                                                                                                                                                                                                                <NA>
## 810                                                                                                                                                                                                                                                                <NA>
## 811                                                                                                                                                                                                                                                                <NA>
## 812                                                                                                                                                                                                                                                                <NA>
## 813                                                                                                                                                                                                                                                                <NA>
## 814                                                                                                                                                                                                                                                                <NA>
## 815                                                                                                                                                                                                                                                                <NA>
## 816                                                                                                                                                                                                                                                                <NA>
## 817                                                                                                                                                                                                                                                                <NA>
## 818                                                                                                                                                                                                                                                                <NA>
## 819                                                                                                                                                                                                                                                                <NA>
## 820                                                                                                                                                                                                                                                                <NA>
## 821                                                                                                                                                                                                                                                                <NA>
## 822                                                                                                                                                                                                                                                                <NA>
## 823                                                                                                                                                                                                                                                                <NA>
## 824                                                                                                                                                                                                                                                                <NA>
## 825                                                                                                                                                                                                                                                                <NA>
## 826                                                                                                                                                                                                                                                                <NA>
## 827                                                                                                                                                                                                                                                                <NA>
## 828                                                                                                                                                                                                                                                                <NA>
## 829                                                                                                                                                                                                                                                                <NA>
## 830                                                                                                                                                                                                                                                                <NA>
## 831                                                                                                                                                                                                                                                                <NA>
## 832                                                                                                                                                                                                                                                                <NA>
## 833                                                                                                                                                                                                                                                                <NA>
## 834                                                                                                                                                                                                                                                                <NA>
## 835                                                                                                                                                                                                                                                                <NA>
## 836                                                                                                                                                                                                                                                                <NA>
## 837                                                                                                                                                                                                                                                                <NA>
## 838                                                                                                                                                                                                                                                                <NA>
## 839                                                                                                                                                                                                                                                                <NA>
## 840                                                                                                                                                                                                                                                                <NA>
## 841                                                                                                                                                                                                                                                                <NA>
## 842                                                                                                                                                                                                                                                                <NA>
## 843                                                                                                                                                                                                                                                                <NA>
## 844                                                                                                                                                                                                                                                                <NA>
## 845                                                                                                                                                                                                                                                                <NA>
## 846                                                                                                                                                                                                                                                                <NA>
## 847                                                                                                                                                                                                                                                                <NA>
## 848                                                                                                                                                                                                                                                                <NA>
## 849                                                                                                                                                                                                                                                                <NA>
## 850                                                                                                                                                                                                                                                                <NA>
## 851                                                                                                                                                                                                                                                                <NA>
## 852                                                                                                                                                                                                                                                                <NA>
## 853                                                                                                                                                                                                                                                                <NA>
## 854                                                                                                                                                                                                                                                                <NA>
## 855                                                                                                                                                                                                                                                                <NA>
## 856                                                                                                                                                                                                                                                                <NA>
## 857                                                                                                                                                                                                                                                                <NA>
## 858                                                                                                                                                                                                                                                                <NA>
## 859                                                                                                                                                                                                                                                                <NA>
## 860                                                                                                                                                                                                                                                                <NA>
## 861                                                                                                                                                                                                                                                                <NA>
## 862                                                                                                                                                                                                                                                                <NA>
## 863                                                                                                                                                                                                                                                                <NA>
## 864                                                                                                                                                                                                                                                                <NA>
## 865                                                                                                                                                                                                                                                                <NA>
## 866                                                                                                                                                                                                                                                                <NA>
## 867                                                                                                                                                                                                                                                                <NA>
## 868                                                                                                                                                                                                                                                                <NA>
## 869                                                                                                                                                                                                                                                                <NA>
## 870                                                                                                                                                                                                                                                                <NA>
## 871                                                                                                                                                                                                                                                                <NA>
## 872                                                                                                                                                                                                                                                                <NA>
## 873                                                                                                                                                                                                                                                                <NA>
## 874                                                                                                                                                                                                                                                                <NA>
## 875                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 876                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 877                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 878                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 879                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 880                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 881                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 882                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 883                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 884                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 885                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 886                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 887                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 888                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 889                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 890                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 891                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 892                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 893                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 894                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 895                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 896                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 897                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 898                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 899                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 900                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 901                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 902                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 903                                                                                                                                                                                  Fire from the army base that jumped into the park, caused by professional hunters.
## 904                                                                                                                                                                                                                                                                <NA>
## 905                                                                                                                                                                                                                                                                <NA>
## 906                                                                                                                                                                                                                                                                <NA>
## 907                                                                                                                                                                                                                                                                <NA>
## 908                                                                                                                                                                                                                                                                <NA>
## 909                                                                                                                                                                                                                                                                <NA>
## 910                                                                                                                                                                                                                                                                <NA>
## 911                                                                                                                                                                                                                                                                <NA>
## 912                                                                                                                                                                                                                                                                <NA>
## 913                                                                                                                                                                                                                                                                <NA>
## 914                                                                                                                                                                                                                                                                <NA>
## 915                                                                                                                                                                                                                                                                <NA>
## 916                                                                                                                                                                                                                                                                <NA>
## 917                                                                                                                                                                                                                                                                <NA>
## 918                                                                                                                                                                                                                                                                <NA>
## 919                                                                                                                                                                                                                                                                <NA>
## 920                                                                                                                                                                                                                                                                <NA>
## 921                                                                                                                                                                                                                                                                <NA>
## 922                                                                                                                                                                                                                                                                <NA>
## 923                                                                                                                                                                                                                                                                <NA>
## 924                                                                                                                                                                                                                                                                <NA>
## 925                                                                                                                                                                                                                                                                <NA>
## 926                                                                                                                                                                                                                                                                <NA>
## 927                                                                                                                                                                                                                                                                <NA>
## 928                                                                                                                                                                                                                                                                <NA>
## 929                                                                                                                                                                                                                                                                <NA>
## 930                                                                                                                                                                                                                                                                <NA>
## 931                                                                                                                                                                                                                                                                <NA>
## 932                                                                                                                                                                                                                                                                <NA>
## 933                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 934                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 935                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 936                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 937                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 938                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 939                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 940                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 941                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 942                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 943                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 944                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 945                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 946                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 947                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 948                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 949                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 950                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 951                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 952                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 953                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 954                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 955                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 956                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 957                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 958                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 959                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 960                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 961                                                                                                                                        The fire was moving very fast because of strong winds, most of the blocks remained unburnt. Damages caused outside the park.
## 962                                                                                                                                                                                                                                                     Nhlangwine camp
## 963                                                                                                                                                                                                                                                     Nhlangwine camp
## 964                                                                                                                                                                                                                                                     Nhlangwine camp
## 965                                                                                                                                                                                                                                                     Nhlangwine camp
## 966                                                                                                                                                                                                                                                     Nhlangwine camp
## 967                                                                                                                                                                                                                                                     Nhlangwine camp
## 968                                                                                                                                                                                                                                                     Nhlangwine camp
## 969                                                                                                                                                                                                                                                     Nhlangwine camp
## 970                                                                                                                                                                                                                                                     Nhlangwine camp
## 971                                                                                                                                                                                                                                                     Nhlangwine camp
## 972                                                                                                                                                                                                                                                     Nhlangwine camp
## 973                                                                                                                                                                                                                                                     Nhlangwine camp
## 974                                                                                                                                                                                                                                                     Nhlangwine camp
## 975                                                                                                                                                                                                                                                     Nhlangwine camp
## 976                                                                                                                                                                                                                                                     Nhlangwine camp
## 977                                                                                                                                                                                                                                                     Nhlangwine camp
## 978                                                                                                                                                                                                                                                     Nhlangwine camp
## 979                                                                                                                                                                                                                                                     Nhlangwine camp
## 980                                                                                                                                                                                                                                                     Nhlangwine camp
## 981                                                                                                                                                                                                                                                     Nhlangwine camp
## 982                                                                                                                                                                                                                                                     Nhlangwine camp
## 983                                                                                                                                                                                                                                                     Nhlangwine camp
## 984                                                                                                                                                                                                                                                     Nhlangwine camp
## 985                                                                                                                                                                                                                                                     Nhlangwine camp
## 986                                                                                                                                                                                                                                                     Nhlangwine camp
## 987                                                                                                                                                                                                                                                     Nhlangwine camp
## 988                                                                                                                                                                                                                                                     Nhlangwine camp
## 989                                                                                                                                                                                                                                                     Nhlangwine camp
## 990                                                                                                                                                                                                                                                     Nhlangwine camp
## 991                                                                                                                                                                                                                                                     Nhlangwine camp
## 992                                                                                                                                                                                                                                                     Nhlangwine camp
## 993                                                                                                                                                                                                                                                                <NA>
## 994                                                                                                                                                                                                                                                                <NA>
## 995                                                                                                                                                                                                                                                                <NA>
## 996                                                                                                                                                                                                                                                                <NA>
## 997                                                                                                                                                                                                                                                                <NA>
## 998                                                                                                                                                                                                                                                                <NA>
## 999                                                                                                                                                                                                                                                                <NA>
## 1000                                                                                                                                                                                                                                                               <NA>
## 1001                                                                                                                                                                                                                                                               <NA>
## 1002                                                                                                                                                                                                                                                               <NA>
## 1003                                                                                                                                                                                                                                                               <NA>
## 1004                                                                                                                                                                                                                                                               <NA>
## 1005                                                                                                                                                                                                                                                               <NA>
## 1006                                                                                                                                                                                                                                                               <NA>
## 1007                                                                                                                                                                                                                                                               <NA>
## 1008                                                                                                                                                                                                                                                               <NA>
## 1009                                                                                                                                                                                                                                                               <NA>
## 1010                                                                                                                                                                                                                                                               <NA>
## 1011                                                                                                                                                                                                                                                               <NA>
## 1012                                                                                                                                                                                                                                                               <NA>
## 1013                                                                                                                                                                                                                                                               <NA>
## 1014                                                                                                                                                                                                                                                               <NA>
## 1015                                                                                                                                                                                                                                                               <NA>
## 1016                                                                                                                                                                                                                                                               <NA>
## 1017                                                                                                                                                                                                                                                               <NA>
## 1018                                                                                                                                                                                                                                                               <NA>
## 1019                                                                                                                                                                                                                                                               <NA>
## 1020                                                                                                                                                                                                                                                               <NA>
## 1021                                                                                                                                                                                                                                                               <NA>
## 1022                                                                                                                                                                                                                                                               <NA>
## 1023                                                                                                                                                                                                                                                               <NA>
## 1024                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1025                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1026                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1027                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1028                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1029                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1030                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1031                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1032                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1033                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1034                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1035                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1036                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1037                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1038                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1039                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1040                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1041                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1042                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1043                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1044                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1045                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1046                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1047                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1048                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1049                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1050                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1051                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1052                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1053                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1054                                                                                                                                                                          A number of plots catched fire the fire that started from block S042 with westerly winds.
## 1055                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1056                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1057                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1058                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1059                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1060                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1061                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1062                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1063                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1064                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1065                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1066                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1067                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1068                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1069                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1070                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1071                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1072                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1073                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1074                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1075                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1076                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1077                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1078                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1079                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1080                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1081                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1082                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1083                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1084                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1085                                                                                                                                                                                            Area was completely blackened. Middle portion only of block S020B burnt
## 1086                                                                                                                                                                                                                              The fire burnt along the railway line
## 1087                                                                                                                                                                                                                              The fire burnt along the railway line
## 1088                                                                                                                                                                                                                              The fire burnt along the railway line
## 1089                                                                                                                                                                                                                              The fire burnt along the railway line
## 1090                                                                                                                                                                                                                              The fire burnt along the railway line
## 1091                                                                                                                                                                                                                              The fire burnt along the railway line
## 1092                                                                                                                                                                                                                              The fire burnt along the railway line
## 1093                                                                                                                                                                                                                              The fire burnt along the railway line
## 1094                                                                                                                                                                                                                              The fire burnt along the railway line
## 1095                                                                                                                                                                                                                              The fire burnt along the railway line
## 1096                                                                                                                                                                                                                              The fire burnt along the railway line
## 1097                                                                                                                                                                                                                              The fire burnt along the railway line
## 1098                                                                                                                                                                                                                              The fire burnt along the railway line
## 1099                                                                                                                                                                                                                              The fire burnt along the railway line
## 1100                                                                                                                                                                                                                              The fire burnt along the railway line
## 1101                                                                                                                                                                                                                              The fire burnt along the railway line
## 1102                                                                                                                                                                                                                              The fire burnt along the railway line
## 1103                                                                                                                                                                                                                              The fire burnt along the railway line
## 1104                                                                                                                                                                                                                              The fire burnt along the railway line
## 1105                                                                                                                                                                                                                              The fire burnt along the railway line
## 1106                                                                                                                                                                                                                              The fire burnt along the railway line
## 1107                                                                                                                                                                                                                              The fire burnt along the railway line
## 1108                                                                                                                                                                                                                              The fire burnt along the railway line
## 1109                                                                                                                                                                                                                              The fire burnt along the railway line
## 1110                                                                                                                                                                                                                              The fire burnt along the railway line
## 1111                                                                                                                                                                                                                              The fire burnt along the railway line
## 1112                                                                                                                                                                                                                              The fire burnt along the railway line
## 1113                                                                                                                                                                                                                              The fire burnt along the railway line
## 1114                                                                                                                                                                                                                              The fire burnt along the railway line
## 1115                                                                                                                                                                                                                              The fire burnt along the railway line
## 1116                                                                                                                                                                                                                              The fire burnt along the railway line
## 1117                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1118                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1119                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1120                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1121                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1122                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1123                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1124                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1125                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1126                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1127                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1128                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1129                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1130                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1131                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1132                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1133                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1134                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1135                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1136                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1137                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1138                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1139                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1140                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1141                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1142                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1143                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1144                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1145                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1146                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1147                                                                                                                                                                                    The Eskom powerline was pushed over by an elephant and then snapped by giraffe.
## 1148                                                                                                                                                                                                                                                               <NA>
## 1149                                                                                                                                                                                                                                                               <NA>
## 1150                                                                                                                                                                                                                                                               <NA>
## 1151                                                                                                                                                                                                                                                               <NA>
## 1152                                                                                                                                                                                                                                                               <NA>
## 1153                                                                                                                                                                                                                                                               <NA>
## 1154                                                                                                                                                                                                                                                               <NA>
## 1155                                                                                                                                                                                                                                                               <NA>
## 1156                                                                                                                                                                                                                                                               <NA>
## 1157                                                                                                                                                                                                                                                               <NA>
## 1158                                                                                                                                                                                                                                                               <NA>
## 1159                                                                                                                                                                                                                                                               <NA>
## 1160                                                                                                                                                                                                                                                               <NA>
## 1161                                                                                                                                                                                                                                                               <NA>
## 1162                                                                                                                                                                                                                                                               <NA>
## 1163                                                                                                                                                                                                                                                               <NA>
## 1164                                                                                                                                                                                                                                                               <NA>
## 1165                                                                                                                                                                                                                                                               <NA>
## 1166                                                                                                                                                                                                                                                               <NA>
## 1167                                                                                                                                                                                                                                                               <NA>
## 1168                                                                                                                                                                                                                                                               <NA>
## 1169                                                                                                                                                                                                                                                               <NA>
## 1170                                                                                                                                                                                                                                                               <NA>
## 1171                                                                                                                                                                                                                                                               <NA>
## 1172                                                                                                                                                                                                                                                               <NA>
## 1173                                                                                                                                                                                                                                                               <NA>
## 1174                                                                                                                                                                                                                                                               <NA>
## 1175                                                                                                                                                                                                                                                               <NA>
## 1176                                                                                                                                                                                                                                                               <NA>
## 1177                                                                                                                                                                                                                  Please see original fire report for further notes
## 1178                                                                                                                                                                                                                  Please see original fire report for further notes
## 1179                                                                                                                                                                                                                  Please see original fire report for further notes
## 1180                                                                                                                                                                                                                  Please see original fire report for further notes
## 1181                                                                                                                                                                                                                  Please see original fire report for further notes
## 1182                                                                                                                                                                                                                  Please see original fire report for further notes
## 1183                                                                                                                                                                                                                  Please see original fire report for further notes
## 1184                                                                                                                                                                                                                  Please see original fire report for further notes
## 1185                                                                                                                                                                                                                  Please see original fire report for further notes
## 1186                                                                                                                                                                                                                  Please see original fire report for further notes
## 1187                                                                                                                                                                                                                  Please see original fire report for further notes
## 1188                                                                                                                                                                                                                  Please see original fire report for further notes
## 1189                                                                                                                                                                                                                  Please see original fire report for further notes
## 1190                                                                                                                                                                                                                  Please see original fire report for further notes
## 1191                                                                                                                                                                                                                  Please see original fire report for further notes
## 1192                                                                                                                                                                                                                  Please see original fire report for further notes
## 1193                                                                                                                                                                                                                  Please see original fire report for further notes
## 1194                                                                                                                                                                                                                  Please see original fire report for further notes
## 1195                                                                                                                                                                                                                  Please see original fire report for further notes
## 1196                                                                                                                                                                                                                  Please see original fire report for further notes
## 1197                                                                                                                                                                                                                  Please see original fire report for further notes
## 1198                                                                                                                                                                                                                  Please see original fire report for further notes
## 1199                                                                                                                                                                                                                  Please see original fire report for further notes
## 1200                                                                                                                                                                                                                  Please see original fire report for further notes
## 1201                                                                                                                                                                                                                  Please see original fire report for further notes
## 1202                                                                                                                                                                                                                  Please see original fire report for further notes
## 1203                                                                                                                                                                                                                  Please see original fire report for further notes
## 1204                                                                                                                                                                                                                  Please see original fire report for further notes
## 1205                                                                                                                                                                                                                  Please see original fire report for further notes
## 1206                                                                                                                                                                                                                                                               <NA>
## 1207                                                                                                                                                                                                                                                               <NA>
## 1208                                                                                                                                                                                                                                                               <NA>
## 1209                                                                                                                                                                                                                                                               <NA>
## 1210                                                                                                                                                                                                                                                               <NA>
## 1211                                                                                                                                                                                                                                                               <NA>
## 1212                                                                                                                                                                                                                                                               <NA>
## 1213                                                                                                                                                                                                                                                               <NA>
## 1214                                                                                                                                                                                                                                                               <NA>
## 1215                                                                                                                                                                                                                                                               <NA>
## 1216                                                                                                                                                                                                                                                               <NA>
## 1217                                                                                                                                                                                                                                                               <NA>
## 1218                                                                                                                                                                                                                                                               <NA>
## 1219                                                                                                                                                                                                                                                               <NA>
## 1220                                                                                                                                                                                                                                                               <NA>
## 1221                                                                                                                                                                                                                                                               <NA>
## 1222                                                                                                                                                                                                                                                               <NA>
## 1223                                                                                                                                                                                                                                                               <NA>
## 1224                                                                                                                                                                                                                                                               <NA>
## 1225                                                                                                                                                                                                                                                               <NA>
## 1226                                                                                                                                                                                                                                                               <NA>
## 1227                                                                                                                                                                                                                                                               <NA>
## 1228                                                                                                                                                                                                                                                               <NA>
## 1229                                                                                                                                                                                                                                                               <NA>
## 1230                                                                                                                                                                                                                                                               <NA>
## 1231                                                                                                                                                                                                                                                               <NA>
## 1232                                                                                                                                                                                                                                                               <NA>
## 1233                                                                                                                                                                                                                                                               <NA>
## 1234                                                                                                                                                                                                                                                               <NA>
## 1235                                                                                                                                                                                                                  Please see original fire report for further notes
## 1236                                                                                                                                                                                                                  Please see original fire report for further notes
## 1237                                                                                                                                                                                                                  Please see original fire report for further notes
## 1238                                                                                                                                                                                                                  Please see original fire report for further notes
## 1239                                                                                                                                                                                                                  Please see original fire report for further notes
## 1240                                                                                                                                                                                                                  Please see original fire report for further notes
## 1241                                                                                                                                                                                                                  Please see original fire report for further notes
## 1242                                                                                                                                                                                                                  Please see original fire report for further notes
## 1243                                                                                                                                                                                                                  Please see original fire report for further notes
## 1244                                                                                                                                                                                                                  Please see original fire report for further notes
## 1245                                                                                                                                                                                                                  Please see original fire report for further notes
## 1246                                                                                                                                                                                                                  Please see original fire report for further notes
## 1247                                                                                                                                                                                                                  Please see original fire report for further notes
## 1248                                                                                                                                                                                                                  Please see original fire report for further notes
## 1249                                                                                                                                                                                                                  Please see original fire report for further notes
## 1250                                                                                                                                                                                                                  Please see original fire report for further notes
## 1251                                                                                                                                                                                                                  Please see original fire report for further notes
## 1252                                                                                                                                                                                                                  Please see original fire report for further notes
## 1253                                                                                                                                                                                                                  Please see original fire report for further notes
## 1254                                                                                                                                                                                                                  Please see original fire report for further notes
## 1255                                                                                                                                                                                                                  Please see original fire report for further notes
## 1256                                                                                                                                                                                                                  Please see original fire report for further notes
## 1257                                                                                                                                                                                                                  Please see original fire report for further notes
## 1258                                                                                                                                                                                                                  Please see original fire report for further notes
## 1259                                                                                                                                                                                                                  Please see original fire report for further notes
## 1260                                                                                                                                                                                                                  Please see original fire report for further notes
## 1261                                                                                                                                                                                                                  Please see original fire report for further notes
## 1262                                                                                                                                                                                                                  Please see original fire report for further notes
## 1263                                                                                                                                                                                                                  Please see original fire report for further notes
## 1264                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1265                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1266                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1267                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1268                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1269                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1270                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1271                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1272                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1273                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1274                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1275                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1276                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1277                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1278                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1279                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1280                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1281                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1282                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1283                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1284                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1285                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1286                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1287                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1288                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1289                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1290                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1291                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1292                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1293                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1294                                                                                                                                               Accidental fire from UCT project into the buffalo excloure.\r\n\r\nPlease see original fire report for further notes
## 1295                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1296                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1297                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1298                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1299                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1300                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1301                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1302                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1303                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1304                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1305                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1306                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1307                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1308                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1309                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1310                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1311                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1312                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1313                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1314                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1315                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1316                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1317                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1318                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1319                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1320                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1321                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1322                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1323                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1324                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1325                                                                                                                                                                                                     All woody and herbaceous vegetation was burnt to its fullests.
## 1326                                                                                                                                                                                                                                                               <NA>
## 1327                                                                                                                                                                                                                                                               <NA>
## 1328                                                                                                                                                                                                                                                               <NA>
## 1329                                                                                                                                                                                                                                                               <NA>
## 1330                                                                                                                                                                                                                                                               <NA>
## 1331                                                                                                                                                                                                                                                               <NA>
## 1332                                                                                                                                                                                                                                                               <NA>
## 1333                                                                                                                                                                                                                                                               <NA>
## 1334                                                                                                                                                                                                                                                               <NA>
## 1335                                                                                                                                                                                                                                                               <NA>
## 1336                                                                                                                                                                                                                                                               <NA>
## 1337                                                                                                                                                                                                                                                               <NA>
## 1338                                                                                                                                                                                                                                                               <NA>
## 1339                                                                                                                                                                                                                                                               <NA>
## 1340                                                                                                                                                                                                                                                               <NA>
## 1341                                                                                                                                                                                                                                                               <NA>
## 1342                                                                                                                                                                                                                                                               <NA>
## 1343                                                                                                                                                                                                                                                               <NA>
## 1344                                                                                                                                                                                                                                                               <NA>
## 1345                                                                                                                                                                                                                                                               <NA>
## 1346                                                                                                                                                                                                                                                               <NA>
## 1347                                                                                                                                                                                                                                                               <NA>
## 1348                                                                                                                                                                                                                                                               <NA>
## 1349                                                                                                                                                                                                                                                               <NA>
## 1350                                                                                                                                                                                                                                                               <NA>
## 1351                                                                                                                                                                                                                                                               <NA>
## 1352                                                                                                                                                                                                                                                               <NA>
## 1353                                                                                                                                                                                                                                                               <NA>
## 1354                                                                                                                                                                                                                                                               <NA>
## 1355                                                                                                                                                                                                                                                               <NA>
## 1356                                                                                                                                                                                                                                                               <NA>
## 1357                                                                                                                                                                                                                                                               <NA>
## 1358                                                                                                                                                                                                                                                               <NA>
## 1359                                                                                                                                                                                                                                                               <NA>
## 1360                                                                                                                                                                                                                                                               <NA>
## 1361                                                                                                                                                                                                                                                               <NA>
## 1362                                                                                                                                                                                                                                                               <NA>
## 1363                                                                                                                                                                                                                                                               <NA>
## 1364                                                                                                                                                                                                                                                               <NA>
## 1365                                                                                                                                                                                                                                                               <NA>
## 1366                                                                                                                                                                                                                                                               <NA>
## 1367                                                                                                                                                                                                                                                               <NA>
## 1368                                                                                                                                                                                                                                                               <NA>
## 1369                                                                                                                                                                                                                                                               <NA>
## 1370                                                                                                                                                                                                                                                               <NA>
## 1371                                                                                                                                                                                                                                                               <NA>
## 1372                                                                                                                                                                                                                                                               <NA>
## 1373                                                                                                                                                                                                                                                               <NA>
## 1374                                                                                                                                                                                                                                                               <NA>
## 1375                                                                                                                                                                                                                                                               <NA>
## 1376                                                                                                                                                                                                                                                               <NA>
## 1377                                                                                                                                                                                                                                                               <NA>
## 1378                                                                                                                                                                                                                                                               <NA>
## 1379                                                                                                                                                                                                                                                               <NA>
## 1380                                                                                                                                                                                                                                                               <NA>
## 1381                                                                                                                                                                                                                                                               <NA>
## 1382                                                                                                                                                                                                                                                               <NA>
## 1383                                                                                                                                                                                                                                                               <NA>
## 1384                                                                                                                                                                                                                                                               <NA>
## 1385                                                                                                                                                                                                                                                               <NA>
## 1386                                                                                                                                                                                                                                                               <NA>
## 1387                                                                                                                                                                                                                                                               <NA>
## 1388                                                                                                                                                                                                                                                               <NA>
## 1389                                                                                                                                                                                                                                                               <NA>
## 1390                                                                                                                                                                                                                                                               <NA>
## 1391                                                                                                                                                                                                                                                               <NA>
## 1392                                                                                                                                                                                                                                                               <NA>
## 1393                                                                                                                                                                                                                                                               <NA>
## 1394                                                                                                                                                                                                                                                               <NA>
## 1395                                                                                                                                                                                                                                                               <NA>
## 1396                                                                                                                                                                                                                                                               <NA>
## 1397                                                                                                                                                                                                                                                               <NA>
## 1398                                                                                                                                                                                                                                                               <NA>
## 1399                                                                                                                                                                                                                                                               <NA>
## 1400                                                                                                                                                                                                                                                               <NA>
## 1401                                                                                                                                                                                                                                                               <NA>
## 1402                                                                                                                                                                                                                                                               <NA>
## 1403                                                                                                                                                                                                                                                               <NA>
## 1404                                                                                                                                                                                                                                                               <NA>
## 1405                                                                                                                                                                                                                                                               <NA>
## 1406                                                                                                                                                                                                                                                               <NA>
## 1407                                                                                                                                                                                                                                                               <NA>
## 1408                                                                                                                                                                                                                                                               <NA>
## 1409                                                                                                                                                                                                                                                               <NA>
## 1410                                                                                                                                                                                                                                                               <NA>
## 1411                                                                                                                                                                                                                                                               <NA>
## 1412                                                                                                                                                                                                                                                               <NA>
## 1413                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1414                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1415                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1416                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1417                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1418                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1419                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1420                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1421                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1422                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1423                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1424                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1425                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1426                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1427                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1428                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1429                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1430                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1431                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1432                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1433                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1434                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1435                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1436                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1437                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1438                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1439                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1440                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1441                                                                                                          The fire was driven by wind and it only touched a small corner of the block. Fire jumped from the back-burns. See original fire report for further notes.
## 1442                                                                                                                                                                                                                                                               <NA>
## 1443                                                                                                                                                                                                                                                               <NA>
## 1444                                                                                                                                                                                                                                                               <NA>
## 1445                                                                                                                                                                                                                                                               <NA>
## 1446                                                                                                                                                                                                                                                               <NA>
## 1447                                                                                                                                                                                                                                                               <NA>
## 1448                                                                                                                                                                                                                                                               <NA>
## 1449                                                                                                                                                                                                                                                               <NA>
## 1450                                                                                                                                                                                                                                                               <NA>
## 1451                                                                                                                                                                                                                                                               <NA>
## 1452                                                                                                                                                                                                                                                               <NA>
## 1453                                                                                                                                                                                                                                                               <NA>
## 1454                                                                                                                                                                                                                                                               <NA>
## 1455                                                                                                                                                                                                                                                               <NA>
## 1456                                                                                                                                                                                                                                                               <NA>
## 1457                                                                                                                                                                                                                                                               <NA>
## 1458                                                                                                                                                                                                                                                               <NA>
## 1459                                                                                                                                                                                                                                                               <NA>
## 1460                                                                                                                                                                                                                                                               <NA>
## 1461                                                                                                                                                                                                                                                               <NA>
## 1462                                                                                                                                                                                                                                                               <NA>
## 1463                                                                                                                                                                                                                                                               <NA>
## 1464                                                                                                                                                                                                                                                               <NA>
## 1465                                                                                                                                                                                                                                                               <NA>
## 1466                                                                                                                                                                                                                                                               <NA>
## 1467                                                                                                                                                                                                                                                               <NA>
## 1468                                                                                                                                                                                                                                                               <NA>
## 1469                                                                                                                                                                                                                                                               <NA>
## 1470                                                                                                                                                                                                                                                               <NA>
## 1471                                                                                                                                                                                                                                                               <NA>
## 1472                                                                                                                                                                                                                                                               <NA>
## 1473                                                                                                                                                                                                                                                               <NA>
## 1474                                                                                                                                                                                                                                                               <NA>
## 1475                                                                                                                                                                                                                                                               <NA>
## 1476                                                                                                                                                                                                                                                               <NA>
## 1477                                                                                                                                                                                                                                                               <NA>
## 1478                                                                                                                                                                                                                                                               <NA>
## 1479                                                                                                                                                                                                                                                               <NA>
## 1480                                                                                                                                                                                                                                                               <NA>
## 1481                                                                                                                                                                                                                                                               <NA>
## 1482                                                                                                                                                                                                                                                               <NA>
## 1483                                                                                                                                                                                                                                                               <NA>
## 1484                                                                                                                                                                                                                                                               <NA>
## 1485                                                                                                                                                                                                                                                               <NA>
## 1486                                                                                                                                                                                                                                                               <NA>
## 1487                                                                                                                                                                                                                                                               <NA>
## 1488                                                                                                                                                                                                                                                               <NA>
## 1489                                                                                                                                                                                                                                                               <NA>
## 1490                                                                                                                                                                                                                                                               <NA>
## 1491                                                                                                                                                                                                                                                               <NA>
## 1492                                                                                                                                                                                                                                                               <NA>
## 1493                                                                                                                                                                                                                                                               <NA>
## 1494                                                                                                                                                                                                                                                               <NA>
## 1495                                                                                                                                                                                                                                                               <NA>
## 1496                                                                                                                                                                                                                                                               <NA>
## 1497                                                                                                                                                                                                                                                               <NA>
## 1498                                                                                                                                                                                                                                                               <NA>
## 1499                                                                                                                                                                                                                                                               <NA>
## 1500                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1501                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1502                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1503                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1504                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1505                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1506                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1507                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1508                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1509                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1510                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1511                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1512                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1513                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1514                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1515                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1516                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1517                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1518                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1519                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1520                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1521                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1522                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1523                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1524                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1525                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1526                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1527                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1528                                                                                                                                                             Although the fire perimeter covered 60% of the block, the fire area was approximately 20% of the block
## 1529                                                                                                                                                                                                                  Please see original fire report for further notes
## 1530                                                                                                                                                                                                                  Please see original fire report for further notes
## 1531                                                                                                                                                                                                                  Please see original fire report for further notes
## 1532                                                                                                                                                                                                                  Please see original fire report for further notes
## 1533                                                                                                                                                                                                                  Please see original fire report for further notes
## 1534                                                                                                                                                                                                                  Please see original fire report for further notes
## 1535                                                                                                                                                                                                                  Please see original fire report for further notes
## 1536                                                                                                                                                                                                                  Please see original fire report for further notes
## 1537                                                                                                                                                                                                                  Please see original fire report for further notes
## 1538                                                                                                                                                                                                                  Please see original fire report for further notes
## 1539                                                                                                                                                                                                                  Please see original fire report for further notes
## 1540                                                                                                                                                                                                                  Please see original fire report for further notes
## 1541                                                                                                                                                                                                                  Please see original fire report for further notes
## 1542                                                                                                                                                                                                                  Please see original fire report for further notes
## 1543                                                                                                                                                                                                                  Please see original fire report for further notes
## 1544                                                                                                                                                                                                                  Please see original fire report for further notes
## 1545                                                                                                                                                                                                                  Please see original fire report for further notes
## 1546                                                                                                                                                                                                                  Please see original fire report for further notes
## 1547                                                                                                                                                                                                                  Please see original fire report for further notes
## 1548                                                                                                                                                                                                                  Please see original fire report for further notes
## 1549                                                                                                                                                                                                                  Please see original fire report for further notes
## 1550                                                                                                                                                                                                                  Please see original fire report for further notes
## 1551                                                                                                                                                                                                                  Please see original fire report for further notes
## 1552                                                                                                                                                                                                                  Please see original fire report for further notes
## 1553                                                                                                                                                                                                                  Please see original fire report for further notes
## 1554                                                                                                                                                                                                                  Please see original fire report for further notes
## 1555                                                                                                                                                                                                                  Please see original fire report for further notes
## 1556                                                                                                                                                                                                                  Please see original fire report for further notes
## 1557                                                                                                                                                                                                                  Please see original fire report for further notes
## 1558                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1559                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1560                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1561                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1562                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1563                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1564                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1565                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1566                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1567                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1568                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1569                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1570                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1571                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1572                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1573                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1574                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1575                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1576                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1577                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1578                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1579                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1580                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1581                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1582                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1583                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1584                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1585                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1586                                                                                                                                            Burnt back from point A-Band C-D, to contain the fire to block C096A and to avoid it spreading to blocks C095A and C096
## 1587                                                                                                                                                                                                                                                               <NA>
## 1588                                                                                                                                                                                                                                                               <NA>
## 1589                                                                                                                                                                                                                                                               <NA>
## 1590                                                                                                                                                                                                                                                               <NA>
## 1591                                                                                                                                                                                                                                                               <NA>
## 1592                                                                                                                                                                                                                                                               <NA>
## 1593                                                                                                                                                                                                                                                               <NA>
## 1594                                                                                                                                                                                                                                                               <NA>
## 1595                                                                                                                                                                                                                                                               <NA>
## 1596                                                                                                                                                                                                                                                               <NA>
## 1597                                                                                                                                                                                                                                                               <NA>
## 1598                                                                                                                                                                                                                                                               <NA>
## 1599                                                                                                                                                                                                                                                               <NA>
## 1600                                                                                                                                                                                                                                                               <NA>
## 1601                                                                                                                                                                                                                                                               <NA>
## 1602                                                                                                                                                                                                                                                               <NA>
## 1603                                                                                                                                                                                                                                                               <NA>
## 1604                                                                                                                                                                                                                                                               <NA>
## 1605                                                                                                                                                                                                                                                               <NA>
## 1606                                                                                                                                                                                                                                                               <NA>
## 1607                                                                                                                                                                                                                                                               <NA>
## 1608                                                                                                                                                                                                                                                               <NA>
## 1609                                                                                                                                                                                                                                                               <NA>
## 1610                                                                                                                                                                                                                                                               <NA>
## 1611                                                                                                                                                                                                                                                               <NA>
## 1612                                                                                                                                                                                                                                                               <NA>
## 1613                                                                                                                                                                                                                                                               <NA>
## 1614                                                                                                                                                                                                                                                               <NA>
## 1615                                                                                                                                                                                                                                                               <NA>
## 1616                                                                                                                                                                                                                                                               <NA>
## 1617                                                                                                                                                                                                                                                               <NA>
## 1618                                                                                                                                                                                                                                                               <NA>
## 1619                                                                                                                                                                                                                                                               <NA>
## 1620                                                                                                                                                                                                                                                               <NA>
## 1621                                                                                                                                                                                                                                                               <NA>
## 1622                                                                                                                                                                                                                                                               <NA>
## 1623                                                                                                                                                                                                                                                               <NA>
## 1624                                                                                                                                                                                                                                                               <NA>
## 1625                                                                                                                                                                                                                                                               <NA>
## 1626                                                                                                                                                                                                                                                               <NA>
## 1627                                                                                                                                                                                                                                                               <NA>
## 1628                                                                                                                                                                                                                                                               <NA>
## 1629                                                                                                                                                                                                                                                               <NA>
## 1630                                                                                                                                                                                                                                                               <NA>
## 1631                                                                                                                                                                                                                                                               <NA>
## 1632                                                                                                                                                                                                                                                               <NA>
## 1633                                                                                                                                                                                                                                                               <NA>
## 1634                                                                                                                                                                                                                                                               <NA>
## 1635                                                                                                                                                                                                                                                               <NA>
## 1636                                                                                                                                                                                                                                                               <NA>
## 1637                                                                                                                                                                                                                                                               <NA>
## 1638                                                                                                                                                                                                                                                               <NA>
## 1639                                                                                                                                                                                                                                                               <NA>
## 1640                                                                                                                                                                                                                                                               <NA>
## 1641                                                                                                                                                                                                                                                               <NA>
## 1642                                                                                                                                                                                                                                                               <NA>
## 1643                                                                                                                                                                                                                                                               <NA>
## 1644                                                                                                                                                                                                                                                               <NA>
## 1645                                                                                                                                                                                                                                                          Roan Camp
## 1646                                                                                                                                                                                                                                                          Roan Camp
## 1647                                                                                                                                                                                                                                                          Roan Camp
## 1648                                                                                                                                                                                                                                                          Roan Camp
## 1649                                                                                                                                                                                                                                                          Roan Camp
## 1650                                                                                                                                                                                                                                                          Roan Camp
## 1651                                                                                                                                                                                                                                                          Roan Camp
## 1652                                                                                                                                                                                                                                                          Roan Camp
## 1653                                                                                                                                                                                                                                                          Roan Camp
## 1654                                                                                                                                                                                                                                                          Roan Camp
## 1655                                                                                                                                                                                                                                                          Roan Camp
## 1656                                                                                                                                                                                                                                                          Roan Camp
## 1657                                                                                                                                                                                                                                                          Roan Camp
## 1658                                                                                                                                                                                                                                                          Roan Camp
## 1659                                                                                                                                                                                                                                                          Roan Camp
## 1660                                                                                                                                                                                                                                                          Roan Camp
## 1661                                                                                                                                                                                                                                                          Roan Camp
## 1662                                                                                                                                                                                                                                                          Roan Camp
## 1663                                                                                                                                                                                                                                                          Roan Camp
## 1664                                                                                                                                                                                                                                                          Roan Camp
## 1665                                                                                                                                                                                                                                                          Roan Camp
## 1666                                                                                                                                                                                                                                                          Roan Camp
## 1667                                                                                                                                                                                                                                                          Roan Camp
## 1668                                                                                                                                                                                                                                                          Roan Camp
## 1669                                                                                                                                                                                                                                                          Roan Camp
## 1670                                                                                                                                                                                                                                                          Roan Camp
## 1671                                                                                                                                                                                                                                                          Roan Camp
## 1672                                                                                                                                                                                                                                                          Roan Camp
## 1673                                                                                                                                                                                                                                                          Roan Camp
## 1674                                                                                                                                                                                                                                                               <NA>
## 1675                                                                                                                                                                                                                                                               <NA>
## 1676                                                                                                                                                                                                                                                               <NA>
## 1677                                                                                                                                                                                                                                                               <NA>
## 1678                                                                                                                                                                                                                                                               <NA>
## 1679                                                                                                                                                                                                                                                               <NA>
## 1680                                                                                                                                                                                                                                                               <NA>
## 1681                                                                                                                                                                                                                                                               <NA>
## 1682                                                                                                                                                                                                                                                               <NA>
## 1683                                                                                                                                                                                                                                                               <NA>
## 1684                                                                                                                                                                                                                                                               <NA>
## 1685                                                                                                                                                                                                                                                               <NA>
## 1686                                                                                                                                                                                                                                                               <NA>
## 1687                                                                                                                                                                                                                                                               <NA>
## 1688                                                                                                                                                                                                                                                               <NA>
## 1689                                                                                                                                                                                                                                                               <NA>
## 1690                                                                                                                                                                                                                                                               <NA>
## 1691                                                                                                                                                                                                                                                               <NA>
## 1692                                                                                                                                                                                                                                                               <NA>
## 1693                                                                                                                                                                                                                                                               <NA>
## 1694                                                                                                                                                                                                                                                               <NA>
## 1695                                                                                                                                                                                                                                                               <NA>
## 1696                                                                                                                                                                                                                                                               <NA>
## 1697                                                                                                                                                                                                                                                               <NA>
## 1698                                                                                                                                                                                                                                                               <NA>
## 1699                                                                                                                                                                                                                                                               <NA>
## 1700                                                                                                                                                                                                                                                               <NA>
## 1701                                                                                                                                                                                                                                                               <NA>
## 1702                                                                                                                                                                                                                                                               <NA>
## 1703                                                                                                                                                                                                                                                      Very hot burn
## 1704                                                                                                                                                                                                                                                      Very hot burn
## 1705                                                                                                                                                                                                                                                      Very hot burn
## 1706                                                                                                                                                                                                                                                      Very hot burn
## 1707                                                                                                                                                                                                                                                      Very hot burn
## 1708                                                                                                                                                                                                                                                      Very hot burn
## 1709                                                                                                                                                                                                                                                      Very hot burn
## 1710                                                                                                                                                                                                                                                      Very hot burn
## 1711                                                                                                                                                                                                                                                      Very hot burn
## 1712                                                                                                                                                                                                                                                      Very hot burn
## 1713                                                                                                                                                                                                                                                      Very hot burn
## 1714                                                                                                                                                                                                                                                      Very hot burn
## 1715                                                                                                                                                                                                                                                      Very hot burn
## 1716                                                                                                                                                                                                                                                      Very hot burn
## 1717                                                                                                                                                                                                                                                      Very hot burn
## 1718                                                                                                                                                                                                                                                      Very hot burn
## 1719                                                                                                                                                                                                                                                      Very hot burn
## 1720                                                                                                                                                                                                                                                      Very hot burn
## 1721                                                                                                                                                                                                                                                      Very hot burn
## 1722                                                                                                                                                                                                                                                      Very hot burn
## 1723                                                                                                                                                                                                                                                      Very hot burn
## 1724                                                                                                                                                                                                                                                      Very hot burn
## 1725                                                                                                                                                                                                                                                      Very hot burn
## 1726                                                                                                                                                                                                                                                      Very hot burn
## 1727                                                                                                                                                                                                                                                      Very hot burn
## 1728                                                                                                                                                                                                                                                      Very hot burn
## 1729                                                                                                                                                                                                                                                      Very hot burn
## 1730                                                                                                                                                                                                                                                      Very hot burn
## 1731                                                                                                                                                                                                                                                      Very hot burn
## 1732                                                                                                                                                                                                                                                               <NA>
## 1733                                                                                                                                                                                                                                                               <NA>
## 1734                                                                                                                                                                                                                                                               <NA>
## 1735                                                                                                                                                                                                                                                               <NA>
## 1736                                                                                                                                                                                                                                                               <NA>
## 1737                                                                                                                                                                                                                                                               <NA>
## 1738                                                                                                                                                                                                                                                               <NA>
## 1739                                                                                                                                                                                                                                                               <NA>
## 1740                                                                                                                                                                                                                                                               <NA>
## 1741                                                                                                                                                                                                                                                               <NA>
## 1742                                                                                                                                                                                                                                                               <NA>
## 1743                                                                                                                                                                                                                                                               <NA>
## 1744                                                                                                                                                                                                                                                               <NA>
## 1745                                                                                                                                                                                                                                                               <NA>
## 1746                                                                                                                                                                                                                                                               <NA>
## 1747                                                                                                                                                                                                                                                               <NA>
## 1748                                                                                                                                                                                                                                                               <NA>
## 1749                                                                                                                                                                                                                                                               <NA>
## 1750                                                                                                                                                                                                                                                               <NA>
## 1751                                                                                                                                                                                                                                                               <NA>
## 1752                                                                                                                                                                                                                                                               <NA>
## 1753                                                                                                                                                                                                                                                               <NA>
## 1754                                                                                                                                                                                                                                                               <NA>
## 1755                                                                                                                                                                                                                                                               <NA>
## 1756                                                                                                                                                                                                                                                               <NA>
## 1757                                                                                                                                                                                                                                                               <NA>
## 1758                                                                                                                                                                                                                                                               <NA>
## 1759                                                                                                                                                                                                                                                               <NA>
## 1760                                                                                                                                                                                                                                                               <NA>
## 1761                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1762                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1763                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1764                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1765                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1766                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1767                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1768                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1769                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1770                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1771                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1772                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1773                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1774                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1775                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1776                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1777                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1778                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1779                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1780                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1781                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1782                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1783                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1784                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1785                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1786                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1787                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1788                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1789                                                                                                                                                                                                        Very high litter load\r\nBlock did not burn for a few years
## 1790                                                                                                                                                                                                                                                               <NA>
## 1791                                                                                                                                                                                                                                                               <NA>
## 1792                                                                                                                                                                                                                                                               <NA>
## 1793                                                                                                                                                                                                                                                               <NA>
## 1794                                                                                                                                                                                                                                                               <NA>
## 1795                                                                                                                                                                                                                                                               <NA>
## 1796                                                                                                                                                                                                                                                               <NA>
## 1797                                                                                                                                                                                                                                                               <NA>
## 1798                                                                                                                                                                                                                                                               <NA>
## 1799                                                                                                                                                                                                                                                               <NA>
## 1800                                                                                                                                                                                                                                                               <NA>
## 1801                                                                                                                                                                                                                                                               <NA>
## 1802                                                                                                                                                                                                                                                               <NA>
## 1803                                                                                                                                                                                                                                                               <NA>
## 1804                                                                                                                                                                                                                                                               <NA>
## 1805                                                                                                                                                                                                                                                               <NA>
## 1806                                                                                                                                                                                                                                                               <NA>
## 1807                                                                                                                                                                                                                                                               <NA>
## 1808                                                                                                                                                                                                                                                               <NA>
## 1809                                                                                                                                                                                                                                                               <NA>
## 1810                                                                                                                                                                                                                                                               <NA>
## 1811                                                                                                                                                                                                                                                               <NA>
## 1812                                                                                                                                                                                                                                                               <NA>
## 1813                                                                                                                                                                                                                                                               <NA>
## 1814                                                                                                                                                                                                                                                               <NA>
## 1815                                                                                                                                                                                                                                                               <NA>
## 1816                                                                                                                                                                                                                                                               <NA>
## 1817                                                                                                                                                                                                                                                               <NA>
## 1818                                                                                                                                                                                                                                                               <NA>
## 1819                                                                                                                                                                                                                                                               <NA>
## 1820                                                                                                                                                                                                                                                               <NA>
## 1821                                                                                                                                                                                                                                                               <NA>
## 1822                                                                                                                                                                                                                                                               <NA>
## 1823                                                                                                                                                                                                                                                               <NA>
## 1824                                                                                                                                                                                                                                                               <NA>
## 1825                                                                                                                                                                                                                                                               <NA>
## 1826                                                                                                                                                                                                                                                               <NA>
## 1827                                                                                                                                                                                                                                                               <NA>
## 1828                                                                                                                                                                                                                                                               <NA>
## 1829                                                                                                                                                                                                                                                               <NA>
## 1830                                                                                                                                                                                                                                                               <NA>
## 1831                                                                                                                                                                                                                                                               <NA>
## 1832                                                                                                                                                                                                                                                               <NA>
## 1833                                                                                                                                                                                                                                                               <NA>
## 1834                                                                                                                                                                                                                                                               <NA>
## 1835                                                                                                                                                                                                                                                               <NA>
## 1836                                                                                                                                                                                                                                                               <NA>
## 1837                                                                                                                                                                                                                                                               <NA>
## 1838                                                                                                                                                                                                                                                               <NA>
## 1839                                                                                                                                                                                                                                                               <NA>
## 1840                                                                                                                                                                                                                                                               <NA>
## 1841                                                                                                                                                                                                                                                               <NA>
## 1842                                                                                                                                                                                                                                                               <NA>
## 1843                                                                                                                                                                                                                                                               <NA>
## 1844                                                                                                                                                                                                                                                               <NA>
## 1845                                                                                                                                                                                                                                                               <NA>
## 1846                                                                                                                                                                                                                                                               <NA>
## 1847                                                                                                                                                                                                                                                               <NA>
##       AREA_HA DATE_START DATE_END FIRE_YEAR WoodyCover   MAP Station Year
## 1      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1987
## 2      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1982
## 3      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1990
## 4      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1991
## 5      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1988
## 6      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1989
## 7      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1981
## 8      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1984
## 9      184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2001
## 10     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1983
## 11     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1998
## 12     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1999
## 13     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1985
## 14     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1986
## 15     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1992
## 16     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1980
## 17     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2000
## 18     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2008
## 19     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2005
## 20     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2006
## 21     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2007
## 22     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1994
## 23     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1997
## 24     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1993
## 25     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2002
## 26     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1995
## 27     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 1996
## 28     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2003
## 29     184.85       <NA>     <NA>      <NA>     38.241 745.7     KRO 2004
## 30   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1987
## 31   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1982
## 32   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1990
## 33   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1991
## 34   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1988
## 35   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1989
## 36   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1981
## 37   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1984
## 38   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2001
## 39   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1983
## 40   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1998
## 41   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1999
## 42   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1985
## 43   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1986
## 44   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1992
## 45   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1980
## 46   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2000
## 47   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2008
## 48   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2005
## 49   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2006
## 50   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2007
## 51   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1994
## 52   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1997
## 53   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1993
## 54   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2002
## 55   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1995
## 56   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 1996
## 57   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2003
## 58   17718.03       <NA>     <NA>      <NA>     31.920 648.1     KRO 2004
## 59    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1987
## 60    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1982
## 61    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1990
## 62    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1991
## 63    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1988
## 64    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1989
## 65    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1981
## 66    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1984
## 67    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2001
## 68    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1983
## 69    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1998
## 70    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1999
## 71    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1985
## 72    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1986
## 73    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1992
## 74    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1980
## 75    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2000
## 76    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2008
## 77    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2005
## 78    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2006
## 79    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2007
## 80    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1994
## 81    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1997
## 82    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1993
## 83    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2002
## 84    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1995
## 85    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 1996
## 86    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2003
## 87    2247.48       <NA>     <NA>      <NA>     35.486 634.9     KRO 2004
## 88    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1987
## 89    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1982
## 90    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1990
## 91    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1991
## 92    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1988
## 93    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1989
## 94    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1981
## 95    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1984
## 96    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2001
## 97    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1983
## 98    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1998
## 99    2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1999
## 100   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1985
## 101   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1986
## 102   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1992
## 103   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1980
## 104   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2000
## 105   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2008
## 106   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2005
## 107   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2006
## 108   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2007
## 109   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1994
## 110   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1997
## 111   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1993
## 112   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2002
## 113   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1995
## 114   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 1996
## 115   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2003
## 116   2149.70       <NA>     <NA>      <NA>     39.474 686.3     KRO 2004
## 117   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1987
## 118   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1982
## 119   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1990
## 120   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1991
## 121   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1988
## 122   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1989
## 123   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1981
## 124   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1984
## 125   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2001
## 126   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1983
## 127   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1998
## 128   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1999
## 129   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1985
## 130   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1986
## 131   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1992
## 132   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1980
## 133   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2000
## 134   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2008
## 135   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2005
## 136   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2006
## 137   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2007
## 138   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1994
## 139   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1997
## 140   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1993
## 141   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2002
## 142   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1995
## 143   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 1996
## 144   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2003
## 145   1188.02       <NA>     <NA>      <NA>     37.143 631.5     KRO 2004
## 146    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2001
## 147    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2000
## 148    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2004
## 149    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2002
## 150    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1986
## 151    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2007
## 152    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1980
## 153    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1999
## 154    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1998
## 155    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1985
## 156    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2003
## 157    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1982
## 158    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1984
## 159    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2006
## 160    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1993
## 161    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1994
## 162    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1996
## 163    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1997
## 164    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1992
## 165    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1983
## 166    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1991
## 167    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2005
## 168    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 2008
## 169    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1987
## 170    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1981
## 171    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1995
## 172    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1988
## 173    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1989
## 174    273.96       <NA>     <NA>      <NA>     28.579 491.8     LET 1990
## 175   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2001
## 176   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2000
## 177   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2004
## 178   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2002
## 179   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1986
## 180   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2007
## 181   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1980
## 182   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1999
## 183   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1998
## 184   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1985
## 185   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2003
## 186   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1982
## 187   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1984
## 188   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2006
## 189   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1993
## 190   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1994
## 191   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1996
## 192   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1997
## 193   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1992
## 194   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1983
## 195   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1991
## 196   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2005
## 197   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 2008
## 198   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1987
## 199   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1981
## 200   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1995
## 201   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1988
## 202   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1989
## 203   7882.86       <NA>     <NA>      <NA>     13.311 489.6     LET 1990
## 204   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1990
## 205   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2008
## 206   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2002
## 207   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2007
## 208   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1999
## 209   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2000
## 210   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1988
## 211   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1998
## 212   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1982
## 213   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1991
## 214   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2005
## 215   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1985
## 216   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1986
## 217   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1987
## 218   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2003
## 219   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1993
## 220   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1997
## 221   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2001
## 222   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2006
## 223   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1981
## 224   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1984
## 225   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1989
## 226   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1996
## 227   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1992
## 228   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 2004
## 229   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1995
## 230   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1994
## 231   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1983
## 232   1247.85       <NA>     <NA>      <NA>     21.006 627.5     OSA 1980
## 233   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1990
## 234   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2008
## 235   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2002
## 236   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2007
## 237   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1999
## 238   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2000
## 239   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1988
## 240   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1998
## 241   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1982
## 242   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1991
## 243   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2005
## 244   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1985
## 245   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1986
## 246   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1987
## 247   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2003
## 248   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1993
## 249   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1997
## 250   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2001
## 251   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2006
## 252   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1981
## 253   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1984
## 254   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1989
## 255   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1996
## 256   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1992
## 257   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 2004
## 258   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1995
## 259   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1994
## 260   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1983
## 261   1427.90       <NA>     <NA>      <NA>     37.553 651.6     OSA 1980
## 262  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1990
## 263  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2008
## 264  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2002
## 265  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2007
## 266  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1999
## 267  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2000
## 268  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1988
## 269  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1998
## 270  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1982
## 271  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1991
## 272  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2005
## 273  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1985
## 274  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1986
## 275  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1987
## 276  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2003
## 277  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1993
## 278  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1997
## 279  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2001
## 280  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2006
## 281  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1981
## 282  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1984
## 283  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1989
## 284  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1996
## 285  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1992
## 286  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 2004
## 287  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1995
## 288  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1994
## 289  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1983
## 290  48722.35       <NA>     <NA>      <NA>     37.131 651.2     OSA 1980
## 291  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1990
## 292  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2008
## 293  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2002
## 294  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2007
## 295  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1999
## 296  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2000
## 297  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1988
## 298  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1998
## 299  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1982
## 300  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1991
## 301  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2005
## 302  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1985
## 303  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1986
## 304  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1987
## 305  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2003
## 306  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1993
## 307  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1997
## 308  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2001
## 309  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2006
## 310  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1981
## 311  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1984
## 312  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1989
## 313  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1996
## 314  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1992
## 315  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 2004
## 316  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1995
## 317  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1994
## 318  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1983
## 319  23473.92       <NA>     <NA>      <NA>     26.039 638.4     OSA 1980
## 320    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1997
## 321    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1993
## 322    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1994
## 323    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1992
## 324    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2008
## 325    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1982
## 326    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1998
## 327    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1991
## 328    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1988
## 329    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2005
## 330    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2004
## 331    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1983
## 332    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1999
## 333    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2006
## 334    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1995
## 335    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2003
## 336    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1989
## 337    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1987
## 338    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2007
## 339    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2002
## 340    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1996
## 341    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1986
## 342    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1980
## 343    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2000
## 344    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 2001
## 345    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1984
## 346    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1985
## 347    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1990
## 348    101.79       <NA>     <NA>      <NA>     34.688 480.1     MAH 1981
## 349    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1997
## 350    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1993
## 351    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1994
## 352    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1992
## 353    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2008
## 354    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1982
## 355    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1998
## 356    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1991
## 357    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1988
## 358    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2005
## 359    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2004
## 360    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1983
## 361    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1999
## 362    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2006
## 363    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1995
## 364    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2003
## 365    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1989
## 366    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1987
## 367    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2007
## 368    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2002
## 369    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1996
## 370    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1986
## 371    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1980
## 372    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2000
## 373    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 2001
## 374    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1984
## 375    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1985
## 376    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1990
## 377    149.52       <NA>     <NA>      <NA>     49.872 503.1     MAH 1981
## 378  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1990
## 379  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1987
## 380  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1985
## 381  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1981
## 382  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1982
## 383  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1993
## 384  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2000
## 385  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2001
## 386  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1988
## 387  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1999
## 388  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2002
## 389  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1994
## 390  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1998
## 391  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1992
## 392  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1984
## 393  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1983
## 394  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1997
## 395  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2007
## 396  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2008
## 397  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1989
## 398  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1986
## 399  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1991
## 400  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1996
## 401  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1980
## 402  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2004
## 403  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 1995
## 404  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2006
## 405  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2003
## 406  17641.20       <NA>     <NA>      <NA>     43.361 757.8     MAL 2005
## 407    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1990
## 408    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1987
## 409    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1985
## 410    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1981
## 411    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1982
## 412    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1993
## 413    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2000
## 414    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2001
## 415    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1988
## 416    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1999
## 417    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2002
## 418    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1994
## 419    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1998
## 420    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1992
## 421    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1984
## 422    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1983
## 423    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1997
## 424    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2007
## 425    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2008
## 426    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1989
## 427    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1986
## 428    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1991
## 429    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1996
## 430    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1980
## 431    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2004
## 432    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 1995
## 433    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2006
## 434    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2003
## 435    135.17       <NA>     <NA>      <NA>     42.116 677.4     MAL 2005
## 436  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1989
## 437  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1984
## 438  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1992
## 439  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1993
## 440  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1987
## 441  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1995
## 442  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1996
## 443  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1997
## 444  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1994
## 445  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1981
## 446  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1986
## 447  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1990
## 448  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1991
## 449  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2003
## 450  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2004
## 451  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2005
## 452  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2006
## 453  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2007
## 454  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1985
## 455  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1999
## 456  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2000
## 457  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1988
## 458  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1980
## 459  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1998
## 460  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2008
## 461  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2002
## 462  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2009
## 463  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2001
## 464  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1983
## 465  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 2010
## 466  10324.73       <NA>     <NA>      <NA>     30.379 489.4     MOO 1982
## 467  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1989
## 468  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1984
## 469  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1992
## 470  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1993
## 471  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1987
## 472  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1995
## 473  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1996
## 474  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1997
## 475  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1994
## 476  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1981
## 477  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1986
## 478  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1990
## 479  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1991
## 480  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2003
## 481  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2004
## 482  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2005
## 483  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2006
## 484  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2007
## 485  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1985
## 486  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1999
## 487  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2000
## 488  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1988
## 489  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1980
## 490  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1998
## 491  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2008
## 492  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2002
## 493  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2009
## 494  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2001
## 495  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1983
## 496  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 2010
## 497  21576.20       <NA>     <NA>      <NA>     17.959 496.3     MOO 1982
## 498   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1987
## 499   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1985
## 500   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2000
## 501   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1990
## 502   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1983
## 503   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1986
## 504   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1988
## 505   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1984
## 506   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2001
## 507   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1993
## 508   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1991
## 509   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1992
## 510   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1982
## 511   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1980
## 512   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1981
## 513   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2005
## 514   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1995
## 515   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1989
## 516   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1999
## 517   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1994
## 518   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2006
## 519   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1996
## 520   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1998
## 521   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2004
## 522   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 1997
## 523   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2003
## 524   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2008
## 525   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2002
## 526   2829.16       <NA>     <NA>      <NA>     48.254 526.5     OLI 2007
## 527   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1987
## 528   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1985
## 529   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2000
## 530   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1990
## 531   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1983
## 532   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1986
## 533   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1988
## 534   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1984
## 535   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2001
## 536   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1993
## 537   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1991
## 538   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1992
## 539   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1982
## 540   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1980
## 541   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1981
## 542   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2005
## 543   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1995
## 544   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1989
## 545   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1999
## 546   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1994
## 547   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2006
## 548   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1996
## 549   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1998
## 550   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2004
## 551   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 1997
## 552   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2003
## 553   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2008
## 554   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2002
## 555   8521.57       <NA>     <NA>      <NA>     21.275 504.9     OLI 2007
## 556    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1987
## 557    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1985
## 558    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2000
## 559    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1990
## 560    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1983
## 561    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1986
## 562    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1988
## 563    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1984
## 564    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2001
## 565    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1993
## 566    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1991
## 567    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1992
## 568    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1982
## 569    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1980
## 570    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1981
## 571    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2005
## 572    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1995
## 573    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1989
## 574    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1999
## 575    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1994
## 576    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2006
## 577    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1996
## 578    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1998
## 579    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2004
## 580    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 1997
## 581    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2003
## 582    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2008
## 583    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2002
## 584    449.25       <NA>     <NA>      <NA>     44.330 504.4     OLI 2007
## 585    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1987
## 586    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1985
## 587    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2000
## 588    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1990
## 589    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1983
## 590    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1986
## 591    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1988
## 592    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1984
## 593    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2001
## 594    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1993
## 595    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1991
## 596    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1992
## 597    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1982
## 598    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1980
## 599    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1981
## 600    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2005
## 601    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1995
## 602    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1989
## 603    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1999
## 604    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1994
## 605    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2006
## 606    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1996
## 607    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1998
## 608    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2004
## 609    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 1997
## 610    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2003
## 611    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2008
## 612    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2002
## 613    120.24       <NA>     <NA>      <NA>     14.407 491.6     OLI 2007
## 614   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1987
## 615   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1985
## 616   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2000
## 617   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1990
## 618   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1983
## 619   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1986
## 620   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1988
## 621   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1984
## 622   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2001
## 623   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1993
## 624   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1991
## 625   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1992
## 626   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1982
## 627   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1980
## 628   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1981
## 629   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2005
## 630   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1995
## 631   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1989
## 632   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1999
## 633   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1994
## 634   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2006
## 635   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1996
## 636   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1998
## 637   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2004
## 638   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 1997
## 639   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2003
## 640   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2008
## 641   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2002
## 642   5657.13       <NA>     <NA>      <NA>     41.817 490.5     OLI 2007
## 643   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1987
## 644   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1985
## 645   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2000
## 646   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1990
## 647   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1983
## 648   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1986
## 649   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1988
## 650   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1984
## 651   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2001
## 652   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1993
## 653   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1991
## 654   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1992
## 655   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1982
## 656   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1980
## 657   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1981
## 658   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2005
## 659   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1995
## 660   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1989
## 661   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1999
## 662   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1994
## 663   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2006
## 664   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1996
## 665   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1998
## 666   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2004
## 667   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 1997
## 668   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2003
## 669   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2008
## 670   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2002
## 671   5654.88       <NA>     <NA>      <NA>     15.931 513.7     OLI 2007
## 672    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2002
## 673    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2008
## 674    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1996
## 675    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1999
## 676    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2001
## 677    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1995
## 678    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1989
## 679    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2003
## 680    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1994
## 681    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1986
## 682    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2007
## 683    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2004
## 684    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1987
## 685    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1992
## 686    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1993
## 687    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2000
## 688    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1988
## 689    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1980
## 690    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1981
## 691    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2006
## 692    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1990
## 693    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1991
## 694    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 2005
## 695    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1998
## 696    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1982
## 697    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1983
## 698    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1997
## 699    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1984
## 700    997.85       <NA>     <NA>      <NA>     10.283 516.7     PAF 1985
## 701    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2002
## 702    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2008
## 703    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1996
## 704    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1999
## 705    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2001
## 706    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1995
## 707    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1989
## 708    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2003
## 709    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1994
## 710    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1986
## 711    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2007
## 712    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2004
## 713    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1987
## 714    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1992
## 715    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1993
## 716    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2000
## 717    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1988
## 718    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1980
## 719    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1981
## 720    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2006
## 721    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1990
## 722    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1991
## 723    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 2005
## 724    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1998
## 725    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1982
## 726    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1983
## 727    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1997
## 728    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1984
## 729    718.01       <NA>     <NA>      <NA>     53.754 528.6     PAF 1985
## 730  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2002
## 731  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2008
## 732  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1996
## 733  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1999
## 734  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2001
## 735  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1995
## 736  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1989
## 737  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2003
## 738  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1994
## 739  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1986
## 740  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2007
## 741  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2004
## 742  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1987
## 743  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1992
## 744  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1993
## 745  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2000
## 746  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1988
## 747  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1980
## 748  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1981
## 749  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2006
## 750  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1990
## 751  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1991
## 752  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 2005
## 753  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1998
## 754  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1982
## 755  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1983
## 756  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1997
## 757  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1984
## 758  11728.28       <NA>     <NA>      <NA>     41.077 493.2     PAF 1985
## 759   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2002
## 760   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2008
## 761   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1996
## 762   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1999
## 763   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2001
## 764   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1995
## 765   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1989
## 766   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2003
## 767   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1994
## 768   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1986
## 769   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2007
## 770   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2004
## 771   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1987
## 772   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1992
## 773   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1993
## 774   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2000
## 775   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1988
## 776   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1980
## 777   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1981
## 778   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2006
## 779   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1990
## 780   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1991
## 781   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 2005
## 782   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1998
## 783   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1982
## 784   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1983
## 785   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1997
## 786   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1984
## 787   2776.70       <NA>     <NA>      <NA>     63.468 481.0     PAF 1985
## 788    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2002
## 789    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2008
## 790    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1996
## 791    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1999
## 792    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2001
## 793    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1995
## 794    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1989
## 795    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2003
## 796    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1994
## 797    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1986
## 798    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2007
## 799    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2004
## 800    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1987
## 801    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1992
## 802    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1993
## 803    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2000
## 804    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1988
## 805    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1980
## 806    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1981
## 807    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2006
## 808    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1990
## 809    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1991
## 810    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 2005
## 811    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1998
## 812    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1982
## 813    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1983
## 814    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1997
## 815    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1984
## 816    149.56       <NA>     <NA>      <NA>     50.619 508.6     PAF 1985
## 817   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2002
## 818   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2008
## 819   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1996
## 820   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1999
## 821   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2001
## 822   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1995
## 823   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1989
## 824   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2003
## 825   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1994
## 826   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1986
## 827   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2007
## 828   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2004
## 829   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1987
## 830   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1992
## 831   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1993
## 832   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2000
## 833   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1988
## 834   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1980
## 835   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1981
## 836   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2006
## 837   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1990
## 838   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1991
## 839   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 2005
## 840   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1998
## 841   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1982
## 842   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1983
## 843   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1997
## 844   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1984
## 845   4192.69       <NA>     <NA>      <NA>     38.370 502.3     PAF 1985
## 846    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2002
## 847    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2008
## 848    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1996
## 849    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1999
## 850    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2001
## 851    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1995
## 852    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1989
## 853    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2003
## 854    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1994
## 855    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1986
## 856    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2007
## 857    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2004
## 858    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1987
## 859    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1992
## 860    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1993
## 861    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2000
## 862    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1988
## 863    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1980
## 864    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1981
## 865    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2006
## 866    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1990
## 867    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1991
## 868    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 2005
## 869    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1998
## 870    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1982
## 871    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1983
## 872    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1997
## 873    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1984
## 874    136.50       <NA>     <NA>      <NA>     47.330 404.8     PAF 1985
## 875    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2002
## 876    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2008
## 877    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1996
## 878    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1999
## 879    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2001
## 880    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1995
## 881    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1989
## 882    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2003
## 883    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1994
## 884    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1986
## 885    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2007
## 886    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2004
## 887    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1987
## 888    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1992
## 889    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1993
## 890    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2000
## 891    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1988
## 892    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1980
## 893    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1981
## 894    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2006
## 895    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1990
## 896    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1991
## 897    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 2005
## 898    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1998
## 899    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1982
## 900    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1983
## 901    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1997
## 902    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1984
## 903    379.01       <NA>     <NA>      <NA>     52.426 437.7     PAF 1985
## 904   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1998
## 905   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1986
## 906   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1991
## 907   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1992
## 908   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1993
## 909   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1981
## 910   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1995
## 911   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1983
## 912   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1984
## 913   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1996
## 914   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1997
## 915   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2000
## 916   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2001
## 917   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2002
## 918   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2003
## 919   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2004
## 920   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2005
## 921   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2006
## 922   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2007
## 923   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 2008
## 924   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1994
## 925   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1999
## 926   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1980
## 927   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1985
## 928   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1982
## 929   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1987
## 930   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1988
## 931   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1989
## 932   7226.38       <NA>     <NA>      <NA>     32.868 516.8     PHA 1990
## 933   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1998
## 934   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1986
## 935   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1991
## 936   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1992
## 937   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1993
## 938   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1981
## 939   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1995
## 940   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1983
## 941   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1984
## 942   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1996
## 943   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1997
## 944   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2000
## 945   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2001
## 946   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2002
## 947   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2003
## 948   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2004
## 949   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2005
## 950   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2006
## 951   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2007
## 952   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 2008
## 953   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1994
## 954   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1999
## 955   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1980
## 956   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1985
## 957   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1982
## 958   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1987
## 959   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1988
## 960   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1989
## 961   3168.02       <NA>     <NA>      <NA>     41.716 531.9     PHA 1990
## 962     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2004
## 963     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2005
## 964     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2007
## 965     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1995
## 966     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1991
## 967     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2006
## 968     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1993
## 969     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1990
## 970     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2002
## 971     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2003
## 972     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1981
## 973     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1989
## 974     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1988
## 975     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1998
## 976     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1996
## 977     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2010
## 978     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1980
## 979     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2000
## 980     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1983
## 981     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2009
## 982     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2008
## 983     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1987
## 984     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 2001
## 985     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1982
## 986     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1986
## 987     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1992
## 988     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1994
## 989     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1984
## 990     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1997
## 991     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1985
## 992     62.00       <NA>     <NA>      <NA>     44.844 735.3     PRE 1999
## 993   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2004
## 994   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2005
## 995   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2007
## 996   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1995
## 997   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1991
## 998   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2006
## 999   2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1993
## 1000  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1990
## 1001  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2002
## 1002  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2003
## 1003  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1981
## 1004  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1989
## 1005  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1988
## 1006  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1998
## 1007  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1996
## 1008  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2010
## 1009  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1980
## 1010  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2000
## 1011  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1983
## 1012  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2009
## 1013  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2008
## 1014  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1987
## 1015  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 2001
## 1016  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1982
## 1017  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1986
## 1018  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1992
## 1019  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1994
## 1020  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1984
## 1021  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1997
## 1022  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1985
## 1023  2084.09       <NA>     <NA>      <NA>     43.944 735.5     PRE 1999
## 1024  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2004
## 1025  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2005
## 1026  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2007
## 1027  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1995
## 1028  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1991
## 1029  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2006
## 1030  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1993
## 1031  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1990
## 1032  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2002
## 1033  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2003
## 1034  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1981
## 1035  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1989
## 1036  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1988
## 1037  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1998
## 1038  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1996
## 1039  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2010
## 1040  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1980
## 1041  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2000
## 1042  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1983
## 1043  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2009
## 1044  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2008
## 1045  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1987
## 1046  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 2001
## 1047  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1982
## 1048  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1986
## 1049  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1992
## 1050  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1994
## 1051  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1984
## 1052  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1997
## 1053  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1985
## 1054  4161.84       <NA>     <NA>      <NA>     44.316 737.6     PRE 1999
## 1055 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2004
## 1056 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2005
## 1057 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2007
## 1058 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1995
## 1059 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1991
## 1060 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2006
## 1061 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1993
## 1062 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1990
## 1063 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2002
## 1064 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2003
## 1065 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1981
## 1066 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1989
## 1067 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1988
## 1068 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1998
## 1069 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1996
## 1070 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2010
## 1071 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1980
## 1072 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2000
## 1073 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1983
## 1074 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2009
## 1075 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2008
## 1076 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1987
## 1077 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 2001
## 1078 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1982
## 1079 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1986
## 1080 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1992
## 1081 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1994
## 1082 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1984
## 1083 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1997
## 1084 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1985
## 1085 41498.94       <NA>     <NA>      <NA>     40.749 697.7     PRE 1999
## 1086   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2004
## 1087   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2005
## 1088   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2007
## 1089   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1995
## 1090   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1991
## 1091   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2006
## 1092   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1993
## 1093   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1990
## 1094   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2002
## 1095   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2003
## 1096   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1981
## 1097   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1989
## 1098   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1988
## 1099   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1998
## 1100   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1996
## 1101   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2010
## 1102   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1980
## 1103   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2000
## 1104   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1983
## 1105   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2009
## 1106   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2008
## 1107   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1987
## 1108   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 2001
## 1109   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1982
## 1110   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1986
## 1111   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1992
## 1112   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1994
## 1113   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1984
## 1114   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1997
## 1115   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1985
## 1116   157.30       <NA>     <NA>      <NA>     38.863 721.5     PRE 1999
## 1117  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2004
## 1118  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2005
## 1119  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2007
## 1120  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1995
## 1121  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1991
## 1122  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2006
## 1123  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1993
## 1124  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1990
## 1125  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2002
## 1126  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2003
## 1127  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1981
## 1128  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1989
## 1129  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1988
## 1130  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1998
## 1131  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1996
## 1132  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2010
## 1133  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1980
## 1134  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2000
## 1135  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1983
## 1136  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2009
## 1137  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2008
## 1138  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1987
## 1139  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 2001
## 1140  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1982
## 1141  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1986
## 1142  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1992
## 1143  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1994
## 1144  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1984
## 1145  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1997
## 1146  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1985
## 1147  2767.95       <NA>     <NA>      <NA>     40.177 705.6     PRE 1999
## 1148    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1997
## 1149    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2004
## 1150    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1981
## 1151    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1996
## 1152    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2001
## 1153    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1986
## 1154    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1980
## 1155    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1995
## 1156    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2007
## 1157    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2008
## 1158    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1998
## 1159    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2000
## 1160    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1999
## 1161    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1983
## 1162    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2005
## 1163    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2006
## 1164    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1990
## 1165    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1992
## 1166    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1984
## 1167    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1985
## 1168    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1994
## 1169    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2003
## 1170    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1993
## 1171    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1982
## 1172    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 2002
## 1173    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1987
## 1174    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1988
## 1175    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1989
## 1176    43.99       <NA>     <NA>      <NA>     12.292 546.3     PUN 1991
## 1177  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1997
## 1178  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2004
## 1179  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1981
## 1180  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1996
## 1181  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2001
## 1182  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1986
## 1183  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1980
## 1184  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1995
## 1185  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2007
## 1186  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2008
## 1187  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1998
## 1188  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2000
## 1189  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1999
## 1190  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1983
## 1191  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2005
## 1192  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2006
## 1193  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1990
## 1194  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1992
## 1195  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1984
## 1196  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1985
## 1197  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1994
## 1198  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2003
## 1199  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1993
## 1200  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1982
## 1201  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 2002
## 1202  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1987
## 1203  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1988
## 1204  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1989
## 1205  6457.60       <NA>     <NA>      <NA>     44.276 505.9     PUN 1991
## 1206  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1997
## 1207  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2004
## 1208  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1981
## 1209  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1996
## 1210  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2001
## 1211  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1986
## 1212  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1980
## 1213  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1995
## 1214  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2007
## 1215  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2008
## 1216  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1998
## 1217  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2000
## 1218  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1999
## 1219  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1983
## 1220  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2005
## 1221  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2006
## 1222  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1990
## 1223  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1992
## 1224  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1984
## 1225  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1985
## 1226  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1994
## 1227  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2003
## 1228  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1993
## 1229  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1982
## 1230  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 2002
## 1231  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1987
## 1232  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1988
## 1233  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1989
## 1234  5272.17       <NA>     <NA>      <NA>     33.262 545.7     PUN 1991
## 1235    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1997
## 1236    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2004
## 1237    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1981
## 1238    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1996
## 1239    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2001
## 1240    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1986
## 1241    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1980
## 1242    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1995
## 1243    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2007
## 1244    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2008
## 1245    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1998
## 1246    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2000
## 1247    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1999
## 1248    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1983
## 1249    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2005
## 1250    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2006
## 1251    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1990
## 1252    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1992
## 1253    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1984
## 1254    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1985
## 1255    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1994
## 1256    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2003
## 1257    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1993
## 1258    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1982
## 1259    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 2002
## 1260    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1987
## 1261    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1988
## 1262    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1989
## 1263    16.41       <NA>     <NA>      <NA>     44.251 516.0     PUN 1991
## 1264   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1981
## 1265   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1983
## 1266   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1980
## 1267   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1985
## 1268   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1982
## 1269   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1987
## 1270   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1984
## 1271   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1993
## 1272   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1994
## 1273   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1995
## 1274   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1996
## 1275   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1997
## 1276   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1998
## 1277   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1986
## 1278   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2000
## 1279   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1988
## 1280   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1989
## 1281   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1990
## 1282   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1991
## 1283   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1992
## 1284   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2005
## 1285   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2006
## 1286   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2007
## 1287   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2008
## 1288   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2009
## 1289   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2010
## 1290   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 1999
## 1291   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2003
## 1292   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2004
## 1293   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2001
## 1294   228.84       <NA>     <NA>      <NA>     20.265 570.8     SAT 2002
## 1295  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1981
## 1296  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1983
## 1297  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1980
## 1298  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1985
## 1299  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1982
## 1300  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1987
## 1301  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1984
## 1302  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1993
## 1303  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1994
## 1304  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1995
## 1305  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1996
## 1306  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1997
## 1307  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1998
## 1308  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1986
## 1309  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2000
## 1310  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1988
## 1311  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1989
## 1312  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1990
## 1313  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1991
## 1314  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1992
## 1315  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2005
## 1316  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2006
## 1317  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2007
## 1318  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2008
## 1319  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2009
## 1320  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2010
## 1321  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 1999
## 1322  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2003
## 1323  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2004
## 1324  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2001
## 1325  1617.36       <NA>     <NA>      <NA>     29.554 568.5     SAT 2002
## 1326  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1996
## 1327  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1997
## 1328  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1998
## 1329  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1999
## 1330  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2000
## 1331  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2001
## 1332  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2002
## 1333  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2003
## 1334  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2004
## 1335  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1990
## 1336  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1991
## 1337  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1988
## 1338  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1989
## 1339  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1992
## 1340  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1993
## 1341  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1982
## 1342  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1983
## 1343  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1985
## 1344  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2005
## 1345  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2006
## 1346  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2007
## 1347  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 2008
## 1348  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1995
## 1349  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1986
## 1350  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1987
## 1351  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1981
## 1352  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1980
## 1353  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1994
## 1354  1249.56       <NA>     <NA>      <NA>     34.256 580.4     SHA 1984
## 1355  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2008
## 1356  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2006
## 1357  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2007
## 1358  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1994
## 1359  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1984
## 1360  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1993
## 1361  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1999
## 1362  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2004
## 1363  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1986
## 1364  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2002
## 1365  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1992
## 1366  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2005
## 1367  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1991
## 1368  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2000
## 1369  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1985
## 1370  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1982
## 1371  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1983
## 1372  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1998
## 1373  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2001
## 1374  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 2003
## 1375  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1997
## 1376  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1990
## 1377  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1987
## 1378  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1988
## 1379  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1989
## 1380  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1996
## 1381  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1981
## 1382  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1995
## 1383  1211.36       <NA>     <NA>      <NA>     32.229 437.3     SHI 1980
## 1384   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2008
## 1385   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2006
## 1386   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2007
## 1387   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1994
## 1388   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1984
## 1389   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1993
## 1390   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1999
## 1391   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2004
## 1392   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1986
## 1393   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2002
## 1394   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1992
## 1395   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2005
## 1396   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1991
## 1397   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2000
## 1398   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1985
## 1399   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1982
## 1400   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1983
## 1401   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1998
## 1402   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2001
## 1403   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 2003
## 1404   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1997
## 1405   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1990
## 1406   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1987
## 1407   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1988
## 1408   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1989
## 1409   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1996
## 1410   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1981
## 1411   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1995
## 1412   806.22       <NA>     <NA>      <NA>     31.280 461.2     SHI 1980
## 1413 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2008
## 1414 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2003
## 1415 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2007
## 1416 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1992
## 1417 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2004
## 1418 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2005
## 1419 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2006
## 1420 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1997
## 1421 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1980
## 1422 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1995
## 1423 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1996
## 1424 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2002
## 1425 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1986
## 1426 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1987
## 1427 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1988
## 1428 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1989
## 1429 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1990
## 1430 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1991
## 1431 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1981
## 1432 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1982
## 1433 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1983
## 1434 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1998
## 1435 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1994
## 1436 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1993
## 1437 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2001
## 1438 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 2000
## 1439 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1984
## 1440 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1985
## 1441 53752.73       <NA>     <NA>      <NA>     41.502 648.0     SKZ 1999
## 1442  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2002
## 1443  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2001
## 1444  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1997
## 1445  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2003
## 1446  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1999
## 1447  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2000
## 1448  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1991
## 1449  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1998
## 1450  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1982
## 1451  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2004
## 1452  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2005
## 1453  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1981
## 1454  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1994
## 1455  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1995
## 1456  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1996
## 1457  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1990
## 1458  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1987
## 1459  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1988
## 1460  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1989
## 1461  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1985
## 1462  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2006
## 1463  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1992
## 1464  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1980
## 1465  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2008
## 1466  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 2007
## 1467  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1993
## 1468  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1984
## 1469  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1986
## 1470  2201.10       <NA>     <NA>      <NA>     44.930 651.3     TSH 1983
## 1471   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2002
## 1472   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2001
## 1473   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1997
## 1474   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2003
## 1475   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1999
## 1476   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2000
## 1477   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1991
## 1478   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1998
## 1479   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1982
## 1480   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2004
## 1481   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2005
## 1482   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1981
## 1483   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1994
## 1484   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1995
## 1485   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1996
## 1486   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1990
## 1487   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1987
## 1488   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1988
## 1489   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1989
## 1490   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1985
## 1491   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2006
## 1492   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1992
## 1493   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1980
## 1494   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2008
## 1495   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 2007
## 1496   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1993
## 1497   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1984
## 1498   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1986
## 1499   143.86       <NA>     <NA>      <NA>     42.453 680.7     TSH 1983
## 1500  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2002
## 1501  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2001
## 1502  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1997
## 1503  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2003
## 1504  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1999
## 1505  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2000
## 1506  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1991
## 1507  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1998
## 1508  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1982
## 1509  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2004
## 1510  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2005
## 1511  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1981
## 1512  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1994
## 1513  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1995
## 1514  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1996
## 1515  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1990
## 1516  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1987
## 1517  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1988
## 1518  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1989
## 1519  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1985
## 1520  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2006
## 1521  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1992
## 1522  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1980
## 1523  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2008
## 1524  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 2007
## 1525  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1993
## 1526  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1984
## 1527  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1986
## 1528  5522.52       <NA>     <NA>      <NA>     47.302 624.2     TSH 1983
## 1529  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2002
## 1530  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2001
## 1531  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1997
## 1532  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2003
## 1533  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1999
## 1534  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2000
## 1535  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1991
## 1536  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1998
## 1537  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1982
## 1538  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2004
## 1539  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2005
## 1540  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1981
## 1541  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1994
## 1542  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1995
## 1543  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1996
## 1544  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1990
## 1545  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1987
## 1546  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1988
## 1547  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1989
## 1548  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1985
## 1549  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2006
## 1550  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1992
## 1551  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1980
## 1552  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2008
## 1553  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 2007
## 1554  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1993
## 1555  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1984
## 1556  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1986
## 1557  2490.33       <NA>     <NA>      <NA>     45.082 667.1     TSH 1983
## 1558  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2002
## 1559  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2001
## 1560  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1997
## 1561  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2003
## 1562  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1999
## 1563  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2000
## 1564  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1991
## 1565  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1998
## 1566  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1982
## 1567  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2004
## 1568  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2005
## 1569  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1981
## 1570  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1994
## 1571  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1995
## 1572  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1996
## 1573  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1990
## 1574  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1987
## 1575  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1988
## 1576  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1989
## 1577  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1985
## 1578  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2006
## 1579  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1992
## 1580  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1980
## 1581  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2008
## 1582  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 2007
## 1583  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1993
## 1584  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1984
## 1585  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1986
## 1586  2775.61       <NA>     <NA>      <NA>     45.836 622.2     TSH 1983
## 1587  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2008
## 1588  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2003
## 1589  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2001
## 1590  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1992
## 1591  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1998
## 1592  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2007
## 1593  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1999
## 1594  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2000
## 1595  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1991
## 1596  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1993
## 1597  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1994
## 1598  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2002
## 1599  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1997
## 1600  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1985
## 1601  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1987
## 1602  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1988
## 1603  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1989
## 1604  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2005
## 1605  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2006
## 1606  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1983
## 1607  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1984
## 1608  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1981
## 1609  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1982
## 1610  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1986
## 1611  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1980
## 1612  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1995
## 1613  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1996
## 1614  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 2004
## 1615  3353.19       <NA>     <NA>      <NA>     29.855 459.3     VLA 1990
## 1616 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2008
## 1617 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2003
## 1618 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2001
## 1619 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1992
## 1620 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1998
## 1621 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2007
## 1622 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1999
## 1623 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2000
## 1624 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1991
## 1625 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1993
## 1626 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1994
## 1627 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2002
## 1628 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1997
## 1629 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1985
## 1630 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1987
## 1631 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1988
## 1632 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1989
## 1633 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2005
## 1634 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2006
## 1635 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1983
## 1636 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1984
## 1637 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1981
## 1638 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1982
## 1639 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1986
## 1640 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1980
## 1641 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1995
## 1642 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1996
## 1643 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 2004
## 1644 17875.37       <NA>     <NA>      <NA>     37.055 484.8     VLA 1990
## 1645    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2008
## 1646    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2003
## 1647    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2001
## 1648    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1992
## 1649    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1998
## 1650    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2007
## 1651    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1999
## 1652    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2000
## 1653    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1991
## 1654    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1993
## 1655    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1994
## 1656    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2002
## 1657    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1997
## 1658    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1985
## 1659    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1987
## 1660    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1988
## 1661    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1989
## 1662    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2005
## 1663    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2006
## 1664    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1983
## 1665    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1984
## 1666    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1981
## 1667    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1982
## 1668    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1986
## 1669    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1980
## 1670    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1995
## 1671    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1996
## 1672    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 2004
## 1673    62.98       <NA>     <NA>      <NA>     15.319 492.4     VLA 1990
## 1674 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2008
## 1675 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2003
## 1676 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2001
## 1677 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1992
## 1678 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1998
## 1679 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2007
## 1680 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1999
## 1681 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2000
## 1682 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1991
## 1683 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1993
## 1684 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1994
## 1685 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2002
## 1686 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1997
## 1687 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1985
## 1688 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1987
## 1689 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1988
## 1690 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1989
## 1691 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2005
## 1692 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2006
## 1693 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1983
## 1694 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1984
## 1695 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1981
## 1696 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1982
## 1697 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1986
## 1698 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1980
## 1699 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1995
## 1700 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1996
## 1701 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 2004
## 1702 11750.99       <NA>     <NA>      <NA>     24.315 448.9     VLA 1990
## 1703  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2008
## 1704  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2003
## 1705  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2001
## 1706  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1992
## 1707  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1998
## 1708  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2007
## 1709  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1999
## 1710  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2000
## 1711  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1991
## 1712  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1993
## 1713  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1994
## 1714  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2002
## 1715  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1997
## 1716  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1985
## 1717  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1987
## 1718  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1988
## 1719  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1989
## 1720  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2005
## 1721  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2006
## 1722  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1983
## 1723  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1984
## 1724  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1981
## 1725  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1982
## 1726  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1986
## 1727  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1980
## 1728  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1995
## 1729  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1996
## 1730  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 2004
## 1731  7715.40       <NA>     <NA>      <NA>      8.825 510.7     VLA 1990
## 1732  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2008
## 1733  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2003
## 1734  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2001
## 1735  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1992
## 1736  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1998
## 1737  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2007
## 1738  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1999
## 1739  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2000
## 1740  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1991
## 1741  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1993
## 1742  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1994
## 1743  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2002
## 1744  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1997
## 1745  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1985
## 1746  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1987
## 1747  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1988
## 1748  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1989
## 1749  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2005
## 1750  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2006
## 1751  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1983
## 1752  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1984
## 1753  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1981
## 1754  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1982
## 1755  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1986
## 1756  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1980
## 1757  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1995
## 1758  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1996
## 1759  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 2004
## 1760  6967.41       <NA>     <NA>      <NA>     19.073 488.5     VLA 1990
## 1761  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2008
## 1762  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2003
## 1763  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2001
## 1764  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1992
## 1765  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1998
## 1766  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2007
## 1767  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1999
## 1768  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2000
## 1769  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1991
## 1770  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1993
## 1771  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1994
## 1772  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2002
## 1773  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1997
## 1774  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1985
## 1775  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1987
## 1776  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1988
## 1777  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1989
## 1778  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2005
## 1779  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2006
## 1780  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1983
## 1781  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1984
## 1782  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1981
## 1783  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1982
## 1784  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1986
## 1785  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1980
## 1786  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1995
## 1787  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1996
## 1788  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 2004
## 1789  8500.34       <NA>     <NA>      <NA>     13.089 481.3     VLA 1990
## 1790  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2008
## 1791  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2003
## 1792  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2001
## 1793  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1992
## 1794  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1998
## 1795  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2007
## 1796  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1999
## 1797  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2000
## 1798  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1991
## 1799  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1993
## 1800  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1994
## 1801  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2002
## 1802  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1997
## 1803  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1985
## 1804  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1987
## 1805  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1988
## 1806  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1989
## 1807  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2005
## 1808  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2006
## 1809  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1983
## 1810  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1984
## 1811  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1981
## 1812  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1982
## 1813  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1986
## 1814  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1980
## 1815  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1995
## 1816  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1996
## 1817  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 2004
## 1818  6446.06       <NA>     <NA>      <NA>      6.644 493.9     VLA 1990
## 1819  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2008
## 1820  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2003
## 1821  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2001
## 1822  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1992
## 1823  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1998
## 1824  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2007
## 1825  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1999
## 1826  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2000
## 1827  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1991
## 1828  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1993
## 1829  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1994
## 1830  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2002
## 1831  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1997
## 1832  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1985
## 1833  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1987
## 1834  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1988
## 1835  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1989
## 1836  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2005
## 1837  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2006
## 1838  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1983
## 1839  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1984
## 1840  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1981
## 1841  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1982
## 1842  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1986
## 1843  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1980
## 1844  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1995
## 1845  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1996
## 1846  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 2004
## 1847  9563.17       <NA>     <NA>      <NA>     11.478 507.5     VLA 1990
##      AnnualPrecip
## 1           596.7
## 2           350.5
## 3           382.6
## 4           497.9
## 5           674.1
## 6           852.6
## 7           575.8
## 8           930.0
## 9           738.3
## 10          515.5
## 11          832.5
## 12          773.6
## 13          569.3
## 14          490.8
## 15          529.9
## 16          516.6
## 17         1190.4
## 18          405.3
## 19          488.1
## 20          703.6
## 21          484.0
## 22          392.9
## 23          713.3
## 24          544.2
## 25          308.6
## 26          691.7
## 27          515.3
## 28          371.5
## 29          716.9
## 30          596.7
## 31          350.5
## 32          382.6
## 33          497.9
## 34          674.1
## 35          852.6
## 36          575.8
## 37          930.0
## 38          738.3
## 39          515.5
## 40          832.5
## 41          773.6
## 42          569.3
## 43          490.8
## 44          529.9
## 45          516.6
## 46         1190.4
## 47          405.3
## 48          488.1
## 49          703.6
## 50          484.0
## 51          392.9
## 52          713.3
## 53          544.2
## 54          308.6
## 55          691.7
## 56          515.3
## 57          371.5
## 58          716.9
## 59          596.7
## 60          350.5
## 61          382.6
## 62          497.9
## 63          674.1
## 64          852.6
## 65          575.8
## 66          930.0
## 67          738.3
## 68          515.5
## 69          832.5
## 70          773.6
## 71          569.3
## 72          490.8
## 73          529.9
## 74          516.6
## 75         1190.4
## 76          405.3
## 77          488.1
## 78          703.6
## 79          484.0
## 80          392.9
## 81          713.3
## 82          544.2
## 83          308.6
## 84          691.7
## 85          515.3
## 86          371.5
## 87          716.9
## 88          596.7
## 89          350.5
## 90          382.6
## 91          497.9
## 92          674.1
## 93          852.6
## 94          575.8
## 95          930.0
## 96          738.3
## 97          515.5
## 98          832.5
## 99          773.6
## 100         569.3
## 101         490.8
## 102         529.9
## 103         516.6
## 104        1190.4
## 105         405.3
## 106         488.1
## 107         703.6
## 108         484.0
## 109         392.9
## 110         713.3
## 111         544.2
## 112         308.6
## 113         691.7
## 114         515.3
## 115         371.5
## 116         716.9
## 117         596.7
## 118         350.5
## 119         382.6
## 120         497.9
## 121         674.1
## 122         852.6
## 123         575.8
## 124         930.0
## 125         738.3
## 126         515.5
## 127         832.5
## 128         773.6
## 129         569.3
## 130         490.8
## 131         529.9
## 132         516.6
## 133        1190.4
## 134         405.3
## 135         488.1
## 136         703.6
## 137         484.0
## 138         392.9
## 139         713.3
## 140         544.2
## 141         308.6
## 142         691.7
## 143         515.3
## 144         371.5
## 145         716.9
## 146         663.1
## 147         952.9
## 148         816.9
## 149         135.3
## 150         279.6
## 151         478.9
## 152            NA
## 153         504.7
## 154         567.1
## 155         427.6
## 156         374.5
## 157         242.9
## 158         504.7
## 159         583.5
## 160         388.0
## 161         292.1
## 162         690.7
## 163         310.7
## 164         236.4
## 165         362.7
## 166         347.7
## 167         240.9
## 168         292.6
## 169         489.7
## 170            NA
## 171         365.6
## 172         380.8
## 173         363.5
## 174         580.5
## 175         663.1
## 176         952.9
## 177         816.9
## 178         135.3
## 179         279.6
## 180         478.9
## 181            NA
## 182         504.7
## 183         567.1
## 184         427.6
## 185         374.5
## 186         242.9
## 187         504.7
## 188         583.5
## 189         388.0
## 190         292.1
## 191         690.7
## 192         310.7
## 193         236.4
## 194         362.7
## 195         347.7
## 196         240.9
## 197         292.6
## 198         489.7
## 199            NA
## 200         365.6
## 201         380.8
## 202         363.5
## 203         580.5
## 204         474.0
## 205         217.5
## 206         436.3
## 207         517.5
## 208         698.4
## 209        1044.5
## 210         388.4
## 211         714.5
## 212         431.0
## 213         232.7
## 214         407.9
## 215         689.4
## 216         609.3
## 217         547.3
## 218         472.4
## 219         614.0
## 220         685.7
## 221         853.8
## 222         690.9
## 223         632.9
## 224         765.2
## 225         648.4
## 226         528.4
## 227         455.2
## 228        1102.1
## 229         554.7
## 230         455.0
## 231         537.2
## 232         575.0
## 233         474.0
## 234         217.5
## 235         436.3
## 236         517.5
## 237         698.4
## 238        1044.5
## 239         388.4
## 240         714.5
## 241         431.0
## 242         232.7
## 243         407.9
## 244         689.4
## 245         609.3
## 246         547.3
## 247         472.4
## 248         614.0
## 249         685.7
## 250         853.8
## 251         690.9
## 252         632.9
## 253         765.2
## 254         648.4
## 255         528.4
## 256         455.2
## 257        1102.1
## 258         554.7
## 259         455.0
## 260         537.2
## 261         575.0
## 262         474.0
## 263         217.5
## 264         436.3
## 265         517.5
## 266         698.4
## 267        1044.5
## 268         388.4
## 269         714.5
## 270         431.0
## 271         232.7
## 272         407.9
## 273         689.4
## 274         609.3
## 275         547.3
## 276         472.4
## 277         614.0
## 278         685.7
## 279         853.8
## 280         690.9
## 281         632.9
## 282         765.2
## 283         648.4
## 284         528.4
## 285         455.2
## 286        1102.1
## 287         554.7
## 288         455.0
## 289         537.2
## 290         575.0
## 291         474.0
## 292         217.5
## 293         436.3
## 294         517.5
## 295         698.4
## 296        1044.5
## 297         388.4
## 298         714.5
## 299         431.0
## 300         232.7
## 301         407.9
## 302         689.4
## 303         609.3
## 304         547.3
## 305         472.4
## 306         614.0
## 307         685.7
## 308         853.8
## 309         690.9
## 310         632.9
## 311         765.2
## 312         648.4
## 313         528.4
## 314         455.2
## 315        1102.1
## 316         554.7
## 317         455.0
## 318         537.2
## 319         575.0
## 320         356.5
## 321         462.6
## 322         227.3
## 323         346.4
## 324         244.5
## 325            NA
## 326         330.5
## 327         283.5
## 328         340.9
## 329         311.9
## 330         755.0
## 331            NA
## 332         755.9
## 333         383.0
## 334         318.8
## 335         302.7
## 336         244.7
## 337         389.8
## 338         499.5
## 339         258.4
## 340         531.6
## 341            NA
## 342            NA
## 343        1081.8
## 344         627.0
## 345            NA
## 346            NA
## 347         308.4
## 348            NA
## 349         356.5
## 350         462.6
## 351         227.3
## 352         346.4
## 353         244.5
## 354            NA
## 355         330.5
## 356         283.5
## 357         340.9
## 358         311.9
## 359         755.0
## 360            NA
## 361         755.9
## 362         383.0
## 363         318.8
## 364         302.7
## 365         244.7
## 366         389.8
## 367         499.5
## 368         258.4
## 369         531.6
## 370            NA
## 371            NA
## 372        1081.8
## 373         627.0
## 374            NA
## 375            NA
## 376         308.4
## 377            NA
## 378         439.6
## 379         578.0
## 380         678.0
## 381         629.5
## 382         203.0
## 383         452.3
## 384        1353.0
## 385         534.3
## 386         501.8
## 387         812.7
## 388         393.4
## 389         274.8
## 390         684.3
## 391         486.7
## 392         793.5
## 393         673.0
## 394         492.0
## 395         460.7
## 396         568.0
## 397         706.4
## 398         595.0
## 399         507.8
## 400         697.5
## 401         718.0
## 402         794.5
## 403         506.3
## 404         866.2
## 405         291.5
## 406         306.2
## 407         439.6
## 408         578.0
## 409         678.0
## 410         629.5
## 411         203.0
## 412         452.3
## 413        1353.0
## 414         534.3
## 415         501.8
## 416         812.7
## 417         393.4
## 418         274.8
## 419         684.3
## 420         486.7
## 421         793.5
## 422         673.0
## 423         492.0
## 424         460.7
## 425         568.0
## 426         706.4
## 427         595.0
## 428         507.8
## 429         697.5
## 430         718.0
## 431         794.5
## 432         506.3
## 433         866.2
## 434         291.5
## 435         306.2
## 436         370.1
## 437         427.8
## 438         331.8
## 439         396.9
## 440         462.1
## 441         314.6
## 442         830.7
## 443         261.0
## 444         291.0
## 445         618.3
## 446         309.9
## 447         284.4
## 448         284.0
## 449         443.3
## 450         733.1
## 451         212.0
## 452         647.5
## 453         434.5
## 454         609.3
## 455         736.8
## 456        1305.8
## 457         467.9
## 458         659.6
## 459         449.3
## 460         355.9
## 461         265.4
## 462         399.8
## 463         583.6
## 464         321.5
## 465         502.3
## 466         194.0
## 467         370.1
## 468         427.8
## 469         331.8
## 470         396.9
## 471         462.1
## 472         314.6
## 473         830.7
## 474         261.0
## 475         291.0
## 476         618.3
## 477         309.9
## 478         284.4
## 479         284.0
## 480         443.3
## 481         733.1
## 482         212.0
## 483         647.5
## 484         434.5
## 485         609.3
## 486         736.8
## 487        1305.8
## 488         467.9
## 489         659.6
## 490         449.3
## 491         355.9
## 492         265.4
## 493         399.8
## 494         583.6
## 495         321.5
## 496         502.3
## 497         194.0
## 498         528.0
## 499         686.5
## 500         774.0
## 501         640.0
## 502         360.7
## 503         284.5
## 504         265.5
## 505         515.9
## 506         632.9
## 507         386.6
## 508         316.3
## 509         362.1
## 510         342.9
## 511         817.6
## 512         682.0
## 513         394.6
## 514         379.8
## 515         190.0
## 516         705.6
## 517         336.6
## 518         574.6
## 519         607.8
## 520         536.6
## 521        1031.5
## 522         400.1
## 523         460.6
## 524         522.9
## 525         291.4
## 526         498.7
## 527         528.0
## 528         686.5
## 529         774.0
## 530         640.0
## 531         360.7
## 532         284.5
## 533         265.5
## 534         515.9
## 535         632.9
## 536         386.6
## 537         316.3
## 538         362.1
## 539         342.9
## 540         817.6
## 541         682.0
## 542         394.6
## 543         379.8
## 544         190.0
## 545         705.6
## 546         336.6
## 547         574.6
## 548         607.8
## 549         536.6
## 550        1031.5
## 551         400.1
## 552         460.6
## 553         522.9
## 554         291.4
## 555         498.7
## 556         528.0
## 557         686.5
## 558         774.0
## 559         640.0
## 560         360.7
## 561         284.5
## 562         265.5
## 563         515.9
## 564         632.9
## 565         386.6
## 566         316.3
## 567         362.1
## 568         342.9
## 569         817.6
## 570         682.0
## 571         394.6
## 572         379.8
## 573         190.0
## 574         705.6
## 575         336.6
## 576         574.6
## 577         607.8
## 578         536.6
## 579        1031.5
## 580         400.1
## 581         460.6
## 582         522.9
## 583         291.4
## 584         498.7
## 585         528.0
## 586         686.5
## 587         774.0
## 588         640.0
## 589         360.7
## 590         284.5
## 591         265.5
## 592         515.9
## 593         632.9
## 594         386.6
## 595         316.3
## 596         362.1
## 597         342.9
## 598         817.6
## 599         682.0
## 600         394.6
## 601         379.8
## 602         190.0
## 603         705.6
## 604         336.6
## 605         574.6
## 606         607.8
## 607         536.6
## 608        1031.5
## 609         400.1
## 610         460.6
## 611         522.9
## 612         291.4
## 613         498.7
## 614         528.0
## 615         686.5
## 616         774.0
## 617         640.0
## 618         360.7
## 619         284.5
## 620         265.5
## 621         515.9
## 622         632.9
## 623         386.6
## 624         316.3
## 625         362.1
## 626         342.9
## 627         817.6
## 628         682.0
## 629         394.6
## 630         379.8
## 631         190.0
## 632         705.6
## 633         336.6
## 634         574.6
## 635         607.8
## 636         536.6
## 637        1031.5
## 638         400.1
## 639         460.6
## 640         522.9
## 641         291.4
## 642         498.7
## 643         528.0
## 644         686.5
## 645         774.0
## 646         640.0
## 647         360.7
## 648         284.5
## 649         265.5
## 650         515.9
## 651         632.9
## 652         386.6
## 653         316.3
## 654         362.1
## 655         342.9
## 656         817.6
## 657         682.0
## 658         394.6
## 659         379.8
## 660         190.0
## 661         705.6
## 662         336.6
## 663         574.6
## 664         607.8
## 665         536.6
## 666        1031.5
## 667         400.1
## 668         460.6
## 669         522.9
## 670         291.4
## 671         498.7
## 672         217.5
## 673            NA
## 674         640.6
## 675         631.3
## 676         631.6
## 677         390.9
## 678         246.0
## 679         283.1
## 680         278.6
## 681         301.3
## 682            NA
## 683         768.8
## 684         357.5
## 685         345.1
## 686         430.2
## 687        1007.6
## 688         483.5
## 689            NA
## 690            NA
## 691            NA
## 692         308.9
## 693         409.4
## 694         135.5
## 695         399.2
## 696         247.5
## 697         189.0
## 698         344.3
## 699         364.6
## 700         528.4
## 701         217.5
## 702            NA
## 703         640.6
## 704         631.3
## 705         631.6
## 706         390.9
## 707         246.0
## 708         283.1
## 709         278.6
## 710         301.3
## 711            NA
## 712         768.8
## 713         357.5
## 714         345.1
## 715         430.2
## 716        1007.6
## 717         483.5
## 718            NA
## 719            NA
## 720            NA
## 721         308.9
## 722         409.4
## 723         135.5
## 724         399.2
## 725         247.5
## 726         189.0
## 727         344.3
## 728         364.6
## 729         528.4
## 730         217.5
## 731            NA
## 732         640.6
## 733         631.3
## 734         631.6
## 735         390.9
## 736         246.0
## 737         283.1
## 738         278.6
## 739         301.3
## 740            NA
## 741         768.8
## 742         357.5
## 743         345.1
## 744         430.2
## 745        1007.6
## 746         483.5
## 747            NA
## 748            NA
## 749            NA
## 750         308.9
## 751         409.4
## 752         135.5
## 753         399.2
## 754         247.5
## 755         189.0
## 756         344.3
## 757         364.6
## 758         528.4
## 759         217.5
## 760            NA
## 761         640.6
## 762         631.3
## 763         631.6
## 764         390.9
## 765         246.0
## 766         283.1
## 767         278.6
## 768         301.3
## 769            NA
## 770         768.8
## 771         357.5
## 772         345.1
## 773         430.2
## 774        1007.6
## 775         483.5
## 776            NA
## 777            NA
## 778            NA
## 779         308.9
## 780         409.4
## 781         135.5
## 782         399.2
## 783         247.5
## 784         189.0
## 785         344.3
## 786         364.6
## 787         528.4
## 788         217.5
## 789            NA
## 790         640.6
## 791         631.3
## 792         631.6
## 793         390.9
## 794         246.0
## 795         283.1
## 796         278.6
## 797         301.3
## 798            NA
## 799         768.8
## 800         357.5
## 801         345.1
## 802         430.2
## 803        1007.6
## 804         483.5
## 805            NA
## 806            NA
## 807            NA
## 808         308.9
## 809         409.4
## 810         135.5
## 811         399.2
## 812         247.5
## 813         189.0
## 814         344.3
## 815         364.6
## 816         528.4
## 817         217.5
## 818            NA
## 819         640.6
## 820         631.3
## 821         631.6
## 822         390.9
## 823         246.0
## 824         283.1
## 825         278.6
## 826         301.3
## 827            NA
## 828         768.8
## 829         357.5
## 830         345.1
## 831         430.2
## 832        1007.6
## 833         483.5
## 834            NA
## 835            NA
## 836            NA
## 837         308.9
## 838         409.4
## 839         135.5
## 840         399.2
## 841         247.5
## 842         189.0
## 843         344.3
## 844         364.6
## 845         528.4
## 846         217.5
## 847            NA
## 848         640.6
## 849         631.3
## 850         631.6
## 851         390.9
## 852         246.0
## 853         283.1
## 854         278.6
## 855         301.3
## 856            NA
## 857         768.8
## 858         357.5
## 859         345.1
## 860         430.2
## 861        1007.6
## 862         483.5
## 863            NA
## 864            NA
## 865            NA
## 866         308.9
## 867         409.4
## 868         135.5
## 869         399.2
## 870         247.5
## 871         189.0
## 872         344.3
## 873         364.6
## 874         528.4
## 875         217.5
## 876            NA
## 877         640.6
## 878         631.3
## 879         631.6
## 880         390.9
## 881         246.0
## 882         283.1
## 883         278.6
## 884         301.3
## 885            NA
## 886         768.8
## 887         357.5
## 888         345.1
## 889         430.2
## 890        1007.6
## 891         483.5
## 892            NA
## 893            NA
## 894            NA
## 895         308.9
## 896         409.4
## 897         135.5
## 898         399.2
## 899         247.5
## 900         189.0
## 901         344.3
## 902         364.6
## 903         528.4
## 904         505.4
## 905         288.7
## 906         319.4
## 907         388.1
## 908         528.6
## 909            NA
## 910         545.8
## 911         375.6
## 912         589.1
## 913         755.2
## 914         374.0
## 915        1067.9
## 916         711.5
## 917         163.3
## 918         306.8
## 919         715.7
## 920         178.6
## 921         547.4
## 922         518.9
## 923         216.6
## 924         315.2
## 925         638.3
## 926            NA
## 927         639.2
## 928            NA
## 929          28.0
## 930         302.0
## 931         546.8
## 932         542.4
## 933         505.4
## 934         288.7
## 935         319.4
## 936         388.1
## 937         528.6
## 938            NA
## 939         545.8
## 940         375.6
## 941         589.1
## 942         755.2
## 943         374.0
## 944        1067.9
## 945         711.5
## 946         163.3
## 947         306.8
## 948         715.7
## 949         178.6
## 950         547.4
## 951         518.9
## 952         216.6
## 953         315.2
## 954         638.3
## 955            NA
## 956         639.2
## 957            NA
## 958          28.0
## 959         302.0
## 960         546.8
## 961         542.4
## 962         912.0
## 963         438.4
## 964         416.8
## 965         857.8
## 966         525.5
## 967         959.4
## 968         532.1
## 969         597.6
## 970         460.9
## 971         385.6
## 972            NA
## 973         728.9
## 974         643.4
## 975         636.5
## 976        1009.2
## 977         878.4
## 978            NA
## 979        1620.8
## 980         633.9
## 981        1045.7
## 982            NA
## 983         811.8
## 984         608.8
## 985            NA
## 986         647.8
## 987         652.0
## 988         281.4
## 989         734.4
## 990         878.9
## 991         749.0
## 992         874.4
## 993         912.0
## 994         438.4
## 995         416.8
## 996         857.8
## 997         525.5
## 998         959.4
## 999         532.1
## 1000        597.6
## 1001        460.9
## 1002        385.6
## 1003           NA
## 1004        728.9
## 1005        643.4
## 1006        636.5
## 1007       1009.2
## 1008        878.4
## 1009           NA
## 1010       1620.8
## 1011        633.9
## 1012       1045.7
## 1013           NA
## 1014        811.8
## 1015        608.8
## 1016           NA
## 1017        647.8
## 1018        652.0
## 1019        281.4
## 1020        734.4
## 1021        878.9
## 1022        749.0
## 1023        874.4
## 1024        912.0
## 1025        438.4
## 1026        416.8
## 1027        857.8
## 1028        525.5
## 1029        959.4
## 1030        532.1
## 1031        597.6
## 1032        460.9
## 1033        385.6
## 1034           NA
## 1035        728.9
## 1036        643.4
## 1037        636.5
## 1038       1009.2
## 1039        878.4
## 1040           NA
## 1041       1620.8
## 1042        633.9
## 1043       1045.7
## 1044           NA
## 1045        811.8
## 1046        608.8
## 1047           NA
## 1048        647.8
## 1049        652.0
## 1050        281.4
## 1051        734.4
## 1052        878.9
## 1053        749.0
## 1054        874.4
## 1055        912.0
## 1056        438.4
## 1057        416.8
## 1058        857.8
## 1059        525.5
## 1060        959.4
## 1061        532.1
## 1062        597.6
## 1063        460.9
## 1064        385.6
## 1065           NA
## 1066        728.9
## 1067        643.4
## 1068        636.5
## 1069       1009.2
## 1070        878.4
## 1071           NA
## 1072       1620.8
## 1073        633.9
## 1074       1045.7
## 1075           NA
## 1076        811.8
## 1077        608.8
## 1078           NA
## 1079        647.8
## 1080        652.0
## 1081        281.4
## 1082        734.4
## 1083        878.9
## 1084        749.0
## 1085        874.4
## 1086        912.0
## 1087        438.4
## 1088        416.8
## 1089        857.8
## 1090        525.5
## 1091        959.4
## 1092        532.1
## 1093        597.6
## 1094        460.9
## 1095        385.6
## 1096           NA
## 1097        728.9
## 1098        643.4
## 1099        636.5
## 1100       1009.2
## 1101        878.4
## 1102           NA
## 1103       1620.8
## 1104        633.9
## 1105       1045.7
## 1106           NA
## 1107        811.8
## 1108        608.8
## 1109           NA
## 1110        647.8
## 1111        652.0
## 1112        281.4
## 1113        734.4
## 1114        878.9
## 1115        749.0
## 1116        874.4
## 1117        912.0
## 1118        438.4
## 1119        416.8
## 1120        857.8
## 1121        525.5
## 1122        959.4
## 1123        532.1
## 1124        597.6
## 1125        460.9
## 1126        385.6
## 1127           NA
## 1128        728.9
## 1129        643.4
## 1130        636.5
## 1131       1009.2
## 1132        878.4
## 1133           NA
## 1134       1620.8
## 1135        633.9
## 1136       1045.7
## 1137           NA
## 1138        811.8
## 1139        608.8
## 1140           NA
## 1141        647.8
## 1142        652.0
## 1143        281.4
## 1144        734.4
## 1145        878.9
## 1146        749.0
## 1147        874.4
## 1148        439.1
## 1149        600.6
## 1150           NA
## 1151        758.1
## 1152        886.5
## 1153        487.3
## 1154           NA
## 1155        491.5
## 1156        685.3
## 1157        358.1
## 1158        578.4
## 1159       1132.9
## 1160        922.1
## 1161        205.1
## 1162        384.0
## 1163        652.0
## 1164        334.3
## 1165        316.2
## 1166        416.4
## 1167        834.3
## 1168        413.1
## 1169        400.1
## 1170        394.2
## 1171        269.1
## 1172        397.2
## 1173        401.7
## 1174        620.2
## 1175        463.7
## 1176        459.4
## 1177        439.1
## 1178        600.6
## 1179           NA
## 1180        758.1
## 1181        886.5
## 1182        487.3
## 1183           NA
## 1184        491.5
## 1185        685.3
## 1186        358.1
## 1187        578.4
## 1188       1132.9
## 1189        922.1
## 1190        205.1
## 1191        384.0
## 1192        652.0
## 1193        334.3
## 1194        316.2
## 1195        416.4
## 1196        834.3
## 1197        413.1
## 1198        400.1
## 1199        394.2
## 1200        269.1
## 1201        397.2
## 1202        401.7
## 1203        620.2
## 1204        463.7
## 1205        459.4
## 1206        439.1
## 1207        600.6
## 1208           NA
## 1209        758.1
## 1210        886.5
## 1211        487.3
## 1212           NA
## 1213        491.5
## 1214        685.3
## 1215        358.1
## 1216        578.4
## 1217       1132.9
## 1218        922.1
## 1219        205.1
## 1220        384.0
## 1221        652.0
## 1222        334.3
## 1223        316.2
## 1224        416.4
## 1225        834.3
## 1226        413.1
## 1227        400.1
## 1228        394.2
## 1229        269.1
## 1230        397.2
## 1231        401.7
## 1232        620.2
## 1233        463.7
## 1234        459.4
## 1235        439.1
## 1236        600.6
## 1237           NA
## 1238        758.1
## 1239        886.5
## 1240        487.3
## 1241           NA
## 1242        491.5
## 1243        685.3
## 1244        358.1
## 1245        578.4
## 1246       1132.9
## 1247        922.1
## 1248        205.1
## 1249        384.0
## 1250        652.0
## 1251        334.3
## 1252        316.2
## 1253        416.4
## 1254        834.3
## 1255        413.1
## 1256        400.1
## 1257        394.2
## 1258        269.1
## 1259        397.2
## 1260        401.7
## 1261        620.2
## 1262        463.7
## 1263        459.4
## 1264        767.3
## 1265        466.0
## 1266        729.7
## 1267        577.4
## 1268        249.9
## 1269        678.2
## 1270        739.5
## 1271        472.0
## 1272        383.6
## 1273        561.2
## 1274        622.7
## 1275        591.0
## 1276        558.5
## 1277        234.4
## 1278        885.5
## 1279        295.1
## 1280        551.6
## 1281        754.8
## 1282        352.9
## 1283        405.0
## 1284        252.1
## 1285        528.3
## 1286        546.5
## 1287        316.1
## 1288        656.8
## 1289        631.9
## 1290        770.2
## 1291        293.1
## 1292        849.1
## 1293        591.2
## 1294        116.5
## 1295        767.3
## 1296        466.0
## 1297        729.7
## 1298        577.4
## 1299        249.9
## 1300        678.2
## 1301        739.5
## 1302        472.0
## 1303        383.6
## 1304        561.2
## 1305        622.7
## 1306        591.0
## 1307        558.5
## 1308        234.4
## 1309        885.5
## 1310        295.1
## 1311        551.6
## 1312        754.8
## 1313        352.9
## 1314        405.0
## 1315        252.1
## 1316        528.3
## 1317        546.5
## 1318        316.1
## 1319        656.8
## 1320        631.9
## 1321        770.2
## 1322        293.1
## 1323        849.1
## 1324        591.2
## 1325        116.5
## 1326        770.7
## 1327        479.4
## 1328        380.1
## 1329        775.4
## 1330       1229.9
## 1331        726.8
## 1332        295.6
## 1333        386.9
## 1334        512.6
## 1335        416.3
## 1336        385.7
## 1337        454.7
## 1338        378.3
## 1339        424.0
## 1340        422.7
## 1341           NA
## 1342           NA
## 1343           NA
## 1344        364.9
## 1345        633.8
## 1346        610.4
## 1347        372.0
## 1348        529.4
## 1349           NA
## 1350        538.7
## 1351           NA
## 1352           NA
## 1353        269.8
## 1354           NA
## 1355        335.2
## 1356        319.4
## 1357        609.2
## 1358        327.2
## 1359        532.4
## 1360        355.4
## 1361        913.1
## 1362        358.2
## 1363        397.5
## 1364        164.1
## 1365        459.4
## 1366        259.7
## 1367        221.2
## 1368       1149.9
## 1369        520.0
## 1370        300.0
## 1371        296.6
## 1372        430.3
## 1373        505.5
## 1374        525.1
## 1375        343.9
## 1376        488.9
## 1377        304.0
## 1378        395.6
## 1379        285.6
## 1380        854.5
## 1381           NA
## 1382        476.7
## 1383           NA
## 1384        335.2
## 1385        319.4
## 1386        609.2
## 1387        327.2
## 1388        532.4
## 1389        355.4
## 1390        913.1
## 1391        358.2
## 1392        397.5
## 1393        164.1
## 1394        459.4
## 1395        259.7
## 1396        221.2
## 1397       1149.9
## 1398        520.0
## 1399        300.0
## 1400        296.6
## 1401        430.3
## 1402        505.5
## 1403        525.1
## 1404        343.9
## 1405        488.9
## 1406        304.0
## 1407        395.6
## 1408        285.6
## 1409        854.5
## 1410           NA
## 1411        476.7
## 1412           NA
## 1413        581.9
## 1414        293.0
## 1415        393.9
## 1416        434.5
## 1417        996.0
## 1418        560.7
## 1419        807.3
## 1420        689.8
## 1421        649.3
## 1422        886.7
## 1423        701.7
## 1424        318.9
## 1425        445.2
## 1426        631.5
## 1427        559.8
## 1428        695.3
## 1429        463.2
## 1430        315.0
## 1431        604.8
## 1432        427.2
## 1433        490.4
## 1434        641.1
## 1435        321.9
## 1436        470.9
## 1437        719.6
## 1438       1909.8
## 1439        727.2
## 1440        607.4
## 1441        980.2
## 1442        203.7
## 1443        720.5
## 1444        601.8
## 1445        230.5
## 1446        767.0
## 1447       1099.6
## 1448        354.5
## 1449        727.9
## 1450        320.6
## 1451        819.5
## 1452        405.1
## 1453        665.3
## 1454        528.6
## 1455        612.1
## 1456        673.4
## 1457        602.1
## 1458        693.7
## 1459        575.4
## 1460        497.8
## 1461        748.0
## 1462        535.7
## 1463        430.2
## 1464        572.1
## 1465        188.4
## 1466        569.6
## 1467        386.9
## 1468        652.3
## 1469        398.4
## 1470        557.5
## 1471        203.7
## 1472        720.5
## 1473        601.8
## 1474        230.5
## 1475        767.0
## 1476       1099.6
## 1477        354.5
## 1478        727.9
## 1479        320.6
## 1480        819.5
## 1481        405.1
## 1482        665.3
## 1483        528.6
## 1484        612.1
## 1485        673.4
## 1486        602.1
## 1487        693.7
## 1488        575.4
## 1489        497.8
## 1490        748.0
## 1491        535.7
## 1492        430.2
## 1493        572.1
## 1494        188.4
## 1495        569.6
## 1496        386.9
## 1497        652.3
## 1498        398.4
## 1499        557.5
## 1500        203.7
## 1501        720.5
## 1502        601.8
## 1503        230.5
## 1504        767.0
## 1505       1099.6
## 1506        354.5
## 1507        727.9
## 1508        320.6
## 1509        819.5
## 1510        405.1
## 1511        665.3
## 1512        528.6
## 1513        612.1
## 1514        673.4
## 1515        602.1
## 1516        693.7
## 1517        575.4
## 1518        497.8
## 1519        748.0
## 1520        535.7
## 1521        430.2
## 1522        572.1
## 1523        188.4
## 1524        569.6
## 1525        386.9
## 1526        652.3
## 1527        398.4
## 1528        557.5
## 1529        203.7
## 1530        720.5
## 1531        601.8
## 1532        230.5
## 1533        767.0
## 1534       1099.6
## 1535        354.5
## 1536        727.9
## 1537        320.6
## 1538        819.5
## 1539        405.1
## 1540        665.3
## 1541        528.6
## 1542        612.1
## 1543        673.4
## 1544        602.1
## 1545        693.7
## 1546        575.4
## 1547        497.8
## 1548        748.0
## 1549        535.7
## 1550        430.2
## 1551        572.1
## 1552        188.4
## 1553        569.6
## 1554        386.9
## 1555        652.3
## 1556        398.4
## 1557        557.5
## 1558        203.7
## 1559        720.5
## 1560        601.8
## 1561        230.5
## 1562        767.0
## 1563       1099.6
## 1564        354.5
## 1565        727.9
## 1566        320.6
## 1567        819.5
## 1568        405.1
## 1569        665.3
## 1570        528.6
## 1571        612.1
## 1572        673.4
## 1573        602.1
## 1574        693.7
## 1575        575.4
## 1576        497.8
## 1577        748.0
## 1578        535.7
## 1579        430.2
## 1580        572.1
## 1581        188.4
## 1582        569.6
## 1583        386.9
## 1584        652.3
## 1585        398.4
## 1586        557.5
## 1587        301.4
## 1588        528.9
## 1589        610.8
## 1590        549.2
## 1591        421.3
## 1592        693.0
## 1593        899.5
## 1594       1044.4
## 1595        330.7
## 1596        463.7
## 1597        360.3
## 1598        218.9
## 1599        278.4
## 1600        768.3
## 1601        450.6
## 1602        469.5
## 1603        461.5
## 1604        287.3
## 1605        640.7
## 1606           NA
## 1607        514.2
## 1608           NA
## 1609           NA
## 1610        309.9
## 1611           NA
## 1612        469.0
## 1613        726.1
## 1614        507.1
## 1615        371.1
## 1616        301.4
## 1617        528.9
## 1618        610.8
## 1619        549.2
## 1620        421.3
## 1621        693.0
## 1622        899.5
## 1623       1044.4
## 1624        330.7
## 1625        463.7
## 1626        360.3
## 1627        218.9
## 1628        278.4
## 1629        768.3
## 1630        450.6
## 1631        469.5
## 1632        461.5
## 1633        287.3
## 1634        640.7
## 1635           NA
## 1636        514.2
## 1637           NA
## 1638           NA
## 1639        309.9
## 1640           NA
## 1641        469.0
## 1642        726.1
## 1643        507.1
## 1644        371.1
## 1645        301.4
## 1646        528.9
## 1647        610.8
## 1648        549.2
## 1649        421.3
## 1650        693.0
## 1651        899.5
## 1652       1044.4
## 1653        330.7
## 1654        463.7
## 1655        360.3
## 1656        218.9
## 1657        278.4
## 1658        768.3
## 1659        450.6
## 1660        469.5
## 1661        461.5
## 1662        287.3
## 1663        640.7
## 1664           NA
## 1665        514.2
## 1666           NA
## 1667           NA
## 1668        309.9
## 1669           NA
## 1670        469.0
## 1671        726.1
## 1672        507.1
## 1673        371.1
## 1674        301.4
## 1675        528.9
## 1676        610.8
## 1677        549.2
## 1678        421.3
## 1679        693.0
## 1680        899.5
## 1681       1044.4
## 1682        330.7
## 1683        463.7
## 1684        360.3
## 1685        218.9
## 1686        278.4
## 1687        768.3
## 1688        450.6
## 1689        469.5
## 1690        461.5
## 1691        287.3
## 1692        640.7
## 1693           NA
## 1694        514.2
## 1695           NA
## 1696           NA
## 1697        309.9
## 1698           NA
## 1699        469.0
## 1700        726.1
## 1701        507.1
## 1702        371.1
## 1703        301.4
## 1704        528.9
## 1705        610.8
## 1706        549.2
## 1707        421.3
## 1708        693.0
## 1709        899.5
## 1710       1044.4
## 1711        330.7
## 1712        463.7
## 1713        360.3
## 1714        218.9
## 1715        278.4
## 1716        768.3
## 1717        450.6
## 1718        469.5
## 1719        461.5
## 1720        287.3
## 1721        640.7
## 1722           NA
## 1723        514.2
## 1724           NA
## 1725           NA
## 1726        309.9
## 1727           NA
## 1728        469.0
## 1729        726.1
## 1730        507.1
## 1731        371.1
## 1732        301.4
## 1733        528.9
## 1734        610.8
## 1735        549.2
## 1736        421.3
## 1737        693.0
## 1738        899.5
## 1739       1044.4
## 1740        330.7
## 1741        463.7
## 1742        360.3
## 1743        218.9
## 1744        278.4
## 1745        768.3
## 1746        450.6
## 1747        469.5
## 1748        461.5
## 1749        287.3
## 1750        640.7
## 1751           NA
## 1752        514.2
## 1753           NA
## 1754           NA
## 1755        309.9
## 1756           NA
## 1757        469.0
## 1758        726.1
## 1759        507.1
## 1760        371.1
## 1761        301.4
## 1762        528.9
## 1763        610.8
## 1764        549.2
## 1765        421.3
## 1766        693.0
## 1767        899.5
## 1768       1044.4
## 1769        330.7
## 1770        463.7
## 1771        360.3
## 1772        218.9
## 1773        278.4
## 1774        768.3
## 1775        450.6
## 1776        469.5
## 1777        461.5
## 1778        287.3
## 1779        640.7
## 1780           NA
## 1781        514.2
## 1782           NA
## 1783           NA
## 1784        309.9
## 1785           NA
## 1786        469.0
## 1787        726.1
## 1788        507.1
## 1789        371.1
## 1790        301.4
## 1791        528.9
## 1792        610.8
## 1793        549.2
## 1794        421.3
## 1795        693.0
## 1796        899.5
## 1797       1044.4
## 1798        330.7
## 1799        463.7
## 1800        360.3
## 1801        218.9
## 1802        278.4
## 1803        768.3
## 1804        450.6
## 1805        469.5
## 1806        461.5
## 1807        287.3
## 1808        640.7
## 1809           NA
## 1810        514.2
## 1811           NA
## 1812           NA
## 1813        309.9
## 1814           NA
## 1815        469.0
## 1816        726.1
## 1817        507.1
## 1818        371.1
## 1819        301.4
## 1820        528.9
## 1821        610.8
## 1822        549.2
## 1823        421.3
## 1824        693.0
## 1825        899.5
## 1826       1044.4
## 1827        330.7
## 1828        463.7
## 1829        360.3
## 1830        218.9
## 1831        278.4
## 1832        768.3
## 1833        450.6
## 1834        469.5
## 1835        461.5
## 1836        287.3
## 1837        640.7
## 1838           NA
## 1839        514.2
## 1840           NA
## 1841           NA
## 1842        309.9
## 1843           NA
## 1844        469.0
## 1845        726.1
## 1846        507.1
## 1847        371.1
## 
## > rm(FireScars_WoodyImpact)
## 
## > rm(krugerFireScarAnalysis_brick)
## 
## > library(ggplot2)
## 
## > library(ggthemes)
## 
## > myTheme_big <- theme_tufte() + theme(text = element_text(family = "sans", 
## +     size = 17), axis.line = element_line(size = 0.3))
## 
## > tileTable <- ddply(.data = firescar_final, .(HERBACEOUS, 
## +     WOODYIMPAC, CAUSE), summarize, Area = sum(AREA_HA), Count = length(WOODYIMPAC))
## 
## > tileTable$CountPercent = tileTable$Count/sum(tileTable$Count)
## 
## > tileTable$AreaPercent = tileTable$Area/sum(tileTable$Area)
## 
## > Tile <- ggplot(aes(x = HERBACEOUS, y = WOODYIMPAC, 
## +     fill = AreaPercent), data = tileTable)
## 
## > Tile + scale_fill_continuous("Percent", high = "red", 
## +     low = "darkblue") + ylab("Impact on Woody Cover") + xlab("\nImpact on Herbaceous Materi ..." ... [TRUNCATED]
```

<figure><img src='figure_rmd/RunAnalysis1.png'  style='display: block'><br><figcaption>Figure 1: plot of chunk RunAnalysis</figcaption></figure><br>

```
## 
## > TileCause <- ggplot(aes(x = CAUSE, y = WOODYIMPAC, 
## +     fill = AreaPercent), data = tileTable)
## 
## > TileCause + scale_fill_continuous("Percent", high = "darkred", 
## +     low = "darkblue") + ylab("Impact on Woody Cover") + xlab("\nCause") + 
## +     g .... [TRUNCATED]
```

<figure><img src='figure_rmd/RunAnalysis2.png'  style='display: block'><br><figcaption>Figure 2: plot of chunk RunAnalysis</figcaption></figure><br>

```
## 
## > yearTable <- ddply(.data = firescar_final, .(HERBACEOUS, 
## +     WOODYIMPAC, Year), summarize, WoodyCount = length(WOODYIMPAC), 
## +     HerbCount = le .... [TRUNCATED]
```
Run the intensity model analysis


**Figures**

<figure><img src='figure_rmd/logicalFramework.png'  style='display: block'><br><figcaption>Figure 5: </figcaption></figure><br>
A. Average woody cover by mean annual precipitation at 50 mm/yr intervals (points); line is significant at P < 0.01, R^2 = 0.5862. B. Relationship between fuel load accumulation (grass biomass) and mean annual precipitation, derived from Figure 1, Govender et al. 2006. C. Relationship between fireline intensity and fuel load given an average rate of spread of .04 m/s. D. Relationship between average fireline intensity and mean annual precipitation, from Figure 2, Govender et al. 2006. E. Probability of top kill as a function of fire intensity by mean annual precipitation, as a function of Higgins et al. 2011 and the relationship shown in D.


<figure><img src='figure_rmd/FRP_by_MAP_Season.png'  style='display: block'><br><figcaption>Figure 6: Fire radiative power by mean annual precipitation, subdivided by (A.) season of burn and (B.) geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FRP_by_WoodyCover.png'  style='display: block'><br><figcaption>Figure 7: Fire radiative power by percent woody cover, subdivided by (A.) season of burn and (B.) geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FRP_by_GLY.png'  style='display: block'><br><figcaption>Figure 8: Fire radiative power by geology and season of burn.</figcaption></figure><br>

<figure><img src='figure_rmd/FRP_avg_by_MAP.png'  style='display: block'><br><figcaption>Figure 9: Fire radiative power, averaged by 50 mm / yr MAP</figcaption></figure><br>
**Kruskal-Wallis Tests**

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  FRP by Geology
## Kruskal-Wallis chi-squared = 13.26, df = 2, p-value = 0.001321
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  FRP by as.factor(Season)
## Kruskal-Wallis chi-squared = 0.0077, df = 1, p-value = 0.9301
```

No significant difference by seasonality, but there is a difference by geologic parent material.



<figure><img src='figure_rmd/ModeledFirelineIntensity.png'  style='display: block'><br><figcaption>Figure 10: Modeled fireline intensity across Kruger National Park</figcaption></figure><br>

<figure><img src='figure_rmd/ModeledProbTopKill.png'  style='display: block'><br><figcaption>Figure 11: Modeled probability of topkill of a one meter sapling across Kruger National Park</figcaption></figure><br>

<figure><img src='figure_rmd/FLI_By_WoodyCover.png'  style='display: block'><br><figcaption>Figure 12: Modeled average fireline intensity by percent woody cover and geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FLI_By_MAP.png'  style='display: block'><br><figcaption>Figure 13: Modeled fireline intensity by MAP and geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FLI_By_NavashniStandards.png'  style='display: block'><br><figcaption>Figure 14: Percentage of probable fires by fireline intensity classes, as defined by Govender et al. 2006</figcaption></figure><br>


<figure><img src='figure_rmd/FLI_By_AvgProbTop.png'  style='display: block'><br><figcaption>Figure 15: Probability of Top Kill by Mean Annual Precipitation, modeled.</figcaption></figure><br>



**Intensity by Woody Cover and MAP: GLM Results**

<table cellspacing="0" style="border: none;">
  <tr>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 1</b></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 2</b></th>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">(Intercept)</td>
    <td style="padding-right: 12px; border: none;">751.11 (2.34)<sup style="vertical-align: 4px;">***</sup></td>
    <td style="padding-right: 12px; border: none;">750.98 (2.34)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">WoodyCover</td>
    <td style="padding-right: 12px; border: none;"></td>
    <td style="padding-right: 12px; border: none;">0.00 (0.00)<sup style="vertical-align: 4px;">**</sup></td>
  </tr>
  <tr>
    <td style="border-top: 1px solid black;">AIC</td>
    <td style="border-top: 1px solid black;">209139.95</td>
    <td style="border-top: 1px solid black;">209134.25</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">BIC</td>
    <td style="padding-right: 12px; border: none;">209155.15</td>
    <td style="padding-right: 12px; border: none;">209157.06</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Log Likelihood</td>
    <td style="padding-right: 12px; border: none;">-104567.97</td>
    <td style="padding-right: 12px; border: none;">-104564.13</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Deviance</td>
    <td style="padding-right: 12px; border: none;">1199390707.30</td>
    <td style="padding-right: 12px; border: none;">1198766849.22</td>
  </tr>
  <tr>
    <td style="border-bottom: 2px solid black;">Num. obs.</td>
    <td style="border-bottom: 2px solid black;">14789</td>
    <td style="border-bottom: 2px solid black;">14789</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;" colspan="3"><span style="font-size:0.8em"><sup style="vertical-align: 4px;">***</sup>p &lt; 0.001, <sup style="vertical-align: 4px;">**</sup>p &lt; 0.01, <sup style="vertical-align: 4px;">*</sup>p &lt; 0.05</span></td>
  </tr>
</table>

<table cellspacing="0" style="border: none;">
  <tr>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 1</b></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 2</b></th>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">(Intercept)</td>
    <td style="padding-right: 12px; border: none;">751.11 (2.34)<sup style="vertical-align: 4px;">***</sup></td>
    <td style="padding-right: 12px; border: none;">-1582.14 (13.62)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">MAP</td>
    <td style="padding-right: 12px; border: none;"></td>
    <td style="padding-right: 12px; border: none;">4.38 (0.03)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="border-top: 1px solid black;">AIC</td>
    <td style="border-top: 1px solid black;">209139.95</td>
    <td style="border-top: 1px solid black;">192866.69</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">BIC</td>
    <td style="padding-right: 12px; border: none;">209155.15</td>
    <td style="padding-right: 12px; border: none;">192889.49</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Log Likelihood</td>
    <td style="padding-right: 12px; border: none;">-104567.97</td>
    <td style="padding-right: 12px; border: none;">-96430.34</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Deviance</td>
    <td style="padding-right: 12px; border: none;">1199390707.30</td>
    <td style="padding-right: 12px; border: none;">399043876.90</td>
  </tr>
  <tr>
    <td style="border-bottom: 2px solid black;">Num. obs.</td>
    <td style="border-bottom: 2px solid black;">14789</td>
    <td style="border-bottom: 2px solid black;">14789</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;" colspan="3"><span style="font-size:0.8em"><sup style="vertical-align: 4px;">***</sup>p &lt; 0.001, <sup style="vertical-align: 4px;">**</sup>p &lt; 0.01, <sup style="vertical-align: 4px;">*</sup>p &lt; 0.05</span></td>
  </tr>
</table>

<figure><img src='figure_rmd/FireScars_by_variables.png'  style='display: block'><br><figcaption>Figure 16: Reported impacts of fire.</figcaption></figure><br>

<figure><img src='figure_rmd/FireScars_by_MAP.png'  style='display: block'><br><figcaption>Figure 17: Impact on woody cover and herbaceous material by mean annual precipitation (A,B). Relationship between percent woody cover and reported impact on woody cover and herbaceous material (C,D).</figcaption></figure><br>


<table cellspacing="0" style="border: none;">
  <tr>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 1</b></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 2</b></th>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">(Intercept)</td>
    <td style="padding-right: 12px; border: none;">3.07 (0.11)<sup style="vertical-align: 4px;">***</sup></td>
    <td style="padding-right: 12px; border: none;">3.96 (0.30)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">WoodyCover</td>
    <td style="padding-right: 12px; border: none;"></td>
    <td style="padding-right: 12px; border: none;">-0.03 (0.01)<sup style="vertical-align: 4px;">**</sup></td>
  </tr>
  <tr>
    <td style="border-top: 1px solid black;">AIC</td>
    <td style="border-top: 1px solid black;">208.01</td>
    <td style="border-top: 1px solid black;">200.20</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">BIC</td>
    <td style="padding-right: 12px; border: none;">212.64</td>
    <td style="padding-right: 12px; border: none;">207.15</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Log Likelihood</td>
    <td style="padding-right: 12px; border: none;">-102.00</td>
    <td style="padding-right: 12px; border: none;">-97.10</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Deviance</td>
    <td style="padding-right: 12px; border: none;">66.67</td>
    <td style="padding-right: 12px; border: none;">58.49</td>
  </tr>
  <tr>
    <td style="border-bottom: 2px solid black;">Num. obs.</td>
    <td style="border-bottom: 2px solid black;">75</td>
    <td style="border-bottom: 2px solid black;">75</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;" colspan="3"><span style="font-size:0.8em"><sup style="vertical-align: 4px;">***</sup>p &lt; 0.001, <sup style="vertical-align: 4px;">**</sup>p &lt; 0.01, <sup style="vertical-align: 4px;">*</sup>p &lt; 0.05</span></td>
  </tr>
</table>

<figure><img src='figure_rmd/PreviousMAP.png'  style='display: block'><br><figcaption>Figure 18: Role of previous two years of mean annual precipitation on reported impacts on woody species (A.) and herbaceous cover (B.)</figcaption></figure><br>

[Analyses of Data on Spatial Drivers of Top Kill](http://github.com/danielg7/SpatialDriversOfTopKill/) by Daniel Godwin is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nc-nd/4.0/).
