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
## > DrySeasonMAP <- extract(krugerBrick, FIRMS_Kruger_DrySeason, 
## +     method = "bilinear", df = TRUE, sp = TRUE)
## 
## > WetSeasonMAP <- extract(krugerBrick, FIRMS_Kruger_WetSeason, 
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
## > krugerFRIRaster <- rasterize(krugerFRI_UTM, krugerFRIRaster, 
## +     field = krugerFRI_UTM$MFRI)
## Found 29349 region(s) and 29487 polygon(s)
## 
## > krugerBrick <- brick(krugerMAP_UTM, krugerFRIRaster)
## 
## > names(krugerBrick) <- c("MAP", "MFRI")
## 
## > krugerFuelLoad <- 382.9 + 3.3 * krugerBrick$MAP + 
## +     979.4 * krugerBrick$MFRI - 0.001 * krugerBrick$MFRI^2 + 0.37 * 
## +     krugerBrick$MAP * kru .... [TRUNCATED] 
## 
## > RelativeHumidity <- 36.6
## 
## > FuelMoisture <- 32.1
## 
## > WindSpeed <- 2.6
## 
## > krugerFirelineIntensity <- 2729 + 0.8684 * krugerFuelLoad - 
## +     530 * sqrt(FuelMoisture) - 0.1907 * RelativeHumidity^2 - 
## +     5961/WindSpeed
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
## +     "FirelineIntensity", "Geology", "MAP")
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
## > extractedFireScars <- extract(x = krugerFireScarAnalysis_brick, 
## +     y = FireScars_WoodyImpact, fun = mean, sp = TRUE)
## 
## > extractedFireScars_df <- as.data.frame(extractedFireScars)
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
## > rm(aggPlied)
## 
## > rm(testPlied)
## 
## > rm(extractedFireScars)
## 
## > rm(extractedFireScars_df)
## 
## > rm(extractedFireScars_wx_df)
## 
## > rm(FireScars_WoodyImpact)
## 
## > rm(krugerFireScarAnalysis_brick)
```
Run the intensity model analysis


**Figures**

<figure><img src='figure_rmd/logicalFramework.png'  style='display: block'><br><figcaption>Figure 3: </figcaption></figure><br>
A. Average woody cover by mean annual precipitation at 50 mm/yr intervals (points); line is significant at P < 0.01, R^2 = 0.5862. B. Relationship between fuel load accumulation (grass biomass) and mean annual precipitation, derived from Figure 1, Govender et al. 2006. C. Relationship between fireline intensity and fuel load given an average rate of spread of .04 m/s. D. Relationship between average fireline intensity and mean annual precipitation, from Figure 2, Govender et al. 2006. E. Probability of top kill as a function of fire intensity by mean annual precipitation, as a function of Higgins et al. 2011 and the relationship shown in D.


<figure><img src='figure_rmd/FRP_by_MAP_Season.png'  style='display: block'><br><figcaption>Figure 4: Fire radiative power by mean annual precipitation, subdivided by (A.) season of burn and (B.) geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FRP_by_WoodyCover.png'  style='display: block'><br><figcaption>Figure 5: Fire radiative power by percent woody cover, subdivided by (A.) season of burn and (B.) geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FRP_by_GLY.png'  style='display: block'><br><figcaption>Figure 6: Fire radiative power by geology and season of burn.</figcaption></figure><br>

<figure><img src='figure_rmd/FRP_avg_by_MAP.png'  style='display: block'><br><figcaption>Figure 7: Fire radiative power, averaged by 50 mm / yr MAP</figcaption></figure><br>
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



<figure><img src='figure_rmd/ModeledFirelineIntensity.png'  style='display: block'><br><figcaption>Figure 8: Modeled fireline intensity across Kruger National Park</figcaption></figure><br>

<figure><img src='figure_rmd/FLI_By_WoodyCover.png'  style='display: block'><br><figcaption>Figure 9: Modeled fireline intensity by percent woody cover and geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FLI_By_MAP.png'  style='display: block'><br><figcaption>Figure 10: Modeled fireline intensity by MAP and geologic parent material.</figcaption></figure><br>

<figure><img src='figure_rmd/FLI_By_NavashniStandards.png'  style='display: block'><br><figcaption>Figure 11: Percentage of probable fires by fireline intensity classes, as defined by Govender et al. 2006</figcaption></figure><br>

**Intensity by Woody Cover and MAP: GLM Results**

<table cellspacing="0" style="border: none;">
  <tr>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 1</b></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 2</b></th>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">(Intercept)</td>
    <td style="padding-right: 12px; border: none;">902.88 (3.15)<sup style="vertical-align: 4px;">***</sup></td>
    <td style="padding-right: 12px; border: none;">902.69 (3.15)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">WoodyCover</td>
    <td style="padding-right: 12px; border: none;"></td>
    <td style="padding-right: 12px; border: none;">0.00 (0.00)<sup style="vertical-align: 4px;">**</sup></td>
  </tr>
  <tr>
    <td style="border-top: 1px solid black;">AIC</td>
    <td style="border-top: 1px solid black;">217862.94</td>
    <td style="border-top: 1px solid black;">217854.87</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">BIC</td>
    <td style="padding-right: 12px; border: none;">217878.14</td>
    <td style="padding-right: 12px; border: none;">217877.68</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Log Likelihood</td>
    <td style="padding-right: 12px; border: none;">-108929.47</td>
    <td style="padding-right: 12px; border: none;">-108924.44</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Deviance</td>
    <td style="padding-right: 12px; border: none;">2165327697.19</td>
    <td style="padding-right: 12px; border: none;">2163854972.63</td>
  </tr>
  <tr>
    <td style="border-bottom: 2px solid black;">Num. obs.</td>
    <td style="border-bottom: 2px solid black;">14788</td>
    <td style="border-bottom: 2px solid black;">14788</td>
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
    <td style="padding-right: 12px; border: none;">902.88 (3.15)<sup style="vertical-align: 4px;">***</sup></td>
    <td style="padding-right: 12px; border: none;">-1427.79 (8.76)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">MAP</td>
    <td style="padding-right: 12px; border: none;"></td>
    <td style="padding-right: 12px; border: none;">4.10 (0.02)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="border-top: 1px solid black;">AIC</td>
    <td style="border-top: 1px solid black;">217862.94</td>
    <td style="border-top: 1px solid black;">191637.42</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">BIC</td>
    <td style="padding-right: 12px; border: none;">217878.14</td>
    <td style="padding-right: 12px; border: none;">191660.23</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Log Likelihood</td>
    <td style="padding-right: 12px; border: none;">-108929.47</td>
    <td style="padding-right: 12px; border: none;">-95815.71</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Deviance</td>
    <td style="padding-right: 12px; border: none;">2165327697.19</td>
    <td style="padding-right: 12px; border: none;">367513431.41</td>
  </tr>
  <tr>
    <td style="border-bottom: 2px solid black;">Num. obs.</td>
    <td style="border-bottom: 2px solid black;">14788</td>
    <td style="border-bottom: 2px solid black;">14788</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;" colspan="3"><span style="font-size:0.8em"><sup style="vertical-align: 4px;">***</sup>p &lt; 0.001, <sup style="vertical-align: 4px;">**</sup>p &lt; 0.01, <sup style="vertical-align: 4px;">*</sup>p &lt; 0.05</span></td>
  </tr>
</table>

<figure><img src='figure_rmd/FireScars_by_variables.png'  style='display: block'><br><figcaption>Figure 12: Reported impacts of fire.</figcaption></figure><br>

<figure><img src='figure_rmd/FireScars_by_MAP.png'  style='display: block'><br><figcaption>Figure 13: Impact on woody cover and herbaceous material by mean annual precipitation (A,B). Relationship between percent woody cover and reported impact on woody cover and herbaceous material (C,D).</figcaption></figure><br>


<table cellspacing="0" style="border: none;">
  <tr>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 1</b></th>
    <th style="text-align: left; border-top: 2px solid black; border-bottom: 1px solid black; padding-right: 12px;"><b>Model 2</b></th>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">(Intercept)</td>
    <td style="padding-right: 12px; border: none;">3.05 (0.11)<sup style="vertical-align: 4px;">***</sup></td>
    <td style="padding-right: 12px; border: none;">3.97 (0.29)<sup style="vertical-align: 4px;">***</sup></td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">WoodyCover</td>
    <td style="padding-right: 12px; border: none;"></td>
    <td style="padding-right: 12px; border: none;">-0.03 (0.01)<sup style="vertical-align: 4px;">**</sup></td>
  </tr>
  <tr>
    <td style="border-top: 1px solid black;">AIC</td>
    <td style="border-top: 1px solid black;">198.98</td>
    <td style="border-top: 1px solid black;">189.90</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">BIC</td>
    <td style="padding-right: 12px; border: none;">203.56</td>
    <td style="padding-right: 12px; border: none;">196.77</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Log Likelihood</td>
    <td style="padding-right: 12px; border: none;">-97.49</td>
    <td style="padding-right: 12px; border: none;">-91.95</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;">Deviance</td>
    <td style="padding-right: 12px; border: none;">61.78</td>
    <td style="padding-right: 12px; border: none;">53.08</td>
  </tr>
  <tr>
    <td style="border-bottom: 2px solid black;">Num. obs.</td>
    <td style="border-bottom: 2px solid black;">73</td>
    <td style="border-bottom: 2px solid black;">73</td>
  </tr>
  <tr>
    <td style="padding-right: 12px; border: none;" colspan="3"><span style="font-size:0.8em"><sup style="vertical-align: 4px;">***</sup>p &lt; 0.001, <sup style="vertical-align: 4px;">**</sup>p &lt; 0.01, <sup style="vertical-align: 4px;">*</sup>p &lt; 0.05</span></td>
  </tr>
</table>

<figure><img src='figure_rmd/PreviousMAP.png'  style='display: block'><br><figcaption>Figure 14: Role of previous two years of mean annual precipitation on reported impacts on woody species (A.) and herbaceous cover (B.)</figcaption></figure><br>
![Creative Commons License](http://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png)

[Analyses of Data on Spatial Drivers of Top Kill](http://github.com/danielg7/SpatialDriversOfTopKill/) by Daniel Godwin is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nc-nd/4.0/).
