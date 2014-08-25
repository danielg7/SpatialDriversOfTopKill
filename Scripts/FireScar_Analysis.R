library(rgdal)
library(maptools)
library(lubridate)

Aggregated_FireScars <- readOGR(dsn="Data/BurnScars/", layer="agg2")
FIRMS <- readOGR(dsn="/Users/danielg7/Documents/FIRMS/",layer="firms138621406915198_MCD14ML")

Aggregated_FireScars$STARTDATE <- ymd(as.character(Aggregated_FireScars$STARTDATE))
Aggregated_FireScars$DATE_START <- ymd(as.character(Aggregated_FireScars$DATE_START))
Aggregated_FireScars$DATE_END <- ymd(as.character(Aggregated_FireScars$DATE_END))

FIRMS$ACQ_DATE <- ymd(as.character(FIRMS$ACQ_DATE))

MergedFires <- merge(na.omit(Aggregated_FireScars),FIRMS,by.x = Aggregated_FireScars$DATE_START,by.y=FIRMS$ACQ_DATE)
