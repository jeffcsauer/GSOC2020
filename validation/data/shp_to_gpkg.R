library(spdep)
library(rgdal)
temp <- st_read("C:/Users/jeffe/Dropbox/GSOC2020/validation/data/baltimore/baltimore_housing.shp")
outname <- "C:/Users/jeffe/Dropbox/GSOC2020/validation/data/commpop.gpkg"
temp_sp <- as_Spatial(temp)
writeOGR(temp_sp, dsn = outname, layer = "commpop", driver = "GPKG")
