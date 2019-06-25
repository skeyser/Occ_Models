#NOAA Precipitation Data 
#test 
pacman::p_load("lubridate", "raster", "rasterVis", "maps", "maptools", "rgdal", "ncdf4", "ncdf4.helpers", "PCICt", "here", "tidyverse")
noaa_prcp_monthly_lt <- nc_open(here("Data_Envi/Precipitation/precip.V1.0.mon.1981-2010.ltm.nc"))
noaa_prcp <- nc_open(here("Data_Envi/Precipitation/precip.V1.0.mon.mean.nc"))

lat <- ncdf4::ncvar_get(noaa_prcp, varid = "lat")
lon <- ncdf4::ncvar_get(noaa_prcp, varid = "lon")
prcp <- ncdf4::ncvar_get(noaa_prcp, varid = "precip")
prcp <- prcp[complete.cases(prcp)]

#NOAA & UoD Temperature Data 
#noaa_temp <- nc_open(here("Data_Envi/Air Temperature/air.mon.mean.nc"))
#noaa_temp_monthly <- nc_open("Data_Envi/Air Temperature/air.mon.1981-2010.ltm.nc")

uod_temp <- nc_open(here("Data_Envi/Air Temperature/uod.air.mon.mean.v401.nc"))
uod_ltm_temp <- nc_open(here("Data_Envi/Air Temperature/uod.air.mon.ltm.v401.nc"))

lat.temp <- ncdf4::ncvar_get(uod_temp, varid = "lat")
lon.temp <- ncdf4::ncvar_get(uod_temp, varid = "lon")

temp <- ncdf4::ncvar_get(uod_temp, varid = "air")
temp <- temp[complete.cases(temp)]

#Precipitation time is being stored as hours since 1700-1-1
time.precip <- ncdf4::ncvar_get(noaa_prcp, varid = "time")

#Need to change time to become readable, uses the lubridate package 
origin <- as.Date("1700-1-1 00:00:0.0", format = "%Y-%m-%d %H:%M:%S")
time.precip <- origin + hours(time.precip)
time.h.precip <- as.Date(time.precip, origin = as.Date("1700-01-01"))

#Time Conversion Temperature UoD
time.temp <- ncdf4::ncvar_get(uod_temp, varid = "time")
time.temp <- as.numeric(time.temp)
origin.temp <- as.Date("1900-1-1 00:00:0.0", format = "%Y-%m-%d %H:%M:%S")
time.temp <- origin_temp + hours(time.temp)
time.h.temp <- as.Date(time.temp, origin_temp = as.Date("1900-01-01"))

#Extract Temperatures for GoM
rts <- read.csv(here("Data_BBS/States_GoM/bbsrts66-98_xy.csv"), header =T)

#Boundaries for GoM rts
#Long -98 to -81, Lat 24 to 31 
lonlat.temp <- as.matrix(expand.grid(lon.temp, lat.temp))
tmp_vec <- as.vector(temp)
temp_df <- data.frame(cbind(lonlat.temp, tmp_vec))
names(temp_df) <- c("lon", "lat", "Deg_C")
temp_df <- temp_df %>% mutate(Deg_C = Deg_K - 273.15)
temp_df <- temp_df[temp_df$lon > 82 & temp_df$lon < 96,]
temp_df <- temp_df[temp_df$lat > 24 & temp_df$lat < 31,]

lonlat.df <- lonlat.temp
lonlat.df <- as.data.frame(lonlat.df)
names(lonlat.df) <- c("lon", "lat")
lonlat.df <- lonlat.df[lonlat.df$lon > 82 & lonlat.df$lon < 96,]
lonlat.df <- lonlat.df[lonlat.df$lat > 24 & lonlat.df$lat < 31,]
ext <- extract(temp, lonlat.df, method = "bilinear")

#Generates of a map of the ncdf data 
map.precip <- raster(here("Data_Envi/Precipitation/precip.V1.0.mon.mean.nc"), header = T)
us.precip <- map("usa", plot = F)
us.sp.precip <- map2SpatialLines(us.precip, proj4string = CRS("+proj=longlat"))
mapTheme <- rasterTheme(states = rev(brewer.pal(10, "RdBu")))
plt.precip <- levelplot(map.precip, margin = F, pretty = T, pare.settings = mapTheme, main = "map")
plt.precip + layer(sp.lines(us.sp.precip, col = "black", lwd = 0.5))

#Map Temp
map.temp <- raster(here("Data_Envi/Air Temperature/air.mon.mean.nc"), header = T)
us.temp <- map("usa", plot = F)
us.sp.temp <- map2SpatialLines(us.temp, proj4string = CRS("+proj=longlat"))
mapTheme <- rasterTheme(states = rev(brewer.pal(10, "RdBu")))
plt.temp <- levelplot(map.temp, margin = F, pretty = T, pare.settings = mapTheme, main = "map")
plt.temp + layer(sp.lines(us.sp.temp, col = "black", lwd = 0.5))



#Extracting Data using Raster and NetCDF 
uod_temp <- here("Data_Envi/Air Temperature/uod.air.mon.mean.v401.nc")

uod_brick <- brick(uod_temp)

plot(uod_brick)

plot(uod_brick[[1]])

ROI <- extent(262,282, 22, 42)

r.crop <- crop(uod_brick, ROI)
plot(r.crop[[1]])

lon.pts <- seq(262, 282, by = 0.5)
lat.pts <- seq(22, 42, by = 0.5)
plot(r.crop)
points(lon.pts, lat.pts, pch = 4, col = "red")

extract.temp <- cbind(lon.pts, lat.pts)
ext <- extract(r.crop, extract.temp, method = "bilinear")
#.rs.unloadPackage("tidyr")
ext.df <- as.data.frame(ext)
