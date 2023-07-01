#WIO VELCOTIY FOR MAINA JM
rm(list = ls())

mainDir <- "E:/10_JMM_CCVA/ClimateriskWIO"
setwd(mainDir)

lsfolder = "E:/LARGE FILE STORAGE (LFS)/CLIM-CCVA/"

library("ncdf4");library("raster");library("sf")
library("rgeos");library("rasterVis");library("gridExtra")
library("doParallel");library("foreach");library("scales");library("data.table")
library("repmis");library("rgdal");library("xts");library("geosphere")

#Import rasters
wio_shp <- st_read("E:/10_JMM_CCVA/ClimateriskWIO/2_Data/shp/country_shape.shp") %>% st_transform(crs = "+proj=longlat")

sst <- c("sst126_ESM25km.nc", "sst245_ESM25km.nc", "sst370_ESM25km.nc", "sst585_ESM25km.nc")
ncfiles <- c(sst)
# loop over the file names, convert netCDFs to raster stack and crop to WIO extent
r <- lapply(1:length(ncfiles), function(x){
  r_rast <-  raster::crop(
    raster::stack(paste0(lsfolder, ncfiles[[x]]), varname = "tos"),
    wio_shp)
  return(r_rast)
})
names(r) <- c("sst126","sst245","sst370","sst585")

#Functions to calculate gradient-based velocity of climate change
# after Garcia Molinos et al.2020 (https://doi.org/10.1111/2041-210X.13295)
source(paste(mainDir, "/1_Codes/v1_spatGrad.R", sep=""))
source(paste(mainDir, "/1_Codes/v2_tempTrend.R", sep=""))
source(paste(mainDir, "/1_Codes/v3_sumSeries.R", sep=""))
source(paste(mainDir, "/1_Codes/v4_angulo.R", sep=""))

#the VoCC function
gVel <- function(r){
  # temporal trend
  vt <- tempTrend(r, th = 10)
  
  # spatial gradient
  rr_mean = calc(r, mean)
  vg <- spatGrad(rr_mean, th = 0.0001, projected = FALSE)
  
  # climate velocity
  VoCC <- vt[[1]]/vg[[1]]
    # velocity angles have opposite direction to the spatial climatic gradient if warming and same direction (cold to warm) if cooling
    ind <- which(getValues(VoCC) > 0)
    VoCCang <- vg[[2]]
    VoCCang[ind] <- vg[[2]][ind] + 180
    VoCCang[] <- ifelse(VoCCang[] >= 360, VoCCang[] - 360, VoCCang[])
    output <- stack(VoCC,VoCCang)
    names(output) <- c("voccMag", "voccAng")
    return(output)
}

#Few parameters to subset required time periods (this works for now but needs improvement!!!!).
nYrs = 30 #30yr climatological period
N = nlayers(r$sst585) #total number of layers in a climate data
s1 = 432-(nYrs*12 - 1) #S1 where to begin the first sub-setting. 2015 to mid-century equals 432 months
s2 = N-(nYrs*12 - 1)

#break files into their respective periods, near future and far future
rr <- lapply(1:length(r), FUN = function(x){
  rr <- list(r[[x]][[s1:432]],r[[x]][[s2:N]])
  names(rr) <- c("near_future", "far_future")
  return(rr)
  })

#Run Velocity for each Scenario and both periods
rr <- unlist(rr)
(gvWIO <- lapply(1:length(rr), FUN = function(x){
  vel <- gVel(rr[[x]])[[1]] 
  return(vel)
  }))
gvWIO <- stack(gvWIO)
names(gvWIO) <- c("vel_ssp126_2050","vel_ssp126_2100","vel_ssp245_2050","vel_ssp245_2100","vel_ssp370_2050","vel_ssp370_2100","vel_ssp585_2050","vel_ssp585_2100")
plot(gvWIO)

#write to drive as NC file
if (require(ncdf4)) {	rnc <- writeRaster(gvWIO, filename=file.path("2_Data/raster/velocity/wio_climate_velocity.nc"), format="CDF", overwrite=TRUE) }

#Visualize
my.at <- seq(0, 100, by = 1)
rasterVis::levelplot(gvWIO[[8]], par.settings = BuRdTheme, at=my.at, main = 'Gradient-based vocc', margin = FALSE)
