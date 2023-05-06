library(ncdf4);library(sp);library(tidyverse);library(sf);library(raster)

rm(list = ls())
lfs.dir = "E:/Datasets/"
#db.dir = "C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/Data/"

#import some important functions
source("1_Codes/2_DataNorms.R")

#Import all data
wio.AOO <- readRDS("2_Data/spreadsheet/2_Ecosystems/wioAOO.wTEV.rds")
wio.AOO.spdf <- st_as_sf(wio.AOO, coords=c('x', 'y'), crs="+proj=longlat")
focusISO3 <- st_read("2_Data/shp/country_shape.shp") %>% st_as_sf() %>% st_transform(crs = "+proj=longlat")

#SEA SURFACE TEMPERATURE
rst.Hist <- stack(paste0(lfs.dir, "CLIM-CCVA/sst_historical_ESM25km.nc")); rst.Hist<- crop(rst.Hist,focusISO3)
rst.sst126 <- stack(paste0(lfs.dir, "CLIM-CCVA/sst126_ESM25km.nc")); rst.sst126<- crop(rst.sst126,focusISO3)
rst.sst245 <- stack(paste0(lfs.dir, "CLIM-CCVA/sst245_ESM25km.nc")); rst.sst245<- crop(rst.sst245,focusISO3)
rst.sst370 <- stack(paste0(lfs.dir, "CLIM-CCVA/sst370_ESM25km.nc")); rst.sst370<- crop(rst.sst370,focusISO3)
rst.sst585 <- stack(paste0(lfs.dir, "CLIM-CCVA/sst585_ESM25km.nc")); rst.sst585<- crop(rst.sst585,focusISO3)

#Plot time series of SST
m1 <- list(rst.sst126, rst.sst245, rst.sst370, rst.sst585)
names(m1) <- c("SSP126", "SSP245", "SSP370", "SSP585")
rr <- lapply(1:length(m1), function(x){raster::cellStats(m1[[x]], stat='mean', na.rm=TRUE)})
dfx <- c()
for (k in 1:length(rr)){
  dfx[[k]] <- rr[[k]] %>% as.data.frame() %>% rownames_to_column()
}

names(dfx) <- c("SSP126", "SSP245", "SSP370", "SSP585")
dfx <- bind_rows(dfx, .id = "column_label");colnames(dfx) <- c("Scenarios", "time", "value")
#dfx <- dfx %>% separate_wider_delim(Mod_Sce_Var, "_", names = c("Model", "Scenario", "Variable")) %>% as.data.frame()
dfx$Date <- format(as.Date(dfx$time, format = "X%Y.%m.%d"), format = "%Y")
ensemble_tas <- dfx %>% group_by(Scenarios,Date) %>% summarise(value = median(value, na.rm=TRUE)) %>% ungroup() %>% data.frame()
#write_csv(ensemble_tas, "outputs/tas_outxx.csv")

#Plot time series for Temperature or Precipitation
(plt2 <- ensemble_tas %>% 
    ggplot(aes(x=as.numeric(Date), y=value, colour=Scenarios)) +
    geom_smooth(method = "loess", se=F)+
    theme_classic(base_size = 10)+
    scale_y_continuous(name = "SST")+scale_colour_viridis_d()+
    theme(axis.text.x = element_text(angle = 0)) + 
    scale_x_continuous(name = "Calendar year", breaks = seq(2015, 2100, 10), expand = c(0, 0)))
ggsave("outputs/SST_series.png", height = 3, width = 5)

#Frequency of extremes
sstx90p.ssp126.2050 <- txpFUN(r=rst.sst126, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90p.ssp126.2050)<-"sstx90p.ssp126.2050"
sstx90p.ssp126.2100 <- txpFUN(r=rst.sst126, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90p.ssp126.2100)<-"sstx90p.ssp126.2100"
sstx90p.ssp245.2050 <- txpFUN(r=rst.sst245, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90p.ssp245.2050)<-"sstx90p.ssp245.2050"
sstx90p.ssp245.2100 <- txpFUN(r=rst.sst245, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90p.ssp245.2100)<-"sstx90p.ssp245.2100"
sstx90p.ssp370.2050 <- txpFUN(r=rst.sst370, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90p.ssp370.2050)<-"sstx90p.ssp370.2050"
sstx90p.ssp370.2100 <- txpFUN(r=rst.sst370, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90p.ssp370.2100)<-"sstx90p.ssp370.2100"
sstx90p.ssp585.2050 <- txpFUN(r=rst.sst585, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90p.ssp585.2050)<-"sstx90p.ssp585.2050"
sstx90p.ssp585.2100 <- txpFUN(r=rst.sst585, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90p.ssp585.2100)<-"sstx90p.ssp585.2100"

#Intensity of extremes
sstx90Int.ssp126.2050 <- txIntFUN(r=rst.sst126, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90Int.ssp126.2050)<-"sstx90Int.ssp126.2050"
sstx90Int.ssp126.2100 <- txIntFUN(r=rst.sst126, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90Int.ssp126.2100)<-"sstx90Int.ssp126.2100"
sstx90Int.ssp245.2050 <- txIntFUN(r=rst.sst245, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90Int.ssp245.2050)<-"sstx90Int.ssp245.2050"
sstx90Int.ssp245.2100 <- txIntFUN(r=rst.sst245, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90Int.ssp245.2100)<-"sstx90Int.ssp245.2100"
sstx90Int.ssp370.2050 <- txIntFUN(r=rst.sst370, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90Int.ssp370.2050)<-"sstx90Int.ssp370.2050"
sstx90Int.ssp370.2100 <- txIntFUN(r=rst.sst370, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90Int.ssp370.2100)<-"sstx90Int.ssp370.2100"
sstx90Int.ssp585.2050 <- txIntFUN(r=rst.sst585, nYrs=30, baseline=rst.Hist, p=0.90, pn="near"); names(sstx90Int.ssp585.2050)<-"sstx90Int.ssp585.2050"
sstx90Int.ssp585.2100 <- txIntFUN(r=rst.sst585, nYrs=30, baseline=rst.Hist, p=0.90, pn="far"); names(sstx90Int.ssp585.2100)<-"sstx90Int.ssp585.2100"

#Trend
sst.trend.ssp126.2050 <- slpFUN(rst.sst126[[73:432]]);names(sst.trend.ssp126.2050)<-"sst.trend.ssp126.2050"
sst.trend.ssp126.2100 <- slpFUN(rst.sst126[[673:1032]]);names(sst.trend.ssp126.2100)<-"sst.trend.ssp126.2100"
sst.trend.ssp245.2050 <- slpFUN(rst.sst245[[73:432]]);names(sst.trend.ssp245.2050)<-"sst.trend.ssp245.2050"
sst.trend.ssp245.2100 <- slpFUN(rst.sst245[[673:1032]]);names(sst.trend.ssp245.2100)<-"sst.trend.ssp245.2100"
sst.trend.ssp370.2050 <- slpFUN(rst.sst370[[73:432]]);names(sst.trend.ssp370.2050)<-"sst.trend.ssp370.2050"
sst.trend.ssp370.2100 <- slpFUN(rst.sst370[[673:1032]]);names(sst.trend.ssp370.2100)<-"sst.trend.ssp370.2100"
sst.trend.ssp585.2050 <- slpFUN(rst.sst585[[73:432]]);names(sst.trend.ssp585.2050)<-"sst.trend.ssp585.2050"
sst.trend.ssp585.2100 <- slpFUN(rst.sst585[[673:1032]]);names(sst.trend.ssp585.2100)<-"sst.trend.ssp585.2100"

hzd.SST <- stack(inormal(sstx90p.ssp126.2050),
                 inormal(sstx90p.ssp126.2100),
                 inormal(sstx90p.ssp245.2050),
                 inormal(sstx90p.ssp245.2100),
                 inormal(sstx90p.ssp370.2050),
                 inormal(sstx90p.ssp370.2100),
                 inormal(sstx90p.ssp585.2050),
                 inormal(sstx90p.ssp585.2100),
  
                 inormal(sstx90Int.ssp126.2050),
                 inormal(sstx90Int.ssp126.2100),
                 inormal(sstx90Int.ssp245.2050),
                 inormal(sstx90Int.ssp245.2100),
                 inormal(sstx90Int.ssp370.2050),
                 inormal(sstx90Int.ssp370.2100),
                 inormal(sstx90Int.ssp585.2050),
                 inormal(sstx90Int.ssp585.2100),
                 
                 inormal(sst.trend.ssp126.2050),
                 inormal(sst.trend.ssp126.2100),
                 inormal(sst.trend.ssp245.2050),
                 inormal(sst.trend.ssp245.2100),
                 inormal(sst.trend.ssp370.2050),
                 inormal(sst.trend.ssp370.2100),
                 inormal(sst.trend.ssp585.2050),
                 inormal(sst.trend.ssp585.2100))
plot(hzd.SST)
summary(hzd.SST)
if (require(ncdf4)) {
  rnc <- raster::writeRaster(stack(sstx90p.ssp126.2050, sstx90p.ssp126.2100,
                                   sstx90p.ssp245.2050, sstx90p.ssp245.2100,
                                   sstx90p.ssp370.2050, sstx90p.ssp370.2100, 
                                   sstx90p.ssp585.2050, sstx90p.ssp585.2100), 
                             filename=file.path("2_Data/raster/wio_esm_sstx90p.nc"), format="CDF", overwrite=TRUE)
}

if (require(ncdf4)) {
  rnc <- raster::writeRaster(stack(sstx90Int.ssp126.2050, sstx90Int.ssp126.2100,
                                   sstx90Int.ssp245.2050, sstx90Int.ssp245.2100,
                                   sstx90Int.ssp370.2050, sstx90Int.ssp370.2100,
                                   sstx90Int.ssp585.2050, sstx90Int.ssp585.2100), 
                             filename=file.path("2_Data/raster/wio_esm_sstx90Int.nc"), format="CDF", overwrite=TRUE)
}

if (require(ncdf4)) {
  rnc <- raster::writeRaster(stack(sst.trend.ssp126.2050,sst.trend.ssp126.2100,
                                   sst.trend.ssp245.2050,sst.trend.ssp245.2100,
                                   sst.trend.ssp370.2050,sst.trend.ssp370.2100,
                                   sst.trend.ssp585.2050,sst.trend.ssp585.2100), 
                             filename=file.path("2_Data/raster/wio_esm_sst_trend.nc"), format="CDF", overwrite=TRUE)
}
#ACIDIFICATION
pH126 <- stack(paste0(lfs.dir, "CLIM-CCVA/ph126_ESM25km.nc")); pH126 <- crop(pH126,focusISO3)
pH245 <- stack(paste0(lfs.dir, "CLIM-CCVA/ph245_ESM25km.nc")); pH245 <- crop(pH245,focusISO3) ##Errors in HDF
pH370 <- stack(paste0(lfs.dir, "CLIM-CCVA/ph370_ESM25km.nc")); pH370 <- crop(pH370,focusISO3)
pH585 <- stack(paste0(lfs.dir, "CLIM-CCVA/ph585_ESM25km.nc")); pH585 <- crop(pH585,focusISO3)

m1 <- list(pH126, pH245, pH370, pH585)
names(m1) <- c("HIST", "SSP126", "SSP245", "SSP370", "SSP585")
rr <- lapply(1:length(m1), function(x){raster::cellStats(m1[[x]], stat='mean', na.rm=TRUE)})
dfx <- c()
for (k in 1:length(rr)){
  dfx[[k]] <- rr[[k]] %>% as.data.frame() %>% rownames_to_column()
}

names(dfx) <- c("SSP126", "SSP245", "SSP370", "SSP585")
dfx <- bind_rows(dfx, .id = "column_label");colnames(dfx) <- c("Mod_Sce_Var", "time", "value")
#dfx <- dfx %>% separate_wider_delim(Mod_Sce_Var, "_", names = c("Model", "Scenario", "Variable")) %>% as.data.frame()
dfx$Date <- format(as.Date(dfx$time, format = "X%Y.%m.%d"), format = "%Y")
ensemble_tas <- dfx %>% group_by(Mod_Sce_Var,Date) %>% summarise(value = median(value, na.rm=TRUE)) %>% ungroup() %>% data.frame()
#write_csv(ensemble_tas, "outputs/tas_outxx.csv")
#Plot time series
(plt2 <- ensemble_tas %>% 
    ggplot(aes(x=as.numeric(Date), y=value, colour=Scenarios)) +
    geom_smooth(method = "loess", se=F)+
    theme_classic(base_size = 10)+
    scale_y_continuous(name = "pH")+scale_colour_viridis_d()+
    theme(axis.text.x = element_text(angle = 0)) + 
    scale_x_continuous(name = "Calendar year", breaks = seq(2015, 2100, 10), expand = c(0, 0)))
ggsave("outputs/ph_series.png", height = 3, width = 5)

pH.ssp126.2050 <- slpFUN(pH126[[73:432]]) %>% mask(., sst.trend.ssp370.2050) ;names(pH.ssp126.2050)<-"pH.ssp126.2050"
pH.ssp126.2100 <- slpFUN(pH126[[673:1032]]) %>% mask(., sst.trend.ssp370.2050) ;names(pH.ssp126.2100)<-"pH.ssp126.2100"
#Remember to remove 126 from the 245 
pH.ssp245.2050 <- slpFUN(pH126[[73:432]]) %>% mask(., sst.trend.ssp370.2050) ;names(pH.ssp245.2050)<-"pH.ssp245.2050"
pH.ssp245.2100 <- slpFUN(pH126[[673:1032]]) %>% mask(., sst.trend.ssp370.2050) ;names(pH.ssp245.2100)<-"pH.ssp245.2100"
pH.ssp370.2050 <- slpFUN(pH370[[73:432]]) %>% mask(., sst.trend.ssp370.2050);names(pH.ssp370.2050)<-"pH.ssp370.2050"
pH.ssp370.2100 <- slpFUN(pH370[[673:1032]]) %>% mask(., sst.trend.ssp370.2050) ;names(pH.ssp370.2100)<-"pH.ssp370.2100"
pH.ssp585.2050 <- slpFUN(pH585[[73:432]]) %>% mask(., sst.trend.ssp370.2050);names(pH.ssp585.2050)<-"pH.ssp585.2050"
pH.ssp585.2100 <- slpFUN(pH585[[673:1032]]) %>% mask(., sst.trend.ssp370.2050) ;names(pH.ssp585.2100)<-"pH.ssp585.2100"
hzd.pH <- stack(inormal(pH.ssp126.2050),
                inormal(pH.ssp126.2100),
                inormal(pH.ssp245.2050),
                inormal(pH.ssp245.2100),
                inormal(pH.ssp370.2050),
                inormal(pH.ssp370.2100),
                inormal(pH.ssp585.2050),
                inormal(pH.ssp585.2100))
plot(hzd.pH)
# if (require(ncdf4)) {
#   rnc <- raster::writeRaster(stack(pH.ssp126.2050,pH.ssp126.2100,
#                                    pH.ssp245.2050,pH.ssp245.2100,
#                                    pH.ssp370.2050,pH.ssp370.2100,
#                                    pH.ssp585.2050,pH.ssp585.2100), 
#                              filename=file.path("2_Data/raster/wio_esm_pH_trend.nc"), format="CDF", overwrite=TRUE)
# }


#LAND SURFACE TEMPERATURE
lst.Hist <- stack(paste0(lfs.dir, "CLIM-CCVA/ts_historical_ESM25km.nc")); lst.Hist<- crop(lst.Hist,focusISO3)
lst.sst126 <- stack(paste0(lfs.dir, "CLIM-CCVA/ts126_ESM25km.nc")); lst.sst126<- crop(lst.sst126,focusISO3)
lst.sst245 <- stack(paste0(lfs.dir, "CLIM-CCVA/ts245_ESM25km.nc")); lst.sst245<- crop(lst.sst245,focusISO3)
lst.sst370 <- stack(paste0(lfs.dir, "CLIM-CCVA/ts370_ESM25km.nc")); lst.sst370<- crop(lst.sst370,focusISO3)
lst.sst585 <- stack(paste0(lfs.dir, "CLIM-CCVA/ts585_ESM25km.nc")); lst.sst585<- crop(lst.sst585,focusISO3)

#Plot timeseries of SST
m1 <- list(lst.sst126, lst.sst245, lst.sst370, lst.sst585)
names(m1) <- c("SSP126", "SSP245", "SSP370", "SSP585")
rr <- lapply(1:length(m1), function(x){raster::cellStats(m1[[x]], stat='mean', na.rm=TRUE)})
dfx <- c()
for (k in 1:length(rr)){
  dfx[[k]] <- rr[[k]] %>% as.data.frame() %>% rownames_to_column()
}

names(dfx) <- c("SSP126", "SSP245", "SSP370", "SSP585")
dfx <- bind_rows(dfx, .id = "column_label");colnames(dfx) <- c("Scenarios", "time", "value")
#dfx <- dfx %>% separate_wider_delim(Mod_Sce_Var, "_", names = c("Model", "Scenario", "Variable")) %>% as.data.frame()
dfx$Date <- format(as.Date(dfx$time, format = "X%Y.%m.%d"), format = "%Y")
ensemble_tas <- dfx %>% group_by(Scenarios,Date) %>% summarise(value = median(value, na.rm=TRUE)) %>% ungroup() %>% data.frame()
#write_csv(ensemble_tas, "outputs/tas_outxx.csv")

#Plot time series
(plt2 <- ensemble_tas %>% 
    ggplot(aes(x=as.numeric(Date), y=value, colour=Scenarios)) +
    geom_smooth(method = "loess", se=F)+
    theme_classic(base_size = 10)+
    scale_y_continuous(name = "LST")+scale_colour_viridis_d()+
    theme(axis.text.x = element_text(angle = 0)) + 
    scale_x_continuous(name = "Calendar year", breaks = seq(2015, 2100, 10), expand = c(0, 0)))
ggsave("outputs/LST_series.png", height = 3, width = 5)

#Frequency of Extremes
lst90p.ssp126.2050 <- txpFUN(r=lst.sst126, nYrs=30, baseline=lst.Hist, p=0.90, pn="near"); names(lst90p.ssp126.2050)<-"lst90p.ssp126.2050"
lst90p.ssp126.2100 <- txpFUN(r=lst.sst126, nYrs=30, baseline=lst.Hist, p=0.90, pn="far"); names(lst90p.ssp126.2100)<-"lst90p.ssp126.2100"
lst90p.ssp245.2050 <- txpFUN(r=lst.sst245, nYrs=30, baseline=lst.Hist, p=0.90, pn="near"); names(lst90p.ssp245.2050)<-"lst90p.ssp245.2050"
lst90p.ssp245.2100 <- txpFUN(r=lst.sst245, nYrs=30, baseline=lst.Hist, p=0.90, pn="far") ; names(lst90p.ssp245.2100)<-"lst90p.ssp245.2100"
lst90p.ssp370.2050 <- txpFUN(r=lst.sst370, nYrs=30, baseline=lst.Hist, p=0.90, pn="near") ; names(lst90p.ssp370.2050)<-"lst90p.ssp370.2050"
lst90p.ssp370.2100 <- txpFUN(r=lst.sst370, nYrs=30, baseline=lst.Hist, p=0.90, pn="far") ; names(lst90p.ssp370.2100)<-"lst90p.ssp370.2100"
lst90p.ssp585.2050 <- txpFUN(r=lst.sst585, nYrs=30, baseline=lst.Hist, p=0.90, pn="near") ; names(lst90p.ssp585.2050)<-"lst90p.ssp585.2050"
lst90p.ssp585.2100 <- txpFUN(r=lst.sst585, nYrs=30, baseline=lst.Hist, p=0.90, pn="far") ; names(lst90p.ssp585.2100)<-"lst90p.ssp585.2100"
#Intensity of extremes
lst90Int.ssp126.2050 <- txIntFUN(r=lst.sst126, nYrs=30, baseline=lst.Hist, p=0.90, pn="near") ; names(lst90Int.ssp126.2050)<-"lst90Int.ssp126.2050"
lst90Int.ssp126.2100 <- txIntFUN(r=lst.sst126, nYrs=30, baseline=lst.Hist, p=0.90, pn="far") ; names(lst90Int.ssp126.2100)<-"lst90Int.ssp126.2100"
lst90Int.ssp245.2050 <- txIntFUN(r=lst.sst245, nYrs=30, baseline=lst.Hist, p=0.90, pn="near") ; names(lst90Int.ssp245.2050)<-"lst90Int.ssp245.2050"
lst90Int.ssp245.2100 <- txIntFUN(r=lst.sst245, nYrs=30, baseline=lst.Hist, p=0.90, pn="far") ; names(lst90Int.ssp245.2100)<-"lst90Int.ssp245.2100"
lst90Int.ssp370.2050 <- txIntFUN(r=lst.sst370, nYrs=30, baseline=lst.Hist, p=0.90, pn="near") ; names(lst90Int.ssp370.2050)<-"lst90Int.ssp370.2050"
lst90Int.ssp370.2100 <- txIntFUN(r=lst.sst370, nYrs=30, baseline=lst.Hist, p=0.90, pn="far") ; names(lst90Int.ssp370.2100)<-"lst90Int.ssp370.2100"
lst90Int.ssp585.2050 <- txIntFUN(r=lst.sst585, nYrs=30, baseline=lst.Hist, p=0.90, pn="near") ; names(lst90Int.ssp585.2050)<-"lst90Int.ssp585.2050"
lst90Int.ssp585.2100 <- txIntFUN(r=lst.sst585, nYrs=30, baseline=lst.Hist, p=0.90, pn="far") ; names(lst90Int.ssp585.2100)<-"lst90Int.ssp585.2100"
#Trend
lst.trend.ssp126.2050 <- slpFUN(lst.sst126[[73:432]]) ;names(lst.trend.ssp126.2050)<-"lst.trend.ssp126.2050"
lst.trend.ssp126.2100 <- slpFUN(lst.sst126[[673:1032]]) ;names(lst.trend.ssp126.2100)<-"lst.trend.ssp126.2100"
lst.trend.ssp245.2050 <- slpFUN(lst.sst245[[73:432]]) ;names(lst.trend.ssp245.2050)<-"lst.trend.ssp245.2050"
lst.trend.ssp245.2100 <- slpFUN(lst.sst245[[673:1032]]);names(lst.trend.ssp245.2100)<-"lst.trend.ssp245.2100"
lst.trend.ssp370.2050 <- slpFUN(lst.sst370[[73:432]]) ;names(lst.trend.ssp370.2050)<-"lst.trend.ssp370.2050"
lst.trend.ssp370.2100 <- slpFUN(lst.sst370[[673:1032]]) ;names(lst.trend.ssp370.2100)<-"lst.trend.ssp370.2100"
lst.trend.ssp585.2050 <- slpFUN(lst.sst585[[73:432]]) ;names(lst.trend.ssp585.2050)<-"lst.trend.ssp585.2050"
lst.trend.ssp585.2100 <- slpFUN(lst.sst585[[673:1032]]) ;names(lst.trend.ssp585.2100)<-"lst.trend.ssp585.2100"

hzd.LST <- stack(inormal(lst90p.ssp126.2050),
                 inormal(lst90p.ssp126.2100),
                 inormal(lst90p.ssp245.2050),
                 inormal(lst90p.ssp245.2100),
                 inormal(lst90p.ssp370.2050),
                 inormal(lst90p.ssp370.2100),
                 inormal(lst90p.ssp585.2050),
                 inormal(lst90p.ssp585.2100),
                 
                 inormal(lst90Int.ssp126.2050),
                 inormal(lst90Int.ssp126.2100),
                 inormal(lst90Int.ssp245.2050),
                 inormal(lst90Int.ssp245.2100),
                 inormal(lst90Int.ssp370.2050),
                 inormal(lst90Int.ssp370.2100),
                 inormal(lst90Int.ssp585.2050),
                 inormal(lst90Int.ssp585.2100),
                 
                 inormal(lst.trend.ssp126.2050),
                 inormal(lst.trend.ssp126.2100),
                 inormal(lst.trend.ssp245.2050),
                 inormal(lst.trend.ssp245.2100),
                 inormal(lst.trend.ssp370.2050),
                 inormal(lst.trend.ssp370.2100),
                 inormal(lst.trend.ssp585.2050),
                 inormal(lst.trend.ssp585.2100))
plot(hzd.LST)
# if (require(ncdf4)) {
#   rnc <- raster::writeRaster(stack((lst90p.ssp126.2050),(lst90p.ssp126.2100),
#                                    (lst90p.ssp245.2050),(lst90p.ssp245.2100),
#                                    (lst90p.ssp370.2050),(lst90p.ssp370.2100),
#                                    (lst90p.ssp585.2050),(lst90p.ssp585.2100)), 
#                              filename=file.path("2_Data/raster/wio_esm_ts90p.nc"), format="CDF", overwrite=TRUE)
# }
# 
# if (require(ncdf4)) {
#   rnc <- raster::writeRaster(stack((lst90Int.ssp126.2050),(lst90Int.ssp126.2100),
#                                    (lst90Int.ssp245.2050),(lst90Int.ssp245.2100),
#                                    (lst90Int.ssp370.2050),(lst90Int.ssp370.2100),
#                                    (lst90Int.ssp585.2050),(lst90Int.ssp585.2100)), 
#                              filename=file.path("2_Data/raster/wio_esm_ts90Int.nc"), format="CDF", overwrite=TRUE)
# }
# 
# if (require(ncdf4)) {
#   rnc <- raster::writeRaster(stack( (lst.trend.ssp126.2050),(lst.trend.ssp126.2100),
#                                     (lst.trend.ssp245.2050),(lst.trend.ssp245.2100),
#                                     (lst.trend.ssp370.2050),(lst.trend.ssp370.2100),
#                                     (lst.trend.ssp585.2050),(lst.trend.ssp585.2100)), 
#                              filename=file.path("2_Data/raster/wio_esm_ts_trend.nc"), format="CDF", overwrite=TRUE)
# }

#RAINFALL
rr.prHist <- stack(paste0(lfs.dir, "CLIM-CCVA/pr_historical_ESM25km.nc")); rr.prHist<- crop(rr.prHist,focusISO3)
rr.pr126 <- stack(paste0(lfs.dir, "CLIM-CCVA/pr126_ESM25km.nc")); rr.pr126<- crop(rr.pr126,focusISO3)
rr.pr245 <- stack(paste0(lfs.dir, "CLIM-CCVA/pr245_ESM25km.nc")); rr.pr245<- crop(rr.pr245,focusISO3)
rr.pr370 <- stack(paste0(lfs.dir, "CLIM-CCVA/pr370_ESM25km.nc")); rr.pr370<- crop(rr.pr370,focusISO3)
rr.pr585 <- stack(paste0(lfs.dir, "CLIM-CCVA/pr585_ESM25km.nc")); rr.pr585<- crop(rr.pr585,focusISO3)

# m1 <- list(rr.prHist, rr.pr126, rr.pr245, rr.pr370, rr.pr585)
# names(m1) <- c("HIST", "SSP126", "SSP245", "SSP370", "SSP585")
# rr <- lapply(1:length(m1), function(x){raster::cellStats(m1[[x]], stat='mean', na.rm=TRUE)})
# dfx <- c()
# for (k in 1:length(rr)){
#   dfx[[k]] <- rr[[k]] %>% as.data.frame() %>% rownames_to_column()
# }
# 
# names(dfx) <- c("HIST", "SSP126", "SSP245", "SSP370", "SSP585")
# dfx <- bind_rows(dfx, .id = "column_label");colnames(dfx) <- c("Mod_Sce_Var", "time", "value")
# #dfx <- dfx %>% separate_wider_delim(Mod_Sce_Var, "_", names = c("Model", "Scenario", "Variable")) %>% as.data.frame()
# dfx$Date <- format(as.Date(dfx$time, format = "X%Y.%m.%d"), format = "%Y")
# ensemble_tas <- dfx %>% group_by(Mod_Sce_Var,Date) %>% summarise(value = median(value, na.rm=TRUE)) %>% ungroup() %>% data.frame()

#Plot time series for Temperature or Precipitation
# (plt2 <- ensemble_tas %>% 
#     ggplot(aes(x=as.numeric(Date), y=value, colour=Mod_Sce_Var)) +
#     geom_line()+
#     theme_classic(base_size = 15)+
#     scale_y_continuous(name = "Ensemble temperature")+
#     theme(axis.text.x = element_text(angle = 0)) + 
#     scale_x_continuous(name = "Calendar year", breaks = seq(1970, 2100, 10), expand = c(0, 0)))
#freq
R10p.ssp126.2050 <- R10Freq(r=rr.pr126, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10p.ssp126.2050)<-"R10p.ssp126.2050"
R10p.ssp126.2100 <- R10Freq(r=rr.pr126, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10p.ssp126.2100)<-"R10p.ssp126.2100"
R10p.ssp245.2050 <- R10Freq(r=rr.pr245, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10p.ssp245.2050)<-"R10p.ssp245.2050"
R10p.ssp245.2100 <- R10Freq(r=rr.pr245, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10p.ssp245.2100)<-"R10p.ssp245.2100"
R10p.ssp370.2050 <- R10Freq(r=rr.pr370, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10p.ssp370.2050)<-"R10p.ssp370.2050"
R10p.ssp370.2100 <- R10Freq(r=rr.pr370, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10p.ssp370.2100)<-"R10p.ssp370.2100"
R10p.ssp585.2050 <- R10Freq(r=rr.pr585, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10p.ssp585.2050)<-"R10p.ssp585.2050"
R10p.ssp585.2100 <- R10Freq(r=rr.pr585, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10p.ssp585.2100)<-"R10p.ssp585.2100"
#Intensity
R10Int.ssp126.2050 <- R10IntFUN(r=rr.pr126, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10Int.ssp126.2050)<-"R10Int.ssp126.2050"
R10Int.ssp126.2100 <- R10IntFUN(r=rr.pr126, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10Int.ssp126.2100)<-"R10Int.ssp126.2100"
R10Int.ssp245.2050 <- R10IntFUN(r=rr.pr245, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10Int.ssp245.2050)<-"R10Int.ssp245.2050"
R10Int.ssp245.2100 <- R10IntFUN(r=rr.pr245, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10Int.ssp245.2100)<-"R10Int.ssp245.2100"
R10Int.ssp370.2050 <- R10IntFUN(r=rr.pr370, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10Int.ssp370.2050)<-"R10Int.ssp370.2050"
R10Int.ssp370.2100 <- R10IntFUN(r=rr.pr370, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10Int.ssp370.2100)<-"R10Int.ssp370.2100"
R10Int.ssp585.2050 <- R10IntFUN(r=rr.pr585, nYears=30, baseR=rr.prHist, p=.1, pn="near") %>% mask(., focusISO3); names(R10Int.ssp585.2050)<-"R10Int.ssp585.2050"
R10Int.ssp585.2100 <- R10IntFUN(r=rr.pr585, nYears=30, baseR=rr.prHist, p=.1, pn="far") %>% mask(., focusISO3); names(R10Int.ssp585.2100)<-"R10Int.ssp585.2100"
#Slope
#To calculate total annual precipitation, sum monthly precipitation values and divide by 30 years. Next multiply by 86400 to convert from 
#kg/m3/s to mm/year
TAP.ssp126.2050 <- (slpFUN(rr.pr126[[73:432]])*86400) %>% mask(., focusISO3); names(TAP.ssp126.2050) <- "dTAP.ssp126.2050"
TAP.ssp126.2100 <- (slpFUN(rr.pr126[[673:1032]])*86400) %>% mask(., focusISO3); names(TAP.ssp126.2100) <- "dTAP.ssp126.2100"
TAP.ssp245.2050 <- (slpFUN(rr.pr245[[73:432]])*86400) %>% mask(., focusISO3); names(TAP.ssp245.2050) <- "dTAP.ssp245.2050"
TAP.ssp245.2100 <- (slpFUN(rr.pr245[[673:1032]])*86400) %>% mask(., focusISO3); names(TAP.ssp245.2100) <- "dTAP.ssp245.2100"
TAP.ssp370.2050 <- (slpFUN(rr.pr370[[73:432]])*86400) %>% mask(., focusISO3); names(TAP.ssp370.2050) <- "dTAP.ssp370.2050"
TAP.ssp370.2100 <- (slpFUN(rr.pr370[[673:1032]])*86400) %>% mask(., focusISO3); names(TAP.ssp370.2100) <- "dTAP.ssp370.2100"
TAP.ssp585.2050 <- (slpFUN(rr.pr585[[73:432]])*86400) %>% mask(., focusISO3); names(TAP.ssp585.2050) <- "dTAP.ssp585.2050"
TAP.ssp585.2100 <- (slpFUN(rr.pr585[[673:1032]])*86400) %>% mask(., focusISO3); names(TAP.ssp585.2100) <- "dTAP.ssp585.2100"

hzd.TAP <- stack(inormal(TAP.ssp126.2050),
                 inormal(TAP.ssp126.2100), 
                 inormal(TAP.ssp245.2050),
                 inormal(TAP.ssp245.2100),
                 inormal(TAP.ssp370.2050), 
                 inormal(TAP.ssp370.2100),
                 inormal(TAP.ssp585.2050), 
                 inormal(TAP.ssp585.2100),

                 inormal(R10p.ssp126.2050),
                 inormal(R10p.ssp126.2100),
                 inormal(R10p.ssp245.2050),
                 inormal(R10p.ssp245.2100),
                 inormal(R10p.ssp370.2050),
                 inormal(R10p.ssp370.2100),
                 inormal(R10p.ssp585.2050),
                 inormal(R10p.ssp585.2100),
                 
                 inormal(R10Int.ssp126.2050),
                 inormal(R10Int.ssp126.2100),
                 inormal(R10Int.ssp245.2050),
                 inormal(R10Int.ssp245.2100),
                 inormal(R10Int.ssp370.2050),
                 inormal(R10Int.ssp370.2100),
                 inormal(R10Int.ssp585.2050),
                 inormal(R10Int.ssp585.2100))

plot(hzd.TAP)
summary(hzd.TAP)
# if (require(ncdf4)) {
#   rnc <- raster::writeRaster(stack((R10p.ssp126.2050),(R10p.ssp126.2100),
#                                    (R10p.ssp245.2050),(R10p.ssp245.2100),
#                                    (R10p.ssp370.2050),(R10p.ssp370.2100),
#                                    (R10p.ssp585.2050),(R10p.ssp585.2100)), 
#                              filename=file.path("2_Data/raster/wio_esm_r10p.nc"), format="CDF", overwrite=TRUE)
# }
# 
# if (require(ncdf4)) {
#   rnc <- raster::writeRaster(stack((R10Int.ssp126.2050),(R10Int.ssp126.2100),
#                                    (R10Int.ssp245.2050),(R10Int.ssp245.2100),
#                                    (R10Int.ssp370.2050),(R10Int.ssp370.2100),
#                                    (R10Int.ssp585.2050),(R10Int.ssp585.2100)), 
#                              filename=file.path("2_Data/raster/wio_esm_r10Int.nc"), format="CDF", overwrite=TRUE)
# }
# 
# if (require(ncdf4)) {
#   rnc <- raster::writeRaster(stack((TAP.ssp126.2050),(TAP.ssp126.2100), 
#                                    (TAP.ssp245.2050),(TAP.ssp245.2100),
#                                    (TAP.ssp370.2050),(TAP.ssp370.2100),
#                                    (TAP.ssp585.2050),(TAP.ssp585.2100)), 
#                              filename=file.path("2_Data/raster/wio_esm_tap_trend.nc"), format="CDF", overwrite=TRUE)
# }

#Net Primary Productivity
#https://www.convertunits.com/from/kg-m/s/to/ton
rr.npp126 <- stack(paste0(lfs.dir, "CLIM-CCVA/npp126_ESM25km.nc")); rr.npp126 <- crop(rr.npp126,focusISO3)
rr.npp245 <- stack(paste0(lfs.dir, "CLIM-CCVA/npp245_ESM25km.nc")); rr.npp245 <- crop(rr.npp245,focusISO3)
rr.npp370 <- stack(paste0(lfs.dir, "CLIM-CCVA/npp370_ESM25km.nc")); rr.npp370 <- crop(rr.npp370,focusISO3)
rr.npp585 <- stack(paste0(lfs.dir, "CLIM-CCVA/npp585_ESM25km.nc")); rr.npp585 <- crop(rr.npp585,focusISO3)

m1 <- list(rr.npp126, rr.npp245, rr.npp370, rr.npp585)
names(m1) <- c("HIST", "SSP126", "SSP245", "SSP370", "SSP585")
rr <- lapply(1:length(m1), function(x){raster::cellStats(m1[[x]], stat='mean', na.rm=TRUE)})
dfx <- c()
for (k in 1:length(rr)){
  dfx[[k]] <- rr[[k]] %>% as.data.frame() %>% rownames_to_column()
}

names(dfx) <- c("SSP126", "SSP245", "SSP370", "SSP585")
dfx <- bind_rows(dfx, .id = "column_label");colnames(dfx) <- c("Mod_Sce_Var", "time", "value")
#dfx <- dfx %>% separate_wider_delim(Mod_Sce_Var, "_", names = c("Model", "Scenario", "Variable")) %>% as.data.frame()
dfx$Date <- format(as.Date(dfx$time, format = "X%Y.%m.%d"), format = "%Y")
ensemble_tas <- dfx %>% group_by(Mod_Sce_Var,Date) %>% summarise(value = median(value, na.rm=TRUE)) %>% ungroup() %>% data.frame()
#write_csv(ensemble_tas, "outputs/tas_outxx.csv")

#Plot time series for Temperature or Precipitation
(plt2 <- ensemble_tas %>% 
    ggplot(aes(x=as.numeric(Date), y=value, colour=Mod_Sce_Var)) +
    geom_line()+
    theme_classic(base_size = 15)+
    scale_y_continuous(name = "Ensemble temperature")+
    theme(axis.text.x = element_text(angle = 0)) + 
    scale_x_continuous(name = "Calendar year", breaks = seq(1970, 2100, 10), expand = c(0, 0)))

NPP.ssp126.2050 <- slpFUN(rr.npp126[[73:432]]) %>% mask(., focusISO3); names(NPP.ssp126.2050) <- "NPP.ssp126.2050"
NPP.ssp126.2100 <- slpFUN(rr.npp126[[673:1032]]) %>% mask(., focusISO3); names(NPP.ssp126.2100) <- "NPP.ssp126.2100"
NPP.ssp245.2050 <- slpFUN(rr.npp245[[73:432]]) %>% mask(., focusISO3); names(NPP.ssp245.2050) <- "NPP.ssp245.2050"
NPP.ssp245.2100 <- slpFUN(rr.npp245[[673:1032]]) %>% mask(., focusISO3); names(NPP.ssp245.2100) <- "NPP.ssp245.2100"
NPP.ssp370.2050 <- slpFUN(rr.npp370[[73:432]]) %>% mask(., focusISO3); names(NPP.ssp370.2050) <- "NPP.ssp370.2050"
NPP.ssp370.2100 <- slpFUN(rr.npp370[[673:1031]]) %>% mask(., focusISO3); names(NPP.ssp370.2100) <- "NPP.ssp370.2100"
NPP.ssp585.2050 <- slpFUN(rr.npp585[[73:432]]) %>% mask(., focusISO3); names(NPP.ssp585.2050) <- "NPP.ssp585.2050"
NPP.ssp585.2100 <- slpFUN(rr.npp585[[673:1031]]) %>% mask(., focusISO3); names(NPP.ssp585.2100) <- "NPP.ssp585.2100"
hzd.NPP <- raster::stack(inormal(NPP.ssp126.2050),
                         inormal(NPP.ssp126.2100),
                         inormal(NPP.ssp245.2050),
                         inormal(NPP.ssp245.2100),
                         inormal(NPP.ssp370.2050),
                         inormal(NPP.ssp370.2100),
                         inormal(NPP.ssp585.2050),
                         inormal(NPP.ssp585.2100))
plot(hzd.NPP)
if (require(ncdf4)) {
  rnc <- raster::writeRaster(stack((NPP.ssp126.2050),(NPP.ssp126.2100),
                                   (NPP.ssp245.2050),(NPP.ssp245.2100),
                                   (NPP.ssp370.2050),(NPP.ssp370.2100),
                                   (NPP.ssp585.2050),(NPP.ssp585.2100)),
                             filename=file.path("2_Data/raster/wio_esm_npp.nc"), format="CDF", overwrite=TRUE)
}



## Consecutive Dry Days
rr.cdd126 <- stack(paste0(lfs.dir, "CLIM-CCVA/cdd126_ESM25km.nc")); rr.cdd126 <- crop(rr.cdd126,focusISO3)
rr.cdd245 <- stack(paste0(lfs.dir, "CLIM-CCVA/cdd245_ESM25km.nc")); rr.cdd245 <- crop(rr.cdd245,focusISO3)
rr.cdd370 <- stack(paste0(lfs.dir, "CLIM-CCVA/cdd370_ESM25km.nc")); rr.cdd370 <- crop(rr.cdd370,focusISO3)
rr.cdd585 <- stack(paste0(lfs.dir, "CLIM-CCVA/cdd585_ESM25km.nc")); rr.cdd585 <- crop(rr.cdd585,focusISO3)

m1 <- list(rr.cdd126, rr.cdd245, rr.cdd370,rr.cdd585)
rr <- lapply(1:length(m1), function(x){raster::cellStats(m1[[x]], stat='mean', na.rm=TRUE)})
dfx <- c()
for (k in 1:length(rr)){
  dfx[[k]] <- rr[[k]] %>% as.data.frame() %>% rownames_to_column()
}

names(dfx) <- c("SSP126", "SSP245", "SSP370", "SSP585")
dfx <- bind_rows(dfx, .id = "column_label");colnames(dfx) <- c("Mod_Sce_Var", "time", "value")
#dfx <- dfx %>% separate_wider_delim(Mod_Sce_Var, "_", names = c("Model", "Scenario", "Variable")) %>% as.data.frame()
dfx$Date <- format(as.Date(dfx$time, format = "X%Y.%m.%d"), format = "%Y")
ensemble_tas <- dfx %>% group_by(Mod_Sce_Var,Date) %>% summarise(value = median(value, na.rm=TRUE)) %>% ungroup() %>% data.frame()
#write_csv(ensemble_tas, "outputs/tas_outxx.csv")

#Plot time series for Temperature or Precipitation
(plt2 <- ensemble_tas %>% 
    ggplot(aes(x=as.numeric(Date), y=value, colour=Mod_Sce_Var)) +
    geom_line()+
    theme_classic(base_size = 15)+
    scale_y_continuous(name = "Ensemble temperature")+
    theme(axis.text.x = element_text(angle = 0)) + 
    scale_x_continuous(name = "Calendar year", breaks = seq(1970, 2100, 10), expand = c(0, 0)))

CDD.ssp126.2050 <- calc(rr.cdd126[[6:35]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp126.2050)<-"CDD.ssp126.2050"
CDD.ssp126.2100 <- calc(rr.cdd126[[57:86]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp126.2100)<-"CDD.ssp126.2100"
CDD.ssp245.2050 <- calc(rr.cdd245[[6:35]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp245.2050)<-"CDD.ssp245.2050"
CDD.ssp245.2100 <- calc(rr.cdd245[[57:86]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp245.2100)<-"CDD.ssp245.2100"
CDD.ssp370.2050 <- calc(rr.cdd370[[6:35]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp370.2050)<-"CDD.ssp370.2050"
CDD.ssp370.2100 <- calc(rr.cdd370[[57:86]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp370.2100)<-"CDD.ssp370.2100"
CDD.ssp585.2050 <- calc(rr.cdd585[[6:35]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp585.2050)<-"CDD.ssp585.2050"
CDD.ssp585.2100 <- calc(rr.cdd585[[57:86]], fun = mean) %>% mask(., focusISO3) ;names(CDD.ssp585.2100)<-"CDD.ssp585.2100"
hzd.CDD <- stack(inormal(CDD.ssp126.2050), 
                 inormal(CDD.ssp126.2100), 
                 inormal(CDD.ssp245.2050),
                 inormal(CDD.ssp245.2100),
                 inormal(CDD.ssp370.2050), 
                 inormal(CDD.ssp370.2100),
                 inormal(CDD.ssp585.2050), 
                 inormal(CDD.ssp585.2100))
plot(hzd.CDD)
if (require(ncdf4)) {
  rnc <- raster::writeRaster(stack((CDD.ssp126.2050),(CDD.ssp126.2100), 
                                   (CDD.ssp245.2050),(CDD.ssp245.2100),
                                   (CDD.ssp370.2050),(CDD.ssp370.2100),
                                   (CDD.ssp585.2050),(CDD.ssp585.2100)),
                             filename=file.path("2_Data/raster/wio_esm_cdd.nc"), format="CDF", overwrite=TRUE)
}


###EVAPOTRANSPIRATION
#https://www.convertunits.com/from/kg-m/s/to/ton
rr.evspsbl126 <- stack(paste0(lfs.dir, "CLIM-CCVA/evspsbl126_ESM25km.nc")); rr.evspsbl126 <- crop(rr.evspsbl126,focusISO3)
rr.evspsbl245 <- stack(paste0(lfs.dir, "CLIM-CCVA/evspsbl245_ESM25km.nc")); rr.evspsbl245 <- crop(rr.evspsbl245,focusISO3)
rr.evspsbl370 <- stack(paste0(lfs.dir, "CLIM-CCVA/evspsbl370_ESM25km.nc")); rr.evspsbl370 <- crop(rr.evspsbl370,focusISO3)
rr.evspsbl585 <- stack(paste0(lfs.dir, "CLIM-CCVA/evspsbl585_ESM25km.nc")); rr.evspsbl585 <- crop(rr.evspsbl585,focusISO3)

m1 <- list(rr.evspsblHist, rr.evspsbl126, rr.evspsbl245, rr.evspsbl370, rr.evspsbl585)
names(m1) <- c("HIST", "SSP126", "SSP245", "SSP370", "SSP585")
rr <- lapply(1:length(m1), function(x){raster::cellStats(m1[[x]], stat='mean', na.rm=TRUE)})
dfx <- c()
for (k in 1:length(rr)){
  dfx[[k]] <- rr[[k]] %>% as.data.frame() %>% rownames_to_column()
}

names(dfx) <- c("SSP126", "SSP245", "SSP370", "SSP585")
dfx <- bind_rows(dfx, .id = "column_label");colnames(dfx) <- c("Mod_Sce_Var", "time", "value")
#dfx <- dfx %>% separate_wider_delim(Mod_Sce_Var, "_", names = c("Model", "Scenario", "Variable")) %>% as.data.frame()
dfx$Date <- format(as.Date(dfx$time, format = "X%Y.%m.%d"), format = "%Y")
ensemble_tas <- dfx %>% group_by(Mod_Sce_Var,Date) %>% summarise(value = median(value, na.rm=TRUE)) %>% ungroup() %>% data.frame()
#write_csv(ensemble_tas, "outputs/tas_outxx.csv")

#Plot time series for Temperature or Precipitation
(plt2 <- ensemble_tas %>% 
    ggplot(aes(x=as.numeric(Date), y=value, colour=Mod_Sce_Var)) +
    geom_line()+
    theme_classic(base_size = 15)+
    scale_y_continuous(name = "Ensemble EVAP")+
    theme(axis.text.x = element_text(angle = 0)) + 
    scale_x_continuous(name = "Calendar year", breaks = seq(1970, 2100, 10), expand = c(0, 0)))

EVSPSBL.ssp126.2050 <- slpFUN(rr.evspsbl126[[73:432]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp126.2050)<-"EVSPSBL.ssp126.2050"
EVSPSBL.ssp126.2100 <- slpFUN(rr.evspsbl126[[673:1032]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp126.2100)<-"EVSPSBL.ssp126.2100"
EVSPSBL.ssp245.2050 <- slpFUN(rr.evspsbl245[[73:432]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp245.2050)<-"EVSPSBL.ssp245.2050"
EVSPSBL.ssp245.2100 <- slpFUN(rr.evspsbl245[[673:1032]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp245.2100)<-"EVSPSBL.ssp245.2100"
EVSPSBL.ssp370.2050 <- slpFUN(rr.evspsbl370[[73:432]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp370.2050)<-"EVSPSBL.ssp370.2050"
EVSPSBL.ssp370.2100 <- slpFUN(rr.evspsbl370[[673:1032]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp370.2100)<-"EVSPSBL.ssp370.2100"
EVSPSBL.ssp585.2050 <- slpFUN(rr.evspsbl585[[73:432]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp585.2050)<-"EVSPSBL.ssp585.2050"
EVSPSBL.ssp585.2100 <- slpFUN(rr.evspsbl585[[673:1032]]) %>% mask(., focusISO3) ;names(EVSPSBL.ssp585.2100)<-"EVSPSBL.ssp585.2100"
hzd.EVAP <- stack(inormal(EVSPSBL.ssp126.2050),
                  inormal(EVSPSBL.ssp126.2100),
                  inormal(EVSPSBL.ssp245.2050), 
                  inormal(EVSPSBL.ssp245.2100),
                  inormal(EVSPSBL.ssp370.2050), 
                  inormal(EVSPSBL.ssp370.2100),
                  inormal(EVSPSBL.ssp585.2050), 
                  inormal(EVSPSBL.ssp585.2100))
plot(hzd.EVAP)
if (require(ncdf4)) {
  rnc <- raster::writeRaster(stack((EVSPSBL.ssp126.2050),(EVSPSBL.ssp126.2100),
                                   (EVSPSBL.ssp245.2050),(EVSPSBL.ssp245.2100),
                                   (EVSPSBL.ssp370.2050),(EVSPSBL.ssp370.2100),
                                   (EVSPSBL.ssp585.2050),(EVSPSBL.ssp585.2100)),
                             filename=file.path("2_Data/raster/wio_esm_evspsbl.nc"), format="CDF", overwrite=TRUE)
}









