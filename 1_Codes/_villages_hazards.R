#PREPARE TO EXTRACT DATA TO AOO
rm(list = ls())

library("tidyverse");library("raster");library("SearchTrees");library("sf");library("sp");library("pracma")
#import some important functions
source("1_Codes/3_KeyFunctions.R")

#Import PU
wio.AOO <- readRDS("2_Data/sheet/2_Ecosystems/wioAOO.wTEV.rds")
wio.AOO.spdf <- st_as_sf(wio.AOO, coords=c('x', 'y'), crs="+proj=longlat")
wio.ISO3 <- st_read("2_Data/shp/country_shape.shp") %>% st_as_sf() %>% st_transform(crs = "+proj=longlat")

#LOAD METRICS
(clim.nc <- list.files("./2_Data/raster", pattern='*.nc',full.names=TRUE))#[-c(1)]
nc <- (as.data.frame(expand.grid(x=c(126,245,370,585), y=c(2050,2100))) %>% arrange(desc(-x)) %>% mutate(cbn = paste(x,y,sep = "_")) %>% dplyr::select(cbn))[,1]
varlst <- c("cdd","evspsbl","npp","pH_trend","r10Int","r10p","sst_trend","sst90Int","sst90p","tap_trend","ts90Int","ts90p")

allRaster <- lapply(1:length(varlst), function(x){
  rr <- raster::brick(clim.nc[grep(varlst[[x]], clim.nc)])
  names(rr) <- paste(varlst[[x]], nc, sep = "_")
  return(rr)}
)

hazard.stk <- stack(allRaster);N <- nlayers(hazard.stk)
df.hzd <- as.data.frame(hazard.stk, xy=TRUE)
#Quantile transform
m1 <- names(hazard.stk)
rr <- lapply(1:length(m1), function(x){quant_transform(df.hzd, m1[[x]])})
rr <- stack(rr);plot(rr)

climdata <- raster::extract(rr, wio.AOO.spdf, sp=TRUE, df=TRUE) %>% as.data.frame()
colnames(climdata)[colnames(climdata) == "coords.x1"] <- "x"
colnames(climdata)[colnames(climdata) == "coords.x2"] <- "y"
climdata <- climdata %>% relocate(x:y, .before = ID)
climdata <- cbind(climdata[1:19], DMwR2::knnImputation(climdata[20:N], k = 3))

namelist <- colnames(climdata)
(climdata <- climdata %>% 
    #Create copies of the Exposed systems. the original versions will be needed later
    mutate(std_corals = CoralExt,
           std_seagrass = seagrassExt,
           std_mangrove = mangroveExt,
           std_cropcover = Cropland,
           std_Nb_sp = Nb_sp,
           std_FRic = FRic,
           std_FDiv = FDiv,
           std_FEve = FEve)%>% 
    mutate(across(std_corals:std_Nb_sp, ~ inormal(.x)))%>% #FRic, FDiv, and FEve are already normalised variables
    dplyr::select(-c(cost.reef,cost.sgrass,cost.mngrov,cost.crop,TEV))%>%
    rowwise() %>%
    mutate(
      #mean across 13 normalised climate metrics
      hzd.mn.ssp126.2050 = median(c_across(namelist[grep("_126_2050", namelist)]), na.rm=TRUE),
      hzd.mn.ssp245.2050 = median(c_across(namelist[grep("_245_2050", namelist)]), na.rm=TRUE),
      hzd.mn.ssp370.2050 = median(c_across(namelist[grep("_370_2050", namelist)]), na.rm=TRUE),
      hzd.mn.ssp585.2050 = median(c_across(namelist[grep("_585_2050", namelist)]), na.rm=TRUE),
      exp.mn.wioo = median(c_across(std_Nb_sp:std_FEve), na.rm=TRUE),
      
      #standard deviations across 13 normalised climate metrics
      hzd.sd.ssp126.2050 = sd(c_across(namelist[grep("_126_2050", namelist)]), na.rm=TRUE),
      hzd.sd.ssp245.2050 = sd(c_across(namelist[grep("_245_2050", namelist)]), na.rm=TRUE),
      hzd.sd.ssp370.2050 = sd(c_across(namelist[grep("_370_2050", namelist)]), na.rm=TRUE),
      hzd.sd.ssp585.2050 = sd(c_across(namelist[grep("_585_2050", namelist)]), na.rm=TRUE),
      exp.sd.wioo = sd(c_across(std_Nb_sp:std_FEve), na.rm=TRUE)) %>% 
    ungroup() %>%
    
    mutate(
      #Deduce inverse variance weights
      wgts.Exposure = (exp.sd.wioo/exp.mn.wioo)^-1,
      wgts.ssp126.2050 = (hzd.sd.ssp126.2050/hzd.mn.ssp126.2050)^-1,
      wgts.ssp245.2050 = (hzd.sd.ssp245.2050/hzd.mn.ssp245.2050)^-1,
      wgts.ssp370.2050 = (hzd.sd.ssp370.2050/hzd.mn.ssp370.2050)^-1,
      wgts.ssp585.2050 = (hzd.sd.ssp585.2050/hzd.mn.ssp585.2050)^-1,
      
      #Estimate potential impacts
      imp.ssp126.2050 = ((hzd.mn.ssp126.2050*wgts.ssp126.2050)*(exp.mn.wioo*wgts.Exposure))/(wgts.ssp126.2050+wgts.Exposure),
      imp.ssp245.2050 = ((hzd.mn.ssp245.2050*wgts.ssp245.2050)*(exp.mn.wioo*wgts.Exposure))/(wgts.ssp245.2050+wgts.Exposure),
      imp.ssp370.2050 = ((hzd.mn.ssp370.2050*wgts.ssp370.2050)*(exp.mn.wioo*wgts.Exposure))/(wgts.ssp370.2050+wgts.Exposure),
      imp.ssp585.2050 = ((hzd.mn.ssp585.2050*wgts.ssp585.2050)*(exp.mn.wioo*wgts.Exposure))/(wgts.ssp585.2050+wgts.Exposure)
    ))
hist(climdata$imp.ssp370.2050, breaks=30)
write_csv(climdata, "2_Data/sheet/2_ClimateMasterAOO.csv")

#########################################################################################################################################
#Merge grid level impacts to the network
#########################################################################################################################################
wio_metrics <- read_csv("2_Data/sheet/2_ClimateMasterAOO.csv")
EucDist <- read_csv("2_Data/sheet/wio_DistMatrix.csv")
villageGridID <- read_csv("2_Data/sheet/VillageWithGridIDs.csv")

EucDist <- EucDist[c("src", "nbr", "EucDist")] %>% filter(EucDist>0) #filter to remove self intersections #
EucDist$inverseDist <- (1/EucDist$EucDist)

colnames(wio_metrics)[colnames(wio_metrics)=="ID"]<-"nbr"
EucDist <- merge(EucDist, wio_metrics, by = "nbr")

#Visual checking of distribution
hist(EucDist$inverseDist, breaks = 30)
(idw.impacts <- EucDist %>% group_by(src) %>%
    summarise(across(CoralExt:imp.ssp585.2050, ~(sum(.x*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE))))%>% mutate(ID = src) %>% dplyr::select(-src)
)
(wio.com.idw.impacts <- merge(villageGridID, idw.impacts, by = "ID"))
plot(wio.com.idw.impacts$Cropland, wio.com.idw.impacts$npp_585_2050)
write_csv(wio.com.idw.impacts, "2_Data/sheet/3_villagesClimateExposureImpacts.csv")
