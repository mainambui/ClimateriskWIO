#PREPARE TO EXTRACT DATA TO AOO
rm(list = ls())

library(tidyverse);library(raster);library(SearchTrees);library(sf);library(sp)
#import some important functions
source("1_Codes/3_KeyFunctions.R")

#Import PU
wio.AOO <- readRDS("2_Data/sheet/2_Ecosystems/wioAOO.wTEV.rds")
wio.AOO.spdf <- st_as_sf(wio.AOO, coords=c('x', 'y'), crs="+proj=longlat")
wio.ISO3 <- st_read("2_Data/shp/country_shape.shp") %>% st_as_sf() %>% st_transform(crs = "+proj=longlat")

#LOAD METRICS
(clim.nc <- list.files("./2_Data/raster", pattern='*.nc',full.names=TRUE))
nc <- (as.data.frame(expand.grid(x=c(126,245,370,585), y=c(2050,2100))) %>% arrange(desc(-x)) %>% mutate(cbn = paste(x,y,sep = "_")) %>% dplyr::select(cbn))[,1]
varlst <- c("cdd", "evspsbl", "npp", "pH_trend", "r10Int","r10p","sst_trend","sst90Int","sst90p","tap_trend","ts90Int","ts90p","ts_trend")
#varlst <- c("sst90Int","sst90p","ts90Int","ts90p")

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
write_rds(climdata, "2_Data/sheet/1_Climate/wioAOO.metrics.rds")

#limit <- max(abs(climdataNew$tx90p_thresh), na.rm = TRUE) * c(-1, 1)
climdata %>% 
  filter(is.finite(imp.ssp370.2050))%>% 
  ggplot()+ geom_sf(data = wio.ISO3, colour = "grey70", fill = "grey95")+
  geom_raster(aes(x,y,fill=imp.ssp370.2050))+
  scale_fill_viridis_c(name = "tx90p[%]", option = "B")+
  #rcartocolor::scale_fill_carto_c(name="", palette = 'Earth', direction = -1)+
  scale_x_continuous(name = "", expand = c(0,0))+
  scale_y_continuous(name = "", expand = c(0,0))+
  #labs(title = "SSTmx greater than 90p baseline threshold [2021-2050, SSP245]")+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(.5, 'cm'))
#ggsave("sst90p2050_ssp245.png", dpi = 1200,width = 11.5, height= 12)

####################################################################################################################################################################
#REPARE TO EXTRACT IMPACTS DATASET TO COMMUNITY LEVEL#
###################################################################################################################################################################

#https://rpubs.com/dieghernan/beautifulmaps_I
library(sf)
library(dplyr)

# Projecting and cleaning
wio_cline <- st_read("C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/Data/2-Exposed systems/shp/wio_coastline_pu.shp")
st_crs(wio_cline)<- "+proj=longlat"
wio_cline <- st_transform(wio_cline, "+proj=moll")

#Read the social datasets
socioecom <- read_csv("2_Data/sheet/3_SocialVulnerability/SocialDataAll.csv")
socioecom <- filter(socioecom, !is.na(x))
#Convert to spatial object
spdf <- st_as_sf(socioecom, coords=c('x', 'y'), crs="+proj=longlat")
spdf <- st_transform(spdf, "+proj=moll")

#Intersect and assign Grid Index IDs to the community level spatial points
spdf <- st_intersection(wio_cline, spdf)
# socID <- data.frame(spdf) %>% mutate(ID = PageName)%>% dplyr::select(ID)
# 
# #Convert the Fishnet to 
# df2 <- cbind.data.frame("ID" = wio_cline$PageName, st_coordinates(st_centroid(st_geometry(wio_cline))))
# colnames(df2)[colnames(df2)=="X"]<-"x"
# colnames(df2)[colnames(df2)=="Y"]<-"y"
# wio.gridNetworks <- expand.grid(src=df2$ID, nbr=df2$ID) %>% filter(src %in% socID$ID)
# 
# colnames(df2)[colnames(df2)=="ID"]<-"src"
# wio.gridNetworks <- inner_join(wio.gridNetworks, df2, by="src")
# colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="x"]<-"src_x"
# colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="y"]<-"src_y"
# 
# colnames(df2)[colnames(df2)=="src"]<-"nbr"
# wio.gridNetworks <- inner_join(wio.gridNetworks, df2, by="nbr")
# colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="x"]<-""
# colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="y"]<-"nbr_y"
# #REMAIN it back to ID
# colnames(df2)[colnames(df2)=="nbr"]<-"ID"
# 
# library(gdistance)
# library(doParallel)
# n.cores <- 6
# #Compute Euclidean distances between RCPs
# system.time({
#   cl <- makeCluster(n.cores)
#   registerDoParallel(cl)
#   clusterExport(cl, c('wio.gridNetworks'))  ## Export the environment variables to each cluster
#   clusterEvalQ(cl, {
#     library(gdistance)
#     library(sf)}
#   )  ## Load the library "gdistance" to each cluster
#   
#   the.edist.function <- function(k) {
#     eucDist <- st_distance(st_as_sf(wio.gridNetworks[k,], coords= c('src_x', 'src_y'), crs="+proj=moll"),
#                            st_as_sf(wio.gridNetworks[k,], coords=c('nbr_x', 'nbr_y'), crs="+proj=moll"))
#     return(eucDist)}
#   eucDist <- clusterApplyLB(cl, x=1:nrow(wio.gridNetworks), fun=the.edist.function)
# })
# 
# stopCluster(cl)
# 
# wio.gridNetworks$EucDist <- as.vector(eucDist %>% do.call("rbind", .))
# hist(wio.gridNetworks$EucDist, breaks=50)
# write_csv(wio.gridNetworks, "wio_DistMatrix.csv")

#########################################################################################################################################
#Merge grid level impacts to the network
#########################################################################################################################################
wio_metrics <- read_rds("2_Data/sheet/1_Climate/wioAOO.metrics.rds")
EucDist <- read_csv("2_Data/sheet/wio_DistMatrix.csv")

EucDist <- EucDist[c("src", "nbr", "EucDist")] %>% filter(EucDist>0) #filter to remove self intersections #
EucDist$inverseDist <- (1/EucDist$EucDist)

colnames(wio_metrics)[colnames(wio_metrics)=="ID"]<-"nbr"
EucDist <- merge(EucDist, wio_metrics, by = "nbr")

#Visual checking of distribution
hist(EucDist$inverseDist, breaks = 30)
names(EucDist)
(idw.impacts <- EucDist %>% group_by(src) %>%
    summarise(
      imp.ssp126.2050 = sum(imp.ssp126.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
      imp.ssp245.2050 = sum(imp.ssp245.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
      imp.ssp370.2050 = sum(imp.ssp370.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
      imp.ssp585.2050 = sum(imp.ssp585.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
      corals = sum(CoralExt*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
      seagrass = sum(seagrassExt*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
      mangrove = sum(mangroveExt*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
      cropcover = sum(Cropland*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE)
    )%>% mutate(ID = src) %>% dplyr::select(-src)
)

(wio.com.idw.impacts <- merge(data.frame(spdf) %>% mutate(ID = PageName) %>% dplyr::select(-c(PageName, PageNumber, geometry)),
                              idw.impacts, by = "ID"))

plot(wio.com.idw.impacts$mangrove, wio.com.idw.impacts$imp.ssp370.2050)
write_csv(wio.com.idw.impacts, "2_Data/sheet/RiskDataFINAL.csv")
