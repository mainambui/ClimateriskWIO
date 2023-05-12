#PREPARE TO EXTRACT DATA TO AOO
library(tidyverse);library(raster);library(SearchTrees);library(sf);library(sp)
#import some important functions
source("1_Codes/2_KeyFunctions.R")

#Import PU
wio.AOO <- readRDS("2_Data/spreadsheet/2_Ecosystems/wioAOO.wTEV.rds")
wio.AOO.spdf <- st_as_sf(wio.AOO, coords=c('x', 'y'), crs="+proj=longlat")
wio.ISO3 <- st_read("2_Data/shp/country_shape.shp") %>% st_as_sf() %>% st_transform(crs = "+proj=longlat")

#LOAD METRICS
clim.nc <- list.files("./2_Data/raster", pattern='*.nc',full.names=TRUE)
hazard.stk <- stack(clim.nc);N <- nlayers(hazard.stk)
hazard.stk <- extract(stack(clim.nc), wio.AOO.spdf, sp=TRUE,df=TRUE) %>% as.data.frame()
colnames(hazard.stk)[colnames(hazard.stk) == "coords.x1"] <- "x"
colnames(hazard.stk)[colnames(hazard.stk) == "coords.x2"] <- "y"
hazard.stk <- DMwR2::knnImputation(hazard.stk[2:123], k = 3)

#Quantile transform
m1 = paste("variable", seq(1,N,1), sep = ".")
rr <- lapply(1:length(m1), function(x){qtrans(hazard.stk, m1[[x]])})
rr <- stack(rr);plot(rr)

# Average across metrics
index <- rep(1:8, times=nlayers(rr)/8)#8=4 ssps * 2 periods
#mean (MN)
mn_stk <- stackApply(rr, indices=index, fun=median, na.rm=TRUE)
names(mn_stk) <- c("hzd.mn.ssp126.2050","hzd.mn.ssp126.2100",
                   "hzd.mn.ssp245.2050","hzd.mn.ssp245.2100",
                   "hzd.mn.ssp370.2050","hzd.mn.ssp370.2100",
                   "hzd.mn.ssp585.2050","hzd.mn.ssp585.2100")
plot(mn_stk)
#Standard Deviation (SD)
sd_stk <- stackApply(rr, indices=index, fun=sd, na.rm = TRUE)
names(sd_stk) <- c("hzd.sd.ssp126.2050","hzd.sd.ssp126.2100","hzd.sd.ssp245.2050",
                   "hzd.sd.ssp245.2100","hzd.sd.ssp370.2050","hzd.sd.ssp370.2100",
                   "hzd.sd.ssp585.2050","hzd.sd.ssp585.2100")
plot(sd_stk)

# climdata <- climdata %>% mutate(dplyr::across(4:94,~ normalize(inormal2(.x))))
climdata <- extract(stack(mn_stk, sd_stk), wio.AOO.spdf, sp=TRUE, df=TRUE) %>% as.data.frame()
colnames(climdata)[colnames(climdata) == "coords.x1"] <- "x"
colnames(climdata)[colnames(climdata) == "coords.x2"] <- "y"

#Calculate Exposure Dimensions of the RISK model.
names(climdata)
grdArea = (25e3)^2
(climdata <- climdata %>% 
    mutate(ExpCoral = normalize(inormal(CoralExt/grdArea)),
           ExpSeagrass = normalize(inormal(seagrassExt/grdArea)), 
           ExpCrop = normalize(inormal(Cropland/grdArea)), 
           ExpMangrove = normalize(inormal(mangroveExt/grdArea)),
           ExpFRic = ((FRic)),
           ExpFDiv = ((FDiv)),
           ExpEve = ((FEve)),
           ExpNObs = normalize((Nb_sp)))%>%
    rowwise() %>%
    mutate(land.sea.mn.exp = mean(c(ExpCoral,ExpSeagrass,ExpCrop,ExpMangrove,ExpFRic,ExpFDiv,ExpEve,ExpNObs), na.rm = TRUE),
           land.sea.sd.exp = sd(c(ExpCoral,ExpSeagrass,ExpCrop,ExpMangrove,ExpFRic,ExpFDiv,ExpEve,ExpNObs), na.rm = TRUE)) %>% ungroup()%>%
    
    mutate(#land.sea.sd.exp = ifelse(is.na(land.sea.sd.exp), 1, land.sea.sd.exp),
           #land.sea.sd.exp = ifelse(is.na(land.sea.sd.exp), agrmt::minnz(na.omit(land.sea.sd.exp)), land.sea.sd.exp),
           wgts.Exposure = (land.sea.sd.exp/land.sea.mn.exp)^-1,
           wgts.ssp126.2050 = (hzd.sd.ssp126.2050/hzd.mn.ssp126.2050)^-1,
           wgts.ssp126.2100 = (hzd.sd.ssp126.2100/hzd.mn.ssp126.2100)^-1,
           wgts.ssp245.2050 = (hzd.sd.ssp245.2050/hzd.mn.ssp245.2050)^-1,
           wgts.ssp245.2100 = (hzd.sd.ssp245.2100/hzd.mn.ssp245.2100)^-1,
           wgts.ssp370.2050 = (hzd.sd.ssp370.2050/hzd.mn.ssp370.2050)^-1,
           wgts.ssp370.2100 = (hzd.sd.ssp370.2100/hzd.mn.ssp370.2100)^-1,
           wgts.ssp585.2050 = (hzd.sd.ssp585.2050/hzd.mn.ssp585.2050)^-1,
           wgts.ssp585.2100 = (hzd.sd.ssp585.2100/hzd.mn.ssp585.2100)^-1) 
  )

(climdata <- climdata %>%
    mutate(imp.ssp126.2050 = ((hzd.mn.ssp126.2050*wgts.ssp126.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp126.2050+wgts.Exposure),
           imp.ssp126.2100 = ((hzd.mn.ssp126.2100*wgts.ssp126.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp126.2100+wgts.Exposure),
           imp.ssp245.2050 = ((hzd.mn.ssp245.2050*wgts.ssp245.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp245.2050+wgts.Exposure),
           imp.ssp245.2100 = ((hzd.mn.ssp245.2100*wgts.ssp245.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp245.2100+wgts.Exposure),
           imp.ssp370.2050 = ((hzd.mn.ssp370.2050*wgts.ssp370.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp370.2050+wgts.Exposure),
           imp.ssp370.2100 = ((hzd.mn.ssp370.2100*wgts.ssp370.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp370.2100+wgts.Exposure),
           imp.ssp585.2050 = ((hzd.mn.ssp585.2050*wgts.ssp585.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp585.2050+wgts.Exposure),
           imp.ssp585.2100 = ((hzd.mn.ssp585.2100*wgts.ssp585.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp585.2100+wgts.Exposure)
    ) %>% ungroup())
hist(climdata$imp.ssp370.2050, breaks=30)
write_rds(climdata, "2_Data/spreadsheet/1_Climate/wioAOO.metrics.rds")

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
