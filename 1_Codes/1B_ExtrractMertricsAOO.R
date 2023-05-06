##PREPARE TO EXTRACT DATA TO AOO
library(SearchTrees)
hzd.CDD <- resample(hzd.CDD, hzd.TAP) #align extent with the risk of the stack
hazard.stk <- stack(hzd.SST,hzd.pH,hzd.CDD,hzd.TAP,hzd.EVAP,hzd.NPP)
index <- rep(1:8, times=nlayers(hazard.stk)/8)#8 = 4 ssps * 2 periods
#mean (MN)
mn_stk <- stackApply(hazard.stk, indices=index, fun=median, na.rm = TRUE)
names(mn_stk) <- c("hzd.mn.ssp126.2050","hzd.mn.ssp126.2100",
                   "hzd.mn.ssp245.2050","hzd.mn.ssp245.2100",
                   "hzd.mn.ssp370.2050","hzd.mn.ssp370.2100",
                   "hzd.mn.ssp585.2050","hzd.mn.ssp585.2100")
plot(mn_stk)
#Standard Deviation (SD)
sd_stk <- stackApply(hazard.stk, indices=index, fun=sd, na.rm = TRUE)
names(sd_stk) <- c("hzd.sd.ssp126.2050","hzd.sd.ssp126.2100","hzd.sd.ssp245.2050",
                   "hzd.sd.ssp245.2100","hzd.sd.ssp370.2050","hzd.sd.ssp370.2100",
                   "hzd.sd.ssp585.2050","hzd.sd.ssp585.2100")
plot(sd_stk)

#EXTRACT DATA TO AOO KEEPING IN MIND THE GAPS AT THE LAND SEA INTERFACE 
hzd.stk <- stack(mn_stk, sd_stk)
hzd.stk <- as.data.frame(hzd.stk, xy = T, na.rm=TRUE)
hzd.stree <- createTree(coordinates(hzd.stk[,c("x","y")]))
inds <- knnLookup(hzd.stree, newdat=coordinates(wio.AOO[, c("x","y")]), k=1)
clm.hzd.pts <- hzd.stk[inds,]

climdata.All <- cbind(wio.AOO, clm.hzd.pts[,3:ncol(clm.hzd.pts)])
# N <- ncol(climdata.All)
# climdata.All <- climdata.All %>% mutate(dplyr::across(4:94,~ normalize(inormal2(.x))))
write_rds(climdata.All, paste0(db.dir, "1-Climate metrics/wioAOO.clim.raw.rds"))

#Import and add Stephs metrics to the data
climdata.All <- readRDS(paste0(db.dir, "1-Climate metrics/wioAOO.clim.raw.rds"))
#Calculate Exposure Dimensions of the RISK model.
names(climdata.All)
grdArea = (25e3)^2
(climdata.All <- climdata.All %>% 
    mutate(ExpCoral = normalize(inormal2(CoralExt)),
           ExpSeagrass = normalize(inormal2(seagrassExt)), 
           ExpCrop = normalize(inormal2(Cropland)), 
           ExpMangrove = normalize(inormal2(mangroveExt)),
           ExpFRic = normalize(inormal2(FRic)),
           ExpFDiv = normalize(inormal2(FDiv)),
           ExpEve = normalize(inormal2(FEve)),
           ExpRich = normalize(inormal2(Nb_sp)))%>%
    rowwise() %>%
    mutate(land.sea.mn.exp = mean(c(ExpCoral,ExpSeagrass,ExpCrop,ExpMangrove,ExpFRic,ExpFDiv,ExpEve), na.rm = TRUE),
           land.sea.sd.exp = sd(c(ExpCoral,ExpSeagrass,ExpCrop,ExpMangrove,ExpFRic,ExpFDiv,ExpEve), na.rm = TRUE),
           #land.sea.sd.exp = ifelse(is.na(land.sea.sd.exp), 1, land.sea.sd.exp),
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
  %>% ungroup())

(climdata.All <- climdata.All %>%
    mutate(imp.ssp126.2050 = ((hzd.mn.ssp126.2050*wgts.ssp126.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp126.2050*wgts.Exposure),
           imp.ssp126.2100 = ((hzd.mn.ssp126.2100*wgts.ssp126.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp126.2100*wgts.Exposure),
           imp.ssp245.2050 = ((hzd.mn.ssp245.2050*wgts.ssp245.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp245.2050*wgts.Exposure),
           imp.ssp245.2100 = ((hzd.mn.ssp245.2100*wgts.ssp245.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp245.2100*wgts.Exposure),
           imp.ssp370.2050 = ((hzd.mn.ssp370.2050*wgts.ssp370.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp370.2050*wgts.Exposure),
           imp.ssp370.2100 = ((hzd.mn.ssp370.2100*wgts.ssp370.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp370.2100*wgts.Exposure),
           imp.ssp585.2050 = ((hzd.mn.ssp585.2050*wgts.ssp585.2050)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp585.2050*wgts.Exposure),
           imp.ssp585.2100 = ((hzd.mn.ssp585.2100*wgts.ssp585.2100)*(land.sea.mn.exp*wgts.Exposure))/(wgts.ssp585.2100*wgts.Exposure)
    ) %>% ungroup())

write_rds(climdata.All, paste0(db.dir, "4-Risk Outputs/wioAOO.All.raw.rds"))
#writexl::write_xlsx(climdata.All, "ClimData.xlsx")
#Import polygons of countries bordering the WIO
#limit <- max(abs(climdata$tx90p_thresh), na.rm = TRUE) * c(-1, 1)
climdata.All %>% 
  filter(is.finite(imp.ssp126.2050))%>% 
  ggplot()+ geom_sf(data = focusISO3, colour = "grey70", fill = "grey95")+
  geom_raster(aes(x,y,fill=imp.ssp126.2050))+
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
