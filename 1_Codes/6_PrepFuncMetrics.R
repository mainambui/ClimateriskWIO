
#Import Steph's corals and fish data
Coral_Data_finalV4 <- readRDS("C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/Data/2-Exposed systems/Steph/Coral_Data_finalV4.rds")
names(Coral_Data_finalV4)
coralDf <- Coral_Data_finalV4[,c("ID","Lat","Lon","Nb_sp","Coral","Mangrove","Seagrass","FRic","FDiv","FEve","FDis","FSpe","FOri","Coral_coun")]
plot(coralDf$Lon, coralDf$Lat)

wio.ccva.pu <-  st_read(paste0(db.dir, "2-Exposed systems/shp/wio_coastline_pu.shp"))
st_crs(wio.ccva.pu)<-"+proj=longlat"
steph.coral.df <- st_as_sf(coralDf, coords=c('Lon', 'Lat'), crs="+proj=longlat")
wio.ccva.pu.corals <- st_intersection(wio.ccva.pu, steph.coral.df)
# st_write(steph.coral.df, "wio.grd.shp")
# st_write(steph.coral.df, "steph.coral.df.shp")
#Write to local disk and use ArcGIS Pro to intersect. Its extremely faster

#Inport intersection outputs from ArcGIS...
wio.ccva.corals <-  read_csv(paste0(db.dir, "2-Exposed systems/StephsMetricsIntesectsWIOext.csv"))
wio.ccva.corals <- as.data.frame(wio.ccva.corals) %>% group_by(ID) %>%
  summarise(Nb_sp = median(Nb_sp, na.rm = TRUE), 
            FRic = median(FRic, na.rm = TRUE),
            FDiv = median(FDiv, na.rm = TRUE),
            FEve = median(FEve, na.rm = TRUE),
            FDis = median(FDis, na.rm = TRUE),
            FSpe = median(FSpe, na.rm = TRUE),
            FOri = median(FOri, na.rm = TRUE))

#Merge fish data to GRID INDEX FEATURE
wio.AOO.all <- readRDS(paste0(db.dir, "2-Exposed systems/wioAOO.rds"))
wio.AOO.all <- merge(wio.AOO.all, wio.ccva.corals, by = "ID", all.x=TRUE)
wio.AOO.all <-wio.AOO.all %>% dplyr::select(-c(PageNumber,coral.presence, seagrass.presence,mangrove.presence, Count))

write_rds(wio.AOO.all, paste0(db.dir, "2-Exposed systems/wioAOO.rds"))
