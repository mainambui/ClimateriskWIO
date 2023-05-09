#https://rpubs.com/dieghernan/beautifulmaps_I
library(sf)
library(dplyr)

# Projecting and cleaning
wio_cline <- st_read("C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/Data/2-Exposed systems/shp/wio_coastline_pu.shp")
st_crs(wio_cline)<- "+proj=longlat"
wio_cline <- st_transform(wio_cline, "+proj=moll")

#Read the social datasets
socioecom <- read_csv("2_Data/spreadsheet/3_SocialVulnerability/SocialDataAll.csv")
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
wio_metrics <- read_rds("2_Data/spreadsheet/1_Climate/wioAOO.metrics.rds")
EucDist <- read_csv("2_Data/spreadsheet/wio_DistMatrix.csv")

EucDist <- EucDist[c("src", "nbr", "EucDist")] %>% filter(EucDist>0,EucDist<1e5) #filter to remove self intersections #
colnames(wio_metrics)[colnames(wio_metrics)=="ID"]<-"nbr"
EucDist <- merge(EucDist, wio_metrics, by = "nbr")
EucDist$inverseDist <- (1/EucDist$EucDist)
#Visual chekcing of distribution
hist(EucDist$inverseDist, breaks = 30)
names(EucDist)
(idw.impacts <- EucDist %>% group_by(src) %>%
  summarise(
    imp.ssp126.2050 = sum(imp.ssp126.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    imp.ssp126.2100 = sum(imp.ssp126.2100*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    imp.ssp245.2050 = sum(imp.ssp245.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    imp.ssp245.2100 = sum(imp.ssp245.2100*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    imp.ssp370.2050 = sum(imp.ssp370.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    imp.ssp370.2100 = sum(imp.ssp370.2100*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    imp.ssp585.2050 = sum(imp.ssp585.2050*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    imp.ssp585.2100 = sum(imp.ssp585.2100*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE),
    
    TEV = sum(TEV*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE)
  )%>% mutate(ID = src) %>% dplyr::select(-src)
)

(wio.com.idw.impacts <- merge(data.frame(spdf) %>% mutate(ID = PageName) %>% dplyr::select(-c(PageName, PageNumber, geometry)),
                              idw.impacts, by = "ID"))

plot(wio.com.idw.impacts$TEV/1e6, wio.com.idw.impacts$imp.ssp370.2050)
write_csv(wio.com.idw.impacts, "2_Data/spreadsheet/RiskDataFINAL.csv")

# initial <- wio_cline
# initial$index_target <- 1:nrow(initial)
# target <- st_geometry(initial)
# 
# grid <- st_make_grid(target,
#                      25*1000,
#                      crs = st_crs(initial),
#                      what = "polygons",
#                      square = TRUE
# )
# 
# # To sf
# grid <- st_sf(index = 1:length(lengths(grid)), grid) # Add index
# 
# # We identify the grids that belongs to a entity by assessing the centroid
# cent_grid <- st_centroid(grid)
# cent_merge <- st_join(cent_grid, initial["index_target"], left = F)
# grid_new <- inner_join(grid, st_drop_geometry(cent_merge))
# 
# # Fishnet
# Fishgeom <- aggregate(grid_new,
#                       by = list(grid_new$index_target),
#                       FUN = min,
#                       do_union = FALSE
# )
# 
# # Lets add the df
# Fishnet <- left_join(
#   Fishgeom %>% dplyr::select(index_target),
#   st_drop_geometry(initial)
# ) %>%
#   dplyr::select(-index_target)
# 
# plot(Fishnet["PageNumber"])


#ExTRACT RASTER VALUES TO POINTS
#Import SE data
# socioecom <- read_csv(paste0(db.dir, "3-Vulnerability/SocialDataAll.csv"))
# socioecom <- filter(socioecom, !is.na(x))
# #Convert to spatial object
# spdf <- st_as_sf(socioecom, coords=c('x', 'y'), crs="+proj=longlat")
# spdf <- st_transform(spdf, "+proj=moll")
# spdf_buff <- st_buffer(spdf, dist = 100000)
# (dfRisk <- raster::extract(cmpd.Hzrd, spdf_buff, weights=TRUE, na.rm=TRUE, df=TRUE) %>% 
#     as.data.frame() %>% group_by(ID) %>% 
#     summarise(mn.spp126.2050 = matrixStats::weightedMedian(climate.imp.spp126.2050, weight, na.rm = TRUE),
#               sd.spp126.2050 = sd(climate.imp.spp126.2050, na.rm = TRUE),
#               
#               mn.spp126.2100 = matrixStats::weightedMedian(climate.imp.ssp126.2100, weight, na.rm = TRUE),
#               sd.spp126.2100 = sd(climate.imp.ssp126.2100, na.rm = TRUE),
#               
#               mn.spp245.2050 = matrixStats::weightedMedian(climate.imp.ssp245.2050, weight, na.rm = TRUE),
#               sd.spp245.2050 = sd(climate.imp.ssp245.2050, na.rm = TRUE),
#               
#               mn.spp245.2100 = matrixStats::weightedMedian(climate.imp.ssp245.2100, weight, na.rm = TRUE),
#               sd.spp245.2100 = sd(climate.imp.ssp245.2100, na.rm = TRUE),
#               
#               mn.spp370.2050 = matrixStats::weightedMedian(climate.imp.ssp370.2050, weight, na.rm = TRUE),
#               sd.spp370.2050 = sd(climate.imp.ssp370.2050, na.rm = TRUE),
#               
#               mn.spp370.2100 = matrixStats::weightedMedian(climate.imp.ssp370.2100, weight, na.rm = TRUE),
#               sd.spp370.2100 = sd(climate.imp.ssp370.2100, na.rm = TRUE)) %>% 
#     cbind(., socioecom))
# 
# writexl::write_xlsx(dfRisk, "WIO_CCVA_COMDATA.xlsx")


