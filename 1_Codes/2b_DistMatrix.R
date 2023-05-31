#https://rpubs.com/dieghernan/beautifulmaps_I
library(tidyverse);library(raster);library(sf);library(sp)

rm(list = ls())
# Projecting and cleaning
wio_cline <- st_read("2_Data/shp/wio_coastline_pu.shp");st_crs(wio_cline)<- "+proj=longlat"
wio_cline <- st_transform(wio_cline, "+proj=moll")

#Read the social datasets
socioecom <- read_csv("2_Data/sheet/3_SocialVulnerability/SocialDataAll.csv")
socioecom <- filter(socioecom, !is.na(x))
#Convert to spatial object
spdf <- st_as_sf(socioecom, coords=c('x', 'y'), crs="+proj=longlat")
spdf <- st_transform(spdf, "+proj=moll")

#Intersect and assign Grid Index IDs to the community level spatial points
spdf <- st_intersection(wio_cline, spdf)
socID <- st_drop_geometry(spdf) %>% mutate(ID = PageName)%>% dplyr::select(ID, Country, Villages)
write_csv(socID, "2_Data/sheet/VillageWithGridIDs.csv")


#Convert the Fishnet to
df2 <- cbind.data.frame("ID" = wio_cline$PageName, st_coordinates(st_centroid(st_geometry(wio_cline))))
colnames(df2)[colnames(df2)=="X"]<-"x"
colnames(df2)[colnames(df2)=="Y"]<-"y"
wio.gridNetworks <- expand.grid(src=df2$ID, nbr=df2$ID) %>% filter(src %in% socID$ID)

colnames(df2)[colnames(df2)=="ID"]<-"src"
wio.gridNetworks <- inner_join(wio.gridNetworks, df2, by="src")
colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="x"]<-"src_x"
colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="y"]<-"src_y"

colnames(df2)[colnames(df2)=="src"]<-"nbr"
wio.gridNetworks <- inner_join(wio.gridNetworks, df2, by="nbr")
colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="x"]<-""
colnames(wio.gridNetworks)[colnames(wio.gridNetworks)=="y"]<-"nbr_y"
#REMAIN it back to ID
colnames(df2)[colnames(df2)=="nbr"]<-"ID"

library(gdistance)
library(doParallel)
n.cores <- 6
#Compute Euclidean distances between RCPs
system.time({
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  clusterExport(cl, c('wio.gridNetworks'))  ## Export the environment variables to each cluster
  clusterEvalQ(cl, {
    library(gdistance)
    library(sf)}
  )  ## Load the library "gdistance" to each cluster

  the.edist.function <- function(k) {
    eucDist <- st_distance(st_as_sf(wio.gridNetworks[k,], coords= c('src_x', 'src_y'), crs="+proj=moll"),
                           st_as_sf(wio.gridNetworks[k,], coords=c('nbr_x', 'nbr_y'), crs="+proj=moll"))
    return(eucDist)}
  eucDist <- clusterApplyLB(cl, x=1:nrow(wio.gridNetworks), fun=the.edist.function)
})

stopCluster(cl)

wio.gridNetworks$EucDist <- as.vector(eucDist %>% do.call("rbind", .))
hist(wio.gridNetworks$EucDist, breaks=50)
write_csv(wio.gridNetworks, "wio_DistMatrix.csv")



