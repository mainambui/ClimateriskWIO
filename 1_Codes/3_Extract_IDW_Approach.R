#https://rpubs.com/dieghernan/beautifulmaps_I
library(sf)
library(dplyr)

# Projecting and cleaning
wio_cline <- st_read("C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/Data/2-Exposed systems/shp/wio_coastline_pu.shp")
st_crs(wio_cline)<- "+proj=longlat"
wio_cline <- st_transform(wio_cline, "+proj=moll")

#Read the social datasets
socioecom <- read_csv(paste0(db.dir, "3-Vulnerability/SocialDataAll.csv"))
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
wio_metrics <- read_rds(paste0(db.dir, "4-Risk Outputs/wioAOO.All.raw.rds"))
EucDist <- read_csv(paste0(db.dir, "4-Risk Outputs/wio_DistMatrix.csv"))

EucDist <- EucDist[c("src", "nbr", "EucDist")] %>% filter(EucDist>0) #filter to remove self intersections #
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
write_csv(wio.com.idw.impacts, paste0(db.dir, "4-Risk Outputs/RiskDataFINAL.csv"))



dfRisk <- wio.com.idw.impacts
#dfRisk <- read_csv("C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/Draft/RiskDataFINAL.csv")
df <- rbind(data.frame(sce = "SSP2-4.5", impact = dfRisk$imp.ssp245.2050, vulnerab = dfRisk$Sensitivity/dfRisk$AdaptiveCapacity, village = dfRisk$Villages, country = dfRisk$Country),
            data.frame(sce = "SSP3-7.0", impact = dfRisk$imp.ssp370.2050, vulnerab = dfRisk$Sensitivity/dfRisk$AdaptiveCapacity, village = dfRisk$Villages, country = dfRisk$Country))

df <- df %>% mutate(country = ifelse(country == "Kenya", "KEN", country),
                    country = ifelse(country == "Tanzania", "TZA", country),
                    country = ifelse(country == "Madagascar", "MDG", country),
                  country = ifelse(country == "Mozambique", "MOZ", country))

yR <- range(df$impact);xR <- range(df$vulnerab)
lgd <- expand.grid(x = seq(xR[1],xR[2], diff(xR)/100), 
                   y = seq(yR[1],yR[2], diff(yR)/100)) %>% 
  mutate(x1 = scales::rescale(x),
         y1 = scales::rescale(y),
         mix = rgb(y1,x1,.5))
# Q1 <- summary(df$impact)[[2]] #first quartile
# Q3 <- summary(df$impact)[[5]] #third quartile
# MdV <- median(df$vulnerab)
(plt1 <- ggplot()+
    
    # geom_rect(aes(fill = "LL", xmin = -Inf, xmax = MdV, ymin = -Inf, ymax = Q1), show.legend = FALSE)+
    # geom_rect(aes(fill = "ML", xmin = -Inf, xmax = MdV, ymin = Q1, ymax = Q3), show.legend = FALSE)+
    # geom_rect(aes(fill = "HL", xmin = -Inf, xmax = MdV, ymin = Q3, ymax = Inf), show.legend = FALSE)+
    # geom_rect(aes(fill = "LH", xmin = MdV, xmax = Inf, ymin = -Inf, ymax = Q1), show.legend = FALSE)+
    # geom_rect(aes(fill = "MH", xmin = MdV, xmax = Inf, ymin = Q1, ymax = Q3), show.legend = FALSE)+
    # geom_rect(aes(fill = "HH", xmin = MdV, xmax = Inf, ymin = Q3, ymax = Inf), show.legend = FALSE)+
    # scale_fill_manual(values = c("LL"="#595757","ML"="#806587","HL"="#A874B8",
    #                             "LH"="#468C38","MH"="#8CA289","HH"="#D3B9DB")) +
    geom_raster(data = lgd, aes(x = x, y = y, fill = mix))+
    scale_fill_identity()+
    geom_point(data = df, aes(x = vulnerab, y = impact, colour = country, shape=sce), size = 1.5) +
    labs(y = "Climate change impacts", x = "", title = "Risk Space")+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50"))+
    theme_bw(base_size = 10)+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
    #geom_vline(xintercept = median(df$vulnerab, na.rm = TRUE), linetype = 2, linewidth = .5, colour = "grey80")+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 2),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.1, 'cm'), 
          
          panel.border = element_blank()))

(plt2 <- ggplot()+
    # geom_rect(aes(fill = "LL", xmin = -Inf, xmax = MdV, ymin = -Inf, ymax = Q1), show.legend = FALSE)+
    # geom_rect(aes(fill = "ML", xmin = -Inf, xmax = MdV, ymin = Q1, ymax = Q3), show.legend = FALSE)+
    # geom_rect(aes(fill = "HL", xmin = -Inf, xmax = MdV, ymin = Q3, ymax = Inf), show.legend = FALSE)+
    # geom_rect(aes(fill = "LH", xmin = MdV, xmax = Inf, ymin = -Inf, ymax = Q1), show.legend = FALSE)+
    # geom_rect(aes(fill = "MH", xmin = MdV, xmax = Inf, ymin = Q1, ymax = Q3), show.legend = FALSE)+
    # geom_rect(aes(fill = "HH", xmin = MdV, xmax = Inf, ymin = Q3, ymax = Inf), show.legend = FALSE)+
    # scale_fill_manual(values = c("LL"="#595757","ML"="#806587","HL"="#A874B8","LH"="#468C38","MH"="#8CA289","HH"="#D3B9DB")) +
    geom_raster(data = lgd, aes(x = x, y = y, fill = mix))+
    scale_fill_identity()+
    
    geom_point(data = df, aes(x = vulnerab, y = impact), alpha = 0) +
    labs(y = "", x = "", title = "Options Space")+
    scale_y_continuous(breaks = c(.3, .43, .55), labels = c("Low", "Medium", "High") , position = "right", expand = c(0,0))+
    scale_x_continuous(breaks = c(1.02, 1.3), labels = c("Low","High"), expand = c(0,0))+
    theme_bw(base_size = 10)+
    theme(axis.text.y = element_text(angle = 90),
          axis.ticks = element_blank(), 
          panel.border = element_blank()))

library(patchwork)
plt1+plt2+plot_layout(ncol = 2)
ggsave("outputs/xx1.png", width = 5.5, height = 4, dpi = 1200)

df <- dfRisk %>% mutate(risk.ssp370.2050 = ((imp.ssp370.2050*Sensitivity)/AdaptiveCapacity),
                        risk.ssp245.2050 = ((imp.ssp245.2050*Sensitivity)/AdaptiveCapacity)) %>% dplyr::select(Country,Villages,risk.ssp370.2050,risk.ssp245.2050)
df <- rbind(data.frame(sce = "SSP2-4.5", risk = (df$risk.ssp245.2050), village = df$Villages, country = df$Country),
            data.frame(sce = "SSP3-7.0", risk = (df$risk.ssp370.2050), village = df$Villages, country = df$Country))
Q1 <- summary(df$risk)[[2]] #first quartile
Q3 <- summary(df$risk)[[5]] #third quartile
(plt2 <- ggplot(data = df)+
    geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Q1)+
    geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = Q1, ymax = Q3)+
    geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = Q3, ymax = Inf)+
    #geom_col(aes(reorder(Villages, risk.ssp370.2050), risk.ssp370.2050), width = .2, fill = "grey90")+
    geom_point(aes(x=reorder(village,-risk), y=risk, colour = country, shape = sce), size = 1.5, position = position_dodge2(width =.5))+
    #geom_text(aes(x=reorder(village,risk), y=risk, label = round(risk,2)), size = 2, colour = "black", fontface = "bold",position = position_dodge(width =.5))+
    labs(x = "Coastal communities", y = "Climate risk index", title = "(A)")+
    #scale_y_continuous(expand = c(0,0), limits = c(.2,.9), breaks = seq(.2,.9, 0.1))+
    theme_classic(base_size = 10)+
    scale_colour_manual(name="", values = c("Kenya"="darkred","Madagascar"="yellow","Mozambique"="dodgerblue4","Tanzania"="grey50"))+
    scale_fill_manual(name="", values = c("Kenya"="darkred","Madagascar"="yellow","Mozambique"="dodgerblue4","Tanzania"="grey50"))+
    theme(legend.position = "", 
          legend.background = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

df <- dfRisk %>% 
  mutate(risk.ssp370.2050 = ((imp.ssp370.2050*Sensitivity)/AdaptiveCapacity),
         risk.ssp245.2050 = ((imp.ssp245.2050*Sensitivity)/AdaptiveCapacity)) %>% 
  group_by(Country) %>% 
  summarise(mn.370 = summary(risk.ssp370.2050)[[3]],
            q1.370 = summary(risk.ssp370.2050)[[2]],
            q3.370 = summary(risk.ssp370.2050)[[5]],
            
            mn.245 = summary(risk.ssp245.2050)[[3]],
            q1.245 = summary(risk.ssp245.2050)[[2]],
            q3.245 = summary(risk.ssp245.2050)[[5]]) %>% ungroup()

df <- rbind(data.frame(sce = "SSP2-4.5", MN = df$mn.245, Q1 = df$q1.245, Q3 = df$q3.245, country = df$Country),
            data.frame(sce = "SSP3-7.0", MN = df$mn.370, Q1 = df$q1.370, Q3 = df$q3.370, country = df$Country)) %>% group_by (country) %>% mutate(sortMag = mean(MN))
(plt3 <- ggplot(data = df)+ 
    geom_pointrange(aes(x = reorder(country, sortMag), y = MN, ymin = Q1, ymax = Q3, colour=country, shape = sce),
                    size=.1, linewidth = .2, position = position_dodge(width =.5))+
    scale_colour_manual(name="", values = c("Kenya"="darkred","Madagascar"="yellow","Mozambique"="dodgerblue4","Tanzania"="grey50"))+
    labs(y = "Climate risk index", x="", title = "(B)")+
    #scale_y_continuous(expand = c(0,0), limits = c(.2,.9), breaks = seq(.2,.9, 0.1), position = "right")+
    theme_classic(base_size = 10)+
    theme(legend.position = "", 
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          panel.border = element_blank(), 
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))
#ggsave("outputs/3_Nations.png", width = 4, height = 4, dpi = 1200)

#plt1 + annotation_custom(ggplotGrob(plt2), xmin = 1, xmax = 14, ymin = .5, ymax = .79)
#((plt1+plt3+plot_layout(ncol = 2, widths = c(2,1)))/plt2)+plot_layout(ncol = 1, heights = c(2,1))
plt2+plt3+plot_layout(ncol = 2, widths = c(5,1))
ggsave("outputs/3_RiskNew.png", dpi = 1200, height = 3, width = 8)




df <- dfRisk %>% 
  mutate(risk.ssp370.2050 = ((imp.ssp370.2050*Sensitivity)/AdaptiveCapacity),
         risk.ssp245.2050 = ((imp.ssp245.2050*Sensitivity)/AdaptiveCapacity))
df %>% ggplot() + geom_point(aes(x = (TEV/1e6), y = risk.ssp245.2050))
df %>% ggplot() + geom_point(aes(x = (TEV/1e6), y = imp.ssp245.2050))


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
socioecom <- read_csv(paste0(db.dir, "3-Vulnerability/SocialDataAll.csv"))
socioecom <- filter(socioecom, !is.na(x))
#Convert to spatial object
spdf <- st_as_sf(socioecom, coords=c('x', 'y'), crs="+proj=longlat")
spdf <- st_transform(spdf, "+proj=moll")
spdf_buff <- st_buffer(spdf, dist = 100000)
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


