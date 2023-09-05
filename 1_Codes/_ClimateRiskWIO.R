#PREPARE TO EXTRACT DATA TO AOO

rm(list = ls())

library("dplyr")
library("tidyterra")
library("terra")
library("sf")
library("pracma")

# #import some important functions
# inormal <- function(x) {
#   qrank <- ((rank(x, na.last = TRUE) - 0.5) / sum(!is.na(x)))
#   z_score <- scales::rescale(sqrt(2)*pracma::erfinv(2*qrank-1))
#   return(z_score)  }
# # 
# #Import PU and convert to a spatial object
# wio.aoo <- readRDS("2_Data/sheet/wio.aoo.rds")
# grdSize <- (25*25*1e4)
# wooFilter <- as.vector(wio.aoo %>% mutate(across(CoralExt:Cropland, ~ .x/grdSize)) %>% rowwise() %>%
#                          mutate(Tot = 100*sum(c(CoralExt,seagrassExt,mangroveExt,Cropland), na.rm = TRUE)) %>% ungroup() %>%
#                          filter(Tot > 0.1)%>% dplyr::select(ID))
# 
# wio.aoo.sub <- filter(wio.aoo, ID %in% wooFilter$ID)
# wio.aoo.spdf <- st_as_sf(wio.aoo.sub, coords=c('x', 'y'), crs="+proj=longlat") %>% vect()
# plot(wio.aoo.spdf)
# 
# #LOAD METRICS
# (clim.nc <- list.files("./2_Data/raster",pattern='*.nc',full.names=TRUE))
# 
# #Note that each NC file contains eight layers which is generally structured as SSP126_2050, SSP126_2100, SSP245_2050, SSP245_2100, SSP370_2050, SSP370_2100, SSP585_2050, SSP585_2100
# nc <- (as.data.frame(expand.grid(x=c(126,245,370,585), y=c(2050,2100))) %>% arrange(desc(-x)) %>% mutate(cbn = paste(x,y,sep = "_")) %>% dplyr::select(cbn))[,1]
# chonic <- c("evspsbl","npp","pH_trend","tap_trend","sst_trend","ts_trend","SLR")
# acute <- c("cdd","r10p","sst90p","ts90p")
# (varlst <- c(chonic,acute))
# 
# hazards <- lapply(1:length(varlst), function(x){
#   rr <- rast(clim.nc[grep(varlst[[x]], clim.nc)])
#   names(rr) <- paste(varlst[[x]],nc, sep = "_")
#   return(rr)}
# ) %>% rast()
# 
# plot(hazards)
# 
# #Put in a function to tidy up the global environment
# QTtransform <- function(r,vlst,scenario,time){
#   df <- as.data.frame(r, xy =TRUE)
#   df$ID <- seq.int(1:nrow(df))
#   nlist <- names(df)
# 
#   xx = c()
#   for(k in seq_along(time)){
#     xx[[k]] <- lapply(1:length(scenario), function(x){
#       d1 <- cbind(df[,c("ID","x","y")], sce = paste("ssp",scenario[[x]],sep=""), df[,nlist[grep(paste(scenario[[x]],time[[k]],sep="_"),nlist)]])
#       colnames(d1) <- c("ID","x","y","Scenario",paste(vlst))
#       return(d1)}) %>% bind_rows()  }
# 
#   names(xx) <- paste(time)
#   xx <- bind_rows(xx, .id = "Period")
# 
#   #Quantile transform
#   inormal <- function(x){
#     qrank <- ((rank(x, na.last = TRUE) - 0.5) / sum(!is.na(x)))
#     z_score <- scales::rescale(sqrt(2)*pracma::erfinv(2*qrank-1))
#     return(z_score)}
# 
#   N <- ncol(xx)
#   dfNmd <- xx %>% as.data.frame() %>%
#     mutate(pH_trend = -1*pH_trend,evspsbl = -1*evspsbl) %>% 
#     #Reserve some key variables before normalising.
#     #reduced primary production vs. reduced in annual rainfall (drought), negative trend in Ph means ocean becoming acidic over time
#     mutate(across(all_of(6:N), ~ scales::rescale(.x)))
# 
#   #Convert to wide format
#   dfWide <- lapply(1:length(vlst), function(x){
#     xv <- reshape2::dcast(dfNmd,ID+x+y~Scenario+Period, value.var= vlst[[x]])
#     scePeriod <- rep(scenario, each = length(time))
#     names(xv) <- c("ID","x","y",paste(rep(vlst[[x]], each = length(scePeriod)), scePeriod, time, sep = "_"))
#     return(xv)})
#   stdHzd <- dfWide %>% purrr::reduce(left_join, by = c("ID","x","y"))
# 
#   #Spatialize
#   vnames <- colnames(stdHzd[,-c(1,2,3)])
#   stdHzdspdf <- lapply(1:length(vnames), function(x){
#     rr <- rast(cbind(stdHzd[c("x","y")], stdHzd[,vnames[[x]]]),
#                type="xyz")
#   })
#   stdHzdspdf <- rast(stdHzdspdf)
#   names(stdHzdspdf) <- paste(vnames)
#   crs(stdHzdspdf) <- "+proj=longlat"
#   return(stdHzdspdf)
# }
# 
# #Inverse normal standardised variables
# (HazardQNormed <- QTtransform(hazards, vlst = varlst, time = c(2050), scenario = c(245,370,585)))
# plot(HazardQNormed[[1:8]])
# 
# #Extract the raw hazards to the WIO's AOO of interest
# # set values below 100 to NA.
# climdata <- terra::extract(HazardQNormed, wio.aoo.spdf, xy=TRUE) %>% as.data.frame()
# climdata <- cbind(as.data.frame(wio.aoo.spdf), climdata[,-1])
# climdata <- climdata %>% relocate(c(x,y), .before = ID)
# N <- ncol(climdata)
# #climdata <- cbind(climdata[1:14], DMwR2::knnImputation(climdata[15:N], k = 3))
# 
# #Import tropical cyclone data.
# slr.tc.data <- rast("C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/3_WIO_SYNTHESIS/Data/1-Climate metrics/TC/TC_count.tif")
# slr.tc.data <- app(slr.tc.data, function(x) (scales::rescale(x)))
# slr.tc.data <- cbind(ID = wio.aoo.spdf$ID, (terra::extract(slr.tc.data, wio.aoo.spdf) %>% as.data.frame())[-1])
# colnames(slr.tc.data)[colnames(slr.tc.data)=="lyr.1"] <- "TC"
# 
# slr.tc.data$tc_245_2050 <- scales::rescale(slr.tc.data$TC)
# slr.tc.data$tc_370_2050 <- scales::rescale(slr.tc.data$TC)
# slr.tc.data$tc_585_2050 <- scales::rescale(slr.tc.data$TC)
# slr.tc.data <- slr.tc.data[,c(1,3:5)] #Select only TC
# 
# climdata <- merge(climdata, slr.tc.data, by = "ID")
# 
# #Normalised exposed systems metrics
# climdata <- climdata %>%
#   #Create copies of the Exposed systems. the original versions will be needed later
#   mutate(std_corals = CoralExt/grdSize,
#          std_seagrass = seagrassExt/grdSize,
#          std_mangrove = mangroveExt/grdSize,
#          std_cropcover = Cropland/grdSize,
#          std_Nb_sp = Nb_sp,
#          std_FRic = FRic,
#          std_FDiv = FDiv,
#          std_FEve = FEve
#   )%>%
#   mutate(across(std_corals:std_Nb_sp,~ scales::rescale((.x))))#FRic, FDiv, and FEve are already normalised variables
# saveRDS(climdata, "2_Data/sheet/all.climdata.rds")

#AGGREGATION NEXT
all.data <- readRDS("2_Data/sheet/all.climdata.rds")
namelist <- colnames(all.data)

# the Choquet Integral of f w.r.t mu
library(kappalab)
mu <- readRDS("2_Data/sheet/mu.rds");summary(mu) #Overall mu
mu.ssp245 <- readRDS("2_Data/sheet/mu.ssp245.rds");summary(mu.ssp245)
mu.ssp370 <- readRDS("2_Data/sheet/mu.ssp370.rds");summary(mu.ssp370)
mu.ssp585 <- readRDS("2_Data/sheet/mu.ssp585.rds");summary(mu.ssp585)

#Choquet.integral(mu,f)
exposure <- c("std_corals","std_seagrass","std_mangrove","std_cropcover","std_Nb_sp","std_FRic","std_FDiv","std_FEve")
(impacts.aoo <- lapply(1:length(exposure),function(i){
  df <- cbind(ID=all.data[,"ID"], all.data[,exposure[[i]]]*all.data[,namelist[grep("_2050",namelist)]])
  n <- ncol(df)
  df <- df %>% 
      mutate(dplyr::across(all_of(2:n),~ ifelse(is.na(.x),0,.x))) %>% rowwise() %>%
      mutate(ssp245.2050 = kappalab::Choquet.integral(mu,c_across(namelist[grep("_245_2050", namelist)])),
             ssp370.2050 = kappalab::Choquet.integral(mu,c_across(namelist[grep("_370_2050", namelist)])),
             ssp585.2050 = kappalab::Choquet.integral(mu,c_across(namelist[grep("_585_2050", namelist)]))) %>% ungroup() %>%
      dplyr::select(ID,ssp245.2050,ssp370.2050,ssp585.2050)
  names(df) <- c("ID",paste(exposure[[i]], c("ssp245.2050","ssp370.2050","ssp585.2050"), sep = "_"))
  return(df)
})
)

impacts.aoo <- impacts.aoo %>% purrr::reduce(left_join, by = "ID")

fuzzyS <- function(x){return(1 - prod((1 - x)))}
(impacts.agg <- impacts.aoo %>% rowwise() %>%
  mutate(imp.ssp245.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp245", colnames(impacts.aoo))])),
         imp.ssp370.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp370", colnames(impacts.aoo))])),
         imp.ssp585.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp585", colnames(impacts.aoo))]))) %>%
  dplyr::select(ID,imp.ssp245.2050,imp.ssp370.2050,imp.ssp585.2050))

impacts.agg <- merge(impacts.agg, all.data[,c("ID","CoralExt","seagrassExt","mangroveExt","Cropland")], by = "ID")

#########################################################################################################################################
#Merge grid level impacts to the network
#########################################################################################################################################
idw.matrix <- read.csv("2_Data/sheet/idw.dist.matrix.csv")
village.aoo.id <- read.csv("2_Data/sheet/village.grid.id.csv")

idw.matrix <- idw.matrix[c("src", "nbr", "EucDist")] %>% filter(EucDist>0) #filter to remove self intersections #
idw.matrix$inverseDist <- (1/(idw.matrix$EucDist))

colnames(impacts.agg)[colnames(impacts.agg)=="ID"]<-"nbr"
idw.matrix <- merge(idw.matrix, impacts.agg, by = "nbr")

#Visual checking of distribution
hist(idw.matrix$inverseDist, breaks = 30)
(idw.impacts <- idw.matrix %>% group_by(src) %>%
    summarise(across(imp.ssp245.2050:Cropland, ~(sum(.x*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE))))%>% mutate(ID = src) %>% dplyr::select(-src)
)
(villageImpacts <- merge(village.aoo.id, idw.impacts, by = "ID"))

######################################################################################################################################################
# NOW ANALYSE VILLAGE LEVEL RISK
#####################################################################################################################################################
#Import socio-economic data

socioecom <- read.csv("2_Data/sheet/3_SocialVulnerability/SocialDataAll.csv")
socioecom <- socioecom %>% mutate(ISO3 = Country,
                                  ISO3 = ifelse(ISO3 == "Kenya", "KEN", ISO3),
                                  ISO3 = ifelse(ISO3 == "Tanzania", "TZA", ISO3),
                                  ISO3 = ifelse(ISO3 == "Madagascar", "MDG", ISO3),
                                  ISO3 = ifelse(ISO3 == "Mozambique", "MOZ", ISO3))
I_Ctrl <- read.csv("2_Data/sheet/impact.control.csv")
socioecom <- merge(socioecom, I_Ctrl, by.x = "ISO3")
plot(socioecom$AdaptiveCapacity, exp(-1*socioecom$ic2020))

library(ggthemes)
library(ggrepel)
socioecom$VillNation <- paste(socioecom$Villages, paste("(",socioecom$ISO3,")", sep = ""))
socioecom$Vulnerable <- scales::rescale((socioecom$Sensitivity/socioecom$AdaptiveCapacity)*exp(-1*socioecom$ic2020), to = c(0.01, 1))
# ggplot(data = socioecom, aes(x="XS", y=Vulnerable,label=VillNation))+
#   geom_boxplot(linewidth = 0.3)+
#   geom_point(aes(colour = ISO3), position=position_jitter(width=.1, height=0))+
#   scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey"))+
#   geom_text_repel(size = 1)+
#   scale_y_continuous(name = "Social vulnerability", expand = c(0,0), limits = c(0,1.02),breaks = seq(0,1,.2))+ labs(x = "")+
#   theme_classic(base_size = 12)+
#   guides(shape="none", colour = "none")+
#   theme(legend.position = "none",
#         legend.title = element_blank(),
#         legend.background = element_rect(fill = NA),
#         legend.key.size = unit(1, 'cm'),
#         axis.text.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.text = element_text(size = 8),
#         legend.key.height = unit(.1, 'cm'),
#         legend.key.width = unit(.2, 'cm'), 
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.line.y = element_line(linewidth = .1),
#         axis.ticks.y = element_line(linewidth = .1)
#         )
#ggsave("3_Outputs/plots/FigS2.png", width = 4, height = 5, dpi = 1200)

##################################################################################################################
#                                  PLOT RISK SPACES
##################################################################################################################
#Import climate data
riskMaster <- merge(socioecom, villageImpacts, by ="Villages")

plot(riskMaster[, c("imp.ssp585.2050","Livelihood","Demography","Cultural","Health","Learning","Assets","Flexibility","Agency","Organisation")])
df <- rbind(data.frame(sce = "SSP2-4.5", 
                       impact = riskMaster$imp.ssp245.2050, 
                       Vulnerability = riskMaster$Vulnerable, 
                       village = riskMaster$Villages, 
                       ISO3 = riskMaster$ISO3),
            data.frame(sce = "SSP5-8.5", 
                       impact = riskMaster$imp.ssp585.2050, 
                       Vulnerability = riskMaster$Vulnerable, 
                       village = riskMaster$Villages, 
                       ISO3 = riskMaster$ISO3))

yR <- range(df$impact);xR <- range(df$Vulnerability)
lgd <- expand.grid(x = seq(0,1,diff(xR)/150),
                   y = seq(0,1,diff(yR)/150)) %>%
  mutate(#mxCol = (1-y)^x, #Yager's complement function
         mxCol = (y*x),
         brks = ntile(mxCol,4),
         brks = ifelse(brks==1,"Low",ifelse(brks==2,"Medium",ifelse(brks==3,"High","Very high")))
  )

lgd$brks <- factor(lgd$brks, levels = c("Low","Medium","High","Very high"))
(riskspace <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+
    scale_fill_viridis_c()+
    scale_fill_manual(values = c("Low"="#d3d3d3", "Medium"="#a88283", "High"="#7e433e", "Very high"="#551601"))+
    geom_point(data = df, aes(x = Vulnerability, y = impact, shape=sce),size = .8, stroke = .2) +
    labs(y = "Climate change impacts", x = "", title = "a. RISK SPACE")+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.2), labels = c("0",".2",".4",".6",".8","1"))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.2), labels = c("0",".2",".4",".6",".8","1"))+ 
    theme_bw(base_size = 8)+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP5-8.5" = 17))+
    guides(shape="none", colour = "none")+
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1,'cm'),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          axis.ticks = element_line(linewidth = .1),
          panel.border = element_blank()))

(optSpace <- ggplot() +
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+
    scale_fill_manual(values = c("Low"="#d3d3d3", "Medium"="#a88283", "High"="#7e433e", "Very high"="#551601"))+
    labs(y = "", x = "", title = "b. OPTION SPACE")+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.1), position = "right")+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.2))+ 
    theme_bw(base_size = 8)+
    guides(shape="none", colour = "none")+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()))

library(patchwork)
# Create grid
grobs <- ggplotGrob(optSpace)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

(rrSpace <- ((riskspace|optSpace+theme(legend.position = "none"))/legend)+plot_layout(heights = c(2,.1)))
ggsave(plot=rrSpace, "3_Outputs/plots/Transition/Fig2N.png", dpi=1200, height=3, width=4)

##############################################################################################################################
# Estimate potential residual risk and plot difference among villages
#############################################################################################################################

riskMaster <- riskMaster %>% mutate(risk585 = ((imp.ssp585.2050*Vulnerable)),
                                    risk370 = ((imp.ssp370.2050*Vulnerable)),
                                    risk245 = ((imp.ssp245.2050*Vulnerable)))

summary(riskMaster$risk585, na.rm=TRUE);sd(riskMaster$risk585, na.rm=TRUE)
summary(riskMaster$risk245, na.rm=TRUE);sd(riskMaster$risk245, na.rm=TRUE)

df <- rbind(data.frame(sce = "SSP2-4.5", risk = (riskMaster$risk245), village = riskMaster$VillNation, ISO3 = riskMaster$ISO3),
            data.frame(sce = "SSP5-8.5", risk = (riskMaster$risk585), village = riskMaster$VillNation, ISO3 = riskMaster$ISO3))
Q1 <- summary(df$risk)[[2]] #first quartile
Q3 <- summary(df$risk)[[5]] #third quartile

(xx <- (lgd %>% group_by(brks) %>% summarise(mx = max(mxCol, na.rm = TRUE)))[2])
df$brks <- ifelse(df$risk < xx[[1]][1], "Low", ifelse(df$risk <  xx[[1]][2], "Medium", ifelse(df$risk < xx[[1]][3],"High","Very high")))
(plt2 <- ggplot(data = df)+
    # geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = .33)+
    # geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = .33, ymax = .66)+
    # geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = .66, ymax = Inf)+
    geom_point(aes(x=reorder(village,risk), y=risk, colour = brks, shape = sce), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    labs(x = "Villages", y = "Climate risk [index]", title = "c. RISK SCORE")+
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,.2), position = "left")+
    scale_shape_manual("", values = c("SSP2-4.5" = 19, "SSP5-8.5" = 17))+
    theme_classic(base_size = 8)+
    scale_colour_manual(values = c("Low" = "#d3d3d3", "Medium" = "#a88283", "High"="#7e433e", "Very high" = "#551601"))+
    guides(colour = "none")+
    theme(legend.position = c(0.1,.95), 
          legend.text = element_text(size = 8),
          legend.spacing.y = element_blank(), 
          legend.background = element_rect(fill = NA),
          panel.grid.major.y = element_line(linewidth = .1),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))
ggsave(plot = plt2, "3_Outputs/plots/Transition/Fig2b.Risk.png", dpi = 1200, height = 3, width = 4)

#Plots bars for each country
# df1 <- riskMaster %>% group_by(ISO3) %>%
#   summarise(mn.585 = mean(risk585),
#             sd.585 = sd(risk585),
#             mn.370 = mean(risk370),
#             sd.370 = sd(risk370),
#             mn.245 = mean(risk245),
#             sd.245 = sd(risk245)) %>% ungroup()
# 
# df2 <- riskMaster %>%
#   summarise(mn.585 = mean(risk585),
#             sd.585 = sd(risk585),
#             mn.370 = mean(risk370),
#             sd.370 = sd(risk370),
#             mn.245 = mean(risk245),
#             sd.245 = sd(risk245))
# df2 <- cbind(ISO3 = "ALL", df2)
# df <- rbind(df1, df2)

df <- rbind(data.frame(sce = "SSP2-4.5", MN = riskMaster$risk245, ISO3 = df$ISO3),
            data.frame(sce = "SSP3-7.0", MN = riskMaster$risk370, ISO3 = df$ISO3),
            data.frame(sce = "SSP5-8.5", MN = riskMaster$risk585, ISO3 = df$ISO3)) %>% group_by (ISO3) %>% mutate(sortMag = mean(MN))
(plt3 <- ggplot(data = df)+
    geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = .33)+
    geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = .33, ymax = .66)+
    geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = .66, ymax = Inf)+
    geom_boxplot(aes(x = reorder(ISO3, sortMag), y = MN, fill = sce),linewidth = .1, position = position_dodge(width =.8))+
    #geom_pointrange(aes(x = reorder(ISO3, sortMag), y = MN, ymin = LL, ymax = UL, colour=ISO3, shape = sce),size=.5, linewidth = .2, position = position_dodge(width =.5))+
    #scale_fill_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50", "ALL"="cyan"))+
    #scale_shape_manual(values = c("SSP2-4.5" = 1, "SSP3-7.0"=0,"SSP5-8.5" = 2))+
    scale_fill_manual(name="", values = c("SSP2-4.5" = "#749B58FF", "SSP3-7.0"="#466983FF","SSP5-8.5" = "#CE3D32FF"))+
    labs(y = "", x="", title = "")+guides(colour = "none")+
    scale_y_continuous(name = "Climate risk [index]", expand = c(0,0), limits = c(0,1), breaks = seq(0,1,.2))+
    theme_classic(base_size = 12)+
    theme(legend.position = c(0.15,.8),
          legend.title = element_blank(),
          legend.text = element_text(size = 5),
          legend.background = element_blank(),
          #panel.background = element_rect(fill = "transparent", colour = NA),
          #plot.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(angle = 0),
          #axis.line.y = element_blank(),
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = .1),
          axis.ticks = element_line(linewidth = .1)))
ggsave(plot = plt3, "3_Outputs/plots/Transition/FigS3b.png", dpi = 1200, height = 4, width = 4)

###################################################################################################################
#                     ECONOMIC VALUATION APPROACHES
##################################################################################################################

#Import value coefficients
econValues <- readxl::read_excel("2_Data/sheet/4_EconomicValuations/EcosystemServiceValueCoefficients.xlsx", sheet = "mini")
riskMaster <- merge(riskMaster, econValues, by.x = "ISO3")

#Find total economic value
riskMaster$TEV = ((riskMaster$CoralExt*riskMaster$CoralsVal)+(riskMaster$seagrassExt*riskMaster$SeagrassVal)+(riskMaster$mangroveExt*riskMaster$MangroveVal)+(riskMaster$Cropland*riskMaster$CropsVal))/1e4 #divide by 10000 to convert from meters to hectares
plot(riskMaster$TEV/1e6, (1-riskMaster$risk585))

riskMaster$fTEV_ssp585 = riskMaster$TEV*(1-riskMaster$risk585)
riskMaster$fTEV_ssp370 = riskMaster$TEV*(1-riskMaster$risk370)
riskMaster$fTEV_ssp245 = riskMaster$TEV*(1-riskMaster$risk245)

#Order of magnitude change
#write_excel_csv(riskMaster, "3_Outputs/sheets/RiskMasterSheet.csv")

#Plot
xmn <- median(riskMaster$TEV/1e6, na.rm=TRUE)
yR <- range(riskMaster$risk585);xR <- range(riskMaster$TEV/1e6)
lgd <- expand.grid(
  x = seq(round(xR[1],0),round(xR[2],2), diff(xR)/150),
  y = seq(0,1, diff(yR)/150)
  ) %>% mutate(
    x1 = ifelse(x<xmn,1,2),
    y1 = ifelse(y<.33,1,ifelse(y<.66,2,3))
    )

custom_pal3 <- c(
  "1-1" = "#d3d3d3",
  "2-1" = "#ba8890",
  "3-1" = "#9e3547",
  "1-2" = "#4279b0",
  "2-2" = "#3a4e78",
  "3-2" = "#311e3b"
)

(plt1 <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = paste(y1,x1,sep="-")))+
    scale_fill_manual(values = custom_pal3)+
    geom_point(data = riskMaster, aes(x = TEV/1e6, y = risk585, shape = "SSP5-8.5"), size = 1, stroke = .2) +
    geom_point(data = riskMaster, aes(x = TEV/1e6, y = risk245, shape = "SSP2-4.5"), size = 1, stroke = .2) +
    labs(y = "Climate risk", x = "Total economic value (Million US$/year)")+
    #scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme_bw(base_size = 8)+
    scale_x_continuous(expand = c(0,0), breaks = seq(round(xR[1],0),round(xR[2],2),1))+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.2))+
    scale_shape_manual(values = c("SSP2-4.5" = 1, "SSP5-8.5" = 2))+
    guides(fill="none", colour = "none")+
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 5),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          panel.border = element_blank()))
ggsave(plot = plt1,"3_Outputs/plots/Transition/TEVBiv.png", width = 4, height = 4, dpi = 1200)

(xx <- rbind(
  data.frame(value=riskMaster$TEV/1e6, sce = "Current", villages=riskMaster$VillNation, ISO3 = riskMaster$ISO3, rankBy = (riskMaster$TEV/1e6)),
  data.frame(value=riskMaster$fTEV_ssp245/1e6, sce = "SSP2-4.5", villages=riskMaster$VillNation, ISO3 = riskMaster$ISO3, rankBy = (riskMaster$TEV/1e6)),
  data.frame(value=riskMaster$fTEV_ssp585/1e6, sce = "SSP5-8.5", villages=riskMaster$VillNation, ISO3 = riskMaster$ISO3, rankBy = (riskMaster$TEV/1e6)))
  )
df <- xx %>% group_by(villages,sce) %>% arrange(desc(villages), .by_group = TRUE) %>%
  ungroup() %>% mutate(paired = rep(1:(n()/3),each=3))
(plt2 <- ggplot()+
    geom_point(data = filter(df, sce == "Current"), aes(x=value, y=reorder(villages,-rankBy), colour = sce, shape=sce), size = 1.2)+
    geom_point(data = filter(df, sce == "SSP2-4.5"), aes(x=value, y=reorder(villages,-rankBy), colour = sce, shape=sce), size = 1.2)+
    geom_point(data = filter(df, sce == "SSP5-8.5"), aes(x=value, y=reorder(villages,-rankBy), colour = sce, shape=sce), size = 1.2)+

    geom_line(data = df, aes(x=value, y=villages, group = paired),color="grey",linewidth=.1, arrow = arrow(ends = "first",type = "closed",length=unit(0.01, "inches")))+
    labs(x = "Total economic value \n(Million US$/year)", y = "")+
    #scale_y_discrete(position = "right")+
    scale_shape_manual("", values = c("Current" = 16, "SSP2-4.5" = 15, "SSP5-8.5" = 17))+
    scale_colour_manual("", values = c("Current" = "grey", "SSP2-4.5" = "cyan", "SSP5-8.5" = "darkred"))+
    theme_classic(base_size = 10)+
    #guides(colour = "none")+
    theme(legend.position = "top",
          #legend.text = element_text(size = 5),
          legend.key.size = unit(.5, 'lines'),
          legend.spacing.y = element_blank(),
          legend.background = element_rect(fill = NA),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          #axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1),
          axis.ticks = element_line(linewidth = .1))
  )

(rr1 <- ggplot()+
    geom_boxplot(data = riskMaster, aes(x="Current", y = TEV/1e6, fill = "Current"), colour = "grey", linewidth =0.1, outlier.size = .01)+
    geom_boxplot(data = riskMaster, aes(x="SSP2-4.5", y = fTEV_ssp245/1e6, fill = "SSP2-4.5"), colour = "grey", linewidth =0.1, outlier.size = .01)+
    geom_boxplot(data = riskMaster, aes(x="SSP5-8.5", y = fTEV_ssp585/1e6, fill = "SSP5-8.5"), colour = "grey", linewidth =0.1, outlier.size = .01)+
    labs(y = "TEV \n(Million US$/year)", x = "")+
    scale_fill_manual(values = c("Current" = "grey", "SSP2-4.5" = "cyan", "SSP5-8.5" = "darkred"))+
    theme_classic(base_size = 10)+    guides(colour = "none")+
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          #legend.text = element_text(size = 5),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_line(linewidth = .1),
          axis.ticks.y = element_line(linewidth = .1)))
(fg3b <- plt2 + inset_element(rr1, 0.6, 0.6, 1, 1))
ggsave(plot = fg3b,"3_Outputs/plots/Fig3b.png", width = 4.5, height = 4.5, dpi = 1200)



