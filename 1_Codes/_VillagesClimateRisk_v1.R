#PREPARE TO EXTRACT DATA TO AOO

rm(list = ls())
library("tidyverse");library("raster");library("SearchTrees");library("sf");library("sp");library("pracma")

#import some important functions
inormal <- function(x) {
  qrank <- ((rank(x, na.last = TRUE, ties.method= "random") - 0.5) / sum(!is.na(x)))
  z_score <- scales::rescale(sqrt(2)*pracma::erfinv(2*qrank-1))
  return(z_score)  }

#Import PU and convert to a spatial object
wio.AOO <- readRDS("2_Data/sheet/2_Ecosystems/wioAOO.rds")
wio.AOO.spdf <- st_as_sf(wio.AOO, coords=c('x', 'y'), crs="+proj=longlat")
#wio.ISO3 <- st_read("2_Data/shp/country_shape.shp") %>% st_as_sf() %>% st_transform(crs = "+proj=longlat")

#LOAD METRICS
(clim.nc <- list.files("./2_Data/raster",pattern='*.nc',full.names=TRUE))

#Note that each NC file contains eight layers which is generally structured as SSP126_2050, SSP126_2100, SSP245_2050, SSP245_2100, SSP370_2050, SSP370_2100, SSP585_2050, SSP585_2100
nc <- (as.data.frame(expand.grid(x=c(126,245,370,585), y=c(2050,2100))) %>% arrange(desc(-x)) %>% mutate(cbn = paste(x,y,sep = "_")) %>% dplyr::select(cbn))[,1]
varlst <- c("cdd","evspsbl","npp","pH_trend","r10p","tap_trend","sst_trend","sst90p","ts_trend","ts90p","ts90Int","r10Int","sst90Int")
allRaster <- lapply(1:length(varlst), function(x){
  rr <- raster::brick(clim.nc[grep(varlst[[x]], clim.nc)])
  names(rr) <- paste(varlst[[x]], nc, sep = "_")
  return(rr)}
)
hazards <- stack(allRaster)

#Extract the raw hazards to the WIO's AOO of interest
climdata <- raster::extract(hazards, wio.AOO.spdf, sp=TRUE, df=TRUE) %>% as.data.frame()
colnames(climdata)[colnames(climdata) == "coords.x1"] <- "x"
colnames(climdata)[colnames(climdata) == "coords.x2"] <- "y"
climdata <- climdata %>% relocate(x:y, .before = ID)
N <- ncol(climdata)
climdata <- cbind(climdata[1:14], DMwR2::knnImputation(climdata[15:N], k = 3))

#Put in a function to tiddy up the Global Environment
QTtransform <- function(df,vlst,scenario,time){
  nlist <- colnames(df)
  #scenario <- match.arg(sce)
  xx = c()
  for(k in seq_along(time)){
    xx[[k]] <- lapply(1:length(scenario), function(x){d1 <- cbind("ID" = df[,"ID"], sce = paste("ssp",scenario[[x]],sep=""), df[,nlist[grep(paste(scenario[[x]],time[[k]],sep="_"),nlist)]])
    colnames(d1) <- c("ID","Scenario",paste(vlst))
    return(d1)}) %>% bind_rows()
  }
  names(xx) <- paste(time)
  xx <- bind_rows(xx, .id = "Period")
  
  #Quantile transform
  inormal <- function(x){
    qrank <- ((rank(x, na.last = TRUE, ties.method= "random") - 0.5) / sum(!is.na(x)))
    z_score <- scales::rescale(sqrt(2)*pracma::erfinv(2*qrank-1))
    return(z_score)}
  dfNmd <- xx %>% as.data.frame() %>% mutate(pH_trend = -1*pH_trend) %>% mutate(across(cdd:sst90Int, ~ inormal(.x)))
  #Convert to wide format
  dfWide <- lapply(1:length(vlst), function(x){
    xv <- reshape2::dcast(dfNmd, ID ~ Scenario+Period, value.var= vlst[[x]])
    scePeriod <- rep(scenario, each = length(time))
    names(xv) <-c("ID", paste(rep(vlst[[x]], each = length(scePeriod)), scePeriod, time, sep = "_"))
    return(xv)})
  
  stdHzd <- dfWide %>% reduce(left_join, by = "ID")
  return(stdHzd)
}

#Inverse Normal Standardised Variables
(HazardQNormed <- QTtransform(climdata, vlst = varlst, time = 2050, scenario = c(245,370,585)))

#Normalised Exposed Systems Metrics
ExposureNormed <- climdata %>% 
  dplyr::select(ID,CoralExt,seagrassExt,mangroveExt,Cropland,Nb_sp,FRic,FDiv,FEve,FDis,FSpe,FOri)%>%
  #Create copies of the Exposed systems. the original versions will be needed later
  mutate(std_corals = CoralExt,
         std_seagrass = seagrassExt,
         std_mangrove = mangroveExt,
         std_cropcover = Cropland,
         std_Nb_sp = Nb_sp,
         std_FRic = FRic,
         std_FDiv = FDiv,
         std_FEve = FEve)%>% 
  mutate(across(std_corals:std_Nb_sp, ~ inormal(.x)))#FRic, FDiv, and FEve are already normalised variables

HazardExposures <- merge(ExposureNormed, HazardQNormed, by = "ID")
namelist <- colnames(HazardExposures)
(ClimImpacts <- HazardExposures %>% 
    rowwise() %>%
    mutate(
      #mean across 13 normalised climate metrics
      hzd.mn.ssp245.2050 = mean(c_across(namelist[grep("_245", namelist)]), na.rm=TRUE),
      hzd.mn.ssp370.2050 = mean(c_across(namelist[grep("_370", namelist)]), na.rm=TRUE),
      hzd.mn.ssp585.2050 = mean(c_across(namelist[grep("_585", namelist)]), na.rm=TRUE),
      exp.mn.wioo = mean(c_across(std_Nb_sp:std_FEve), na.rm=TRUE),
      
      #standard deviations across 13 normalised climate metrics
      hzd.sd.ssp245.2050 = sd(c_across(namelist[grep("_245", namelist)]), na.rm=TRUE),
      hzd.sd.ssp370.2050 = sd(c_across(namelist[grep("_370", namelist)]), na.rm=TRUE),
      hzd.sd.ssp585.2050 = sd(c_across(namelist[grep("_585", namelist)]), na.rm=TRUE),
      exp.sd.wioo = sd(c_across(std_Nb_sp:std_FEve), na.rm=TRUE)) %>% ungroup() %>%
    
    mutate(
      #Deduce inverse variance weights
      wgts.Exposure = (exp.sd.wioo/exp.mn.wioo)^-1,
      wgts.ssp245.2050 = (hzd.sd.ssp245.2050/hzd.mn.ssp245.2050)^-1,
      wgts.ssp370.2050 = (hzd.sd.ssp370.2050/hzd.mn.ssp370.2050)^-1,
      wgts.ssp585.2050 = (hzd.sd.ssp585.2050/hzd.mn.ssp585.2050)^-1,
      
      #Estimate potential impacts
      imp.ssp245.2050 = ((hzd.mn.ssp245.2050*wgts.ssp245.2050)*(exp.mn.wioo*wgts.Exposure))/(wgts.ssp245.2050+wgts.Exposure),
      imp.ssp370.2050 = ((hzd.mn.ssp370.2050*wgts.ssp370.2050)*(exp.mn.wioo*wgts.Exposure))/(wgts.ssp370.2050+wgts.Exposure),
      imp.ssp585.2050 = ((hzd.mn.ssp585.2050*wgts.ssp585.2050)*(exp.mn.wioo*wgts.Exposure))/(wgts.ssp585.2050+wgts.Exposure)
    ))
hist(ClimImpacts$imp.ssp585.2050, breaks=30)
write_csv(ClimImpacts, "2_Data/sheet/2_ClimateImpactsGRIDs.csv")

#########################################################################################################################################
#Merge grid level impacts to the network
#########################################################################################################################################

ImpactsGRIDs <- read_csv("2_Data/sheet/2_ClimateImpactsGRIDs.csv")
EucDist <- read_csv("2_Data/sheet/DistMatrix.csv")
villageGridID <- read_csv("2_Data/sheet/VillageGridIDs.csv")

EucDist <- EucDist[c("src", "nbr", "EucDist")] %>% filter(EucDist>0) #filter to remove self intersections #
EucDist$inverseDist <- (1/EucDist$EucDist)

colnames(ImpactsGRIDs)[colnames(ImpactsGRIDs)=="ID"]<-"nbr"
EucDist <- merge(EucDist, ImpactsGRIDs, by = "nbr")

#Visual checking of distribution
hist(EucDist$inverseDist, breaks = 30)
(idw.impacts <- EucDist %>% group_by(src) %>%
    summarise(across(CoralExt:imp.ssp585.2050, ~(sum(.x*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE))))%>% mutate(ID = src) %>% dplyr::select(-src)
)
(wio.com.idw.impacts <- merge(villageGridID, idw.impacts, by = "ID"))
plot(wio.com.idw.impacts$Cropland, wio.com.idw.impacts$sst90Int_ssp585)
write_csv(wio.com.idw.impacts, "2_Data/sheet/3_VillageImpacts.csv")


######################################################################################################################################################
# NOW ANALYSE VILLAGE LEVEL RISK
#####################################################################################################################################################
#Import socio-economic data

socioecom <- read_csv("2_Data/sheet/3_SocialVulnerability/SocialDataAll.csv")
socioecom <- socioecom %>% mutate(ISO3 = Country,
                                  ISO3 = ifelse(ISO3 == "Kenya", "KEN", ISO3),
                                  ISO3 = ifelse(ISO3 == "Tanzania", "TZA", ISO3),
                                  ISO3 = ifelse(ISO3 == "Madagascar", "MDG", ISO3),
                                  ISO3 = ifelse(ISO3 == "Mozambique", "MOZ", ISO3))
I_Ctrl <- read_csv("2_Data/sheet/4_ImpactControl.csv")
socioecom <- merge(socioecom, I_Ctrl, by.x = "ISO3")
plot(socioecom$AdaptiveCapacity, exp(-socioecom$ic2020))

library(ggthemes)
library(ggrepel)
socioecom$VillNation <- paste(socioecom$Villages, paste("(",socioecom$ISO3,")", sep = ""))
socioecom$Vulnerable <- socioecom$Sensitivity*(1-socioecom$AdaptiveCapacity)

ggplot(data = socioecom, aes(x="XS", y=Vulnerable,label=VillNation))+
  geom_boxplot(linewidth = 0.3)+
  geom_point(aes(colour = ISO3), position=position_jitter(width=.1, height=0))+
  scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey"))+
  geom_text_repel(size = 1)+
  labs(y = "Social vulnerability", x = "")+
  theme_bw(base_size = 15)+
  guides(shape="none", colour = "none")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(1, 'cm'),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.height = unit(.1, 'cm'),
        legend.key.width = unit(.2, 'cm'), 
        panel.border = element_blank())
#ggsave("3_Outputs/plots/FigS3.png", width = 3, height = 4, dpi = 1200)

##################################################################################################################
#                                  PLOT RISK SPACES
##################################################################################################################
#Import climate data
villageImpacts <- read_csv("2_Data/sheet/3_VillageImpacts.csv")
riskMaster <- merge(socioecom, villageImpacts, by ="Villages")

df <- rbind(data.frame(sce = "SSP2-4.5", impact = riskMaster$imp.ssp245.2050, Vulnerability = riskMaster$Vulnerable, village = riskMaster$Villages, ISO3 = riskMaster$ISO3),
            data.frame(sce = "SSP5-8.5", impact = riskMaster$imp.ssp585.2050, Vulnerability = riskMaster$Vulnerable, village = riskMaster$Villages, ISO3 = riskMaster$ISO3))
yR <- range(df$impact);xR <- range(df$Vulnerability)
lgd <- expand.grid(x = seq(0,.5, diff(xR)/1000),
                   y = seq(0,1, diff(yR)/1000)) %>%
  mutate(x1 = scales::rescale(x),
         y1 = scales::rescale(y),
         mxCol = x1^2*y1^2,
         brks = ntile(mxCol,4),
         brks = ifelse(brks==1,"Low",ifelse(brks==2,"Medium",ifelse(brks==3,"High","Very high")))
  )

lgd$brks <- factor(lgd$brks, levels = c("Low","Medium","High","Very high"))
(riskspace <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+
    #scale_fill_manual(values = c("Low" = "grey90", "Medium" = "grey80", "High"="grey70", "Very high" = "grey60"))+
    scale_fill_manual(values = c("Low" = "#d3d3d3", "Medium" = "#a88283", "High"="#7e433e", "Very high" = "#551601"))+
    geom_point(data = df, aes(x = Vulnerability, y = impact, shape=sce)) +
    labs(y = "Climate change impacts [Index]", x = "", title = "a. RISK SPACE")+
    scale_y_continuous(expand = c(0,0), breaks = seq(.3,1,.1))+
    theme_bw(base_size = 8)+
    scale_x_continuous(expand = c(0,0))+
    scale_shape_manual(values = c("SSP2-4.5" = 1, "SSP5-8.5" = 2))+
    guides(shape="none", colour = "none")+
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          panel.border = element_blank()))

(optSpace <- ggplot() +
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+
    #scale_fill_manual(values = c("Low" = "grey90", "Medium" = "grey80", "High"="grey70", "Very high" = "grey60"))+
    scale_fill_manual(values = c("Low" = "#d3d3d3", "Medium" = "#a88283", "High"="#7e433e", "Very high" = "#551601"))+
    labs(y = "", x = "", title = "b. OPTION SPACE")+
    scale_y_continuous(expand = c(0,0), breaks = seq(.3,1, .1), position = "right")+
    scale_x_continuous(expand = c(0,0))+ 
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
(rrSpace <- ((riskspace|optSpace+theme(legend.position = "none"))/legend)+plot_layout(heights = c(2, .1)))
#ggsave(plot = rrSpace, "3_Outputs/plots/Fig2a.png", dpi = 1200, height = 3, width = 4)

##############################################################################################################################
#                                 Estimate potential residual risk and plot difference among villages
#############################################################################################################################

riskMaster <- riskMaster %>% 
  mutate(risk585 = (imp.ssp585.2050*Vulnerable),
         risk370 = (imp.ssp370.2050*Vulnerable),
         risk245 = (imp.ssp245.2050*Vulnerable),
         brk1 = (scales::rescale(imp.ssp585.2050)^2*scales::rescale(Vulnerable)^2),
         brk2 = (scales::rescale(imp.ssp245.2050)^2*scales::rescale(Vulnerable)^2),
         rr.ssp585 = risk585*exp(-ic2020),
         rr.ssp370 = risk370*exp(-ic2020),
         rr.ssp245 = risk245*exp(-ic2020))
summary(riskMaster$risk585, na.rm=TRUE);sd(riskMaster$risk585, na.rm=TRUE)
summary(riskMaster$risk245, na.rm=TRUE);sd(riskMaster$risk245, na.rm=TRUE)

df <- rbind(data.frame(sce = "SSP2-4.5", risk = (riskMaster$risk245), village = riskMaster$VillNation, ISO3 = riskMaster$ISO3, col = riskMaster$brk2),
            data.frame(sce = "SSP5-8.5", risk = (riskMaster$risk585), village = riskMaster$VillNation, ISO3 = riskMaster$ISO3, col = riskMaster$brk1))
Q1 <- summary(df$risk)[[2]] #first quartile
Q3 <- summary(df$risk)[[5]] #third quartile

(xx <- (lgd %>% group_by(brks) %>% summarise(mx = max(mxCol, na.rm = TRUE)))[2])
df$brks <- ifelse(df$col < 0.00457, "Low", ifelse(df$col <  0.0348, "Medium", ifelse(df$col < 0.146,"High","Very high")))
(plt2 <- ggplot(data = df)+
    #geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = .33)+http://127.0.0.1:13755/graphics/plot_zoom_png?width=3072&height=1223
    #geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = .33, ymax = .66)+
    #geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = .66, ymax = Inf)+
    geom_point(aes(x=reorder(village,risk), y=risk, colour = brks, shape = sce), size = 1, position = position_dodge2(width =.5))+
    labs(x = "Villages", y = "Climate risk", title = "c. RISK SCORE")+
    scale_y_continuous(expand = c(0,0), limits = c(.3,1), breaks = seq(.3,1, .1))+
    scale_shape_manual("", values = c("SSP2-4.5" = 1, "SSP5-8.5" = 2))+
    theme_classic(base_size = 8)+
    scale_colour_manual(values = c("Low" = "#d3d3d3", "Medium" = "#a88283", "High"="#7e433e", "Very high" = "#551601"))+
    #scale_colour_manual(values = c("Low" = "grey90", "Medium" = "grey80", "High"="grey70", "Very high" = "grey60"))+
    guides(colour = "none")+
    theme(legend.position = c(0.1,.9), 
          legend.text = element_text(size = 5),
          legend.spacing.y = element_blank(), 
          legend.background = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

ggsave(plot = plt2, "3_Outputs/plots/Fig2b.png", dpi = 1200, height = 3, width = 4)

#Plots bars for each country
# df1 <- dfRisk %>% group_by(ISO3) %>% 
#   summarise(mn.370 = mean(rr.ssp370),
#             sd.370 = sd(rr.ssp370),
#             mn.245 = mean(rr.ssp245),
#             sd.245 = sd(rr.ssp245)) %>% ungroup()
# 
# df2 <- dfRisk %>% 
#   summarise(mn.370 = mean(rr.ssp370),
#             sd.370 = sd(rr.ssp370),
#             mn.245 = mean(rr.ssp245),
#             sd.245 = sd(rr.ssp245)) 
# df2 <- cbind(ISO3 = "ALL", df2)
# df <- rbind(df1, df2)
# 
# df <- rbind(data.frame(sce = "SSP2-4.5", MN = df$mn.245, LL = df$mn.245-df$sd.245, UL = df$mn.245+df$sd.245, ISO3 = df$ISO3),
#             data.frame(sce = "SSP3-7.0", MN = df$mn.370, LL = df$mn.370-df$sd.370, UL = df$mn.370+df$sd.370, ISO3 = df$ISO3)) %>% group_by (ISO3) %>% mutate(sortMag = mean(MN))
# (plt3 <- ggplot(data = df)+
#     geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Q1)+
#     geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = Q1, ymax = Q3)+
#     geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = Q3, ymax = Inf)+
#     geom_pointrange(aes(x = reorder(ISO3, sortMag), y = MN, ymin = LL, ymax = UL, colour=ISO3, shape = sce),
#                     size=.1, linewidth = .2, position = position_dodge(width =.5))+
#     scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50", "ALL"="cyan"))+
#     scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
#     labs(y = "", x="", title = "(B)")+
#     #scale_y_continuous(expand = c(0,0), limits = c(.25,.55), breaks = seq(.25,.55, 0.05),position = "right")+
#     theme_classic(base_size = 10)+
#     theme(legend.position = "", 
#           panel.background = element_rect(fill = "transparent", colour = NA),
#           plot.background = element_rect(fill = "transparent"),
#           axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#           axis.line.y = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           panel.border = element_blank(), 
#           axis.line = element_line(linewidth = .1), 
#           axis.ticks.x = element_line(linewidth = .1)))
#library(patchwork)


###################################################################################################################
#                     ECONOMIC VALUATION APPROACHES
##################################################################################################################

#Import value coefficients
econValues <- readxl::read_excel("2_Data/sheet/4_EconomicValuations/EcosystemServiceValueCoefficients.xlsx", sheet = "mini")
riskMaster <- merge(riskMaster, econValues, by.x = "ISO3")
#Find total economic value
riskMaster$TEV = ((riskMaster$CoralExt*riskMaster$CoralsVal)+(riskMaster$seagrassExt*riskMaster$SeagrassVal)+(riskMaster$mangroveExt*riskMaster$MangroveVal)+(riskMaster$Cropland*riskMaster$CropsVal))/1e4 #divide by 10000 to convert from meters to hectares
plot(riskMaster$TEV/1e6, riskMaster$rr.ssp585)
riskMaster$ld_ssp585 = riskMaster$TEV*riskMaster$rr.ssp585
riskMaster$ld_ssp370 = riskMaster$TEV*riskMaster$rr.ssp370
riskMaster$ld_ssp245 = riskMaster$TEV*riskMaster$rr.ssp245
sum(riskMaster$ld_ssp585)

dfs <- riskMaster[,c("ISO3","Villages","Sensitivity","AdaptiveCapacity","imp.ssp245.2050","imp.ssp370.2050","imp.ssp585.2050",
                     "rr.ssp245","rr.ssp370","rr.ssp585","TEV","ld_ssp245","ld_ssp370","ld_ssp585")]
write_excel_csv(dfs, "2_Data/sheet/TableS3.csv")

library(ggthemes)
library(ggrepel)
ggplot()+
  geom_boxplot(data = dfs, aes(x="SSP2-4.5", y=ld_ssp245/1e6, fill = "SSP2-4.5"), linewidth =0.5)+
  geom_boxplot(data = dfs, aes(x="SSP5-8.5", y=ld_ssp585/1e6, fill = "SSP5-8.5"), linewidth =0.5)+
  scale_y_continuous("Potential loss & damages (Million US$/year)")+
  scale_x_discrete("")+
  scale_fill_manual(values = c("SSP2-4.5" = "darkgreen", "SSP5-8.5" = "darkred"))+
  #theme_minimal(base_size = 14)+
  theme(legend.position = "none", 
        legend.title = element_blank())
ggsave("3_Outputs/plots/potentialloss.png", dpi = 1200, width = 2.39, height = 5.84)

(rr1 <- ggplot()+
    geom_boxplot(data = dfs, aes(x="SSP2-4.5", y=ld_ssp245/1e6, fill = "SSP2-4.5"), linewidth =0.5)+
    geom_boxplot(data = dfs, aes(x="SSP5-8.5", y=ld_ssp585/1e6, fill = "SSP5-8.5"), linewidth =0.5)+
    labs(y = "", x = "")+
    scale_fill_manual(values = c("SSP2-4.5" = "darkgreen", "SSP5-8.5" = "darkred"))+
    theme_bw(base_size = 10)+coord_flip()+
    guides(colour = "none")+
    theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(1, 'cm'),
        #axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.height = unit(.1, 'cm'),
        legend.key.width = unit(.2, 'cm'), 
        panel.border = element_blank(),plot.margin = unit(c(0,0,0,-5.5), "mm")))

yq2<- summary(dfs$rr.ssp585)[2]
xq2<- summary(dfs$TEV/1e6)[3]
#Plot
yR <- range(dfs$rr.ssp585);xR <- range(dfs$TEV/1e6)
lgd <- expand.grid(x = seq(1,7.5, diff(xR)/150),
                   y = seq(0,1, diff(yR)/150)) %>% 
  mutate(x1 = ifelse(x<4.25,1,2),
         y1 = ifelse(y<0.5,1,2)
  )

custom_pal3 <- c(
  "1-1" = "#d3d3d3", # low x, low y
  "1-2" = "#9e3547", # high x, low y
  "2-1" = "#4279b0", # low x, high y
  "2-2" = "#311e3b" # high x, high y
)
(plt1 <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = paste(x1,y1,sep="-")))+
    scale_fill_manual(values = custom_pal3)+
    geom_point(data = dfs, aes(x = TEV/1e6, y = rr.ssp585, shape = "SSP5-8.5"), size = 1.5) +
    geom_point(data = dfs, aes(x = TEV/1e6, y = rr.ssp245, shape = "SSP2-4.5"), size = 1.5) +
    labs(y = "Potential residual risk", x = "Total Economic Value (Million US$/year)")+
    #scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme_bw(base_size = 10)+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
    scale_shape_manual(values = c("SSP2-4.5" = 1, "SSP5-8.5" = 2))+
    guides(fill="none", colour = "none")+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'), 
          panel.border = element_blank()))


#legNtake<-ggExtra::ggMarginal(legNtake, type="histogram", col="grey30", fill="white", size = 3)
densLUI <- ggplot(data = dataFig3PAs[!is.na(dataFig3PAs$iucn),],aes(x = resLUI)) + 
  geom_histogram(color = "grey10", fill = "white", size = .5) + 
  scale_x_continuous(expand = c(0,0))+theme_void()+theme(plot.margin = unit(c(0,0,0,-5.5), "mm"))

library(patchwork)
(legxx <- rr1/plt1 + plot_layout(ncol =1,heights = c(.5, 4)))
ggsave("3_Outputs/plots/L&D22.png", width = 4, height = 5, dpi = 1200)
