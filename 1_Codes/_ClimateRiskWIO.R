############################################################################################################################
### Supporting Information to ###

# Title: Rising climate risk and loss and damage to coastal small-scale fisheries livelihoods
# Authors: Maina, Asamoah, et al. 202x
# Journal: Nature Sustainability
# School of Natural Sciences, Macquarie University, Sydney, Australia.
# Codes by: asamoahfrt@gmail.com
# Last updated: 

############################################################################################################################
# SUPPORTING SCRIPT 1: ESTIMATING CLIMATE RISK FOR VILLAGES

# Steps in this script:
#  1. Load, understand and prepare the dataset for analysis.
#  2. Extract climate data and ecosystem data to AOO.
#  3. Deduce and transform data for use in the integrative framework.
#  3. Estimate climate impacts based on Choquet integral
#  5. Extract data to
# .....

# First, clean workspace:
rm(list = ls()); gc()

# Install and load R packages needed to run the analysis:
needed_packages <- c("dplyr", "tidyterra", "terra", "sf", "ggplot2")
new.packages <- needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages,new.packages)

# Import some important functions
inormal <- function(x) {
  qrank <- ((rank(x, na.last = TRUE) - 0.5) / sum(!is.na(x)))
  z_score <- scales::rescale(sqrt(2)*pracma::erfinv(2*qrank-1))
  return(z_score)  }

#Import PU and convert to a spatial object
wio.aoo <- readRDS("2_Data/sheet/wio.aoo.rds")
grdSize <- (25*25*1e4)
wooFilter <- as.vector(wio.aoo %>% mutate(across(CoralExt:Cropland, ~ .x/grdSize)) %>% rowwise() %>%
                         mutate(Tot = 100*sum(c(CoralExt,seagrassExt,mangroveExt,Cropland), na.rm = TRUE)) %>% ungroup() %>%
                         filter(Tot > 0.1)%>% dplyr::select(ID))

wio.aoo.sub <- filter(wio.aoo, ID %in% wooFilter$ID)
wio.aoo.spdf <- st_as_sf(wio.aoo.sub, coords=c('x', 'y'), crs="+proj=longlat") %>% vect()
plot(wio.aoo.spdf)

#LOAD METRICS
(clim.nc <- list.files("./2_Data/raster",pattern='*.nc',full.names=TRUE))

#Note that each NC file contains eight layers which is generally structured as SSP126_2050, SSP126_2100, SSP245_2050, SSP245_2100, SSP370_2050, SSP370_2100, SSP585_2050, SSP585_2100
nc <- (as.data.frame(expand.grid(x=c(126,245,370,585), y=c(2050,2100))) %>% arrange(desc(-x)) %>% mutate(cbn = paste(x,y,sep = "_")) %>% dplyr::select(cbn))[,1]
chonic <- c("evspsbl","npp","pH_trend","tap_trend","sst_trend","ts_trend","SLR")
acute <- c("cdd","r10p","sst90p","ts90p")
(varlst <- c(chonic,acute))

hazards <- lapply(1:length(varlst), function(x){
  rr <- rast(clim.nc[grep(varlst[[x]], clim.nc)])
  names(rr) <- paste(varlst[[x]],nc, sep = "_")
  return(rr)}
) %>% rast()

plot(hazards)

#Extract the raw hazards to the WIO's AOO of interest
# set values below 100 to NA.
climdata <- terra::extract(hazards, wio.aoo.spdf, xy=TRUE) %>% as.data.frame()
climdata <- cbind(as.data.frame(wio.aoo.spdf), climdata[,-1])
climdata <- climdata %>% relocate(c(x,y), .before = ID)
N <- ncol(climdata)
#climdata <- cbind(climdata[1:14], DMwR2::knnImputation(climdata[15:N], k = 3))

#Inverse normal standardised variables
source("1_Codes/qTransform.R")
(HazardQNormed <- qTranform(climdata, vlst = varlst, time = c(2050), scenario = c(245,370,585)))

# #Import tropical cyclone data.
slr.tc.data <- rast("C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/3_WIO_SYNTHESIS/Data/1-Climate metrics/TC/TC_count.tif")
#slr.tc.data <- app(slr.tc.data, function(x) (scales::rescale(x)))
slr.tc.data <- cbind(ID = wio.aoo.spdf$ID, (terra::extract(slr.tc.data, wio.aoo.spdf) %>% as.data.frame())[-1])
colnames(slr.tc.data)[colnames(slr.tc.data)=="layer"] <- "TC"

#Note here tropical cyclone does not change with scenario and its available for the recent past period
slr.tc.data$tc_245_2050 <- scales::rescale(slr.tc.data$TC)
slr.tc.data$tc_370_2050 <- scales::rescale(slr.tc.data$TC)
slr.tc.data$tc_585_2050 <- scales::rescale(slr.tc.data$TC)
slr.tc.data <- slr.tc.data[,c(1,3:5)] #Select only TC

all.climdata <- merge(HazardQNormed, slr.tc.data, by = "ID")

# #Normalised exposed systems metrics
all.climdata <- merge(all.climdata, climdata[,1:14], by = "ID") %>%
   #Create copies of the Exposed systems. the original versions will be needed later
   mutate(std_corals = log10(CoralExt/grdSize),
          std_seagrass = log10(seagrassExt/grdSize),
          std_mangrove = log10(mangroveExt/grdSize),
          std_cropcover = log10(Cropland/grdSize),
          std_FDiv = FDiv, 
          std_Nb_sp = FRic,
          std_FRic = FRic,
          std_FEve = FEve )%>%
 mutate(across(std_corals:std_FEve,~ scales::rescale(.x)))#FRic, FDiv, and FEve are already normalised variables
saveRDS(all.climdata, "2_Data/sheet/all.climdata.rds")






#AGGREGATION NEXT
all.data <- readRDS("2_Data/sheet/all.climdata.rds")
namelist <- colnames(all.data)
fuzzyS <- function(x){return(1 - prod((1 - x)))}

#Combine Corals and fish data as one variable
coralSys <- all.data %>% 
  dplyr::select(std_Nb_sp,std_FRic,std_FDiv,std_FEve) %>%
  mutate(dplyr::across(all_of(1:4),~ ifelse(is.na(.x),0,.x))) %>%
  rowwise() %>% 
  mutate(std_fish.div = mean(c(std_Nb_sp,std_FRic,std_FDiv,std_FEve))) %>% ungroup()
all.data <- cbind(all.data, coralSys[,"std_fish.div"])

# the Choquet Integral of f w.r.t mu
library(kappalab)
mu <- readRDS("2_Data/sheet/mu.rds");summary(mu) #Overall mu

#Choquet.integral(mu,f)
exposure <- c("std_seagrass","std_mangrove","std_cropcover","std_corals","std_fish.div")
#exposure <- c("std_seagrass","std_mangrove","std_cropcover","std_corals.all")
(impacts.aoo <- lapply(1:length(exposure),function(i){
  df <- cbind(ID=all.data[,"ID"], all.data[,exposure[[i]]]*all.data[,namelist[grep("_2050",namelist)]])
  n <- ncol(df)
  df <- df %>% 
    mutate(dplyr::across(all_of(2:n),~ ifelse(is.na(.x),0,.x))) %>% 
    rowwise() %>%
    mutate(ssp245.2050 = kappalab::Choquet.integral(mu,c_across(namelist[grep("_245_2050", namelist)])),
           ssp370.2050 = kappalab::Choquet.integral(mu,c_across(namelist[grep("_370_2050", namelist)])),
           ssp585.2050 = kappalab::Choquet.integral(mu,c_across(namelist[grep("_585_2050", namelist)]))) %>% ungroup() %>%
      dplyr::select(ID,ssp245.2050,ssp370.2050,ssp585.2050)
  names(df) <- c("ID",paste(exposure[[i]], c("ssp245.2050","ssp370.2050","ssp585.2050"), sep = "_"))
  return(df)
})
)

impacts.aoo <- impacts.aoo %>% purrr::reduce(left_join, by = "ID")
(impacts.agg <- impacts.aoo %>% rowwise() %>%
  mutate(imp.ssp245.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp245", colnames(impacts.aoo))])),
         imp.ssp370.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp370", colnames(impacts.aoo))])),
         imp.ssp585.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp585", colnames(impacts.aoo))]))))

impacts.agg <- merge(impacts.agg, all.data[,c("ID","CoralExt","seagrassExt","mangroveExt","Cropland")], by = "ID")

#########################################################################################################################################
#Merge grid level impacts to the network
#########################################################################################################################################
idw.matrix <- data.table::fread("2_Data/sheet/idw.dist.matrix.csv", stringsAsFactors=TRUE, encoding="UTF-8")
village.aoo.id <- data.table::fread("2_Data/sheet/village.grid.id.csv", stringsAsFactors=TRUE, encoding="UTF-8")

idw.matrix <- idw.matrix[,c("src", "nbr", "EucDist")] %>% filter(EucDist>0) #filter to remove self intersections #
idw.matrix$inverseDist <- (1/(idw.matrix$EucDist/1e6)^2)
hist(idw.matrix$inverseDist , breaks = 50)
idw.matrix <- subset(idw.matrix, inverseDist >.1)

colnames(impacts.agg)[colnames(impacts.agg)=="ID"]<-"nbr"
idw.matrix <- merge(idw.matrix, impacts.agg, by = "nbr")

(idw.impacts <- idw.matrix %>% group_by(src) %>% 
    summarise(across(where(is.numeric), ~(sum(.x*inverseDist, na.rm = TRUE)/sum(inverseDist, na.rm = TRUE))))%>% mutate(ID = src) %>% dplyr::select(-c(src,EucDist,inverseDist))
)
(villageImpacts <- merge(village.aoo.id, idw.impacts, by = "ID"))

##############################
#FIGURE 3: PLOT IMPACTS for EACH SYSTEM
##############################
data <- as.data.frame(villageImpacts)[,c("Country","Villages",colnames(villageImpacts)[grep("std_",colnames(villageImpacts))])]
data <- data.frame(
  mapply(c,
         cbind(data[,c("Country","Villages",colnames(data)[grep("ssp245",colnames(data))])], sce = "SSP2-4.5"),
         cbind(data[,c("Country","Villages",colnames(data)[grep("ssp585",colnames(data))])], sce = "SSP5-8.5")))
colnames(data) <- c("ISO3", "Villages", "Seagrass","Mangrove","Crop","Corals","Functional diversity","Scenario")
data <- data %>% mutate(
  ISO3 = ifelse(ISO3=="Kenya", "KEN", ISO3),
  ISO3 = ifelse(ISO3=="Mozambique", "MOZ", ISO3),
  ISO3 = ifelse(ISO3=="Madagascar", "MDG", ISO3),
  ISO3 = ifelse(ISO3=="Tanzania", "TZA", ISO3)
  )

data <- filter(data, Scenario == "SSP5-8.5")
# convert data from wide to long format
data_long <- tidyr::gather(data, key = "ecosystems", value = "impacts", -ISO3, -Villages, -Scenario)
data_long$ISO3 <- as.factor(data_long$ISO3)
data_long$impacts <- as.numeric(data_long$impacts)

data_long <- data_long %>% group_by(Villages) %>% mutate(totImpact = sum(impacts, na.rm = TRUE)) %>% ungroup()
data_long$impacts <- 100*(data_long$impacts/data_long$totImpact)

data_long <- data_long %>% group_by(ISO3, Villages) %>% mutate(tot=sum(impacts)) %>% ungroup()
nc_uni <- unique(data$ISO3)

data_long <- data_long %>% arrange(ISO3, tot) %>%
  add_row(ISO3 = rep(nc_uni, 15)) %>% # add a couple of empty rows to separate countries
  arrange(ISO3,Villages) %>% ungroup()
data_long$id <- rep(seq(1, nrow(data_long)/nlevels(as.factor(data_long$ecosystems))), each=nlevels(as.factor(data_long$ecosystems)))

# Get the name and the y position of each label
label_data <- data_long %>% group_by(id, Villages) %>% summarize(tot=sum(impacts))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar# I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

brks <- data_long %>%
  group_by(ISO3) %>%
  summarize(start = min(id), end = max(id) - 2) %>%
  rowwise() %>%
  mutate(title = mean(c(start, end))) %>%
  ungroup() %>%
  mutate(end = data.table::shift(end + 1, n = 1, type = "shift", fill = max(end) + 1), start = start -1)

brks$start[brks$ISO3 == "KEN"] <- -1
brks$end[brks$ISO3 == "TZA"] <- 0

max_value <- max(data_long$tot, na.rm = T); y_max <- max_value+.1; v <- c(0, 25, 50, 75, 100)
brks <- brks %>% mutate(v = list(v)) %>% tidyr::unnest()

# Prepare a data frame for grid (scales)
grid_data <- brks
grid_data$end <- grid_data$end[c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

(p <- data_long %>%
    ggplot() + geom_bar(aes(x = as.factor(id), y = impacts, fill = ecosystems), position="stack", stat="identity", width = 0.5) +
    annotate("text", x = rep(max(data_long$id[data_long$ISO3=="TZA"]),length(v)), y = v-.1,label = paste(head(v), "%"),color = "grey", size = 2.5, angle = 0, fontface = "bold", hjust = 0.7)+
    ylim(-50, y_max+10)+
    scale_fill_manual(
      values = c("Seagrass" = "darkred",
                 "Mangrove" = "yellow",
                 "Functional diversity" = "dodgerblue4",
                 "Corals" = "cyan",
                 "Crop" = "#F2BA49"))+
    #paletteer::scale_fill_paletteer_d(palette = "rcartocolor::Pastel")+
    theme_minimal(base_size = 15) +
    guides(fill = guide_legend(ncol = 1))+
    coord_polar()+
    theme(
      legend.position = c(.52,.5),
      legend.direction = "vertical",
      legend.key.height = unit(.1, 'cm'),
      legend.key.width = unit(.2, 'cm'),
      legend.title = element_blank(),
      legend.text = element_text(size = 5),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      #plot.margin = unit(c(1,-1,-1,-1), "cm")
    )+
    geom_text(data = label_data, aes(x = id, y = tot+5, label = Villages, hjust = hjust), color="black",  size = 3.5, angle = label_data$angle, inherit.aes = FALSE)+
    geom_text(data = brks, aes(x = title, y = -10, label = ISO3),  colour = "black", alpha = 0.8, size = 3.5, fontface = "bold", inherit.aes = FALSE)
)
ggsave(plot=p, filename="3_Outputs/plots/Fig_n.png", dpi = 1200, height = 5, width = 5)
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
I_Ctrl <- read.csv("2_Data/sheet/3_SocialVulnerability/impact.control.csv")
socioecom <- merge(socioecom, I_Ctrl, by.x = "ISO3")
plot(socioecom$AdaptiveCapacity, exp(-1*socioecom$ic2020))

library(ggthemes);library(ggrepel)
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
    geom_point(data = df, aes(x = Vulnerability, y = impact, shape=sce, colour = ISO3),size = .8, stroke = .2) +
    labs(y = "Climate change impacts", x = "", #title = "a. RISK SPACE"
         )+
    scale_y_continuous(expand = c(0,0),breaks = seq(0,1,.2), labels = c("0",".2",".4",".6",".8","1"))+
    scale_x_continuous(expand = c(0,0),breaks = seq(0,1,.2), labels = c("0",".2",".4",".6",".8","1"))+ 
    theme_bw(base_size = 8)+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP5-8.5" = 17))+
    scale_colour_manual(values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50"))+
    guides(shape="none")+
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1,'cm'),
          legend.text = element_text(size = 5),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          axis.ticks = element_line(linewidth = .1),
          panel.border = element_blank()))

(optSpace <- ggplot() +
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+
    scale_fill_manual(values = c("Low"="#d3d3d3", "Medium"="#a88283", "High"="#7e433e", "Very high"="#551601"))+
    labs(y = "", x = "",# title = "b. OPTION SPACE"
         )+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.1), position = "right")+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.2))+ 
    theme_bw(base_size = 8)+
    guides(shape="none", colour = "none")+
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
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
grobs1 <- ggplotGrob(riskspace)$grobs
legend1 <- grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]

(rrSpace <- ((riskspace+theme(legend.position = "none")|optSpace+theme(legend.position = "none"))/legend1)+plot_layout(heights = c(2,.1)))
ggsave(plot=rrSpace, "3_Outputs/plots/Fig2a.png", dpi=1200, height=3, width=4)

##############################################################################################################################
# Estimate potential residual risk and plot difference among villages
#############################################################################################################################

riskMaster <- riskMaster %>% 
  mutate(risk585 = (imp.ssp585.2050*Vulnerable),
         risk370 = (imp.ssp370.2050*Vulnerable),
         risk245 = (imp.ssp245.2050*Vulnerable))

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
ggsave(plot = plt2, "3_Outputs/plots/Fig2b.png", dpi = 1200, height = 3, width = 4)

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
ggsave(plot = plt3, "3_Outputs/plots/FigS3.png", dpi = 1200, height = 4, width = 4)

###################################################################################################################
#                     ECONOMIC VALUATION APPROACHES
##################################################################################################################
r = 1+0.024 #https://www.cbo.gov/publication/58957

#Import value coefficients
#econValues <- readxl::read_excel("2_Data/sheet/4_EconomicValuations/EcosystemServiceValueCoefficients.xlsx", sheet = "mini")
econValues <- readr::read_csv("2_Data/sheet/4_EconomicValuations/Unit_Service_Values_Villages.csv")[,c(2,11,12,13,14)]
colnames(econValues) <- c("Villages","CoralUnitValue","CropUnitValue","MangroveUnitValue","SeagrassUnitValue")
tev.data <- merge(riskMaster, econValues, by = "Villages")

#Find total economic value
tev.data$TEV = ((tev.data[,"CoralExt"]*tev.data[,"CoralUnitValue"])+(tev.data[,"seagrassExt"]*tev.data[,"SeagrassUnitValue"])+(tev.data[,"mangroveExt"]*tev.data[,"MangroveUnitValue"])+(tev.data[,"Cropland"]*tev.data[,"CropUnitValue"]))/1e4 #divide by 10000 to convert from meters to hectares
plot(tev.data$TEV/1e6, (1-tev.data$risk585))

sum(tev.data$TEV, na.rm = TRUE)

tev.data$fTEV_ssp585 = r*(tev.data$TEV*(1-tev.data$risk585))
sum(tev.data$fTEV_ssp585);mean(tev.data$fTEV_ssp585);sd(tev.data$fTEV_ssp585)
median(tev.data$fTEV_ssp585, na.rm = TRUE)/median(tev.data$TEV, na.rm = TRUE)

tev.data$fTEV_ssp370 = r*(tev.data$TEV*(1-tev.data$risk370))
sum(tev.data$fTEV_ssp370);mean(tev.data$fTEV_ssp370);sd(tev.data$fTEV_ssp370)
sum(tev.data$fTEV_ssp370, na.rm = TRUE)/sum(tev.data$TEV, na.rm = TRUE)

tev.data$fTEV_ssp245 = r*(tev.data$TEV*(1-tev.data$risk245))
sum(tev.data$fTEV_ssp245);mean(tev.data$fTEV_ssp245);sd(tev.data$fTEV_ssp245)
median(tev.data$fTEV_ssp245, na.rm = TRUE)/median(tev.data$TEV, na.rm = TRUE)

tev.data$perc_tev_585 <- (tev.data$TEV-tev.data$fTEV_ssp585)/tev.data$TEV
tev.data$perc_tev_245 <- (tev.data$TEV-tev.data$fTEV_ssp245)/tev.data$TEV


#Order of magnitude change
#write_excel_csv(tev.data, "3_Outputs/sheets/RiskMasterSheet.csv")

#import bivariate codes
source("E:/4A_FLII_RISKS/codes/ColMatrix.R")

# Define the number of breaks
nBreaks <- 50
# Create the colour matrix
col.matVul <- colmat(nbreaks = nBreaks, breakstyle = "quantile",
                     xlab = "X", ylab = "Y",
                     upperleft = "#4279b0",
                     upperright = "#311e3b",
                     bottomleft =  "#d3d3d3",
                     bottomright = "#9e3547",
                     saveLeg = FALSE, plotLeg = TRUE)
# Retrieve bivariate colour pallet data
lgdBiv <- BivLegend$data; names(lgdBiv) <- c("binY", "binX", "BivCol", "UID")

#Plot
ymn <- median(tev.data$TEV/1e6, na.rm=TRUE)
xR <- range(tev.data$risk585);yR <- range(tev.data$TEV/1e6)
lgd <- expand.grid(
  y = seq(0,30, diff(yR)/150),
  x = seq(0,1, diff(xR)/150)
  ) %>% mutate(
    binY = ntile(y,50),
    binX = ntile(x,50)
    )%>% 
  inner_join(y = lgdBiv, by = c("binY", "binX"))

# custom_pal3 <- c(
#   "1-1" = "#d3d3d3", # low x, low y
#   "2-1" = "#ba8890",
#   "3-1" = "#9e3547", # high x, low y
#   "1-2" = "#8aa6c2",
#   "2-2" = "#7a6b84", # medium x, medium y
#   "3-2" = "#682a41",
#   "1-3" = "#4279b0", # low x, high y
#   "2-3" = "#3a4e78",
#   "3-3" = "#311e3b" # high x, high y
# )
(plt1 <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = BivCol))+
    #scale_fill_manual(values = custom_pal3)+
    scale_fill_identity()+
    geom_point(data = tev.data, aes(x = risk585, y = TEV/1e6, shape = "SSP5-8.5"), size = 2, stroke = .2) +
    #geom_point(data = tev.data, aes(x = risk245, y = TEV/1e6, shape = "SSP2-4.5"), size = 1, stroke = .2) +
    labs(x = "Climate risk", y = "Total economic value \n(Million US$/year)")+
    #scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme_bw(base_size = 15)+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,30,5), limits = c(0,30))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.2))+
    scale_shape_manual(values = c("SSP2-4.5" = 1, "SSP5-8.5" = 2))+
    guides(fill="none", colour = "none")+
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          #legend.text = element_text(size = 5),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          panel.border = element_blank()))
ggsave(plot = plt1,"3_Outputs/plots/Fig3b.png", width = 4, height = 3.5, dpi = 1200)

(xx <- rbind(
  data.frame(value=tev.data$TEV/1e6, sce = "Current", villages=tev.data$VillNation, ISO3 = tev.data$ISO3, rankBy = (tev.data$TEV/1e6)),
  data.frame(value=tev.data$fTEV_ssp245/1e6, sce = "SSP2-4.5", villages=tev.data$VillNation, ISO3 = tev.data$ISO3, rankBy = (tev.data$TEV/1e6)),
  data.frame(value=tev.data$fTEV_ssp585/1e6, sce = "SSP5-8.5", villages=tev.data$VillNation, ISO3 = tev.data$ISO3, rankBy = (tev.data$TEV/1e6)))
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
    scale_colour_manual("", values = c("Current" = "grey", "SSP2-4.5" = "green4", "SSP5-8.5" = "darkred"))+
    theme_classic(base_size = 10)+
    #guides(colour = "none")+
    theme(legend.position = "bottom",
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
    geom_boxplot(data = tev.data, aes(x="Current", y = TEV/1e6, fill = "Current"), colour = "grey", linewidth =0.1, outlier.size = .01)+
    geom_boxplot(data = tev.data, aes(x="SSP2-4.5", y = fTEV_ssp245/1e6, fill = "SSP2-4.5"), colour = "grey", linewidth =0.1, outlier.size = .01)+
    geom_boxplot(data = tev.data, aes(x="SSP5-8.5", y = fTEV_ssp585/1e6, fill = "SSP5-8.5"), colour = "grey", linewidth =0.1, outlier.size = .01)+
    labs(y = "TEV \n(Million US$/year)", x = "")+
    scale_fill_manual(values = c("Current" = "grey", "SSP2-4.5" = "green4", "SSP5-8.5" = "darkred"))+
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
ggsave(plot = fg3b,"3_Outputs/plots/Fig3b.png", width = 4.5, height = 6, dpi = 1200)


################################
#Adaptation gap plot as radial bar chart 
###############################
# Libraries
library(dplyr)
library(tidyverse)
library(hrbrthemes)
bar_width <- .9 # default width of bars in geom_bar
#plot
tev.data %>%
  filter(!is.na(risk585)) %>%
  arrange(risk585) %>%
  #tail(6) %>%
  mutate(Villages=factor(Villages, Villages)) %>%
  mutate(adapt.gap585 = ifelse(risk585 >0.25, risk585-0.25, risk585))%>%
  ggplot( aes(x=Villages, y=risk585) ) +
  geom_bar(fill="#69b3a2", stat="identity") +
  geom_text(hjust = 1, size = 3, aes( y = 0, label = paste(Villages,""))) +
  theme_ipsum() +   
  geom_rect(aes(
    xmin = as.numeric(Villages) - bar_width / 2,
    xmax = as.numeric(Villages) + bar_width / 2,
    ymin = 0.25,
    ymax = risk585,
  ), fill = "blue") +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.title=element_blank(),
    axis.text.y=element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y") +
  ylim(0,1) 




