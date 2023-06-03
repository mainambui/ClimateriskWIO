library(tidyverse);library(sf);library(sp)

socioecom <- read_csv("2_Data/sheet/3_SocialVulnerability/SocialDataAll.csv")
socioecom <- socioecom %>% mutate(ISO3 = Country,
                            ISO3 = ifelse(ISO3 == "Kenya", "KEN", ISO3),
                            ISO3 = ifelse(ISO3 == "Tanzania", "TZA", ISO3),
                            ISO3 = ifelse(ISO3 == "Madagascar", "MDG", ISO3),
                            ISO3 = ifelse(ISO3 == "Mozambique", "MOZ", ISO3))
I_Ctrl <- read_csv("2_Data/sheet/4_ImpactControl.csv")
socioecom <- merge(socioecom, I_Ctrl, by.x = "ISO3")
plot(socioecom$AdaptiveCapacity, exp(-socioecom$ic2020))


#Import climate data
villageImpacts <- read_csv("2_Data/sheet/3_villagesClimateExposureImpacts.csv")
riskMaster <- merge(socioecom, villageImpacts, by ="Villages")

##################################################################################################################
#                                  PLOT RISK SPACES
##################################################################################################################

riskMaster$Vulnerable <- riskMaster$Sensitivity/riskMaster$AdaptiveCapacity
df <- rbind(data.frame(sce = "SSP2-4.5", impact = riskMaster$imp.ssp245.2050, Vulnerability = riskMaster$Vulnerable, village = riskMaster$Villages, ISO3 = riskMaster$ISO3),
            data.frame(sce = "SSP3-7.0", impact = riskMaster$imp.ssp370.2050, Vulnerability = riskMaster$Vulnerable, village = riskMaster$Villages, ISO3 = riskMaster$ISO3))
# yR <- range(df$impact);xR <- range(df$Vulnerability)
# lgd <- expand.grid(x = seq(xR[1],xR[2], diff(xR)/1500),
#                    y = seq(yR[1],yR[2], diff(yR)/1500)) %>% 
#   mutate(x1 = scales::rescale(x),
#          y1 = scales::rescale(y),
#          mxCol = x1^2+y1^2,
#          brks = ntile(mxCol,3),
#          brks = ifelse(brks==1,"Acceptable",ifelse(brks==2,"Tolerable","Intolerable"))
#   )
# lgd$brks <- factor(lgd$brks, levels = c("Acceptable","Tolerable","Intolerable"))
# (plt1 <- ggplot()+
#     geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+
#     #scale_fill_viridis_d(option = "B", direction = -1, begin = .25,end = .75)+
#     scale_fill_manual(values = c("Acceptable" = "grey90", "Tolerable" = "grey80", "Intolerable"="grey70"))+
#     geom_point(data = df, aes(x = Vulnerability, y = impact, colour = ISO3, shape=sce), size = 1.5) +
#     labs(y = "Climate Change Impacts [Index]", x = "Vulnerability [Index]")+
#     scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
#     theme_bw(base_size = 10)+
#     scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
#     scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
#     guides(shape="none", colour = "none")+
#     theme(legend.position = "bottom",
#           legend.title = element_blank(),
#           legend.background = element_rect(fill = NA),
#           legend.key.size = unit(1, 'cm'),
#           legend.text = element_text(size = 8),
#           legend.key.height = unit(.1, 'cm'),
#           legend.key.width = unit(.2, 'cm'), 
#           panel.border = element_rect(linewidth = 1)))
yR <- range(df$impact);xR <- range(df$Vulnerability)
yq1<- summary(df$impact)[2]
yq3<- summary(df$impact)[5]
xq2<- summary(df$Vulnerability)[3]

lgd <- expand.grid(x = seq(xR[1],xR[2], diff(xR)/1500),
                   y = seq(yR[1],yR[2], diff(yR)/1500)) %>% 
   mutate(x1 = ifelse(x<xq2,1,2),
          y1 = ifelse(y<yq1,1,ifelse(y<yq3,2,3)), 
          bicol = paste(y1,x1,sep = "-"))

custom_bicol <- c(
  "1-1" = "#d3d3d3", # low x, low y
  "2-1" = "#ba8890",
  "3-1" = "#9e3547", # high x, low y
  "1-2" = "#4279b0", # low x, high y
  "2-2" = "#3a4e78",
  "3-2" = "#311e3b" # high x, high y
)
(riskspace <- ggplot() +
    geom_raster(data = lgd, aes(x = x, y = y, fill = bicol))+
    scale_fill_manual(values = custom_bicol)+
    guides(fill = "none")+
    geom_point(data = df, aes(x = Vulnerability, y = impact, colour = ISO3, shape=sce), size = 1.5) +
    labs(y = "Climate Change Impacts [Index]", x = "", title = "RISK SPACE")+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme_bw(base_size = 10)+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
    guides(shape="none", colour = "none")+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'), 
          panel.border = element_blank()))
(optSpace <- ggplot() +
    geom_raster(data = lgd, aes(x = x, y = y, fill = bicol))+
    scale_fill_manual(values = custom_bicol)+
    labs(y = "", x = "", title = "OPTIONS SPACE")+
    theme_bw(base_size = 10)+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0), position = "right")+
    guides(fill = "none",shape="none", colour = "none")+
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
riskspace+optSpace
#ggsave("3_Outputs/plots/RS2.png", width = 5, height = 4, dpi = 1200)
##############################################################################################################################
#                                 Estimate potential residual risk and plot difference among villages
#############################################################################################################################
riskMaster <- riskMaster %>% 
  mutate(risk370 = ((imp.ssp370.2050*Sensitivity)/AdaptiveCapacity),
         risk245 = ((imp.ssp245.2050*Sensitivity)/AdaptiveCapacity),
         rr.ssp370 = ((imp.ssp370.2050*Sensitivity)/AdaptiveCapacity)*exp(-ic2020),
         rr.ssp245 = ((imp.ssp245.2050*Sensitivity)/AdaptiveCapacity)*exp(-ic2020))

df <- rbind(data.frame(sce = "SSP2-4.5", risk = (riskMaster$rr.ssp245), village = riskMaster$Villages, ISO3 = riskMaster$ISO3),
            data.frame(sce = "SSP3-7.0", risk = (riskMaster$rr.ssp370), village = riskMaster$Villages, ISO3 = riskMaster$ISO3))
Q1 <- summary(df$risk)[[2]] #first quartile
Q3 <- summary(df$risk)[[5]] #third quartile
(plt2 <- ggplot(data = df)+
    geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Q1)+
    geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = Q1, ymax = Q3)+
    geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = Q3, ymax = Inf)+
    #geom_col(aes(reorder(Villages, risk.ssp370.2050), risk.ssp370.2050), width = .2, fill = "grey90")+
    geom_point(aes(x=reorder(village,-risk), y=risk, colour = ISO3, shape = sce), size = 1.5, position = position_dodge2(width =.5))+
    #geom_text(aes(x=reorder(village,risk), y=risk, label = round(risk,2)), size = 2, colour = "black", fontface = "bold",position = position_dodge(width =.5))+
    labs(x = "Villages", y = "Residual risk [index]")+
    scale_y_continuous(expand = c(0,0), limits = c(.3,.7), breaks = seq(.3,.7, 0.05))+
    #scale_x_discrete(expand = c(0,0), position = "top")+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
    theme_classic(base_size = 14)+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    scale_fill_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme(legend.position = "", 
          legend.background = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))
ggsave(plot=plt2, "3_Outputs/plots/3_RiskResidual.png", dpi = 1200, height = 5.37, width = 8.6)


wio.ISO3 <- st_read("2_Data/shp/country_shape.shp") %>% st_as_sf() %>% st_transform(crs = "+proj=longlat")
(alWIO = ggplot() +
    geom_sf(data = wio.ISO3, colour = "white", fill = "grey70", linewidth=0.05)+
    geom_point(data = riskMaster, aes(x = x, y = y, colour = ISO3))+
    theme_bw(base_size = 14) +
    guides(colour = guide_legend(title.position = "left",title = "",title.hjust = 0.5))+
    scale_x_continuous(name = "", expand = c(0,0))+
    scale_y_continuous(name = "", expand = c(0,0))+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    coord_sf(xlim = c(20, 60), ylim = c(12, -35), expand = TRUE)+
    theme(legend.position = c(0.85,0.2),
        legend.direction = "vertical",
        legend.title = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        panel.grid = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(.3, 'cm'),
        legend.key.width = unit(.5, 'cm')))

library(ggthemes)
library(ggrepel)
mdg_data <- riskMaster %>% filter(ISO3 == "MDG")
(MDG <- ggplot() +
    geom_sf(data = wio.ISO3, colour = "white", fill = "grey70", linewidth=0.05)+
    geom_point(data = mdg_data,aes(x = x, y = y, colour = ISO3), size = 1)+
    geom_text_repel(data = mdg_data,aes(x = x, y = y, label=Villages), size = 1)+
    theme_bw(base_size = 4) +
    scale_x_continuous(name = "", expand = c(0,0), limits = c(46.5, 48))+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-16, -14.19710))+
    scale_colour_manual(name="", values = c("MDG"="yellow"))+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()))

moz_data <- riskMaster %>% filter(ISO3 == "MOZ")
(MOZ <- ggplot() +
    geom_sf(data = wio.ISO3, colour = "white", fill = "grey70", linewidth=0.05)+
    geom_point(data = moz_data,aes(x = x, y = y, colour = ISO3), size = 1)+
    geom_text_repel(data = moz_data,aes(x = x, y = y, label=Villages), size = 1)+
    theme_bw(base_size = 4) +
    scale_x_continuous(name = "", expand = c(0,0), limits = c(32.5, 33.74740))+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-26.10780, -25.1))+
    scale_colour_manual(name="", values = c("MOZ"="dodgerblue4"))+
    theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()))

ken_tza_data <- riskMaster %>% filter(ISO3 %in% c("KEN", "TZA"))
(KENTZA <- ggplot() +
    geom_sf(data = wio.ISO3, colour = "white", fill = "grey70", linewidth=0.05)+
    geom_point(data = ken_tza_data,aes(x = x, y = y, colour = ISO3), size = 1)+
    geom_text_repel(data = ken_tza_data,aes(x = x, y = y, label=Villages), size=1)+
    theme_bw(base_size = 4) +
    scale_x_continuous(name = "", expand = c(0,0), limits = c(39, 41.3))+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-5.2, -1.8))+
    scale_colour_manual(name="", values = c("KEN"="darkred","TZA"="cyan"))+
    theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()))

library(patchwork)
myPatch <- (KENTZA+MDG+MOZ+plot_layout(ncol = 3))
xx <- alWIO/myPatch + plot_layout(heights = c(3,1))
ggsave("3_outputs/etx10p.png", height = 8.39, width = 5.55, dpi = 1200)

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
plot(riskMaster$TEV/1e6, riskMaster$rr.ssp370)
riskMaster$ld_ssp370 = riskMaster$TEV*riskMaster$rr.ssp370
riskMaster$ld_ssp245 = riskMaster$TEV*riskMaster$rr.ssp245
summary(riskMaster$ld_ssp370)
sum(riskMaster$ld_ssp370)
sd(riskMaster$ld_ssp370)

dfs <- riskMaster[,c("ISO3","Villages","Sensitivity","AdaptiveCapacity","imp.ssp245.2050","imp.ssp370.2050","rr.ssp245", "rr.ssp370","TEV","ld_ssp245","ld_ssp370")]
#write_excel_csv(dfs, "2_Data/sheet/TableS3.csv")


dfs
ggplot()+
  geom_boxplot(data = dfs, aes(x = "SSP2-4.5", y = ld_ssp245/1e6, fill = "SSP2-4.5"))+
  geom_boxplot(data = dfs, aes(x = "SSP3-7.0", y = ld_ssp370/1e6, fill = "SSP3-7.0"))+
  scale_y_continuous("Potential loss & damages (Million US$/year)")+
  scale_x_discrete("")+
  scale_fill_manual(values = c("SSP2-4.5" = "darkgreen", "SSP3-7.0" = "darkred"))+
  theme_minimal(base_size = 14)+
  theme(legend.position = "none", 
        legend.title = element_blank())
ggsave("3_Outputs/plots/potentialloss.png", dpi = 1200, wdith = 2.39, height = 5.84)











yq2<- summary(dfs$rr.ssp370)[2]
xq2<- summary(dfs$TEV/1e6)[3]
#Plot
yR <- range(dfs$rr.ssp370);xR <- range(dfs$TEV/1e6)
lgd <- expand.grid(x = seq(xR[1],xR[2], diff(xR)/1500),
                   y = seq(yR[1],yR[2], diff(yR)/1500)) %>% 
  mutate(x1 = ifelse(x<xq2,1,2),
         y1 = ifelse(y<yq2,1,2)
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
    geom_point(data = dfs, aes(x = TEV/1e6, y = rr.ssp370, shape = ISO3), size = 1.5) +
    labs(y = "Potential residual risk", x = "Total Economic Value (Million US$/year)")+
    #scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme_bw(base_size = 10)+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
    #scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
    guides(shape="none", colour = "none")+
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'), 
          panel.border = element_blank()))
ggsave("3_Outputs/plots/L&D.png", width = 4, height = 4, dpi = 1200)
