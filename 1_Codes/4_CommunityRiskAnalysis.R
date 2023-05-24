library(tidyverse)
dfRisk <- read_csv("2_Data/sheet/RiskDataFINAL.csv")
dfRisk <- dfRisk %>% mutate(ISO3 = Country,
                            ISO3 = ifelse(ISO3 == "Kenya", "KEN", ISO3),
                            ISO3 = ifelse(ISO3 == "Tanzania", "TZA", ISO3),
                            ISO3 = ifelse(ISO3 == "Madagascar", "MDG", ISO3),
                            ISO3 = ifelse(ISO3 == "Mozambique", "MOZ", ISO3))
I_Ctrl <- read_csv("2_Data/sheet/ImpactControl.csv")
dfRisk <- merge(dfRisk, I_Ctrl, by.x = "ISO3")
plot(dfRisk$AdaptiveCapacity, exp(-dfRisk$ic2020))


##################################################################################################################
#                                  PLOT RISK SPACES
##################################################################################################################
dfRisk$Vulnerable <- dfRisk$Sensitivity/dfRisk$AdaptiveCapacity
df <- rbind(data.frame(sce = "SSP2-4.5", impact = dfRisk$imp.ssp245.2050, Vulnerability = dfRisk$Vulnerable, village = dfRisk$Villages, ISO3 = dfRisk$ISO3),
            data.frame(sce = "SSP3-7.0", impact = dfRisk$imp.ssp370.2050, Vulnerability = dfRisk$Vulnerable, village = dfRisk$Villages, ISO3 = dfRisk$ISO3))
yR <- range(df$impact);xR <- range(df$Vulnerability)
lgd <- expand.grid(x = seq(xR[1],xR[2], diff(xR)/1500),
                   y = seq(yR[1],yR[2], diff(yR)/1500)) %>% 
  mutate(x1 = scales::rescale(x),
         y1 = scales::rescale(y),
         mxCol = x1^2+y1^2,
         brks = ntile(mxCol,3),
         brks = ifelse(brks==1,"Acceptable",ifelse(brks==2,"Tolerable","Intolerable"))
  )

lgd$brks <- factor(lgd$brks, levels = c("Acceptable","Tolerable","Intolerable"))
(plt1 <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+
    #scale_fill_viridis_d(option = "B", direction = -1, begin = .25,end = .75)+
    scale_fill_manual(values = c("Acceptable" = "grey90", "Tolerable" = "grey80", "Intolerable"="grey70"))+
    geom_point(data = df, aes(x = Vulnerability, y = impact, colour = ISO3, shape=sce), size = 1.5) +
    labs(y = "Climate Change Impacts [Index]", x = "Vulnerability [Index]", title = "Risk Space")+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="Cyan"))+
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
          panel.border = element_rect(linewidth = 1)))
ggsave("3_Outputs/plots/Risk Space.png", width = 4, height = 4, dpi = 1200)

##############################################################################################################################
#                                 Estimate potential residual risk and plot difference among villages
#############################################################################################################################
dfRisk <- dfRisk %>% 
  mutate(rr.ssp585 = ((imp.ssp585.2050*Sensitivity)/AdaptiveCapacity)*exp(-ic2020),
         rr.ssp370 = ((imp.ssp370.2050*Sensitivity)/AdaptiveCapacity)*exp(-ic2020),
         rr.ssp245 = ((imp.ssp245.2050*Sensitivity)/AdaptiveCapacity)*exp(-ic2020))

df <- rbind(data.frame(sce = "SSP2-4.5", risk = (dfRisk$rr.ssp245), village = dfRisk$Villages, ISO3 = dfRisk$ISO3),
            data.frame(sce = "SSP3-7.0", risk = (dfRisk$rr.ssp370), village = dfRisk$Villages, ISO3 = dfRisk$ISO3))
Q1 <- summary(df$risk)[[2]] #first quartile
Q3 <- summary(df$risk)[[5]] #third quartile
(plt2 <- ggplot(data = df)+
    geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Q1)+
    geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = Q1, ymax = Q3)+
    geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = Q3, ymax = Inf)+
    #geom_col(aes(reorder(Villages, risk.ssp370.2050), risk.ssp370.2050), width = .2, fill = "grey90")+
    geom_point(aes(x=reorder(village,-risk), y=risk, colour = ISO3, shape = sce), size = 1.5, position = position_dodge2(width =.5))+
    #geom_text(aes(x=reorder(village,risk), y=risk, label = round(risk,2)), size = 2, colour = "black", fontface = "bold",position = position_dodge(width =.5))+
    labs(x = "Coastal communities", y = "Residual risk index")+
    #scale_y_continuous(expand = c(0,0), limits = c(.25,.55), breaks = seq(.25,.55, 0.05))+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
    theme_classic(base_size = 10)+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    scale_fill_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme(legend.position = "", 
          legend.background = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

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
ggsave(plot=plt2, "3_Outputs/plots/3_RiskResidual.png", dpi = 1200, height = 4, width = 6)

###################################################################################################################
#                     ECONOMIC VALUATION APPROACHES
##################################################################################################################
#Import value coefficients
econValues <- readxl::read_excel("2_Data/sheet/4_EconomicValuations/EcosystemServiceValueCoefficients.xlsx", sheet = "mini")
dfRisk <- merge(dfRisk, econValues, by.x = "ISO3")
#Find total economic value
dfRisk$TEV = ((dfRisk$corals*dfRisk$CoralsVal)+(dfRisk$seagrass*dfRisk$SeagrassVal)+(dfRisk$mangrove*dfRisk$MangroveVal)+(dfRisk$cropcover*dfRisk$CropsVal))/1e4 #divide by 10000 to convert from meters to hectares
plot(dfRisk$TEV/1e6, dfRisk$rr.ssp370)
dfRisk$ld_ssp370 = dfRisk$TEV*dfRisk$rr.ssp370
dfRisk$ld_ssp245 = dfRisk$TEV*dfRisk$rr.ssp245
summary(dfRisk$ld_ssp370)
sum(dfRisk$ld_ssp370)
sd(dfRisk$ld_ssp370)

dfs <- dfRisk[,c("Country","Villages","Sensitivity","AdaptiveCapacity","imp.ssp245.2050","imp.ssp370.2050","rr.ssp245", "rr.ssp370","TEV","ld_ssp245","ld_ssp370")]
write_excel_csv(dfs, "2_Data/sheet/TableS3.csv")
