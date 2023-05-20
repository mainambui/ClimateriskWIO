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

df <- rbind(data.frame(sce = "SSP2-4.5", impact = dfRisk$imp.ssp245.2050, vulnerab = dfRisk$Sensitivity/dfRisk$AdaptiveCapacity, village = dfRisk$Villages, ISO3 = dfRisk$ISO3),
            data.frame(sce = "SSP3-7.0", impact = dfRisk$imp.ssp370.2050, vulnerab = dfRisk$Sensitivity/dfRisk$AdaptiveCapacity, village = dfRisk$Villages, ISO3 = dfRisk$ISO3))
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
    geom_point(data = df, aes(x = vulnerab, y = impact, colour = ISO3, shape=sce), size = 1.5) +
    labs(y = "Climate change impacts", x = "", title = "Risk Space")+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50"))+
    theme_bw(base_size = 10)+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
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
#ggsave("outputs/xx1.png", width = 5.5, height = 4, dpi = 1200)

dfRisk<- dfRisk %>% mutate(rr.ssp585 = ((imp.ssp370.2050*Sensitivity)/AdaptiveCapacity)*exp(-ic2020),
                           rr.ssp245 = ((imp.ssp245.2050*Sensitivity)/AdaptiveCapacity)*exp(-ic2020))
df <- rbind(data.frame(sce = "SSP2-4.5", risk = (dfRisk$rr.ssp245), village = dfRisk$Villages, ISO3 = dfRisk$ISO3),
            data.frame(sce = "SSP5-8.5", risk = (dfRisk$rr.ssp585), village = dfRisk$Villages, ISO3 = dfRisk$ISO3))
Q1 <- summary(df$risk)[[2]] #first quartile
Q3 <- summary(df$risk)[[5]] #third quartile
(plt2 <- ggplot(data = df)+
    geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Q1)+
    geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = Q1, ymax = Q3)+
    geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = Q3, ymax = Inf)+
    #geom_col(aes(reorder(Villages, risk.ssp370.2050), risk.ssp370.2050), width = .2, fill = "grey90")+
    geom_point(aes(x=reorder(village,-risk), y=risk, colour = ISO3, shape = sce), size = 1.5, position = position_dodge2(width =.5))+
    #geom_text(aes(x=reorder(village,risk), y=risk, label = round(risk,2)), size = 2, colour = "black", fontface = "bold",position = position_dodge(width =.5))+
    labs(x = "Coastal communities", y = "Residual risk index", title = "(A)")+
    #scale_y_continuous(expand = c(0,0), limits = c(.25,.55), breaks = seq(.25,.55, 0.05))+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP5-8.5" = 17))+
    theme_classic(base_size = 10)+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50"))+
    scale_fill_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50"))+
    theme(legend.position = "", 
          legend.background = element_rect(fill = NA),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

#Plots bars for each country
df1 <- dfRisk %>% group_by(ISO3) %>% 
  summarise(mn.585 = mean(rr.ssp585),
            sd.585 = sd(rr.ssp585),
            mn.245 = mean(rr.ssp245),
            sd.245 = sd(rr.ssp245)) %>% ungroup()

df2 <- dfRisk %>% 
  summarise(mn.585 = mean(rr.ssp585),
            sd.585 = sd(rr.ssp585),
            mn.245 = mean(rr.ssp245),
            sd.245 = sd(rr.ssp245)) 
df2 <- cbind(ISO3 = "ALL", df2)
df <- rbind(df1, df2)

df <- rbind(data.frame(sce = "SSP2-4.5", MN = df$mn.245, LL = df$mn.245-df$sd.245, UL = df$mn.245+df$sd.245, ISO3 = df$ISO3),
            data.frame(sce = "SSP5-8.5", MN = df$mn.585, LL = df$mn.585-df$sd.585, UL = df$mn.585+df$sd.585, ISO3 = df$ISO3)) %>% group_by (ISO3) %>% mutate(sortMag = mean(MN))
(plt3 <- ggplot(data = df)+
    geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Q1)+
    geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = Q1, ymax = Q3)+
    geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = Q3, ymax = Inf)+
    geom_pointrange(aes(x = reorder(ISO3, sortMag), y = MN, ymin = LL, ymax = UL, colour=ISO3, shape = sce),
                    size=.1, linewidth = .2, position = position_dodge(width =.5))+
    scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50", "ALL"="cyan"))+
    scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP5-8.5" = 17))+
    labs(y = "", x="", title = "(B)")+
    #scale_y_continuous(expand = c(0,0), limits = c(.25,.55), breaks = seq(.25,.55, 0.05),position = "right")+
    theme_classic(base_size = 10)+
    theme(legend.position = "", 
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_blank(), 
          axis.line = element_line(linewidth = .1), 
          axis.ticks.x = element_line(linewidth = .1)))
#plt1 + annotation_custom(ggplotGrob(plt2), xmin = 1, xmax = 14, ymin = .5, ymax = .79)
plt2+plt3+plot_layout(ncol = 2, widths = c(5,1))
ggsave("3_Outputs/plots/3_RiskResidual.png", dpi = 1200, height = 4, width = 8)


dfRisk %>% ggplot() + geom_point(aes(x = (TEV/1e6), y = rr.ssp585))
dfRisk %>% ggplot() + geom_point(aes(x = (TEV/1e6), y = imp.ssp370.2050))
