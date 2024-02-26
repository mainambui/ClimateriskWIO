library(dplyr);library(ggplot2);library(readxl)
data_ <- read_excel("2_Data/blue_ventures_data/20240214_CCVA_Analysis.xlsx")
colnames(data_) <- c("HH","Livelihood","Demography","Cultural","Health","Learning","Assets","Flexibility","Agency","Organisation","SS","AC")

data_$villages <- stringr::str_extract(data_$HH, "^.{3}")

data_gg <- aggregate( .~ villages, data_[-1], mean)
(ss_bar <- ggplot(data = data_gg)+geom_point(aes(x=reorder(villages,SS), y=SS), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    labs(x = "", y = "Social sensitivity")+ scale_y_continuous(expand = c(0,0), limits = c(0,.2))+theme_classic(base_size = 15)+
    theme(legend.text = element_text(size = 8),
          legend.spacing.y = element_blank(), 
          legend.background = element_rect(fill = NA),
          panel.grid.major.y = element_line(linewidth = .1),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

(ac_bar <- ggplot(data = data_gg)+geom_point(aes(x=reorder(villages,AC), y=AC), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    labs(x = "Villages", y = "Adaptive capacity")+ scale_y_continuous(expand = c(0,0), limits = c(0.4,.5))+theme_classic(base_size = 15)+
    theme(legend.text = element_text(size = 8),
          legend.spacing.y = element_blank(), 
          legend.background = element_rect(fill = NA),
          panel.grid.major.y = element_line(linewidth = .1),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

library(patchwork)
ss_bar/ac_bar
ggsave("3_Outputs/plots/mdg_files/bar.png", dpi = 1200, height = 5, width = 4)


library(ggthemes);library(ggrepel)
yR <- range(data_gg$SS);xR <- range(data_gg$AC)
lgd <- expand.grid(x = seq(0.4,.5,diff(xR)/500), y = seq(0,.2,diff(yR)/500)) %>% mutate(brk1 = ntile(x,2), brk2 = ntile(y,2),brks = paste(brk1, brk2))
(biplot <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+scale_fill_viridis_c()+
    scale_fill_manual(values = c("1 1"="#d3d3d3", "1 2"="dodgerblue4", "2 1"="yellow4", "2 2"="red4"))+
    geom_point(data = data_gg, aes(x = AC, y = SS),size = 3, stroke = .2) +
    geom_text_repel(data = data_gg, aes(x = AC, y = SS, label = villages),size = 3, colour = "white")+
    scale_y_continuous(name = "Social sensitivity",expand = c(0,0))+
    scale_x_continuous(name = "Adaptive capacity",expand = c(0,0))+ 
    theme_bw(base_size = 15)+
    guides(shape="none")+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1,'cm'),
          legend.text = element_text(size = 5),
          legend.key.height = unit(.1, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          axis.ticks = element_line(linewidth = .1),
          panel.border = element_blank()))
ggsave("3_Outputs/plots/mdg_/mdg_biplot.png", dpi = 1200, height = 4, width = 4)

