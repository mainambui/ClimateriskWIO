#sumamry stats
#non-parametric test of sensitivity and adaptove capacity between villages











##plotting

library(dplyr);library(ggplot2);library(readxl)
data_ <- read_excel("2_Data/blue_ventures_data/20240214_CCVA_Analysis_9Mar.xlsx", skip=1)
#colnames(data_) <- c("HH","Livelihood","Demography","Cultural","Health","Learning","Assets","Flexibility","Agency","Organisation","SS","AC")

#data_$villages <- stringr::str_extract(data_$HH, "^.{3}")

data_gg <- aggregate( .~ Village, data_[-c(1,3)], mean)
stderror <- function(x) sd(x)/sqrt(length(x))
data_se <- aggregate( .~ Village, data_[-c(1,3)], stderror)
data_gg$ss.se<-data_se$Sensitivity
data_gg$ac.se<-data_se$`Adaptive capacity`
(ss_bar <- ggplot(data = data_gg)+
    geom_point(aes(x=reorder(Village,Sensitivity), y=Sensitivity), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    #geom_errorbar(aes(ymin=Sensitivity-ss.se, ymax=Sensitivity+ss.se), colour="black", width=.1)+
    labs(x = "", y = "Social sensitivity")+ scale_y_continuous(expand = c(0,0), limits = c(0,.8))+theme_classic(base_size = 15)+
    theme(legend.text = element_text(size = 8),
          legend.spacing.y = element_blank(), 
          legend.background = element_rect(fill = NA),
          panel.grid.major.y = element_line(linewidth = .1),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

(ac_bar <- ggplot(data = data_gg)+geom_point(aes(x=reorder(Village,`Adaptive capacity`), y=`Adaptive capacity`), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    labs(x = "Villages", y = "Adaptive capacity")+ scale_y_continuous(expand = c(0,0), limits = c(0.2,.8))+theme_classic(base_size = 15)+
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
ggsave("3_Outputs/plots/mdg_/bar.png", dpi = 1200, height = 5, width = 4)


###Biplot woth gradient color
source("1_Codes/make_color_gradient.R")
#create gradient
g <- make_gradient(
  deg = 45, n = 500, cols = brewer.pal(9, "Spectral")
)


library(ggthemes);library(ggrepel)
#yR <- range(data_$Sensitivity);xR <- range(data_$`Adaptive capacity`)
#lgd <- expand.grid(x = seq(0.2,.8,diff(xR)/500), y = seq(0.2,.8,diff(yR)/500)) %>% mutate(brk1 = ntile(x,2), brk2 = ntile(y,2),brks = paste(brk1, brk2))
(biplot <- ggplot()+
    #scale_fill_manual(values = c("1 1"="#d3d3d3", "1 2"="dodgerblue4", "2 1"="yellow4", "2 2"="red4"))+
    annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
    geom_point(data = data_, aes(x = `Adaptive capacity`, y = Sensitivity,    shape=factor(Village)),size = 3, stroke = .2) +
    #geom_text_repel(data = data_, aes(x = `Adaptive capacity`, y = Sensitivity),size = 3, colour = "white")+
    scale_shape_manual(values=seq(0,13))+
    scale_y_continuous(name = "Social sensitivity",expand = c(0,0),limits = c(0.2,0.8))+
    scale_x_continuous(name = "Adaptive capacity",expand = c(0,0),limits = c(0.2,0.8))+ 
    geom_hline(yintercept = quantile(data_$Sensitivity, .5), lty=2)+
    geom_vline(xintercept = quantile(data_$Sensitivity, .5), lty=2)+
    theme_bw(base_size = 15)+
    guides(shape="legend")+
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.key.size = unit(1,'cm'),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.2, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          axis.ticks = element_line(linewidth = .1),
          panel.border = element_blank()))
ggsave("3_Outputs/plots/mdg_/mdg_biplot.png", dpi = 1200, height = 6, width = 6)



# other plot
require(reshape2);library(ggsignif)
domains<-ggplot(data = melt(data_[,-c(1,2,3)]), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))+ theme(legend.position="none", axis.text.x=element_text(angle = -90, hjust = 0))
ggsave("3_Outputs/plots/mdg_/mdg_domains.png", dpi = 1200, height = 4, width = 4) +  geom_jitter() 



##boxplots with stat summary - this function calculates staus sumary for geon stat summary
get_box_stats <- function(y, upper_limit = max(data_$Sensitivity) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(length(y), "\n",
                  round(mean(y), 2), "\n",
                  round(median(y), 2), "\n"
    )
  ))
}

#Sensitivity
SS<-ggplot(data_, aes(x = reorder(Village,Sensitivity), y = `Sensitivity`, fill = `Sensitivity`)) +
  geom_boxplot( outlier.colour = "red",outlier.fill = "red") +
  #stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 1) +
  theme_classic() + geom_jitter() +
  #ylim(0.2, 8)+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Villages") 
  #geom_signif(test="wilcox.test", comparisons = combn(levels(as.factor(data_$Village)),2, simplify = F)[-4],step_increase = 0.2)
  #stat_compare_means(comparisons=as.factor(data_$Village), method="wilcox.test", label="p.signif", color="red")

ggsave("3_Outputs/plots/mdg_/mdg_SSplot.png", dpi = 1200, height = 4, width = 8)

#adaptive capacity
AC<-ggplot(data_, aes(x = reorder(Village,`Adaptive capacity`), y = `Adaptive capacity`, fill = `Adaptive capacity`)) +
  geom_boxplot( outlier.colour = "red",outlier.fill = "red") +
  #stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 1) +
  theme_classic() + 
  geom_jitter() +
  #ylim(0.2, 8)+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Villages") 
#geom_signif(test="wilcox.test", comparisons = combn(levels(as.factor(data_$Village)),2, simplify = F)[-4],step_increase = 0.2)
#stat_compare_means(comparisons=as.factor(data_$Village), method="wilcox.test", label="p.signif", color="red")

ggsave("3_Outputs/plots/mdg_/mdg_ACplot.png", dpi = 1200, height = 4, width = 8)

library(patchwork)
SS/AC
ggsave("3_Outputs/plots/mdg_/cc.ac.boxplots.png", dpi = 1200, height = 5, width = 4)
