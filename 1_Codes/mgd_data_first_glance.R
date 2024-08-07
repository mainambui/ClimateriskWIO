<<<<<<< HEAD
library(dplyr);library(ggplot2);library(readr)
data_ <- read_csv("2_Data/blue_ventures_data/ccva_data_bv.csv")
colnames(data_) <- c("HH","SS","AC","Livelihood","Demography","Cultural","Health","Learning","Assets","Flexibility","Agency","Organisation")
=======
>>>>>>> ec1b050b332c8ef38ba1ef2c92d53049783b4537

library(dplyr);library(ggplot2);library(readxl); 
library(rstatix);library(ggstatsplot);library(plotly)
data_ <- read_excel("2_Data/blue_ventures_data/20240214_CCVA_Analysis_9Mar.xlsx", skip=1)
#colnames(data_) <- c("HH","Livelihood","Demography","Cultural","Health","Learning","Assets","Flexibility","Agency","Organisation","SS","AC")
data_sub<-data_[-c(1,3)]

<<<<<<< HEAD
data_gg <- aggregate( .~ villages, data_[-1], mean)
(ss_bar <- ggplot(data = data_gg)+geom_point(aes(x=reorder(villages,SS), y=SS), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    labs(x = "", y = "Social sensitivity")+ scale_y_continuous(expand = c(0,0), limits = c(0,.6))+theme_classic(base_size = 15)+
=======
#summary stats
#non-parametric test of sensitivity and adaptive capacity between villages
#summary stats
#https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/

data_ %>% 
  group_by(Village) %>%
  get_summary_stats(Sensitivity, type = "common")

ggboxplot(data_, x = "Village", y = "Sensitivity")

#comparisons
ss.kruskal <- data_ %>% kruskal_test(Sensitivity ~ Village)
ss.kruskal

#y.             n     statistic  df      p        method        
#1 Sensitivity 244    59.9       13 0.0000000555 Kruskal-Wallis
#we see significant differences

#effect size
#The interpretation values commonly in published literature are:
#0.01- < 0.06 (small effect), 0.06 - < 0.14 (moderate effect) and >= 0.14 (large effect).

data_ %>% kruskal_effsize(Sensitivity ~ Village)

#.y.             n effsize method  magnitude
#1 Sensitivity   244   0.204 eta2[H] large e    

#From the output of the Kruskal-Wallis test, we know that there is a significant difference between groups,
#but we don’t know which pairs of groups are different.
 
#Pairwise comparisons using Dunn’s test:x
 # Pairwise comparisons
 pwc_ss <- data_ %>% 
   dunn_test(Sensitivity ~ Village, p.adjust.method = "bonferroni") 
 pwc_ss
 
 #Pairwise comparisons using Wilcoxon’s test:
 #confidence level modif9ed tp 0.90 i.e sognofcance level of p<=0.1
 pwc2_ss <- data_ %>% 
 wilcox_test(Sensitivity ~ Village, p.adjust.method = "bonferroni")
 pwc2_ss
 ##The pairwise comparison shows that, only Ampan and Nosy are significantly different .
 
 #Visualise
 pwc2_ss <- pwc2_ss %>% add_xy_position(x = "Village")
 ggboxplot(data_, x = "Village", y = "Sensitivity") +
   stat_pvalue_manual(pwc2_ss, hide.ns = TRUE) +
   labs(subtitle = get_test_label(ss.kruskal, detailed = TRUE),
     caption = get_pwc_label(pwc2_ss))
   #scale_y_continuous(name = "Social sensitivity",expand = c(0,0),limits = c(0.2,0.8))
 


##plotting
 ##correlate the domains
 #https://indrajeetpatil.github.io/ggstatsplot/
 corr_domains<-ggcorrmat(
   data     = data_[,4:12],
  #colors   = c("#B2182B", "white", "#4D4D4D"),
   title    = "Correlalogram for SS and AC domains"
 )
ggsave("3_Outputs/plots/mdg_/domaina_corr.png", dpi = 1200, height = 6, width = 6) 
 
#b
ss.pair<-ggbetweenstats(
  data  = data_,
  x     = Village,
  y     = Sensitivity,
  type = "nonparametric",
  pairwise.display = "significant",
  title = "Sensitivity across villages"
)
ss.pair  + theme(
  axis.text.x=element_text(size = 12, angle = -90, hjust = 0.5),
  axis.line = element_line(linewidth = .1, color="black"), 
  axis.ticks = element_line(linewidth = .1, color="black"),
  axis.title.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 12)
)

ggsave("3_Outputs/plots/mdg_/ss_stats_plots.png", dpi = 1200, height = 8, width = 8)


ac.pair<-ggbetweenstats(
  data  = data_,
  x     = Village,
  y     = `Adaptive capacity`,
  type = "nonparametric",
  pairwise.display = "significant",
  title = "Adaptive Capacity across villages"
)
ac.pair<-ac.pair  + theme(
  axis.text.x=element_text(size = 12, angle = -90, hjust = 0.5),
  axis.line = element_line(linewidth = .1, color="black"), 
  axis.ticks = element_line(linewidth = .1, color="black"),
  axis.title.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 12)
)

ggsave("3_Outputs/plots/mdg_/ac_stats_plots.png", dpi = 1200, height = 8, width = 8)

#Circular plots of domains

ddata_long<-melt(data_[,c(2,4:7)])



##Other plots
#plotting sensitivity and adaptive capacity means 
data_gg <- aggregate( .~ Village, data_[-c(1,3)], mean)
stderror <- function(x) sd(x)/sqrt(length(x))
data_se <- aggregate( .~ Village, data_[-c(1,3)], stderror)
data_gg$ss.se<-data_se$Sensitivity
data_gg$ac.se<-data_se$`Adaptive capacity`
(ss_bar <- ggplot(data = data_gg)+
    geom_point(aes(x=reorder(Village,Sensitivity), y=Sensitivity), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    #geom_errorbar(aes(ymin=Sensitivity-ss.se, ymax=Sensitivity+ss.se), colour="black", width=.1)+
    labs(x = "", y = "Social sensitivity")+ scale_y_continuous(expand = c(0,0), limits = c(0,.8))+theme_classic(base_size = 15)+
>>>>>>> ec1b050b332c8ef38ba1ef2c92d53049783b4537
    theme(legend.text = element_text(size = 8),
          legend.spacing.y = element_blank(), 
          legend.background = element_rect(fill = NA),
          panel.grid.major.y = element_line(linewidth = .1),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
          axis.line = element_line(linewidth = .1), 
          axis.ticks = element_line(linewidth = .1)))

<<<<<<< HEAD
(ac_bar <- ggplot(data = data_gg)+geom_point(aes(x=reorder(villages,AC), y=AC), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    labs(x = "Villages", y = "Adaptive capacity")+ scale_y_continuous(expand = c(0,0), limits = c(0.35,.5))+theme_classic(base_size = 15)+
=======
(ac_bar <- ggplot(data = data_gg)+geom_point(aes(x=reorder(Village,`Adaptive capacity`), y=`Adaptive capacity`), size = 1.5, position = position_dodge2(width =.5), stroke = .2)+
    labs(x = "Villages", y = "Adaptive capacity")+ scale_y_continuous(expand = c(0,0), limits = c(0.2,.8))+theme_classic(base_size = 15)+
>>>>>>> ec1b050b332c8ef38ba1ef2c92d53049783b4537
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
<<<<<<< HEAD
yR <- range(data_$SS);xR <- range(data_$AC)
lgd <- expand.grid(x = seq(0,1,diff(xR)/500), y = seq(0,1,diff(yR)/500)) %>% mutate(brk1 = ntile(x,2), brk2 = ntile(y,2),brks = paste(brk1, brk2))
(biplot <- ggplot()+
    geom_raster(data = lgd, aes(x = x, y = y, fill = brks))+scale_fill_viridis_c()+
    scale_fill_manual(values = c("1 1"="#d3d3d3", "1 2"="dodgerblue4", "2 1"="yellow4", "2 2"="red4"))+
    geom_point(data = data_, aes(x = AC, y = SS),size = 3, stroke = .2) +
    #geom_text_repel(data = data_, aes(x = AC, y = SS, label = villages),size = 3, colour = "white")+
    scale_y_continuous(name = "Social sensitivity",expand = c(0,0))+
    scale_x_continuous(name = "Adaptive capacity",expand = c(0,0))+ 
=======
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
>>>>>>> ec1b050b332c8ef38ba1ef2c92d53049783b4537
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
                  round(mean(y,na.r=T), 2), "\n",
                  round(median(y,na.r=T), 2), "\n"
    )
  ))
}

<<<<<<< HEAD


(SS<-ggplot(na.omit(data_), aes(x = reorder(villages,SS), y = SS, fill = SS)) +
  geom_boxplot() +
  stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9) +
  theme_classic() + geom_jitter() + theme(legend.position="none") + xlab("Villages") + scale_y_continuous(limit= c(0, 1), expand = c(0,0)))
=======
#Sensitivity
SS<-ggplot(data_, aes(x = reorder(Village,Sensitivity), y = `Sensitivity`, fill = `Sensitivity`)) +
  geom_boxplot( outlier.colour = "red",outlier.fill = "red") +
  #stat_summary(fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 1) +
  theme_classic() + geom_jitter() +
  #ylim(0.2, 8)+
  theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Villages") 
  #geom_signif(test="wilcox.test", comparisons = combn(levels(as.factor(data_$Village)),2, simplify = F)[-4],step_increase = 0.2)
  #stat_compare_means(comparisons=as.factor(data_$Village), method="wilcox.test", label="p.signif", color="red")
>>>>>>> ec1b050b332c8ef38ba1ef2c92d53049783b4537

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
