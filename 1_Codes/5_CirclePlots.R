# Climate Risks Analysis for WIO
# Written by Ernest Frimpong Asamoah (asamfrt@gmail.com)
library(ncdf4);library(sp);library(tidyverse);library(sf)
library(raster) #some functions conflict with dplyr

rm(list = ls())
# = "C:/Users/MQ45019738/Dropbox/6_WIO_CCVA/WIOProjects/Data/"
db.dir = "~/Documents/Mygitprojects/ClimateriskWIO/2_Data/sheet/" #jm local git dor


#Import SE data
#socioecom <- read_csv(paste0(db.dir, "3-Vulnerability//SocialDataAll.csv"))
socioecom <- read_csv(paste0(db.dir, "3_SocialVulnerability/SocialDataAll.csv"))
socioecom <- filter(socioecom, !is.na(x))

##############################
#FIGURE 3: PLOT CLIMATE POLAR
##############################

CircData <- socioecom %>% group_by(Country) %>%
  summarise(Livelihood = mean(Livelihood, na.rm=TRUE),
            Demography = mean(Demography, na.rm=TRUE),
            Cultural = mean(Cultural, na.rm=TRUE),
            Health = mean(Health, na.rm=TRUE),
            Learning = mean(Learning, na.rm=TRUE),
            Assets = mean(Assets, na.rm=TRUE),
            Flexibility = mean(Flexibility, na.rm=TRUE),
            Agency = mean(Agency, na.rm=TRUE),
            Organisation = mean(Organisation, na.rm=TRUE)) %>% ungroup() %>%
  mutate(Country = ifelse(Country=="Kenya", "KEN", Country),
         Country = ifelse(Country=="Mozambique", "MOZ", Country),
         Country = ifelse(Country=="Madagascar", "MDG", Country),
         Country = ifelse(Country=="Tanzania", "TZA", Country))

CircData <- rbind(
  data.frame(ISO3=CircData$Country, value=CircData$Livelihood, varr="Livelihood", cat="AC"),
  data.frame(ISO3=CircData$Country, value=CircData$Demography, varr="Demography", cat="AC"),
  data.frame(ISO3=CircData$Country, value=CircData$Cultural, varr="Cultural", cat="AC"),
  data.frame(ISO3=CircData$Country, value=CircData$Health, varr="Health", cat="AC"),
  data.frame(ISO3=CircData$Country, value=CircData$Learning, varr="Learning", cat="SS"),
  data.frame(ISO3=CircData$Country, value=CircData$Assets, varr="Assets", cat="SS"),
  data.frame(ISO3=CircData$Country, value=CircData$Flexibility, varr="Flexibility", cat="SS"),
  data.frame(ISO3=CircData$Country, value=CircData$Agency, varr="Agency", cat="SS"),
  data.frame(ISO3=CircData$Country, value=CircData$Organisation, varr="Organisation", cat="SS"))

CircDataSub <- filter(CircData, cat=="AC")
rm(data)
data <- CircDataSub

data$varr <- as.factor(data$varr)
data$ISO3 <- as.factor(data$ISO3)

data <- data %>% group_by(varr, ISO3) %>% mutate(tot=sum(value)) %>% ungroup()
data <- data %>% arrange(ISO3, tot) %>% 
  add_row(ISO3 = rep(unique(data$ISO3), 1)) %>% # add a couple of empty rows to separate countries
  arrange(ISO3) %>% ungroup()
data$id <- rep(seq(1, nrow(data)/nlevels(as.factor(data$cat))), each=nlevels(as.factor(data$cat)))

# Get the name and the y position of each label
label_data <- data %>% group_by(id, varr) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar# I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

brks <- data %>% 
  group_by(ISO3) %>% 
  summarize(start = min(id), end = max(id) - 2) %>% 
  rowwise() %>% 
  mutate(title = mean(c(start, end))) %>% 
  ungroup() %>% 
  mutate(end = data.table::shift(end + 1, n = 1, type = "shift", fill = max(end) + 1), start = start -1) 

brks$start[brks$ISO3 == "KEN"] <- -1
brks$end[brks$ISO3 == "TZA"] <- 0

max_value <- max(data$tot, na.rm = T); y_max <- max_value+.1; v <- c(0,.1,.2,.3,.4,.50)
brks <- brks %>% mutate(v = list(v)) %>% unnest()

(p <- data %>% 
    ggplot() + geom_bar(aes(x = as.factor(id), y = value, fill = ISO3), position="stack", stat="identity", width = 0.5) +
    annotate("text", x = rep(max(data$id[data$ISO3 == "TZA"]), length(v)), y = v-.1, 
             label = paste0(head(v), "%"), color = "grey", size = 2.5, angle = 0, fontface = "bold", hjust = 0.7)+
    ylim(-.3, y_max)+
    scale_fill_manual(values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="grey50"))+
    #paletteer::scale_fill_paletteer_d(palette = "rcartocolor::Pastel")+
    theme_minimal(base_size = 15) +
    guides(fill = guide_legend(nrow = 1))+
    coord_polar()+
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(c(1,-1,-1,-1), "cm")
    )+
    geom_text(data = label_data, aes(x = id, y = tot+0.02, label = varr, hjust = hjust), color="black",  size = 3.5, angle = label_data$angle, inherit.aes = FALSE)+
    geom_text(data = brks, aes(x = title, y = -.08, label = ISO3),  colour = "black", alpha = 0.8, size = 3.5, fontface = "bold", inherit.aes = FALSE)
)
# cairo_pdf(filename = here::here("outputs", "SS.pdf"), width = 8.2, height = 8, fallback_resolution = 600)
# print(p)
# dev.off()
#ggsave(plot = p,"./outputs/AC.png", dpi = 1200, height = 4, width = 4)