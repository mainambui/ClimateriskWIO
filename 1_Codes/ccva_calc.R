#############################################################################################
# WIO CLIMATE CHANGE VULNERABIOLTIY ASSESSEMENTS
# author Ernest Frimpong Asamoah
# date: July 2023
############################################################################################

# library
library(ggplot2)
library(dplyr)
library(PCAmixdata) # for PCA mix
library(eatATA) # for  dummiesToFactor()
library(readxl)

#Clear global environment
rm(list=ls())

# data prep to calculate MSL
# Codes modified from Stephanie D'agata (update March 2022)
# outputs: MSL HH scores between 0 and 1
dat.n <- read_excel("data/q45_MSL.xlsx",sheet = "mdg_q45_with_count") # with counts
head(dat.n,1)
colnames(dat.n)
summary(dat.n)
class(dat.n)
names(dat.n)

# dummies to factor variables
# on MD data, number
dat.n <- dummiesToFactor(dat.n, dummies = c("electricity_solar","electricity_generator","electricity_grid","electricity_none","electricity_other"),facVar="electricity")
dat.n <- dummiesToFactor(dat.n, dummies = c("floor_material_dirt_soil","floor_material_wood","floor_material_concrete","floor_material_tile"),facVar="floor")
dat.n <- dummiesToFactor(dat.n, dummies = c("wall_material_bamboo_thatch","wall_material_wood","wall_material_metal","wall_material_cement","wall_material_other"),facVar="wall")
dat.n <- dummiesToFactor(dat.n, dummies = c("roof_material_bamboo_thatch","roof_material_wood","roof_material_metal","roof_material_tile","roof_material_other"),facVar="roof")
summary(dat.n)
names(dat.n)

# select quantity + vector variables for gower

rm(dat.n.gower)
dat.n.gower <- dat.n %>% 
  dplyr::select(survey_no,cooking_pots_nb,radio_nb,dvd_nb,mattresses_nb,mobile_phone_nb,smart_phone_nb,
                flushing_toilet_nb,indoor_piped_nb,computers_nb,cow_nb,sheep_nb,goat_nb,pig_nb,television_nb)

# hh.id to rownames
rownames(dat.n.gower) <- dat.n.gower$survey_no
dat.n.gower$survey_no <- NULL
summary(dat.n.gower)
dim(dat.n.gower)

# select numeric variables > 0
dat.n.gower <- as.data.frame(dat.n.gower) %>% dplyr::select_if(colSums(., na.rm=TRUE) != 0)

# scale numeric
dat.n.gower.sc <- scale(dat.n.gower)

# add factor variables
dat.n.gower.sc <- cbind(dat.n.gower.sc, dat.n[,c("electricity","floor","wall","roof")])
summary(dat.n.gower.sc)

# order factor variable
dat.n.gower.sc$electricity <- ordered(dat.n.gower.sc$electricity, levels = c("_none_","electricity_none","electricity_solar","electricity_grid"))
dat.n.gower.sc$floor <- ordered(dat.n.gower.sc$floor, levels = c("_none_","floor_material_dirt_soil","floor_material_wood","floor_material_concrete","floor_material_tile"))
dat.n.gower.sc$wall <- ordered(dat.n.gower.sc$wall, levels = c("wall_material_other","wall_material_bamboo_thatch","wall_material_wood","wall_material_metal","wall_material_cement"))
dat.n.gower.sc$roof <- ordered(dat.n.gower.sc$roof, levels = c("roof_material_bamboo_thatch","roof_material_metal"))

summary(dat.n.gower.sc)
dim(dat.n.gower.sc)

# TEST PCA MIX
X.quanti <- dat.n.gower.sc[,1:12]
X.quali <- dat.n.gower.sc[,13:16]
pca <- PCAmix(X.quanti,X.quali,ndim=20,rename.level=TRUE)

## output PCA.mix
pca$eig
pca$ind$coord

plot(pca, choice="sqload",main="Squared correlations")
plot(pca, choice="cor",main="Correlation circle")

pca$quali
pca$quanti$coord
pca$quanti$cos2

###### Create composite index with the first 10 axis (eig >1, cum variance > 65%)
material_style_of_life <- (pca$ind$coord[,1]*pca$eig[1,2]+
                             pca$ind$coord[,2]*pca$eig[2,2]+
                             pca$ind$coord[,3]*pca$eig[3,2]+
                             pca$ind$coord[,4]*pca$eig[4,2]+
                             pca$ind$coord[,5]*pca$eig[5,2]+
                             pca$ind$coord[,6]*pca$eig[6,2]+
                             pca$ind$coord[,7]*pca$eig[7,2]+
                             pca$ind$coord[,8]*pca$eig[8,2]+
                             pca$ind$coord[,9]*pca$eig[9,2]+
                             pca$ind$coord[,10]*pca$eig[10,2])/pca$eig[10,3]

material_style_of_life <- cbind.data.frame( survey_no = row.names(dat.n.gower.sc),
                               material_style_of_life)
material_style_of_life$material_style_of_life <- scales::rescale(material_style_of_life$material_style_of_life,to = c(0, 1))

#Check the distributions
hist(material_style_of_life$material_style_of_life, breaks = 30)


#############################################################################################################################################
#Calculate Climate Change Vulnerability
# by ERNEST FRIMPONG ASAMOAH
############################################################################################################################################

#import 
dat.m <- read_excel("data/mdg_data.xlsx", sheet = "scored")

#Merge MSL data back to the scored datasets
dat.m <- merge(dat.m, material_style_of_life, by = "survey_no")
class(dat.m)
str(dat.m)

names(dat.m)

#Sensitivy domains
livelihoods <- c("employment_status","percentage_of_catch_from_fishing_sold","percentage_of_income_from_the_main_activity","time_conducting_the_activity")
demography <- c("gender","years_living_in_the_village","percentage_of_children_in_the_family_members","family_dependency")
health <- c("age","food_security_and_wellbeing","sense_of_place" )
cultural <- c("appreciation_of_biodiversity","identity_and_pride","appreciation_of_lifestyle")

#adapative domains
learning <- c("education","knowledge_of_rules","access_to_information")
organisation <- c("community_cohesion","trust_in_organization","linking_social_capital")
agency <- c("level_of_participation","perceived_capacity_to_change","recognition_of_causality")
asset <- c("access_to_credit","community_infrastructures","material_style_of_life")
flexibility <- c("adapt_to_live_without_fishing","gear","livelihood_multiplicity","spatial_mobility")

#import
wght.ind <- read_excel("data/weights.xlsx", sheet = "indicator")
wght.dom <- read_excel("data/weights.xlsx", sheet = "domain")

#INDICATOR LEVEL WEIGTHING
indicators <- c(livelihoods,demography,health,cultural,learning,organisation,agency,asset,flexibility)
weighted.indicators = c()
for(var in indicators){
  weighted.indicators[var] <-  (dat.m[var] * (filter(wght.ind, indicator %in% var))[[2]])
}

#cbind all rows of the resultant list of vectors
weighted.indicators <- do.call(cbind, weighted.indicators) %>% as.data.frame()
#colnames(weighted.indicators) <- paste(indicators, "weighted", sep = "_")

#Now add all 
domains <- list(livelihoods,demography,health,cultural,learning,organisation,agency,asset,flexibility)

(domain.sum <- lapply(1:length(domains), function(x){
  mysub <- domains[[x]]
  xx <- rowSums(weighted.indicators[mysub], na.rm = TRUE)
  return(xx)
}))

domain.sum <- do.call(cbind, domain.sum) %>% as.data.frame()
colnames(domain.sum) <- c("livelihoods","demography","health","cultural","learning","organisation","agency","asset","flexibility")

## DOMAIN LEVEL WEIGHING
sensitivity <- c("livelihoods","demography","health","cultural")
adaptive.capacity <- c("learning","organisation","agency","asset","flexibility")

#loop over a list of domains and their respective weights
d.list <- c(sensitivity,adaptive.capacity)
weighted.domains = c()
for(d in d.list){
  weighted.domains[d] <-  (domain.sum[d] * (filter(wght.dom, domain %in% d))[[2]])
}
weighted.domains <- do.call(cbind, weighted.domains) %>% as.data.frame()

#Calculate sensitivity and adaptive capacity
weighted.domains$sensitivity <- rowSums(weighted.domains[sensitivity], na.rm = TRUE)
weighted.domains$adaptive.capacity <- rowSums(weighted.domains[adaptive.capacity], na.rm = TRUE)

weighted.domains <- cbind.data.frame(survey_no = dat.m["survey_no"], weighted.domains)

#Visualise the estimates
ggplot(data = weighted.domains)+
  geom_point(aes(sensitivity, adaptive.capacity), size = 3)

write.csv(weighted.domains, "data/ccva_data_out.csv")
