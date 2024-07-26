############################################################################################################################

# Supporting workflow for "Rising climate risk and loss and damage to coastal small-scale fisheries livelihoods"
# Authors: Maina, Asamoah, et al. 202x
# Journal: Nature Sustainability
# School of Natural Sciences, Macquarie University, Sydney, Australia.
# Codes by: asamoahfrt@gmail.com
# Last updated: 9/01/2024

############################################################################################################################

# SUPPORTING SCRIPT 1: ESTIMATING CLIMATE RISK FOR VILLAGES
# Steps in this script:
#  1. Load, understand and prepare the dataset for analysis.
#  2. Extract climate data and ecosystem data to AOO.
#  3. Deduce and transform data for use in the integrative framework.
#  3. Estimate climate impacts based on Choquet integral
#  5. Extract data to
# .....

# Clear the workspace and install/load necessary packages
rm(list = ls()); gc()

needed_packages <- c("dplyr", "tidyterra", "terra", "sf", "ggplot2")
new.packages <- needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages, new.packages)

# Import some important functions
inormal <- function(x) {
  qrank <- ((rank(x, na.last = TRUE) - 0.5) / sum(!is.na(x)))
  z_score <- scales::rescale(sqrt(2) * pracma::erfinv(2 * qrank - 1))
  return(z_score)
}

# Import PU and convert to a spatial object
wio.aoo <- readRDS("2_Data/sheet/wio.aoo.rds")
grdSize <- (25 * 25 * 1e4)

# Filter data
wooFilter <- wio.aoo %>%
  mutate(across(CoralExt:Cropland, ~ .x / grdSize)) %>%
  rowwise() %>%
  mutate(Tot = 100 * sum(c(CoralExt, seagrassExt, mangroveExt, Cropland), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Tot > 0.1) %>%
  select(ID)

wio.aoo.sub <- filter(wio.aoo, ID %in% wooFilter$ID)
wio.aoo.spdf <- st_as_sf(wio.aoo.sub, coords = c('x', 'y'), crs = "+proj=longlat") %>% vect()
plot(wio.aoo.spdf)

# Load climate metrics
clim.nc <- list.files("./2_Data/raster", pattern = '*.nc', full.names = TRUE)

# Extract variable names and scenarios
nc <- expand.grid(x = c(126, 245, 370, 585), y = c(2050, 2100)) %>%
  arrange(desc(-x)) %>%
  mutate(cbn = paste(x, y, sep = "_")) %>%
  select(cbn) %>%
  pull()

chronic <- c("evspsbl", "npp", "pH_trend", "tap_trend", "sst_trend", "ts_trend", "SLR")
acute <- c("cdd", "r10p", "sst90p", "ts90p")
varlst <- c(chronic, acute)

# Load and plot hazard data
hazards <- lapply(varlst, function(var) {
  rr <- rast(clim.nc[grep(var, clim.nc)])
  names(rr) <- paste(var, nc, sep = "_")
  return(rr)
}) %>% rast()

plot(hazards)

# Extract raw hazards
climdata <- terra::extract(hazards, wio.aoo.spdf, xy = TRUE) %>% as.data.frame()
climdata <- cbind(as.data.frame(wio.aoo.spdf), climdata[,-1])
climdata <- climdata %>%
  relocate(c(x, y), .before = ID) #%>% cbind(DMwR2::knnImputation(climdata[15:ncol(climdata)], k = 1))

# Inverse normal standardization
source("1_Codes/qTransform.R")
HazardQNormed <- qTranform(climdata, vlst = varlst, time = c(2050), scenario = c(245, 370, 585))

# Load tropical cyclone data
slr.tc.data <- rast("./2_Data/raster/TC_count.tif")
slr.tc.data <- cbind(ID = wio.aoo.spdf$ID, (terra::extract(slr.tc.data, wio.aoo.spdf) %>% as.data.frame())[-1])
colnames(slr.tc.data)[colnames(slr.tc.data) == "layer"] <- "TC"

# Standardize and merge tropical cyclone data
slr.tc.data <- slr.tc.data %>%
  mutate(across(starts_with("TC"), ~ scales::rescale(.x))) %>%
  select(ID, tc_245_2050 = TC, tc_370_2050 = TC, tc_585_2050 = TC)

all.climdata <- merge(HazardQNormed, slr.tc.data, by = "ID")
all.climdata <- merge(climdata[, 1:3], all.climdata, by = "ID") %>% as.data.frame()

# Convert to spatial object and plot
# all.climdata <- all.climdata %>%
#   filter(!is.na(y), !is.na(x))
# 
# normed_hazards_ <- st_as_sf(all.climdata, coords = c("x", "y"), crs = st_crs(4326))
# plot(normed_hazards_[15])

# Normalize exposed systems metrics and save
ecosystem_exposure <- climdata[, 3:14] %>%
  mutate(across(c(CoralExt:Cropland), log10),
         across(c(CoralExt:Cropland), scales::rescale),
         across(c(Nb_sp:FOri), scales::rescale))
colnames(ecosystem_exposure) <- c("ID", paste0("std_", colnames(ecosystem_exposure[-1])))

all.climdata <- merge(all.climdata, ecosystem_exposure, by = "ID")
saveRDS(all.climdata, "2_Data/sheet/all.climdata.rds")

# Aggregate data
all.data <- readRDS("2_Data/sheet/all.climdata.rds")
namelist <- colnames(all.data)
fuzzyS <- function(x) {return(1 - prod((1 - x)))}

# Combine corals and fish data
coralSys <- all.data %>%
  select(std_Nb_sp, std_FRic, std_FDiv, std_FEve) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), 0, .x))) %>%
  rowwise() %>%
  mutate(std_fish.div = mean(c(std_Nb_sp, std_FRic, std_FDiv, std_FEve))) %>%
  ungroup()

all.data <- cbind(all.data, coralSys["std_fish.div"])

# Compute impacts
mu <- readRDS("2_Data/sheet/mu.rds")
exposure <- c("std_seagrassExt", "std_mangroveExt", "std_Cropland", "std_CoralExt", "std_fish.div")

impacts.aoo <- lapply(exposure, function(exp) {
  # Calculate ecosystem impact (exposure * hazard) for year 2050
  df <- cbind(ID = all.data["ID"], all.data[,exp] * all.data[grep("_2050", colnames(all.data))])
  df <- df %>% mutate(across(-ID, ~ ifelse(is.na(.x), 0, .x)))
  
  # Compound impact across for each climate scenario
  df$impact_ssp245.2050 <- apply(df[grep("_245_2050", colnames(df))], 1, function(x) kappalab::Choquet.integral(mu,x))
  df$impact_ssp370.2050 <- apply(df[grep("_370_2050", colnames(df))], 1, function(x) kappalab::Choquet.integral(mu,x))
  df$impact_ssp585.2050 <- apply(df[grep("_585_2050", colnames(df))], 1, function(x) kappalab::Choquet.integral(mu,x))
  
  df <- df %>% dplyr::select(ID, impact_ssp245.2050,impact_ssp370.2050, impact_ssp585.2050)
  names(df) <- c("ID", paste(exp, colnames(df[-1]), sep = "_"))
  return(df) 
  }  ) %>% purrr::reduce(left_join, by = "ID")

(impacts.agg <- impacts.aoo %>% rowwise() %>%
               mutate(imp.ssp245.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp245", colnames(impacts.aoo))])),
                      imp.ssp370.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp370", colnames(impacts.aoo))])),
                      imp.ssp585.2050 = fuzzyS(c_across(colnames(impacts.aoo)[grep("_ssp585", colnames(impacts.aoo))]))))

impacts.agg <- merge(impacts.agg, climdata[, c("ID", "CoralExt", "seagrassExt", "mangroveExt", "Cropland")], by = "ID")




# Merge grid level impacts to the network
idw.matrix <- data.table::fread("2_Data/sheet/idw.dist.matrix.csv", stringsAsFactors = TRUE, encoding = "UTF-8") %>%
  filter(EucDist > 0) %>%
  mutate(inverseDist = 1 / (EucDist / 1e6)^2) %>%
  filter(inverseDist > 0.1)

idw.matrix <- merge(idw.matrix, impacts.agg, by.x = "nbr", by.y = "ID")

idw.impacts <- idw.matrix %>%
  group_by(src) %>%
  summarise(across(where(is.numeric), ~ sum(.x * inverseDist, na.rm = TRUE) / sum(inverseDist, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(ID = src) %>%
  select(-src, -EucDist, -inverseDist)


villageImpacts <- merge(data.table::fread("2_Data/sheet/village.grid.id.csv", stringsAsFactors = TRUE, encoding = "UTF-8"), 
                        idw.impacts, by = "ID")


# Plot impacts
data <- villageImpacts %>%
  select(Country, Villages, starts_with("std_")) %>%
  pivot_longer(-c(Country, Villages), names_to = "ecosystems", values_to = "impacts") %>%
  mutate(ISO3 = case_when(
    Country == "Kenya" ~ "KEN",
    Country == "Mozambique" ~ "MOZ",
    Country == "Madagascar" ~ "MDG",
    Country == "Tanzania" ~ "TZA"
  ))

data <- filter(data, Scenario == "SSP5-8.5")

data_long <- data %>%
  group_by(Villages) %>%
  mutate(totImpact = sum(impacts, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_longer(-c(Country, Villages, ISO3, Scenario), names_to = "ecosystem", values_to = "impacts") %>%
  group_by(Country, ecosystem) %>%
  summarise(avg_impact = mean(impacts, na.rm = TRUE), sd_impact = sd(impacts, na.rm = TRUE)) %>%
  ungroup()

ggplot(data_long, aes(x = reorder(ecosystem, avg_impact), y = avg_impact, fill = Country)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = avg_impact - sd_impact, ymax = avg_impact + sd_impact), width = 0.2, position = position_dodge(0.9)) +
  labs(title = "Impacts of Climate Change on Different Ecosystems by Country", x = "Ecosystem", y = "Average Impact") +
  theme_minimal()



# Load socio-economic data
socioecom <- read.csv("2_Data/sheet/3_SocialVulnerability/SocialDataAll.csv")
socioecom <- socioecom %>%
  mutate(ISO3 = case_when(
    Country == "Kenya" ~ "KEN",
    Country == "Tanzania" ~ "TZA",
    Country == "Madagascar" ~ "MDG",
    Country == "Mozambique" ~ "MOZ",
    TRUE ~ Country
  ))

I_Ctrl <- read.csv("2_Data/sheet/3_SocialVulnerability/impact.control.csv")
socioecom <- merge(socioecom, I_Ctrl, by = "ISO3")

# Plot Adaptive Capacity vs. Impact Control (2020)
plot(socioecom$AdaptiveCapacity, exp(-1 * socioecom$ic2020))

# Prepare socio-economic data for plotting
socioecom$VillNation <- paste(socioecom$Villages, paste("(", socioecom$ISO3, ")", sep = ""))
socioecom$Vulnerable <- scales::rescale((socioecom$Sensitivity / socioecom$AdaptiveCapacity) * exp(-1 * socioecom$ic2020), to = c(0.01, 1))

# Plot social vulnerability
ggplot(data = socioecom, aes(x = "XS", y = Vulnerable, label = VillNation)) +
  geom_boxplot(linewidth = 0.3) +
  geom_point(aes(colour = ISO3), position = position_jitter(width = .1, height = 0)) +
  scale_colour_manual(name = "", values = c("KEN" = "darkred", "MDG" = "yellow", "MOZ" = "dodgerblue4", "TZA" = "grey")) +
  geom_text_repel(size = 1) +
  scale_y_continuous(name = "Social vulnerability", expand = c(0, 0), limits = c(0, 1.02), breaks = seq(0, 1, .2)) +
  labs(x = "") +
  theme_classic(base_size = 12) +
  guides(shape = "none", colour = "none") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(1, 'cm'),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.height = unit(.1, 'cm'),
    legend.key.width = unit(.2, 'cm'),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(linewidth = .1),
    axis.ticks.y = element_line(linewidth = .1)
  )
ggsave("3_Outputs/plots/FigS2.png", width = 4, height = 5, dpi = 1200)

# Import and merge climate impact data
riskMaster <- merge(socioecom, villageImpacts, by = "Villages")

# Plot impact vs. socio-economic factors
plot(riskMaster[, c("imp.ssp585.2050", "Livelihood", "Demography", "Cultural", "Health", "Learning", "Assets", "Flexibility", "Agency", "Organisation")])

# Prepare data for risk space plot
df <- rbind(
  data.frame(sce = "SSP2-4.5", impact = riskMaster$imp.ssp245.2050, Vulnerability = riskMaster$Vulnerable, village = riskMaster$Villages, ISO3 = riskMaster$ISO3),
  data.frame(sce = "SSP5-8.5", impact = riskMaster$imp.ssp585.2050, Vulnerability = riskMaster$Vulnerable, village = riskMaster$Villages, ISO3 = riskMaster$ISO3)
)

# Define risk space legend
yR <- range(df$impact)
xR <- range(df$Vulnerability)
lgd <- expand.grid(x = seq(0, 1, diff(xR) / 150), y = seq(0, 1, diff(yR) / 150)) %>%
  mutate(
    mxCol = y * x,
    brks = ntile(mxCol, 4),
    brks = case_when(
      brks == 1 ~ "Low",
      brks == 2 ~ "Medium",
      brks == 3 ~ "High",
      brks == 4 ~ "Very high"
    )
  )
lgd$brks <- factor(lgd$brks, levels = c("Low", "Medium", "High", "Very high"))

# Plot risk space
riskspace <- ggplot() +
  geom_raster(data = lgd, aes(x = x, y = y, fill = brks)) +
  scale_fill_manual(values = c("Low" = "#d3d3d3", "Medium" = "#a88283", "High" = "#7e433e", "Very high" = "#551601")) +
  geom_point(data = df, aes(x = Vulnerability, y = impact, shape = sce, colour = ISO3), size = .8, stroke = .2) +
  labs(y = "Climate change impacts", x = "") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, .2), labels = c("0", ".2", ".4", ".6", ".8", "1")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, .2), labels = c("0", ".2", ".4", ".6", ".8", "1")) +
  theme_bw(base_size = 8) +
  scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP5-8.5" = 17)) +
  scale_colour_manual(values = c("KEN" = "darkred", "MDG" = "yellow", "MOZ" = "dodgerblue4", "TZA" = "grey50")) +
  guides(shape = "none") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(1, 'cm'),
    legend.text = element_text(size = 5),
    legend.key.height = unit(.1, 'cm'),
    legend.key.width = unit(.2, 'cm'),
    axis.ticks = element_line(linewidth = .1),
    panel.border = element_blank()
  )

# Plot option space
optSpace <- ggplot() +
  geom_raster(data = lgd, aes(x = x, y = y, fill = brks)) +
  scale_fill_manual(values = c("Low" = "#d3d3d3", "Medium" = "#a88283", "High" = "#7e433e", "Very high" = "#551601")) +
  labs(y = "", x = "") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, .1), position = "right") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, .2)) +
  theme_bw(base_size = 8) +
  guides(shape = "none", colour = "none") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.background = element_rect(fill = NA),
    legend.key.size = unit(1, 'cm'),
    legend.text = element_text(size = 8),
    legend.key.height = unit(.1, 'cm'),
    legend.key.width = unit(.2, 'cm'),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Combine risk space and option space plots
grobs1 <- ggplotGrob(riskspace)$grobs
legend1 <- grobs1[[which(sapply(grobs1, function(x) x$name) == "guide-box")]]
rrSpace <- ((riskspace + theme(legend.position = "none") | optSpace + theme(legend.position = "none")) / legend1) + plot_layout(heights = c(2, .1))
ggsave(plot = rrSpace, "3_Outputs/plots/Fig2a.png", dpi = 1200, height = 3, width = 4)



# Estimation of Residual Climate Risk
# First, the script calculates residual climate risk for different scenarios (SSP2-4.5, SSP3-7.0, SSP5-8.5) using vulnerability and impact data.
riskMaster <- riskMaster %>% 
  mutate(risk585 = (imp.ssp585.2050*Vulnerable),
         risk370 = (imp.ssp370.2050*Vulnerable),
         risk245 = (imp.ssp245.2050*Vulnerable))

# Visualization of Risk Scores
# The script then visualizes the risk scores for villages under different scenarios using a scatter plot.
# The plot is created using ggplot2, displaying risk scores for each village, colored by risk level (Low, Medium, High, Very high).
df <- rbind(data.frame(sce = "SSP2-4.5", risk = (riskMaster$risk245), village = riskMaster$VillNation, ISO3 = riskMaster$ISO3),
            data.frame(sce = "SSP5-8.5", risk = (riskMaster$risk585), village = riskMaster$VillNation, ISO3 = riskMaster$ISO3))

plt2 <- ggplot(data = df)+
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
        axis.ticks = element_line(linewidth = .1))
ggsave(plot = plt2, "3_Outputs/plots/Fig2b.png", dpi = 1200, height = 3, width = 4)


# 3. Box Plot of Climate Risk by Country
# Next, the script visualizes the climate risk for each country using a box plot.
# The plot shows risk scores across different scenarios.

df <- rbind(data.frame(sce = "SSP2-4.5", MN = riskMaster$risk245, ISO3 = df$ISO3),
            data.frame(sce = "SSP3-7.0", MN = riskMaster$risk370, ISO3 = df$ISO3),
            data.frame(sce = "SSP5-8.5", MN = riskMaster$risk585, ISO3 = df$ISO3)) %>% group_by (ISO3) %>% mutate(sortMag = mean(MN))

plt3 <- ggplot(data = df)+
  geom_rect(fill = "grey90", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = .33)+
  geom_rect(fill = "grey80", xmin = -Inf, xmax = Inf, ymin = .33, ymax = .66)+
  geom_rect(fill = "grey70", xmin = -Inf, xmax = Inf, ymin = .66, ymax = Inf)+
  geom_boxplot(aes(x = reorder(ISO3, sortMag), y = MN, fill = sce),linewidth = .1, position = position_dodge(width =.8))+
  scale_fill_manual(name="", values = c("SSP2-4.5" = "#749B58FF", "SSP3-7.0"="#466983FF","SSP5-8.5" = "#CE3D32FF"))+
  labs(y = "", x="", title = "")+
  scale_y_continuous(name = "Climate risk [index]", expand = c(0,0), limits = c(0,1), breaks = seq(0,1,.2))+
  theme_classic(base_size = 12)+
  theme(legend.position = c(0.15,.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 0),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = .1),
        axis.ticks = element_line(linewidth = .1))
ggsave(plot = plt3, "3_Outputs/plots/FigS3.png", dpi = 1200, height = 4, width = 4)



# 4. Economic Valuation
# The script then reads in economic valuation data and calculates the total economic value (TEV) and its residual value after accounting for climate risk.
econValues <- readr::read_csv("2_Data/sheet/4_EconomicValuations/Unit_Service_Values_Villages.csv")[,c(2,11,12,13,14)]
colnames(econValues) <- c("Villages","CoralUnitValue","CropUnitValue","MangroveUnitValue","SeagrassUnitValue")
tev.data <- merge(riskMaster, econValues, by = "Villages")

tev.data$TEV = ((tev.data[,"CoralExt"]*tev.data[,"CoralUnitValue"])+
                  (tev.data[,"seagrassExt"]*tev.data[,"SeagrassUnitValue"])+
                  (tev.data[,"mangroveExt"]*tev.data[,"MangroveUnitValue"])+
                  (tev.data[,"Cropland"]*tev.data[,"CropUnitValue"])
                )/1e4

r = 1+0.024 #https://www.cbo.gov/publication/58957
plot(tev.data$TEV/1e6, (1-tev.data$risk585))

sum(tev.data$TEV, na.rm = TRUE)

tev.data$fTEV_ssp585 = r*(tev.data$TEV*(1-exp(-tev.data$risk585)))
sum(tev.data$fTEV_ssp585);mean(tev.data$fTEV_ssp585);sd(tev.data$fTEV_ssp585)
log10(sum(tev.data$fTEV_ssp585, na.rm = TRUE)/sum(tev.data$TEV, na.rm = TRUE))

tev.data$fTEV_ssp370 = r*(tev.data$TEV*(1-exp(-tev.data$risk370)))
sum(tev.data$fTEV_ssp370);mean(tev.data$fTEV_ssp370);sd(tev.data$fTEV_ssp370)
log10(sum(tev.data$fTEV_ssp370, na.rm = TRUE)/sum(tev.data$TEV, na.rm = TRUE))

tev.data$fTEV_ssp245 = r*(tev.data$TEV*(1-tev.data$risk245))
sum(tev.data$fTEV_ssp245);mean(tev.data$fTEV_ssp245);sd(tev.data$fTEV_ssp245)
log10(sum(tev.data$fTEV_ssp245, na.rm = TRUE)/sum(tev.data$TEV, na.rm = TRUE))

tev.data$perc_tev_585 <- (tev.data$TEV-tev.data$fTEV_ssp585)/tev.data$TEV
tev.data$perc_tev_245 <- (tev.data$TEV-tev.data$fTEV_ssp245)/tev.data$TEV

#Order of magnitude change
lnd_dat <- tev.data[,c("ISO3", "Villages","Sensitivity","AdaptiveCapacity","imp.ssp245.2050","imp.ssp370.2050","imp.ssp585.2050","risk585","risk370","risk245","TEV","fTEV_ssp585","fTEV_ssp370","fTEV_ssp245")]
write.csv(lnd_dat, "3_Outputs/sheets/SupplementalTable4.csv", row.names = F)

#import bivariate codes
source("1_Codes/ColMatrix.R")

# Define the number of breaks
nBreaks <- 20
# Create the colour matrix
col.matVul <- colmat(nbreaks = nBreaks, breakstyle = "quantile",
                     xlab = "X", ylab = "Y",
                     upperleft = "#4279b0",
                     upperright = "#9e3547",
                     bottomleft =  "#d3d3d3",
                     bottomright = "#311e3b",
                     saveLeg = FALSE, plotLeg = TRUE)
# Retrieve bivariate colour pallet data
lgdBiv <- BivLegend$data; names(lgdBiv) <- c("binY", "binX", "BivCol", "UID")

#Plot
ymn <- median(tev.data$TEV/1e6, na.rm=TRUE)
xR <- range(tev.data$risk585);yR <- range(tev.data$TEV/1e6)
lgd <- expand.grid(
  y = seq(0,50, diff(yR)/150),
  x = seq(0,1, diff(xR)/150)
) %>% mutate(
  binY = ntile(y,20),
  binX = ntile(x,20)
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
    labs(x = "Climate risk index", y = "Total economic value \n(Million US$/year)")+
    #scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))+
    theme_bw(base_size = 15)+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,50,5), limits = c(0,50))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.2), limits = c(-0.01,1))+
    scale_shape_manual(values = c("SSP5-8.5" = 2))+
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
# Adaptation gap plot as radial bar chart 
###############################
# Libraries
library(dplyr)
library(tidyverse)
library(hrbrthemes)
bar_width <- .9 # default width of bars in geom_bar

# Plot
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




