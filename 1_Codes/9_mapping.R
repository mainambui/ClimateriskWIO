kwale = st_read('vilages/Kwale.shp')
africa = rnaturalearth::ne_countries(continent = 'africa', returnclass = 'sf')
africa2 = vect(africa)

# create a list to store output data frames
cdd_list <- list()

# Define the file names
cdd <- c("cdd126_ESM25km.nc", "cdd245_ESM25km.nc", "cdd370_ESM25km.nc", "cdd585_ESM25km.nc")
evs <- c("evspsbl126_ESM25km.nc", "evspsbl245_ESM25km.nc", "evspsbl370_ESM25km.nc", "evspsbl585_ESM25km.nc")
npp <- c("npp126_ESM25km.nc", "npp245_ESM25km.nc", "npp370_ESM25km.nc", "npp585_ESM25km.nc")
ph <- c("ph126_ESM25km.nc", "ph245_ESM25km.nc", "ph370_ESM25km.nc", "ph585_ESM25km.nc")
pr <- c("pr126_ESM25km.nc", "pr245_ESM25km.nc", "pr370_ESM25km.nc", "pr585_ESM25km.nc")
sst <- c("sst126_ESM25km.nc", "sst245_ESM25km.nc", "sst370_ESM25km.nc", "sst585_ESM25km.nc")
ts <- c("ts126_ESM25km.nc", "ts245_ESM25km.nc", "ts370_ESM25km.nc", "ts585_ESM25km.nc")
tx <- c("tx90126_ESM25km.nc", "tx90245_ESM25km.nc", "tx90370_ESM25km.nc", "tx90585_ESM25km.nc")

files <- c(tx)


# loop over the file names
for (file in files) {
  # read in raster object
  cdd_rast <- raster(file)
  # mask raster object
  cdd_masked <- raster::mask(cdd_rast, kwale)
  # add masked raster object to list
  cdd_list[[file]] <- cdd_masked
}

# combine list of raster objects into a single stack
tx_stack <- stack(cdd_list)


pr_stack <- cdd_stack*86400
evs_stack = (cdd_stack*31536000)/1000000

# convert stack to data frame
tx_df <- raster::as.data.frame(cdd_stack, xy = TRUE, na.rm = TRUE)
colnames(tx_df) <- c('Lon','Lat', 'TX90 SSP126','TX90 SSP245','TX90 SSP370','TX90 SSP585')
colnames(cdd_df) <- c('CDD SSP126','CDD SSP245','CDD SSP370','CDD SSP585')




library(ggplot2)
library(tidyr)

cdd.breaks <- c(10, 20, 30, 40, 50, 70, 90, 120, 150)
tx.breaks <- c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
ts.breaks <- c(10, 13, 16, 19, 22, 25, 28, 31, 34)
pr.breaks <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)


splashdown_cols <- c("#A56723", "#FF9933", "#FFC847", "#FFF67A", "#EDFBCD", "#B8E6D9", "#84D1D4", "#5398C4", "#13378E")
#splashdown_cols <- c("#FF9933", "#ffc02e", "#FFC847", "#FFF67A", "#EDFBCD", "#B8E6D9", "#84D1D4", "#5398C4", "#13378E")

# convert data from wide to long format
tx_long <- gather(tx_df, key = "variable", value = "value", -Lon, -Lat)

# plot the data
ggplot(data = tx_long, aes(x = Lon, y = Lat, fill = value)) + 
  geom_raster() +
  geom_sf(data = kwale, inherit.aes = FALSE, fill = NA) +
  scale_fill_stepsn(colours = splashdown_cols, breaks= tx.breaks, limits=c(20,50), name= 'mm/yr') +
  #scale_fill_gradient(name = "CDD", low = "white", high = "red") +
  facet_wrap(~variable, ncol = 2) + 
  theme(text=element_text(size=20)) + theme(legend.key.height =unit(1, "cm")) + 
  labs(title = "", subtitle = "", x = "Longitude", y = "Latitude") + 
  guides(fill = guide_legend(title.position = "right",direction = "vertical",
                             title.theme = element_text(angle = 90, size = 14, colour = "black"),
                             barheight = .3, barwidth = .95,
                             title.hjust = 0.5, raster = FALSE,
                             title = expression('mm/day'))) +
  theme_minimal()

ggsave("tx.png", bg = 'white', width = 25, height = 25, units = "cm")
