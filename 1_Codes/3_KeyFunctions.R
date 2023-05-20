#Some key functions
norm_scale = scales::rescale

#https://gist.github.com/variani/d6a42ac64f05f8ed17e6c9812df5492b
#inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
inormal <- function(x) {
  qrank <- ((rank(x, na.last = TRUE, ties.method= "random") - 0.5) / sum(!is.na(x)))
  z_score <- norm_scale(sqrt(2)*pracma::erfinv(2*qrank-1))
  return(z_score)
}

quant_transform <- function(df, var){
  mn <- min(df[[var]], na.rm = TRUE)
  if (mn<0) {
    df[[var]] <- df[[var]]+abs(mn)
  }else{df[[var]] <- df[[var]]}
  #Quantile normalisation
  df$score <- ((rank(df[[var]], na.last = TRUE, ties.method= "random") - 0.5) / sum(!is.na(df[[var]])))
  df$score <- norm_scale(sqrt(2)*pracma::erfinv(2*df[,"score"]-1))
  #Spatialize
  rr <- rasterFromXYZ(df[,c("x","y","score")])
  names(rr) <- paste0(var)
  crs(rr) <- "+proj=longlat"
  return(rr)
}


rs_min_max <- function(df, var){
  mn <- min(df[[var]], na.rm = TRUE)
  if (mn<0) {
    df[[var]] <- df[[var]]+abs(mn)
  }else{df[[var]] <- df[[var]]}

  df$score <- norm_scale(df[[var]])
  #Spatialise
  rr <- rasterFromXYZ(df[,c("x","y","score")])
  names(rr) <- paste0(var)
  crs(rr) <- "+proj=longlat"
  return(rr)
}

#Slope/Trends
slpFUN = function(k) {
  time <- 1:nlayers(k) 
  fun <- function(x) {if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}
  mn <- raster::calc(k, fun = fun)
  return(mn <- mn*(length(time)/12))
}

# normRaster <- function(x) {
#   x = x + minValue(x)
#   return((x-minValue(x)) /(maxValue(x)-minValue(x)))
# }

##Extreme temperature above the top 90th percentile of the baseline (1981-2010)
txpFUN <- function(r, nYrs, baseline, p, pn=c("near","far")){
  
  thresh <- raster::calc(baseline, fun = function(x){ quantile(x, probs = p, na.rm=TRUE)})
  names(thresh) <- "threshold"
  
  clim.years.thresh = nYrs
  period = match.arg(pn)
  if(period=="far"){
    end <- nlayers(r)
  } else{
    end <- 540
  }
  start = end-(clim.years.thresh*12 - 1)
  rr <- stack(thresh, r[[start:end]])
  rr <- as.data.frame(rr, xy = TRUE, na.rm = TRUE)
  N = ncol(rr)
  
  spdf <- rr %>% 
    mutate(dplyr::across(starts_with("X2"),~ ifelse(.x >= threshold, 1, 0)))%>%
    mutate(value = 100*(rowSums(.[3:N], na.rm = TRUE)/length(3:N))) %>% 
    dplyr::select(x,y,value)
  
  coordinates(spdf) <- ~ x + y
  gridded(spdf) <- TRUE
  spdf <- raster(spdf)
  
  crs(spdf) <- crs(r)
  return(spdf)
}


txIntFUN <- function(r, nYrs, baseline, p, pn=c("near","far")){
  
  thresh <- raster::calc(baseline, fun = function(x){ quantile(x, probs = p, na.rm=TRUE)})
  names(thresh) <- "threshold"
  
  clim.years.thresh = nYrs
  period = match.arg(pn)
  if(period=="far"){
    end <- nlayers(r)
  } else{
    end <- 540
  }
  start = end-(clim.years.thresh*12 - 1)
  rr <- stack(thresh, r[[start:end]])
  rr <- as.data.frame(rr, xy = TRUE, na.rm = TRUE)
  N = ncol(rr)
  
  spdf <- rr %>% 
    mutate(dplyr::across(starts_with("X2"),~ ifelse(.x >= threshold, 100*((.x-threshold)/threshold), NA)))%>%
    mutate(value = rowMeans(.[3:N], na.rm = TRUE)) %>% 
    dplyr::select(x,y,value)
  
  coordinates(spdf) <- ~ x + y
  gridded(spdf) <- TRUE
  spdf <- raster(spdf)
  
  crs(spdf) <- crs(r)
  return(spdf)
}


#Extreme rainfall below the 10th percentile of the baseline (1981-2010)
r95Freq <- function(r, nYears, baseR, p, pn=c("near","far")){
  
  thresh <- raster::calc(baseR, fun = function(x){ quantile(x, probs = p, na.rm=TRUE)})
  names(thresh) <- "threshold"
  
  clim.years.thresh = nYears
  period = match.arg(pn)
  if(period=="far"){
    end <- nlayers(r)
  } else{
    end <- 432 #the projected periods starts from 2015 and ends at 2100. So by 2050, 12*36yrs = 432 months
  }
  start = end-(clim.years.thresh*12 - 1)
  rr <- stack(thresh, r[[start:end]])
  rr <- as.data.frame(rr, xy = TRUE, na.rm = TRUE)
  n = ncol(rr)
  
  spdf <- rr %>% 
    mutate(dplyr::across(starts_with("X2"),~ ifelse(.x >= threshold, 1, 0)))%>%
    mutate(value = 100*(rowSums(.[3:n], na.rm = TRUE)/length(3:n))) %>% 
    dplyr::select(x,y,value)
  
  coordinates(spdf) <- ~ x + y
  gridded(spdf) <- TRUE
  spdf <- raster(spdf)
  
  crs(spdf) <- crs(r)
  return(spdf)
}

#Extreme rainfall above 95th percentile of the baseline (1981-2010)
r95IntFUN <- function(r, nYears, baseR, p, pn=c("near","far")){
  
  thresh <- raster::calc(baseR, fun = function(x){ quantile(x, probs = p, na.rm=TRUE)})
  names(thresh) <- "threshold"
  
  clim.years.thresh = nYears
  period = match.arg(pn)
  if(period=="far"){
    end <- nlayers(r)
  } else{
    end <- 432 #the projected periods starts from 2015 and ends at 2100. So by 2050, 12*36yrs = 432 months
  }
  start = end-(clim.years.thresh*12 - 1)
  rr <- stack(thresh, r[[start:end]])
  rr <- as.data.frame(rr, xy = TRUE, na.rm = TRUE)
  n = ncol(rr)
  
  spdf <- rr %>% 
    mutate(dplyr::across(starts_with("X2"),~ ifelse(.x <= threshold, NA, 100*((.x-threshold)/threshold))))%>%
    mutate(value = rowMeans(.[3:n], na.rm = TRUE)) %>% 
    dplyr::select(x,y,value)
  
  coordinates(spdf) <- ~ x + y
  gridded(spdf) <- TRUE
  spdf <- raster(spdf)
  
  crs(spdf) <- crs(r)
  return(spdf)
}

#Here impacts is conceptualised a correlative measure
# library("lavaan")
# c_stress <- ' chronic =~ std_id.1 + std_id.2 + std_id.3 '
# fit <- lavaan::sem(c_stress, data = stdIMpact1, std.lv=T)
# summary(fit)
# stdIMpact1$PredImp <- scales::rescale(mahalanobis(stdIMpact1[,3:5], colMeans(stdIMpact1[,3:5]), cov(stdIMpact1[,3:5])))

#Import the exposed system and regrid
# regrdFun <- function(df,var,shp){
#   spdf <- df %>% dplyr::select(x,y,var) %>% na.omit()
#   coordinates(spdf) <- ~ x + y
#   gridded(spdf) <- TRUE
#   spdf <- crop(raster(spdf), shp)
#   
#   crs(spdf) <- "+proj=longlat"
#   spdf <- raster::projectRaster(spdf, res = .25, method="bilinear", crs = "+proj=longlat")
#   #name(spdf)<-paste0(var)
# }


#Estimate potential impacts as function of hazard and exposure
#Normalise the datasets using quantile transformation
# rescale_QN_rr <- function(rr, nm){
#    #
#    df <- as.data.frame(rr[[nm]], xy=TRUE, na.rm=TRUE)
# 
#    #Quantile normalisation
#    df <- df[order(df[,nm], decreasing = FALSE),]
#    df$id <- 1:nrow(df)
#    df$score <- df$id/nrow(df)
#    df$z_score <- sqrt(2)*pracma::erfinv(2*df[,"score"]-1)
# 
#    df$std_id <- scales::rescale(df[, "z_score"])
#    df <- df %>% dplyr::select(x,y,std_id) %>% filter(is.finite(std_id))
# 
#   #Convert to raster
#    spdf <- df
#    coordinates(spdf) <- ~ x + y
#    gridded(spdf) <- TRUE
#    spdf <- raster(spdf)
# 
#    crs(spdf) <- "+proj=longlat"
# 
#    spdf <- raster::resample(spdf, tx90p_stack, method="bilinear")
#    names(spdf) <- paste0("std_", nm)
#    return(spdf)
#  }



# decayNormed <- function(x, rt, type = c("decay", "growth")){
#   
# }
