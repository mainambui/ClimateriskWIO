
#Put in a function to tidy up the global environment
qTranform <- function(df,vlst,scenario,time){
  nlist <- names(df)
  xx = c()
  for(k in seq_along(time)){
    xx[[k]] <- lapply(1:length(scenario), function(x){
      d1 <- cbind(df[,"ID"], sce = paste("ssp",scenario[[x]],sep=""), df[,nlist[grep(paste(scenario[[x]],time[[k]],sep="_"),nlist)]])
      colnames(d1) <- c("ID","Scenario",paste(vlst))
      return(d1)}) %>% bind_rows()  }

  names(xx) <- paste(time)
  xx <- bind_rows(xx, .id = "Period")

  #Quantile transform
  inormal <- function(x){
    qrank <- ((rank(x, na.last = TRUE) - 0.5) / sum(!is.na(x)))
    z_score <- scales::rescale(sqrt(2)*pracma::erfinv(2*qrank-1))
    return(z_score)}

  N <- ncol(xx)
  dfNmd <- xx %>% as.data.frame() %>%
    mutate(pH_trend = -1*pH_trend, evspsbl = -1*evspsbl) %>%
    #Reserve some key variables before normalising.
    #reduced primary production vs. reduced in annual rainfall (drought), negative trend in Ph means ocean becoming acidic over time
    mutate(across(all_of(4:N), ~ scales::rescale(.x)))

  #Convert to wide format
  dfWide <- lapply(1:length(vlst), function(x){
    xv <- reshape2::dcast(dfNmd,ID~Scenario+Period, value.var= vlst[[x]])
    scePeriod <- rep(scenario, each = length(time))
    names(xv) <- c("ID",paste(rep(vlst[[x]], each = length(scePeriod)), scePeriod, time, sep = "_"))
    return(xv)})
  stdHzd <- dfWide %>% purrr::reduce(left_join, by = "ID")
  return(stdHzd)
  }
