
library(raster)
rr <- stack(list.dirs(xx))
time <- 1:nlayers(rr)

#Linear trend
fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[2] }}
rr.slope=calc(rr, fun)
plot(rr.slope)

#estimate the pVaues
fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[8] }}
p <- calc(rr, fun=fun)
plot(p, main="p-Value")

m = c(0, 0.05, 1, 0.05, 1, 0)
rclmat = matrix(m, ncol=3, byrow=TRUE)
p.mask = reclassify(p, rclmat)
fun=function(x) { x[x<1] <- NA; return(x)}
p.mask.NA = calc(p.mask, fun)


#Mask all insignificant values in the trend map, so we only get change significant at the 95% level:
trend.sig = mask(gimms.slope, p.mask.NA)
plot(trend.sig, main="significant change")