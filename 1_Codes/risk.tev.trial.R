# Constructing a contour of circles
xColor <- seq(0,1,length.out=1000) # scale of color for x and y, 
yColor <- seq(0,1,length.out=1000)
x <- seq(0,7.5,length.out=1000) # on respective x and y axis
y <- seq(0,1,length.out=1000)
df <- cbind(expand.grid(x=xColor, y=yColor), expand.grid(x=x, y=y)) #grid for colors
colnames(df)<-c("xColor","yColor","x","y")
df$zColor <- (df$xColor^2+df$yColor^2) # the color factor for radius
df$zColor <- ntile(df$zColor, 10)

(plt1 <- ggplot()+
  geom_raster(data = df, aes(x = x, y = y, fill = zColor))+
  scale_fill_distiller(type = "div", palette = "Spectral", direction=-1)+
  geom_point(data = dfs, aes(x = TEV/1e6, y = rr.ssp370, shape = ISO3), size = 1.5) +
  labs(y = "Potential residual risk", x = "Total Economic Value (Million US$/year)")+
  #scale_colour_manual(name="", values = c("KEN"="darkred","MDG"="yellow","MOZ"="dodgerblue4","TZA"="cyan"))
  theme_bw(base_size = 10)+
  scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
  #scale_shape_manual(values = c("SSP2-4.5" = 19, "SSP3-7.0" = 17))+
  guides(shape="none", colour = "none")+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size = 8),
        legend.key.height = unit(.1, 'cm'),
        legend.key.width = unit(.2, 'cm'), 
        panel.border = element_blank()))
