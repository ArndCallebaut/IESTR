

opt_each_site = apply(indices_pres2,FUN=max,2)

diag_cms = matrix(0,N_cycles,length(gsc[,1]))
for (i in 1:N_cycles){
  cm_loc = cm[[i]]
  diag_cms[i,] = diag(cm_loc)
}

maxi_diags = apply(diag_cms,FUN=max,2)

plot(maxi_diags,opt_each_site)

### ### ###
efficienty3 = maxi_diags / cost[gsc[,c(1,2)]+1]
eff3 = efficienty3 / max(efficienty3)
#eff = (log(efficienty)-min(log(efficienty)))/(max(log(efficienty))-min(log(efficienty)))

eff_tensor3 = array(0.5,c(nrow,ncol,3))

for (i in 1:length(efficienty3)){
  #rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=1 #R
  eff_tensor3[gsc[i,1]+1,gsc[i,2]+1,2]= eff3[i]  #G
  eff_tensor3[gsc[i,1]+1,gsc[i,2]+1,3]= 1-eff3[i]#B
}

eff_tensor3 [is.na(map)]=1
eff_tensor3 = round(eff_tensor3*255)
raster_RGB = stack(raster(eff_tensor3[,,1]),raster(eff_tensor3[,,2]),raster(eff_tensor3[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))
legend("topright",c("High effect ratio", "Low effect ratio"), cex=1.0, bty="y",
       fill=c("green","purple"),inset=.04)








library(ggplot2)
theme_set(theme_minimal())

a = data.frame(0,c(2,2))
a$X = c(0,1)
a$Colonization.probability = c(0,1)
sp <- ggplot(a, aes(X, Colonization.probability))+
  theme(legend.text = element_text(size=10))+
  geom_point(aes(color = Colonization.probability))
sp + scale_color_gradient(low = "darkgreen", high = "green3")


library(ggplot2)
theme_set(theme_minimal())

a = data.frame(0,c(2,2))
a$X = c(1,20)
a$Introduction.time = c(1,20)
sp <- ggplot(a, aes(X, Introduction.time))+
  theme(legend.text = element_text(size=10))+
  geom_point(aes(color = Introduction.time))
sp + scale_color_gradient(low = "yellow", high = "brown",breaks=c(1,10,20))









colfunc <- colorRampPalette(c("green", "darkgreen"))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

layout(matrix(1:2,ncol=2), width = c(1,1),height = c(1,1))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', 
     main = 'Colonization probability', cex.main=2.5)
text(x=1.3, y = seq(0,1,l=5), labels = seq(0,1,l=5),cex = 2)
rasterImage(legend_image, 0, 0, 1,1)




colfunc <- colorRampPalette(c("yellow", "brown"))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', 
     main = 'Colonization probability', cex.main=2)
text(x=1.3, y = seq(0,1,l=5), labels = seq(0,1,l=5),cex = 2)
rasterImage(legend_image, 0, 0, 1,1)



plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1% - 10% choosen site", '11% - 50% choosen site', '50% - 100% choosen site',
                            'Surviving sites', 'Not surviving sites'), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c('grey60', 'grey40', 'black', 'blue', 'orange'))
mtext("Frequency of choices", at=0.2, cex=2)



plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend=c('Suitable conditions lost',
                           'Suitable conditions maintained', 
                           'New suitable conditions',
                           'Maintained range',
                           'Lost range'),
       pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c('grey40', 'grey50', 'gray90', 'blue', 'red'))



plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Suitable conditions lost',
                            'Suitable conditions maintained', 
                            'New suitable conditions',
                            'Initially present sites'),
       pch=16, 
       pt.cex=3, 
       cex=1.5, 
       bty='n',
       col = c('grey40', 'grey50', 'gray60', 'green3'))

mtext("Species", at=0.2, cex=2)

nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)
data <- expand.grid(X=nr, Y=nc)
cccc = 600-cost3
data$IntroductionTime <- c(as.matrix( (0 + 30*(cccc-min(cccc,na.rm=TRUE)) / (max(cccc,na.rm=TRUE)-min(cccc,na.rm=TRUE)))))
ggplot() + 
  coord_fixed()+
  scale_fill_continuous(na.value = "white",low="yellow", high="brown")+
  geom_raster(data = data , aes(x = X, y=Y,fill = IntroductionTime )) + 
  theme(panel.background = element_blank(),
        strip.background = element_blank() ,
        line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


countmap = array(0,c(nrow,ncol))


nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)
data <- expand.grid(X=nr, Y=nc)
cccc = 600-cost3
data$IntroductionTime <- c(as.matrix( (0 + 30*(cccc-min(cccc,na.rm=TRUE)) / (max(cccc,na.rm=TRUE)-min(cccc,na.rm=TRUE)))))
ggplot() + 
  coord_fixed()+
  scale_fill_continuous(na.value = "white",low="yellow", high="brown")+
  geom_raster(data = data , aes(x = X, y=Y,fill = IntroductionTime )) + 
  theme(panel.background = element_blank(),
        strip.background = element_blank() ,
        line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

