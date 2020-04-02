
#Load the simulation results
# FSkq indicate the final size for eta k and scenario q

#eta 0.25
FS11<-FinSize
CP11<-PeakInc
EXT11<-Ext.prob

FS12<-FinSize
CP12<-PeakInc
EXT12<-Ext.prob

FS13<-FinSize
CP13<-PeakInc
EXT13<-Ext.prob

#eta 0.5
FS21<-FinSize
CP21<-PeakInc
EXT21<-Ext.prob

FS22<-FinSize
CP22<-PeakInc
EXT22<-Ext.prob

FS23<-FinSize
CP23<-PeakInc
EXT23<-Ext.prob

#eta 0.75
FS31<-FinSize
CP31<-PeakInc
EXT31<-Ext.prob

FS32<-FinSize
CP32<-PeakInc
EXT32<-Ext.prob

FS33<-FinSize
CP33<-PeakInc
EXT33<-Ext.prob

#save(FS11,FS12,FS13,FS21,FS22,FS23,FS31,FS32,FS33, CP11,CP12,CP13,CP21,CP22,CP23,CP31,CP32,CP33, EXT11,EXT12,EXT13,EXT21,EXT22,EXT23,EXT31,EXT32,EXT33, file = "PlotAsy31Ra25Rs25r025.RData")

library(latex2exp)
par(mar=c(3,5,1,3))
boxplot(FS11,FS12,FS13,FS21,FS22,FS23,FS31,FS32,FS33, main="", at=c(1,2,3,5,6,7,9,10,11), border = "brown", col = c("goldenrod2", "darkgreen","blue"), las=1, notch= FALSE, xaxt="n", xlim=c(0.5,14), ylim=c(50,1000), ylab="Final Size")
par(new=T)
plot(y=c(EXT11,EXT12,EXT13,EXT21,EXT22,EXT23,EXT31,EXT32,EXT33), x=c(1,2,3,5,6,7,9,10,11), pch=8, col="purple", lwd=1, xlim=c(0.5,14), ylim=c(0,1), xaxt="n", yaxt="n", xlab = " ", ylab = " " )
axis(4,)
axis(1, at=c(2,6,10), labels=c(TeX("$\\eta =0.25$"),TeX("$\\eta =0.5$"), TeX("$\\eta =0.75$")))
legend(12,0.9, legend = c("IAS","IBS","IBTBS" ), col=c("goldenrod2", "darkgreen", "blue"), pch=c(15,15,15), cex = 0.7, box.lty = 0)

CP11m<-mean(CP11)
CP11q1<-quantile(CP11, probs = 0.025)
CP11q2<-quantile(CP11, probs = 0.975)

CP12m<-mean(CP12)
CP12q1<-quantile(CP12, probs = 0.025)
CP12q2<-quantile(CP12, probs = 0.975)

CP13m<-mean(CP13)
CP13q1<-quantile(CP13, probs = 0.025)
CP13q2<-quantile(CP13, probs = 0.975)

CP21m<-mean(CP21)
CP21q1<-quantile(CP21, probs = 0.025)
CP21q2<-quantile(CP21, probs = 0.975)

CP22m<-mean(CP22)
CP22q1<-quantile(CP22, probs = 0.025)
CP22q2<-quantile(CP22, probs = 0.975)

CP23m<-mean(CP23)
CP23q1<-quantile(CP23, probs = 0.025)
CP23q2<-quantile(CP23, probs = 0.975)

CP31m<-mean(CP31)
CP31q1<-quantile(CP31, probs = 0.025)
CP31q2<-quantile(CP31, probs = 0.975)

CP32m<-mean(CP32)
CP32q1<-quantile(CP32, probs = 0.025)
CP32q2<-quantile(CP32, probs = 0.975)

CP33m<-mean(CP33)
CP33q1<-quantile(CP33, probs = 0.025)
CP33q2<-quantile(CP33, probs = 0.975)


par(mar=c(5,5,1,1))
plot(x=c(CP31m,CP32m,CP33m,CP21m,CP22m,CP23m,CP11m,CP12m,CP13m), y=c(1,1.25,1.5,2.75,3,3.25,4.5,4.75,5),col=c("goldenrod2","darkgreen","blue"), ylim = c(0.5,6), pch=16, yaxt="n", ylab = " ", xlab = "Peak Incidence ", xlim = c(0,450))
axis(2, at=c(4.75,3,1.25), labels=c(TeX("$\\eta =0.25$"),TeX("$\\eta =0.5$"), TeX("$\\eta =0.75$")), las=2)

arrows(x0=CP31m,y0=1, y1=1, x1=CP31q1, length = 0, col = "goldenrod2")
arrows(x0=CP31m,y0=1, y1=1, x1=CP31q2, length = 0, col="goldenrod2")
arrows(x0=CP32m,y0=1.25, y1=1.25, x1=CP32q1, length = 0, col="darkgreen")
arrows(x0=CP32m,y0=1.25, y1=1.25, x1=CP32q2, length = 0,col="darkgreen")
#arrows(x0=CP33m,y0=1.5, y1=1.5, x1=CP33q1, length = 0,col="red")
#arrows(x0=CP33m,y0=1.5, y1=1.5, x1=CP33q2, length = 0,col="red")
arrows(x0=CP33m,y0=1.5, y1=1.5, x1=CP33q1, length = 0,col="blue")
arrows(x0=CP33m,y0=1.5, y1=1.5, x1=CP33q2, length = 0,col="blue")

arrows(x0=CP21m,y0=2.75, y1=2.75, x1=CP21q1, length = 0, col = "goldenrod2")
arrows(x0=CP21m,y0=2.75, y1=2.75, x1=CP21q2, length = 0, col="goldenrod2")
arrows(x0=CP22m,y0=3, y1=3, x1=CP22q1, length = 0, col="darkgreen")
arrows(x0=CP22m,y0=3, y1=3, x1=CP22q2, length = 0,col="darkgreen")
#arrows(x0=CP23m,y0=3.5, y1=3.5, x1=CP23q1, length = 0,col="red")
#arrows(x0=CP23m,y0=3.5, y1=3.5, x1=CP23q2, length = 0,col="red")
arrows(x0=CP23m,y0=3.25, y1=3.25, x1=CP23q1, length = 0,col="blue")
arrows(x0=CP23m,y0=3.25, y1=3.25, x1=CP23q2, length = 0,col="blue")

arrows(x0=CP11m,y0=4.5, y1=4.5, x1=CP11q1, length = 0, col = "goldenrod2")
arrows(x0=CP11m,y0=4.5, y1=4.5, x1=CP11q2, length = 0, col="goldenrod2")
arrows(x0=CP12m,y0=4.75, y1=4.75, x1=CP12q1, length = 0, col="darkgreen")
arrows(x0=CP12m,y0=4.75, y1=4.75, x1=CP12q2, length = 0,col="darkgreen")
#arrows(x0=CP13m,y0=5.5, y1=5.5, x1=CP13q1, length = 0,col="red")
#arrows(x0=CP13m,y0=5.5, y1=5.5, x1=CP13q2, length = 0,col="red")
arrows(x0=CP13m,y0=5, y1=5, x1=CP13q1, length = 0,col="blue")
arrows(x0=CP13m,y0=5, y1=5, x1=CP13q2, length = 0,col="blue")

legend(380,5.5, legend = c("IAS","IBS","IBTBS" ), col=c("goldenrod2", "darkgreen", "blue" ), pch=c(16,16,16),lty=c(1,1,1), cex = 0.7, box.lty = 0)


