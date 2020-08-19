######################################################################################################################################################
########## The CBP KIX domain regulates long-term memory and circadian activity: RNA-Seq Analysis
########## Ethan Bahl
########## May 28, 2020

######################################################################################################################################################
#################################################################################################### SETUP.
library(EDASeq) # EDASeq_2.20.0
library(edgeR) # edgeR_3.28.0
library(RUVSeq) # RUVSeq_1.20.0
library(RRHO)
library(RColorBrewer) # RColorBrewer_1.1-2
library(scatterplot3d) # scatterplot3d_0.3-41
# library(org.Mm.eg.db) # org.Mm.eg.db_3.10.0
# library(biomaRt) # biomaRt_2.42.0
library(ComplexHeatmap) # ComplexHeatmap_2.2.0
# library(gridBase) # gridBase_0.4-7
library(circlize) # circlize_0.4.8
library(svglite) # svglite_1.2.3

######################################################################################################################################################
#################################################################################################### LEARNING ANALYSIS: EDASeq normalization.
########## set up objects and parameters for plotting.
set = fig2.data$set
dataWithin = fig2.data$dataWithin
dataNorm = fig2.data$dataNorm
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
legend.cols = unique(pData(set)$group.color)
# layout matrix.
m = matrix(
    c(
        1,2,3,
        4,5,6,
        7,8,9,
        10,11,12
    ),
    byrow=TRUE,
    nrow=4
)

########## save plot as SVG.
# svglite("figure_S5.svg", height=16, width=18, system_fonts = list(sans = "Arial"))

layout(mat = m, heights = c(0.04, rep(0.24, 3)))
par(mar=c(0,0,0,0))
# column 1 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("raw"), cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 2 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC normalization"), cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 3 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC / depth normalization"), cex = 3, font=2, col = "black")  
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# update parameters for the main plots in the figure.
par(mar=c(5,5,4,2))

########## GC content bias plots.
# graphical parameters.
ylab = "log(gene counts)"
xlab = "GC content"
ylim=c(0,8.5)
xlim=c(min(fData(set)$gc),max(fData(set)$gc))
c1 = alpha("#333333", 0.8)
c2 = alpha("#333333", 0.2)
# GC content bias plots (raw).
biasPlot(set, "gc", log=T, col=pData(set)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(set)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(set)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)
# GC content bias plots (GC normalization).
biasPlot(dataWithin, "gc", log=T, col=pData(dataWithin)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(dataWithin)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(dataWithin)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)
# GC content bias plots (GC, depth normalization).
biasPlot(dataNorm, "gc", log=T, col=pData(dataNorm)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(dataNorm)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(dataNorm)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

########## RLE plots.
# graphical parameters.
ylim=c(-2, 2)
# RLE plots (raw).
plotRLE(set, col=pData(set)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plots (GC normalization).
plotRLE(dataWithin, col=pData(dataWithin)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plots (GC, depth normalization).
plotRLE(dataNorm, col=pData(dataNorm)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, cex=1.1, bty="n")

########## PCA plots.
# graphical parameters.
xlim=c(-0.6,0.6)
ylim=c(-0.6, 0.6)
# PCA plots (raw).
plotPCA(set, col=pData(set)$group.color, xlim=xlim, ylim=xlim, labels=F, pch=16, main="PCA plot", cex=1.3, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, bty="n", cex=1.1)
# PCA plots (GC normalization).
plotPCA(dataWithin, col=pData(dataWithin)$group.color, xlim=xlim, ylim=ylim, labels=F, pch=16, main="PCA plot", cex=1.3, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, bty="n", cex=1.1)
# PCA plots (GC, depth normalization).
plotPCA(dataNorm, col=pData(dataNorm)$group.color, xlim=xlim, ylim=ylim, labels=F, pch=16, main="PCA plot", cex=1.3, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, bty="n", cex=1.1)
# save plot.
# dev.off()

######################################################################################################################################################
#################################################################################################### WILDTYPE ANALYSIS: EDASeq normalization.
########## set up objects and parameters for plotting.
set = fig3.data$set
dataWithin = fig3.data$dataWithin
dataNorm = fig3.data$dataNorm
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
legend.cols = unique(pData(set)$group.color)
# layout matrix.
m = matrix(
    c(
        1,2,3,
        4,5,6,
        7,8,9,
        10,11,12
    ),
    byrow=TRUE,
    nrow=4
)

########## save plot as SVG.
#svglite("figure_S4.svg", height=16, width=18, system_fonts = list(sans = "Arial"))
layout(mat = m, heights = c(0.04, rep(0.24, 3)))
par(mar=c(0,0,0,0))
# column 1 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("raw"), cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 2 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC normalization"), cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 3 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC / depth normalization"), cex = 3, font=2, col = "black")  
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# update parameters for the main plots in the figure.
par(mar=c(5,5,4,2))

########## GC content bias plots.
# graphical parameters.
ylab = "log(gene counts)"
xlab = "GC content"
ylim=c(0,8.5)
xlim=c(min(fData(set)$gc),max(fData(set)$gc))
c1 = alpha("#333333", 0.8)
c2 = alpha("#333333", 0.2)
# GC content bias plots (raw).
biasPlot(set, "gc", log=T, col=pData(set)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(set)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(set)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)
# GC content bias plots (GC normalization).
biasPlot(dataWithin, "gc", log=T, col=pData(dataWithin)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(dataWithin)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(dataWithin)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)
# GC content bias plots (GC, depth normalization).
biasPlot(dataNorm, "gc", log=T, col=pData(dataNorm)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(dataNorm)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(dataNorm)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

########## RLE plots.
# graphical parameters.
ylim=c(-2, 2)
# RLE plots (raw).
plotRLE(set, col=pData(set)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plots (GC normalization).
plotRLE(dataWithin, col=pData(dataWithin)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plots (GC, depth normalization).
plotRLE(dataNorm, col=pData(dataNorm)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, cex=1.1, bty="n")

########## PCA plots.
# graphical parameters.
xlim=c(-0.6,0.6)
ylim=c(-0.6, 0.6)
# PCA plots (raw).
plotPCA(set, col=pData(set)$group.color, xlim=xlim, ylim=xlim, labels=F, pch=16, main="PCA plot", cex=1.3, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, bty="n", cex=1.1)
# PCA plots (GC normalization).
plotPCA(dataWithin, col=pData(dataWithin)$group.color, xlim=xlim, ylim=ylim, labels=F, pch=16, main="PCA plot", cex=1.3, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, bty="n", cex=1.1)
# PCA plots (GC, depth normalization).
plotPCA(dataNorm, col=pData(dataNorm)$group.color, xlim=xlim, ylim=ylim, labels=F, pch=16, main="PCA plot", cex=1.3, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, bty="n", cex=1.1)
# save plot.
#dev.off()

######################################################################################################################################################
#################################################################################################### HOMECAGE ANALYSIS: EDASeq normalization.
########## set up objects and parameters for plotting.
set = homecage.data$set
dataWithin = homecage.data$dataWithin
dataNorm = homecage.data$dataNorm
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
legend.cols = rev(unique(pData(set)$group.color))
# layout matrix.
m = matrix(
    c(
        1,2,3,
        4,5,6,
        7,8,9,
        10,11,12
    ),
    byrow=TRUE,
    nrow=4
)

########## save plot as SVG.
# svglite("figures/homecage_sfig1.svg", height=16, width=18, system_fonts = list(sans = "Arial"))
layout(mat = m, heights = c(0.04, rep(0.24, 3)))
par(mar=c(0,0,0,0))
# column 1 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("raw"), cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 2 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC normalization"), cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 3 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC / depth normalization"), cex = 3, font=2, col = "black")  
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")
# update parameters for the main plots in the figure.
par(mar=c(5,5,4,2))

########## GC content bias plots.
# graphical parameters.
ylab = "log(gene counts)"
xlab = "GC content"
ylim=c(0,8.5)
xlim=c(min(fData(set)$gc),max(fData(set)$gc))
c1 = alpha("#333333", 0.8)
c2 = alpha("#333333", 0.2)
# GC content bias plots (raw).
biasPlot(set, "gc", log=T, col=pData(set)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(set)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(set)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)
# GC content bias plots (GC normalization).
biasPlot(dataWithin, "gc", log=T, col=pData(dataWithin)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(dataWithin)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(dataWithin)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)
# GC content bias plots (GC, depth normalization).
biasPlot(dataNorm, "gc", log=T, col=pData(dataNorm)$group.color, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(dataNorm)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c1, whiskcol=c1, staplecol=c1, boxcol=c1, outcol=c1, pch=23, cex=0, xaxt="n")
text(x=median(fData(dataNorm)$gc), y=0.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

########## RLE plots.
# graphical parameters.
ylim=c(-2, 2)
# RLE plots (raw).
plotRLE(set, col=pData(set)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plots (GC normalization).
plotRLE(dataWithin, col=pData(dataWithin)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plots (GC, depth normalization).
plotRLE(dataNorm, col=pData(dataNorm)$group.color, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, cex=1.1, bty="n")

########## PCA plots.
# graphical parameters.
xlim=c(-0.8,0.8)
ylim=c(-0.8, 0.8)
# PCA plots (raw).
plotPCA(set, col=pData(set)$group.color, xlim=xlim, ylim=xlim, labels=F, pch=16, main="PCA plot", cex=1.5, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, bty="n", cex=1.1)
# PCA plots (GC normalization).
plotPCA(dataWithin, col=pData(dataWithin)$group.color, xlim=xlim, ylim=ylim, labels=F, pch=16, main="PCA plot", cex=1.5, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataWithin)$group), fill=legend.cols, bty="n", cex=1.1)
# PCA plots (GC, depth normalization).
plotPCA(dataNorm, col=pData(dataNorm)$group.color, xlim=xlim, ylim=ylim, labels=F, pch=16, main="PCA plot", cex=1.5, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(pData(dataNorm)$group), fill=legend.cols, bty="n", cex=1.1)
# save plot.
#dev.off()

