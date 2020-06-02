######################################################################################################################################################
########## The CBP KIX domain regulates long-term memory and circadian activity: RNA-Seq Analysis
########## Ethan Bahl
########## May 28, 2020

######################################################################################################################################################
#################################################################################################### SECTION 01: SETUP.


library(EDASeq) # EDASeq_2.20.0
library(edgeR) # edgeR_3.28.0
library(RUVSeq) # RUVSeq_1.20.0
library(scatterplot3d) # scatterplot3d_0.3-41
# library(org.Mm.eg.db) # org.Mm.eg.db_3.10.0
# library(biomaRt) # biomaRt_2.42.0
library(RColorBrewer) # RColorBrewer_1.1-2
library(ComplexHeatmap) # ComplexHeatmap_2.2.0
# library(gridBase) # gridBase_0.4-7
library(circlize) # circlize_0.4.8
library(svglite) # svglite_1.2.3

########## load data.
# read count data.
count.data = as.matrix(read.csv("data/counts_all.csv", row.names=1, header=TRUE))
# read sample information, set up factors/levels.
sample.info = read.csv("data/sample_info.csv", row.names=1, header=TRUE, stringsAsFactors=FALSE)
sample.info$experiment = factor(sample.info$experiment, levels=c("experiment1", "experiment2"))
sample.info$genotype = factor(sample.info$genotype, levels=c("WT", "CBP KIX/KIX"))
sample.info$exposure = factor(sample.info$exposure, levels=c("homecage", "learning"))
sample.info$group = factor(sample.info$group, levels=c("EXP1.WT.HOME", "EXP1.WT.MWM", "EXP1.KIX.MWM", "EXP2.WT.HOME", "EXP2.WT.CFC", "EXP2.KIX.CFC"))
# read gene information.
gene.info = read.csv("data/gene_info.csv", row.names=1, header=TRUE, stringsAsFactors=FALSE)

########## create seqExpressionSet.
set.all = newSeqExpressionSet(counts=count.data, phenoData=sample.info, featureData=gene.info)

######################################################################################################################################################
#################################################################################################### SECTION 02: Figure 2 (differential expression).

########## select samples from learning groups, create DGEList for edgeR.
set = set.all[, pData(set.all)$exposure == "learning"]
pData(set) = droplevels(pData(set))
y = DGEList(counts(set), samples=pData(set), genes=fData(set))

########## filter undetected / lowly expressed genes.
design_filter = model.matrix(~ experiment + genotype, data=pData(set))
keep = filterByExpr(y, design=design_filter)
y = y[keep, , keep.lib.sizes=FALSE]
set = set[keep,]

########## EDASeq normalization.
# normalize for GC content differences.
dataWithin <- withinLaneNormalization(set, "gc", which="full", offset=TRUE)
# normalize for sequencing depth.
dataNorm <- betweenLaneNormalization(dataWithin, which="upper", offset=TRUE)

########## create lists to hold data and results.
fig2.data = list(set=set, dataWithin=dataWithin, dataNorm=dataNorm)
fig2.results = list()

########## first pass of edgeR without RUV normalization.
# set up objects, attach EDASeq offsets to DGEList, create design matrix.
set = fig2.data$dataNorm
y = DGEList(counts=counts(set), samples=pData(set), genes=fData(set))
y$offset = -offst(set)
design = model.matrix(~experiment + genotype, data=pData(set))
# edgeR glmQLF (robust) workflow.
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
test = glmQLFTest(fit)
top = topTags(test, Inf)$table
# store results from first pass.
fig2.results[[1]] = list(set=set, y=y, fit=fit, test=test, top=top)

########## second pass of edgeR with RUV normalization.
# RUVr normalization.
set.ruv = RUVr(set, rownames(set), k=3, residuals(fit, type="deviance"))
# reattach gene information.
fData(set.ruv) = fData(set)[rownames(set.ruv),]
# build new design matrix with latent factors of unwanted variation included.
design = model.matrix(~ experiment + genotype + W_1 + W_2 + W_3, data=pData(set.ruv))
# second pass of edgeR.
y = DGEList(counts=counts(set), samples=pData(set), genes=fData(set))
y$offset = -offst(set)
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
test = glmQLFTest(fit, coef=3)
top = topTags(test, Inf)$table
# store results from second pass.
fig2.results[[2]] = list(set=set.ruv, y=y, fit=fit, test=test, top=top)
names(fig2.results) = c("W0", "W3")

######################################################################################################################################################
#################################################################################################### SECTION 03: Figure 3 (differential expression).

########## select wildtype samples, create DGEList for edgeR.
set = set.all[, pData(set.all)$genotype == "WT"]
pData(set) = droplevels(pData(set))
y = DGEList(counts(set), samples=pData(set), genes=fData(set))

########## filter undetected / lowly expressed genes.
design_filter = model.matrix(~ experiment + exposure, data=pData(set))
keep = filterByExpr(y, design=design_filter)
y = y[keep, , keep.lib.sizes=FALSE]
set = set[keep,]

########## EDASeq normalization.
# normalize for GC content differences.
dataWithin <- withinLaneNormalization(set, "gc", which="full", offset=TRUE)
# normalize for sequencing depth.
dataNorm <- betweenLaneNormalization(dataWithin, which="upper", offset=TRUE)

########## create lists to hold data and results.
fig3.data = list(set=set, dataWithin=dataWithin, dataNorm=dataNorm)
fig3.results = list()

########## first pass of edgeR without RUV normalization.
# set up objects, attach EDASeq offsets to DGEList, create design matrix.
set = fig3.data$dataNorm
y = DGEList(counts=counts(set), samples=pData(set), genes=fData(set))
y$offset = -offst(set)
design = model.matrix(~experiment + exposure, data=pData(set))
# edgeR glmQLF (robust) workflow.
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
test = glmQLFTest(fit)
top = topTags(test, Inf)$table
# store results from first pass.
fig3.results[[1]] = list(set=set, y=y, fit=fit, test=test, top=top)

########## second pass of edgeR with RUV normalization.
# RUVr normalization.
set.ruv = RUVr(set, rownames(set), k=4, residuals(fit, type="deviance"))
# reattach gene information.
fData(set.ruv) = fData(set)[rownames(set.ruv),]
# build new design matrix with latent factors of unwanted variation included.
design = model.matrix(~ experiment + exposure + W_1 + W_2 + W_3 + W_4, data=pData(set.ruv))
# second pass of edgeR.
y = DGEList(counts=counts(set), samples=pData(set), genes=fData(set))
y$offset = -offst(set)
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
test = glmQLFTest(fit, coef=3)
top = topTags(test, Inf)$table
# store results from second pass.
fig3.results[[2]] = list(set=set.ruv, y=y, fit=fit, test=test, top=top)
names(fig3.results) = c("W0", "W4")

######################################################################################################################################################
#################################################################################################### SECTION 04: Enrichment Analysis.

########## enrichment analysis.
# testing DEGs (DE between WT-learning and WT-homecage) for enrichment of positive controls.
fisher.test(table(fig3.results$W4$top$FDR<=0.05, fig3.results$W4$top$control=="positive"))
# testing genes regulated by KIX after learning for enrichment of spatial learning regulated genes.
fisher.test(table(fig2.results$W3$top$FDR<=0.05, fig2.results$W3$top$ensembl_gene_id %in% subset(fig3.results$W4$top, FDR<=0.05)$ensembl_gene_id))

######################################################################################################################################################
#################################################################################################### SECTION 05: Figure S5 (EDASeq normalization, corresponding to Figure 2 in main text).

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
svglite("figure_S5.svg", height=16, width=18, system_fonts = list(sans = "Arial"))
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
dev.off()

######################################################################################################################################################
#################################################################################################### SECTION 06: Figure S4 (EDASeq, corresponding to Figure 3 in main text).

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
svglite("figure_S4.svg", height=16, width=18, system_fonts = list(sans = "Arial"))
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
dev.off()

######################################################################################################################################################
#################################################################################################### SECTION 07: Figure S7 (RUVSeq, corresponding to Figure 2 in main text).

########## prepare data for plotting.
# get differential expression results, with and without RUV normalization.
tt = fig2.results$W0$top
tt.ruv = fig2.results$W3$top
set = fig2.results$W0$set
set.ruv = fig2.results$W3$set
# get PCA plot data using adapted source code from EDASeq's plotPCA().
object = normCounts(set)
Y <- apply(log(object+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
S <- svd(Y)
object = normCounts(set.ruv)
Y <- apply(log(object+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
S.ruv <- svd(Y)

########## plot set up.
# graphical parameters.
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
legend.cols = unique(pData(set)$group.color)
# layout matrix.
m = matrix(
    c(
        1,2,
        3,4,
        5,6,
        7,8,
        9,10
    ),
    byrow=TRUE,
    nrow=5
)

########## save plot as SVG.
svglite("figure_S7.svg", height=20, width=12, system_fonts = list(sans = "Arial"))
layout(mat = m, heights = c(0.04, rep(0.24, 4)))
par(mar=c(0,0,0,0))
# column 1 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("GC / depth normalization"), cex = 3, font=2, col = "black")
segments(x0=0.05, x1=0.95, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 2 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("GC / depth / RUV normalization"), cex = 3, font=2, col = "black")  
segments(x0=0.05, x1=0.95, y0=0.2, y1=0.2, lwd=2.5, col="black")
# update parameters for the main plots in the figure.
par(mar=c(5,5,3,2))

########## PValue distribution plots.
# PValue distribution of differential expression results (without RUV normalization).
hist(tt$PValue, xlab="PValue", ylab="number of genes", main="distribution of P-values", cex.main=cex.main, ylim=c(0,3000), cex.lab=cex.lab, cex.axis=cex.axis, col="#d3d3d3")
# PValue distribution of differential expression results (with RUV normalization).
hist(tt.ruv$PValue, xlab="PValue", ylab="number of genes", main="distribution of P-values", cex.main=cex.main, ylim=c(0,3000), cex.lab=cex.lab, cex.axis=cex.axis, col="#d3d3d3")

########## mean-difference plots.
# mean-difference plot (without RUV normalization).
sig = tt[,"FDR"] <= 0.05
plot(tt[-which(sig),"logCPM"], tt[-which(sig),"logFC"],
    pch=16,cex=0.3,
    ylab="logFC (relative to control)",
    xlab="average log(CPM)",
    xlim=c(min(tt$logCPM), max(tt$logCPM)),
    ylim=c(min(tt$logFC), max(tt$logFC)),
    main="mean-difference plot",
    col=scales::alpha("#333333", 0.6),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)
points(tt[sig, "logCPM"],tt[sig, "logFC"],
    pch=16, cex=1,
    col=scales::alpha(ifelse(tt[sig,"logFC"] > 0, "red", "blue"), 0.8)
)
abline(h=0, lty=2, col="#AFAFAF")
legend("topright", legend=c("upregulated", "downregulated", "not significant"), fill=c("red", "blue", "#333333"), bty="n", cex=1.2)
# mean-difference plot (with RUV normalization).
sig = tt.ruv[,"FDR"] <= 0.05
plot(tt.ruv[-which(sig),"logCPM"], tt.ruv[-which(sig),"logFC"],
    pch=16,cex=0.3,
    ylab="logFC (relative to control)",
    xlab="average log(CPM)",
    xlim=c(min(tt.ruv$logCPM), max(tt.ruv$logCPM)),
    ylim=c(min(tt.ruv$logFC), max(tt.ruv$logFC)),
    main="mean-difference plot",
    col=scales::alpha("#333333", 0.6),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)
points(tt.ruv[sig, "logCPM"],tt.ruv[sig, "logFC"],
    pch=16, cex=1,
    col=scales::alpha(ifelse(tt.ruv[sig,"logFC"] > 0, "red", "blue"), 0.8)
)
abline(h=0, lty=2, col="#AFAFAF")
legend("topright", legend=c("upregulated", "downregulated", "not significant"), fill=c("red", "blue", "#333333"), bty="n", cex=1.2)

########## RLE plots.
# RLE plot (without RUV normalization).
plotRLE(set, col=pData(set)$group.color, ylim=c(-2,2), main="RLE plot", cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plot (with RUV normalization).
plotRLE(set.ruv, col=pData(set.ruv)$group.color, ylim=c(-2,2), main="RLE plot", cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(set.ruv)$group), fill=legend.cols, cex=1.1, bty="n")

########## 3D PCA plots.
# 3D PCA plot (without RUV normalization).
s = S
# since the first component is correlated with experimental batch (cor > 0.99), and we account for experimental batch directly in the DE model, we remove the first component for RUV assessment and visualization purposes.
cor(as.numeric(pData(set)$experiment), s$u[,1])
s$u = s$u[,-1]
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=pData(set)$group.color, main="3D PCA plot", cex.symbols=2, angle=40, type="h", xlab="PC2", ylab="", zlab="PC4", xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1), cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.6, y = -1.6, labels="PC3", srt = 35, cex=cex.lab-0.2)
# 3D PCA plot (with RUV normalization).
s = S.ruv
# since the first component is correlated with experimental batch (cor > 0.99), and we account for experimental batch directly in the DE model, we remove the first component for RUV assessment and visualization purposes.
cor(as.numeric(pData(set.ruv)$experiment), s$u[,1])
s$u = s$u[,-1]
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=pData(set)$group.color, main="3D PCA plot", cex.symbols=2, angle=40, type="h", xlab="PC2", ylab="", zlab="PC4", xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1), cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.6, y = -1.6, labels="PC3", srt = 35, cex=cex.lab-0.2)
# save plot.
dev.off()

######################################################################################################################################################
#################################################################################################### SECTION 08: Figure S6 (RUVSeq, corresponding to Figure 3 in main text).

########## prepare data for plotting.
# get differential expression results, with and without RUV normalization.
tt = fig3.results$W0$top
tt.ruv = fig3.results$W4$top
set = fig3.results$W0$set
set.ruv = fig3.results$W4$set
# get PCA plot data using adapted source code from EDASeq's plotPCA().
object = normCounts(set)
Y <- apply(log(object+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
S <- svd(Y)
object = normCounts(set.ruv)
Y <- apply(log(object+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
S.ruv <- svd(Y)

########## plot set up.
# graphical parameters.
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
legend.cols = unique(pData(set)$group.color)
# layout matrix.
m = matrix(
    c(
        1,2,
        3,4,
        5,6,
        7,8,
        9,10
    ),
    byrow=TRUE,
    nrow=5
)

########## save plot as SVG.
svglite("figure_S6.svg", height=20, width=12, system_fonts = list(sans = "Arial"))
layout(mat = m, heights = c(0.04, rep(0.24, 4)))
par(mar=c(0,0,0,0))
# column 1 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("GC / depth normalization"), cex = 3, font=2, col = "black")
segments(x0=0.05, x1=0.95, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 2 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("GC / depth / RUV normalization"), cex = 3, font=2, col = "black")  
segments(x0=0.05, x1=0.95, y0=0.2, y1=0.2, lwd=2.5, col="black")
# update parameters for the main plots in the figure.
par(mar=c(5,5,3,2))

########## PValue distribution plots.
# PValue distribution of differential expression results (without RUV normalization).
hist(tt$PValue, xlab="PValue", ylab="number of genes", main="distribution of P-values", cex.main=cex.main, ylim=c(0,3000), cex.lab=cex.lab, cex.axis=cex.axis, col="#d3d3d3")
# PValue distribution of differential expression results (with RUV normalization).
hist(tt.ruv$PValue, xlab="PValue", ylab="number of genes", main="distribution of P-values", cex.main=cex.main, ylim=c(0,3000), cex.lab=cex.lab, cex.axis=cex.axis, col="#d3d3d3")

########## mean-difference plots.
# mean-difference plot (without RUV normalization).
sig = tt[,"FDR"] <= 0.05
plot(tt[-which(sig),"logCPM"], tt[-which(sig),"logFC"],
    pch=16,cex=0.3,
    ylab="logFC (relative to control)",
    xlab="average log(CPM)",
    xlim=c(min(tt$logCPM), max(tt$logCPM)),
    ylim=c(min(tt$logFC), max(tt$logFC)),
    main="mean-difference plot",
    col=scales::alpha("#333333", 0.6),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)
points(tt[sig, "logCPM"],tt[sig, "logFC"],
    pch=16, cex=1,
    col=scales::alpha(ifelse(tt[sig,"logFC"] > 0, "red", "blue"), 0.8)
)
abline(h=0, lty=2, col="#AFAFAF")
legend("topright", legend=c("upregulated", "downregulated", "not significant"), fill=c("red", "blue", "#333333"), bty="n", cex=1.2)
# mean-difference plot (with RUV normalization).
sig = tt.ruv[,"FDR"] <= 0.05
plot(tt.ruv[-which(sig),"logCPM"], tt.ruv[-which(sig),"logFC"],
    pch=16,cex=0.3,
    ylab="logFC (relative to control)",
    xlab="average log(CPM)",
    xlim=c(min(tt.ruv$logCPM), max(tt.ruv$logCPM)),
    ylim=c(min(tt.ruv$logFC), max(tt.ruv$logFC)),
    main="mean-difference plot",
    col=scales::alpha("#333333", 0.6),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)
points(tt.ruv[sig, "logCPM"],tt.ruv[sig, "logFC"],
    pch=16, cex=1,
    col=scales::alpha(ifelse(tt.ruv[sig,"logFC"] > 0, "red", "blue"), 0.8)
)
abline(h=0, lty=2, col="#AFAFAF")
legend("topright", legend=c("upregulated", "downregulated", "not significant"), fill=c("red", "blue", "#333333"), bty="n", cex=1.2)

########## RLE plots.
# RLE plot (without RUV normalization).
plotRLE(set, col=pData(set)$group.color, ylim=c(-2,2), main="RLE plot", cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(set)$group), fill=legend.cols, cex=1.1, bty="n")
# RLE plot (with RUV normalization).
plotRLE(set.ruv, col=pData(set.ruv)$group.color, ylim=c(-2,2), main="RLE plot", cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(pData(set.ruv)$group), fill=legend.cols, cex=1.1, bty="n")

########## 3D PCA plots.
# 3D PCA plot (without RUV normalization).
s = S
# since the first component is correlated with experimental batch (cor > 0.99), and we account for experimental batch directly in the DE model, we remove the first component for RUV assessment and visualization purposes.
cor(as.numeric(pData(set)$experiment), s$u[,1])
s$u = s$u[,-1]
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=pData(set)$group.color, main="3D PCA plot", cex.symbols=2, angle=40, type="h", xlab="PC2", ylab="", zlab="PC4", xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1), cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.6, y = -1.6, labels="PC3", srt = 35, cex=cex.lab-0.2)
# 3D PCA plot (with RUV normalization).
s = S.ruv
# since the first component is correlated with experimental batch (cor > 0.99), and we account for experimental batch directly in the DE model, we remove the first component for RUV assessment and visualization purposes.
cor(as.numeric(pData(set.ruv)$experiment), s$u[,1])
s$u = s$u[,-1]
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=pData(set)$group.color, main="3D PCA plot", cex.symbols=2, angle=40, type="h", xlab="PC2", ylab="", zlab="PC4", xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1), cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.6, y = -1.6, labels="PC3", srt = 35, cex=cex.lab-0.2)
# save plot.
dev.off()

######################################################################################################################################################
#################################################################################################### SECTION 09: Figure 4 (analysis of Figure 2,3 overlap).

########## prep data.
# get differential expression results and merge them together.
res2 = fig2.results$W3$top
res3 = fig3.results$W4$top
int.genes = intersect(rownames(res2), rownames(res3))
res = cbind(res2[int.genes,], res3[int.genes, c("logFC", "logCPM", "F", "PValue", "FDR")])
colnames(res)[5:9] = paste0("fig2.", colnames(res)[5:9])
colnames(res)[10:14] = paste0("fig3.", colnames(res)[10:14])
# compute maximum FDR value for the two analyses.
res$max.fdr = pmax(res$fig2.FDR, res$fig3.FDR)
# find genes that were both regulated by learning, and dysregulated after learning by the CBP KIX/KIX genotype.
int.genes = rownames(subset(res, max.fdr <= 0.05))
## set up row order.
tmp = subset(res, max.fdr<= 0.05)
tmp = tmp[order(tmp$max.fdr),]

########## (for figure 2 data) get RUV-normalized counts, find overlapping genes, scale for visualization.
# get RUV-normalized counts.
de.heat2 = normCounts(fig2.results$W3$set)
# compute counts-per-million.
de.heat2 = cpm(de.heat2, prior.count=2, log=TRUE)
# select genes significant from fig2 and fig3 analysis.
de.heat2 = de.heat2[int.genes,]
de.heat2 = de.heat2[rownames(tmp),]
# change rownames to symbol.
rownames(de.heat2) = res[rownames(de.heat2), "gene_name"]
# because we account for experimental batch directly in the DE model, here we split data by experimental batch for scaling.
de.heat2 = list(EXP1 = de.heat2[, c(1:6)], EXP2 = de.heat2[, c(7:ncol(de.heat2))])
# scaling.
de.heat2 = lapply(de.heat2, function(x) t(scale(t(x))))
# merge back together.
analysis2 = do.call('cbind', de.heat2)
# set up column title splits and colors.
column.title.color2 = unique(pData(fig2.data$set)$group.color)
col_split2 = data.frame(group=paste0(ifelse(pData(fig2.data$set)$genotype=="WT", "WT", "KIX"), " ", ifelse(pData(fig2.data$set)$paradigm == "homecage", "HOME", ifelse(pData(fig2.data$set)$paradigm=="CFC", "CFC", "MWM"))))
rownames(col_split2) = rownames(pData(fig2.data$set))
col_split2 = col_split2[colnames(analysis2),,drop=F]
col_split2$group = factor(as.character(col_split2$group), levels=c("WT MWM", "WT CFC", "KIX MWM", "KIX CFC"))
# set up row split for up v down genes.
row_split2 = data.frame(direction=ifelse(tmp$fig2.logFC>0, "down in KIX after learning\n(relative to WT after learning)", "up in KIX after learning\n(relative to WT after learning)" ))
row_split2$direction = relevel(row_split2$direction, ref=levels(row_split2$direction)[2])
# update names to show replicate number.
colnames(analysis2) = pData(fig2.data$set)[colnames(analysis2), "replicate"]


########## (for figure 3 data) get RUV-normalized counts, find overlapping genes, scale for visualization.
# get RUV-normalized counts.
de.heat3 = normCounts(fig3.results$W4$set)
# compute counts-per-million.
de.heat3 = cpm(de.heat3, prior.count=2, log=TRUE)
# select genes significant from fig2 and fig3 analysis.
de.heat3 = de.heat3[int.genes,]
de.heat3 = de.heat3[rownames(tmp),]
# change rownames to symbol.
rownames(de.heat3) = res[rownames(de.heat3), "gene_name"]
# because we account for experimental batch directly in the DE model, here we split data by experimental batch for scaling.
de.heat3 = list(EXP1 = de.heat3[, c(1:6)], EXP2 = de.heat3[, c(7:ncol(de.heat3))])
# scaling.
de.heat3 = lapply(de.heat3, function(x) t(scale(t(x))))
# merge back together.
analysis3 = do.call('cbind', de.heat3)
# set up column title splits and colors.
column.title.color3 = unique(pData(fig3.data$set)$group.color)
col_split3 = data.frame(group=paste0(ifelse(pData(fig3.data$set)$genotype=="WT", "WT", "KIX"), " ", ifelse(pData(fig3.data$set)$paradigm == "homecage", "HOME", ifelse(pData(fig3.data$set)$paradigm=="CFC", "CFC", "MWM")), c(rep(" (EXP1)", 3), rep("", 3), rep(" (EXP2)", 3), rep("", 4))))
rownames(col_split3) = rownames(pData(fig3.data$set))
col_split3 = col_split3[colnames(analysis3),,drop=F]
col_split3$group = factor(as.character(col_split3$group), levels=c("WT HOME (EXP1)", "WT HOME (EXP2)", "WT MWM", "WT CFC"))
# set up row split for up v down genes.
row_split3 = data.frame(direction=ifelse(tmp$fig3.logFC>0, "up in WT after learning\n(relative to WT homecage)", "down in WT after learning\n(relative to WT homecage)" ))
row_split3$direction = relevel(row_split3$direction, ref=levels(row_split3$direction)[2])
# update names to show replicate number.
colnames(analysis3) = pData(fig3.data$set)[colnames(analysis3), "replicate"]

########## get colors for heatmap cells.
colz = colorRamp2(seq(c(-1)* max(abs(c(analysis2, analysis3))), max(abs(c(analysis2, analysis3))), length = 101), colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(101))

####################
# create the left side of figure 4.
set.seed(777)
hleft = Heatmap(
    analysis3,
    col=colz, 
    cluster_columns=FALSE,
    column_split = col_split3,
    column_gap = unit(c(1,4, 1), "mm"),
    column_title_gp = gpar(col="black", fill = column.title.color3[c(1,3,2,4)], font = 2, fontsize=10),
    column_names_rot = 0,
    column_names_gp = gpar(col="#333333", fontsize=10),
    column_names_centered = TRUE,
    row_title_gp = gpar(col="#333333", font =2, fontsize=13),
    row_names_gp = gpar(col="#333333", fontsize=14, font=2, hjust="center", vjust="center"),
    row_split = row_split3,
    row_gap = unit(1, "mm"),
    row_names_side = "right",
    row_title_side = "left",
    row_names_centered = TRUE,
    cluster_rows=FALSE,
    show_row_names=TRUE,
    width= unit(6, "in"),
    height=unit(6, "in"),
    show_heatmap_legend=FALSE,
    border=TRUE,
    rect_gp = gpar(col= "#FFFFFF", lwd=2),
    show_parent_dend_line= TRUE,
    name="scaled logCPM\n(RUV normalized)"
)
# create the right side of figure 4.
set.seed(777)
hright = Heatmap(
    analysis2,
    col=colz, 
    cluster_columns=FALSE,
    column_split = col_split2,
    column_gap = unit(c(1,4, 1), "mm"),
    column_title_gp = gpar(col="black", fill = column.title.color2[c(1,3,2,4)], font = 2, fontsize=10),
    column_names_rot = 0,
    column_names_gp = gpar(col="#333333", fontsize=10),
    column_names_centered = TRUE,
    row_title_gp = gpar(col="#333333", font =2, fontsize=13),
    row_names_gp = gpar(col="#333333", fontsize=14, font=2, hjust="center", vjust="center"),
    row_split = row_split2,
    row_gap = unit(1, "mm"),
    row_names_side = "right",
    row_dend_side = "right",
    row_title_side = "right",
    cluster_rows=FALSE,
    show_row_names=FALSE,
    width= unit(6, "in"),
    height=unit(6, "in"),
    show_heatmap_legend=TRUE,
    border=TRUE,
    rect_gp = gpar(col= "#FFFFFF", lwd=2),
    name="scaled logCPM\n(RUV normalized)"
)

# save as SVG using svglite().
# left.
svglite("figures/figure_4_left.svg", system_fonts = list(sans = "Arial"))
hleft
dev.off()
# right.
svglite("figures/figure_4_right.svg", system_fonts = list(sans = "Arial"))
hright
dev.off()