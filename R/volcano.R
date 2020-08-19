######################################################################################################################################################
########## The CBP KIX domain regulates long-term memory and circadian activity: RNA-Seq Analysis
########## Ethan Bahl
########## May 28, 2020

########## VOLCANO PLOTS.

######################################################################################################################################################
#################################################################################################### SETUP.
library(EDASeq) # EDASeq_2.20.0
library(edgeR) # edgeR_3.28.0
library(RUVSeq) # RUVSeq_1.20.0
library(RRHO) # RRHO_1.26.0
library(RColorBrewer) # RColorBrewer_1.1-2
library(scatterplot3d) # scatterplot3d_0.3-41
library(org.Mm.eg.db) # org.Mm.eg.db_3.10.0
library(biomaRt) # biomaRt_2.42.0
library(ComplexHeatmap) # ComplexHeatmap_2.2.0
library(circlize) # circlize_0.4.8
library(svglite) # svglite_1.2.3
library(clusterProfiler) # clusterProfiler_3.14.0
library(clusterProfiler.dplyr) # clusterProfiler.dplyr_0.0.2
library(ggplot2) # ggplot2_3.3.0

######################################################################################################################################################
#################################################################################################### VOLCANO PLOT 1
### effect of KIX genotype relative to wild type (after learning).

# plot setup.
tt = fig2.results$W3$top
ylim = c(0, max(-log10(tt[,"FDR"])))
xlim=c(-2.2,1.9)
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
sig = tt[,"FDR"] <= 0.05

### save plot to SVG.
#svglite("figures/volcano.svg", height=8, width=8, system_fonts = list(sans = "Arial"))

par(mar=c(5,5,4,2))

# plot insignificant genes.
plot(tt[-which(sig),"logFC"], -log10(tt[-which(sig), "FDR"]),
    pch=16, cex=0.7,
    xlab="logFC (relative to control)",
    ylab=expression('-log'[10]*'(FDR)'),
    main="volcano plot",
    xlim=xlim,
    ylim=ylim,
    col=scales::alpha("#333333", 0.3),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)

# add colored points for significant genes.
points(tt[sig, "logFC"],-log10(tt[sig, "FDR"]),
    pch=16,
    col=scales::alpha(ifelse(tt[sig,"logFC"] > 0, "red", "dodgerblue"), 0.7)
)

# add gridlines.
abline(v=0, lty=2)
abline(h= -log10(0.05), lty=2)

# dev.off()

######################################################################################################################################################
#################################################################################################### VOLCANO PLOT 2
### effect of learning relative to homecage (wild type).

# plot setup.
tt = fig3.results$W4$top
ylim = c(0, max(-log10(tt[,"FDR"])))
xlim=c(-2.5,3.7)
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
sig = tt[,"FDR"] <= 0.05

### save plot to SVG.
#svglite("figures/volcano.svg", height=8, width=8, system_fonts = list(sans = "Arial"))

par(mar=c(5,5,4,2))

# plot insignificant genes.
plot(tt[-which(sig),"logFC"], -log10(tt[-which(sig), "FDR"]),
    pch=16, cex=0.7,
    xlab="logFC (relative to control)",
    ylab=expression('-log'[10]*'(FDR)'),
    main="volcano plot",
    xlim=xlim,
    ylim=ylim,
    col=scales::alpha("#333333", 0.3),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)

# add colored points for significant genes.
points(tt[sig, "logFC"],-log10(tt[sig, "FDR"]),
    pch=16,
    col=scales::alpha(ifelse(tt[sig,"logFC"] > 0, "red", "dodgerblue"), 0.7)
)

# add gridlines.
abline(v=0, lty=2)
abline(h= -log10(0.05), lty=2)

# dev.off()

######################################################################################################################################################
#################################################################################################### VOLCANO PLOT 3
### effect of KIX genotype relative to wild type (at baseline).

# plot setup.
tt = homecage.results$W0$top
ylim = c(0, 4)
xlim=c(-3.2,3.9)
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
sig = tt[,"FDR"] <= 0.05

### save plot to SVG.
svglite("figures/volcano_homecage.svg", height=8, width=8, system_fonts = list(sans = "Arial"))

par(mar=c(5,5,4,2))

# plot insignificant genes.
plot(tt[-which(sig),"logFC"], -log10(tt[-which(sig), "FDR"]),
    pch=16, cex=0.7,
    xlab="logFC (relative to control)",
    ylab=expression('-log'[10]*'(FDR)'),
    main="volcano plot",
    xlim=xlim,
    ylim=ylim,
    col=scales::alpha("#333333", 0.3),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)

# add colored points for significant genes.
points(tt[sig, "logFC"],-log10(tt[sig, "FDR"]),
    pch=16,
    col=scales::alpha(ifelse(tt[sig,"logFC"] > 0, "red", "dodgerblue"), 0.7)
)

# add gridlines.
abline(v=0, lty=2)
abline(h= -log10(0.05), lty=2)

dev.off()
