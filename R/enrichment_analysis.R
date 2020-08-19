######################################################################################################################################################
########## The CBP KIX domain regulates long-term memory and circadian activity: RNA-Seq Analysis
########## Ethan Bahl
########## May 28, 2020

########## ENRICHMENT / RRHO.

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
#################################################################################################### FISHER TEST.

# testing DEGs (DE between WT-learning and WT-homecage) for enrichment of positive controls.
fisher.test(table(fig3.results$W4$top$FDR<=0.05, fig3.results$W4$top$control=="positive"))

# testing genes regulated by KIX after learning for enrichment of spatial learning regulated genes.
fisher.test(table(fig2.results$W3$top$FDR<=0.05, fig2.results$W3$top$ensembl_gene_id %in% subset(fig3.results$W4$top, FDR<=0.05)$ensembl_gene_id))

######################################################################################################################################################
#################################################################################################### RRHO.
########## RRHO analysis.

fig2.top = fig2.results$W3$top
fig2.top$signedF = sign(fig2.top$logFC) * fig2.top$F 
fig3.top = fig3.results$W4$top
fig3.top$signedF = sign(fig3.top$logFC) * fig3.top$F 

ov.top = cbind(
    fig2.top[intersect(rownames(fig2.top), rownames(fig3.top)), ],
    fig3.top[intersect(rownames(fig2.top), rownames(fig3.top)), c(5:10)]
)

colnames(ov.top)[5:16] = paste0(colnames(ov.top)[5:16],rep(c(".fig2",".fig3"), each=6))

gene.list1 = cbind(gene=ov.top[,"ensembl_gene_id", drop=F], rank=rank(ov.top$F.fig2))
gene.list2 = cbind(gene=ov.top[,"ensembl_gene_id", drop=F], rank=rank(ov.top$F.fig3))

# set step size, run RRHO.
ss = 100
rrho.result = RRHO(gene.list1, gene.list2, stepsize=ss, BY=TRUE, alternative="enrichment", labels=c("figure 2 results", "figure 3 results"), log10.ind=TRUE)

# Pval testing.
set.seed(777)
pval.testing <- pvalRRHO(rrho.result, 50)
pval.testing$pval

# get BY-adjusted signficance matrix.
m = rrho.result$hypermat.by
# identify bins corresponding to max enrichment.
max.ind = which(m == max(m), arr.ind=TRUE)[1,]

# svglite("figures/rrho.svg", height=8, width=8, system_fonts = list(sans = "Arial"))
par(mar=c(6,5,6,7) , mgp=c(4,1,0), font.axis=1)
image(m, col=colorRampPalette(rev(brewer.pal(11, "Spectral")))(301),
   # xlab="genes regulated by KIX genotype after learning\n(P-value rank)",
   # ylab="genes regulated by learning\n(P-value rank)",
    xlab="", ylab="",
    xaxt="n", yaxt="n",
    main="RRHO: effect of KIX after learning\ncompared to effect of learning",
    cex.lab=1.4, cex.main=1.6
)

axis(1, at=c(0), labels=c("differentially\nexpressed"), tick=FALSE, hadj=0, cex=0.8)
axis(2, at=c(0), labels=c("differentially\nexpressed"), tick=FALSE, hadj=0, padj=0.5, cex=0.8)
axis(1, at=c(1), labels=c("not significant"), tick=FALSE, hadj=1, padj=-1.5, cex=0.8)
axis(2, at=c(1), labels=c("not significant"), tick=FALSE, hadj=1, padj=1, cex=0.8)
points(max.ind["row"] / nrow(m), max.ind["col"]/ncol(m), col="white", lwd=3, cex=1.5)

title(xlab="genes regulated by KIX genotype after learning\n(P-value rank)", cex.lab=1.2, font.lab=2)

par(mgp=c(2,1,0))
title(ylab="genes regulated by learning\n(P-value rank)", cex.lab=1.2, font.lab=2)

par(xpd=TRUE)
legend(1.02, 0.8, c("maximum\nenrichment"), pch = 1, col="#afafaf", pt.cex=1.5, pt.lwd=3, bty="n", cex=1)

legend(1.02, 0.6, rev(c("0", rep("", 9), round(max(m), 2))), fill=(brewer.pal(11, "Spectral")), bty="n", cex=1, y.intersp=0.6, title=expression('-log'[10]*'(FDR)'))

#dev.off()

######################################################################################################################################################
#################################################################################################### RRHO PART2
### Comparing biological signal using traditional significance vs RRHO-identified thresholds.

library(clusterProfiler)
library(org.Mm.eg.db)
library(clusterProfiler.dplyr)


ngene.row = ss * max.ind["row"]
ngene.col = ss * max.ind["col"]
fdrcut.row = ov.top[order(ov.top$F.fig2, decreasing=TRUE),"FDR.fig2"][c(ngene.row)]
fdrcut.col = ov.top[order(ov.top$F.fig3, decreasing=TRUE),"FDR.fig3"][c(ngene.col)]

# create gene groups: one for traiditional threshold, one for RRHO threshold.
gene.groups = list(FDR_strict=subset(ov.top, FDR.fig2<=0.05 & FDR.fig3<=0.05)$ensembl_gene_id, FDR_relaxed=subset(ov.top, FDR.fig2<=fdrcut.row & FDR.fig3<=fdrcut.col)$ensembl_gene_id)

# clusterProfiler.
ck = compareCluster(geneCluster=gene.groups,
    fun="enrichGO",
    OrgDb=org.Mm.eg.db,
    keyType="ENSEMBL",
    ont="ALL",
    pAdjustMethod="fdr", pvalueCutoff=1, qvalueCutoff=1,
    universe=ov.top$ensembl_gene_id,
    readable=TRUE
)

ck = filter(ck, p.adjust<=0.05)

myPalette <- colorRampPalette(rev(brewer.pal(9, "BuPu")[3:9]))
sc <- scale_colour_gradientn(colours = myPalette(101), limits=c(0, 0.05))

#svglite("figures/rrho_go_comp.svg", height=8, width=12, system_fonts = list(sans = "Arial"))
dotplot(ck, showCategory=Inf) + sc + ggtitle("gene ontology enrichment of overlapping genes", subtitle="comparison between traditional and RRHO FDR thresholds") + 
    theme(plot.title = element_text(size = 18, face = "bold"), plot.subtitle=element_text(size=14, face="bold"), axis.text.y.left=element_text(size=10)) +
    scale_x_discrete(labels= c("FDR <= 0.05", "RRHO thresholds"))
#dev.off()