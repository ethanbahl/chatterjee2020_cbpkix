######################################################################################################################################################
########## The CBP KIX domain regulates long-term memory and circadian activity: RNA-Seq Analysis
########## Ethan Bahl
########## May 28, 2020

########## HEATMAP.

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
#################################################################################################### HEATMAP
### genes regulated by learning & KIX genotype in response to learning

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
#svglite("figures/figure_4_left.svg", system_fonts = list(sans = "Arial"))
hleft
#dev.off()
# right.
#svglite("figures/figure_4_right.svg", system_fonts = list(sans = "Arial"))
hright
#dev.off()
