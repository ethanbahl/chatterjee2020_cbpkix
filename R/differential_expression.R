######################################################################################################################################################
########## The CBP KIX domain regulates long-term memory and circadian activity: RNA-Seq Analysis
########## Ethan Bahl
########## May 28, 2020

########## DIFFERENTIAL EXPRESSION.

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

########## load data.
# read count data.
count.data = as.matrix(read.csv("data/counts_all.csv", row.names=1, header=TRUE))
# read sample information, set up factors/levels.
sample.info = read.csv("data/sample_info.csv", row.names=1, header=TRUE, stringsAsFactors=FALSE)
sample.info$experiment = factor(sample.info$experiment, levels=c("experiment1", "experiment2"))
sample.info$genotype = factor(sample.info$genotype, levels=c("WT", "CBP KIX/KIX"))
sample.info$exposure = factor(sample.info$exposure, levels=c("homecage", "learning"))
sample.info$group = factor(sample.info$group, levels=c("EXP1.WT.HOME", "EXP1.WT.MWM", "EXP1.KIX.MWM", "EXP2.WT.HOME", "EXP2.WT.CFC", "EXP2.KIX.HOME", "EXP2.KIX.CFC"))
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
#################################################################################################### SECTION 03: Figure 3 (differential expression).
########## select homecage samples, create DGEList for edgeR.
set = set.all[, pData(set.all)$exposure == "homecage" & pData(set.all)$experiment=="experiment2"]
pData(set) = droplevels(pData(set))
y = DGEList(counts(set), samples=pData(set), genes=fData(set))

########## filter undetected / lowly expressed genes.
design_filter = model.matrix(~ genotype, data=pData(set))
keep = filterByExpr(y, design=design_filter)
y = y[keep, , keep.lib.sizes=FALSE]
set = set[keep,]

########## EDASeq normalization.
# normalize for GC content differences.
dataWithin <- withinLaneNormalization(set, "gc", which="full", offset=TRUE)
# normalize for sequencing depth.
dataNorm <- betweenLaneNormalization(dataWithin, which="upper", offset=TRUE)

########## create lists to hold data and results.
homecage.data = list(set=set, dataWithin=dataWithin, dataNorm=dataNorm)
homecage.results = list()

########## first pass of edgeR without RUV normalization.
# set up objects, attach EDASeq offsets to DGEList, create design matrix.
set = homecage.data$dataNorm
y = DGEList(counts=counts(set), samples=pData(set), genes=fData(set))
y$offset = -offst(set)
design = model.matrix(~genotype, data=pData(set))
# edgeR glmQLF (robust) workflow.
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
test = glmQLFTest(fit)
top = topTags(test, Inf)$table
# store results from first pass.
homecage.results[[1]] = list(set=set, y=y, fit=fit, test=test, top=top)

# ########## second pass of edgeR with RUV normalization.
# # RUVr normalization.
# set.ruv = RUVr(set, rownames(set), k=1, residuals(fit, type="deviance"))
# # reattach gene information.
# fData(set.ruv) = fData(set)[rownames(set.ruv),]
# # build new design matrix with latent factors of unwanted variation included.
# design = model.matrix(~ genotype + W_1, data=pData(set.ruv))
# # second pass of edgeR.
# y = DGEList(counts=counts(set), samples=pData(set), genes=fData(set))
# y$offset = -offst(set)
# y = estimateDisp(y, design=design, robust=TRUE)
# fit = glmQLFit(y, design, robust=TRUE)
# test = glmQLFTest(fit, coef=2)
# top = topTags(test, Inf)$table
# # store results from second pass.
# homecage.results[[2]] = list(set=set.ruv, y=y, fit=fit, test=test, top=top)
#ames(homecage.results) = c("W0", "W1")
names(homecage.results) = c("W0")
