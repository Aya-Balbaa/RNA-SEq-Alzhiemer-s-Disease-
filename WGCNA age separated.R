nine= cbind(wt9m, nlf9, nlgf9)
nine2=round(nine)
row.names(nine2)=row.names(cts3)
conditionnin = factor((c(rep("nwt",12 ), rep("nlf", 12), rep("nlgf",12 ))))
sampleTable <- data.frame(conditionnin = factor((c(rep("wt",12 ), rep("nlf", 12), rep("nlgf",12 )))))
ddsage<- DESeqDataSetFromMatrix(countData=nine2, colData=sampleTable, design = ~ conditionnin)
dds <- estimateSizeFactors(ddsage)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ conditionnin, colData(ddsage))
mod0 <- model.matrix(~   1, colData(ddsage))
svseq <- svaseq(dat, mod, mod0, n.sv = 2 )
svseq$sv

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]


write.csv(dat, "nlf-aging- aging after sv=3 counts.csv")
design(ddssva) <- ~  SV1 + SV2 + conditionnin
dds9s=DESeq(ddssva)
vsti=vst(dds9s, blind=FALSE)
vsti5=assay(vsti)
design3= model.matrix(~0+conditionnin)
library(limma)
svs=model.matrix(~0+ddssva$SV1+ddssva$SV2)
rbi <- limma::removeBatchEffect(vsti5, covariate= svs, design=design3)
##for PCA
vari=apply(rbi, 1, var)
var1=as.matrix(vari)
varb2=var1[order(var1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
varb3=varb2[1:10000,]
varb3=as.matrix(varb3)
i2=c(row.names(varb3))
rlif=rbi[i2,]

##for WGCNA
madi=apply(rbi, 1, mad)
mad1=as.matrix(madi)
varb2=mad1[order(mad1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
varb3=varb2[1:15000,]
varb3=as.matrix(varb3)
i2=c(row.names(varb3))
rlif=rbi[i2,]
trlif=t(rlif)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(trlif, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net= blockwiseModules(trlif, corType="bicor", power=6, TOMType= "signed", minModuleSize=30, deepsplit= 2, pamRespectDendro= FALSE)
v=net[["colors"]]
write.csv(v, "9m-15k-bicor-P6-min 30-vst blind FALSE-SV1-2 removed.csv")
netLabels = matchLabels(net$colors, moduleLabels);
moduleLabels = net$colors
netLabels = matchLabels(net$colors, moduleLabels);
moduleLabels = net$colors
write.csv(netLabels, "color-confim.csv")
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
##relating modules to traits 
dim(trlif)
nGenes = ncol(trlif);
nSamples = nrow(trlif);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(trlif, moduleLabels)$eigengenes
MEs0 = moduleEigengenes(trlif, moduleColors)$eigengenes
write.csv(MEs0, "9m-MES0 bicor pwer6-min 30-modules-vst FALSE SV=2 removed.csv")
MEs = orderMEs(MEs0)
write.csv(MEs, "9m-MES-ordered bicor pwer6-min 30-modules-vst FALSE SV=2 removed.csv")

moduleTraitCor = cor(MEs, design3, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.csv(moduleTraitPvalue, "9m-15k-bicor pwer6-min 30-modules-vst FALSE SV removed-pvale.csv")
write.csv(moduleTraitCor, "9m-15k-bicor pwer8-min 30-vst FALSE SV=2 removed-moduleTraitCor.csv")
##then
sizeGrWindow(10,8)
sizeGrWindow(12,9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = paste(signif(moduleTraitCor, 2));
textMatrix = paste(signif(moduleTraitPvalue, 1));

par(mar = c(8, 9, 2, 2));
design3=model.matrix(~0+c)
design2=as.data.frame(design3)
dim(design2)
dim(moduleTraitCor)
class(design2)
class(moduleTraitCor)
sizeGrWindow(9,7)
par(mfrow = c(2,2))
par(mar = c(4, 5, 4, 6));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1))






##with batch
ddsage2<- DESeqDataSetFromMatrix(countData=nine2, colData=sampleTable, design = ~ batch9+conditionnin)
sampleTable <- data.frame(conditionnin = factor(c(rep("nwt",12 ), rep("nlf", 12), rep("nlgf",12 ))), batch9=batch9)
batch9=factor(c(2,2,2,2,2,2,2,2,1,1,1,1, 1,1,2,2,2,1,1,1,2,2,2,1, 2,2,2,2,2,2,2,2,2,1,1,1))
dds9=DESeq(ddsage2)
vsti=vst(dds9, blind=FALSE)
vsti5=assay(vsti)
design3= model.matrix(~0+conditionnin)
library(limma)
rbi <- limma::removeBatchEffect(vsti5, dds9$batch9, design=design3)
vari=apply(rbi, 1, var)
var1=as.matrix(vari)
varb2=var1[order(var1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
varb3=varb2[1:10000,]
varb3=as.matrix(varb3)
i2=c(row.names(varb3))
rlif=rbi[i2,]
##PCA didn't work


##18 month
eighteen=cbind(wt18m, nlf18, nlgf18)
eighteen2=round(eighteen)
row.names(eighteen2)=row.names(cts3)
condition18 = factor((c(rep("wt",12 ), rep("nlf", 11), rep("nlgf",11 ))))
sampleTable <- data.frame(condition18 = factor((c(rep("wt",12 ), rep("nlf", 11), rep("nlgf",11 )))))
ddsage<- DESeqDataSetFromMatrix(countData=eighteen2, colData=sampleTable, design = ~ condition18)
dds <- estimateSizeFactors(ddsage)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition18, colData(ddsage))
mod0 <- model.matrix(~   1, colData(ddsage))
svseq <- svaseq(dat, mod, mod0, n.sv = 2 )
svseq$sv

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

design(ddssva) <- ~  SV1 + SV2 +condition18
dds18s=DESeq(ddssva)
vsti=vst(dds18s, blind=FALSE)
vsti5=assay(vsti)
design3= model.matrix(~0+condition18)
library(limma)
svs=model.matrix(~0+ddssva$SV1+ddssva$SV2)
rbi <- limma::removeBatchEffect(vsti5, covariate= svs, design=design3)
##for PCA
vari=apply(rbi, 1, var)
var1=as.matrix(vari)
varb2=var1[order(var1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
varb3=varb2[1:10000,]
varb3=as.matrix(varb3)
i2=c(row.names(varb3))
rlif=rbi[i2,]

##for WGCNA
madi=apply(rbi, 1, mad)
mad1=as.matrix(madi)
varb2=mad1[order(mad1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
varb3=varb2[1:15000,]
varb3=as.matrix(varb3)
i2=c(row.names(varb3))
rlif=rbi[i2,]
trlif=t(rlif)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(trlif, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net= blockwiseModules(trlif, corType="bicor", power=6, TOMType= "signed", minModuleSize=30, deepsplit= 2, pamRespectDendro= FALSE)
v=net[["colors"]]
write.csv(v, "18m-15k-bicor-P6-min 30-vst blind FALSE-SV1-2 removed-repeat.csv")
netLabels = matchLabels(net$colors, moduleLabels);
moduleLabels = net$colors
netLabels = matchLabels(net$colors, moduleLabels);
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
write.csv(moduleColors, "18m-15k-bicor-P6-min 30-vst blind FALSE-SV1-2 removed-2.csv")

MEs = net$MEs;
geneTree = net$dendrograms[[1]];
##relating modules to traits 
dim(trlif)
nGenes = ncol(trlif);
nSamples = nrow(trlif);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(trlif, moduleLabels)$eigengenes
MEs0 = moduleEigengenes(trlif, moduleColors)$eigengenes
write.csv(MEs0, "18m-MES0 bicor pwer6-min 30-modules-vst FALSE SV=2 removed.csv")
MEs = orderMEs(MEs0)
write.csv(MEs, "18m-MES-ordered bicor pwer6-min 30-modules-vst FALSE SV=2 removed.csv")

moduleTraitCor = cor(MEs, design3, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.csv(moduleTraitPvalue, "18m-15k-bicor pwer6-min 30-modules-vst FALSE SV removed-pvale.csv")
write.csv(moduleTraitCor, "18m-15k-bicor pwer8-min 30-vst FALSE SV=2 removed-moduleTraitCor.csv")
##then
sizeGrWindow(10,8)
sizeGrWindow(12,9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = paste(signif(moduleTraitCor, 2));
textMatrix = paste(signif(moduleTraitPvalue, 1));

par(mar = c(8, 9, 2, 2));
design3=model.matrix(~0+c)
design2=as.data.frame(design3)
dim(design2)
dim(moduleTraitCor)
class(design2)
class(moduleTraitCor)
sizeGrWindow(9,7)
par(mfrow = c(2,2))
par(mar = c(4, 5, 4, 6));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1))


##24 months
twfm=cbind(wt24m, nlf24, nlgf24)
twfm2=round(twfm)
row.names(twfm2)=row.names(cts3)
condition24 = factor((c(rep("wt",10 ), rep("nlf", 9), rep("nlgf",6 ))))
sampleTable <- data.frame(condition24 = factor((c(rep("wt",10 ), rep("nlf", 9), rep("nlgf",6 )))))
ddsage<- DESeqDataSetFromMatrix(countData=twfm2, colData=sampleTable, design = ~ condition24)
dds <- estimateSizeFactors(ddsage)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ condition24, colData(ddsage))
mod0 <- model.matrix(~   1, colData(ddsage))
svseq <- svaseq(dat, mod, mod0, n.sv = 3 )
svseq$sv

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

design(ddssva) <- ~  SV1 + SV2 + SV3+ condition24
dds24s=DESeq(ddssva)
vsti=vst(dds24s, blind=FALSE)
vsti5=assay(vsti)
design3= model.matrix(~0+condition24)
library(limma)
svs=model.matrix(~0+ddssva$SV1+ddssva$SV2+ddssva$SV3)
rbi <- limma::removeBatchEffect(vsti5, covariate= svs, design=design3)
##for PCA
vari=apply(rbi, 1, var)
var1=as.matrix(vari)
varb2=var1[order(var1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
varb3=varb2[1:10000,]
varb3=as.matrix(varb3)
i2=c(row.names(varb3))
rlif=rbi[i2,]

##for WGCNA
madi=apply(rbi, 1, mad)
mad1=as.matrix(madi)
varb2=mad1[order(mad1[,1],decreasing=TRUE),]
varb2=as.matrix(varb2)
varb3=varb2[1:15000,]
varb3=as.matrix(varb3)
i2=c(row.names(varb3))
rlif=rbi[i2,]
trlif=t(rlif)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(trlif, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net= blockwiseModules(trlif, corType="bicor", power=6, TOMType= "signed", minModuleSize=30, deepsplit= 2, pamRespectDendro= FALSE)
v=net[["colors"]]
write.csv(v, "24m-15k-bicor-P6-min 30-vst blind FALSE-SV1-2 removed.csv")
netLabels = matchLabels(net$colors, moduleLabels);
moduleLabels = net$colors
netLabels = matchLabels(net$colors, moduleLabels);
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
##relating modules to traits 
dim(trlif)
nGenes = ncol(trlif);
nSamples = nrow(trlif);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(trlif, moduleLabels)$eigengenes
MEs0 = moduleEigengenes(trlif, moduleColors)$eigengenes
write.csv(MEs0, "24m-MES0 bicor pwer6-min 30-modules-vst FALSE SV=2 removed.csv")
MEs = orderMEs(MEs0)
write.csv(MEs, "24m-MES-ordered bicor pwer6-min 30-modules-vst FALSE SV=2 removed.csv")

moduleTraitCor = cor(MEs, design3, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.csv(moduleTraitPvalue, "24m-15k-bicor pwer6-min 30-modules-vst FALSE SV removed-pvale.csv")
write.csv(moduleTraitCor, "24m-15k-bicor pwer8-min 30-vst FALSE SV=2 removed-moduleTraitCor.csv")
##then
sizeGrWindow(10,8)
sizeGrWindow(12,9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = paste(signif(moduleTraitCor, 2));
textMatrix = paste(signif(moduleTraitPvalue, 1));

par(mar = c(8, 8, 1, 1));
design3=model.matrix(~0+c)
design2=as.data.frame(design3)
dim(design2)
dim(moduleTraitCor)
class(design2)
class(moduleTraitCor)
sizeGrWindow(9,7)
par(mfrow = c(2,2))
par(mar = c(4, 5, 4, 6));
table(moduleColors)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               cex.lab.y = 0.7,
               font.lab.y = 1,
               zlim = c(-1,1))
merge = mergeCloseModules(trlif, moduleColors, cutHeight = 0.25, verbose = 3)
