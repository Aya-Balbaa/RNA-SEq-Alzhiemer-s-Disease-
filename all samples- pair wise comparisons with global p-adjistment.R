alls=cbind(allwt,allnlf,allnlgf,htauage,allnfht,nghtage)
dim(alls)
alls=round(alls)
s=alls[rowSums(alls)>0,]
dim(s)
batchalls=factor(c(1,2,2,1,1,2,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,2,2,2,2,2,1,1,1,1,2,2,2,2,2,3,2,2,2,1,1,1,1,2,2,2,1,1,1,1,1,1,2,2,3,4,4,4,2,3,2,2,3,3,1,2,2,2,3,1,1,1,1,2,2,2,1,1,2,2,2,1,1,1,2,2,2,1,2,2,3,3,3,2,2,2,2,3,1,1,2,2,2,3,4,2,1,1,1,1,2,4,4,4,4,3,2,3,2,2,2,2,2,1,1,2,2,2,2,1,1,2,2,2,1,1,1,1,3,2,2,1,1,2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,3,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,1,1,1,3,3,3,1,1,1,3,3,3,3,3,3,3,3,3,3,3,4,4,3,3,3,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4))
b=as.matrix(batchalls)
conditionalls=factor(c(rep("wt2",12), rep("wt4",12), rep("wt9",12), rep("wt14",12), rep("wt18",12), rep("wt24",10), rep("nlf4",12), rep("nlf9", 12), rep("nlf14", 12), rep("nlf18", 11), rep("nlf24",9), rep("nlgf2",11), rep("nlgf4",12), rep("nlgf9",12), rep("nlgf18",11), rep("nlgf24", 6), rep("htau9",10), rep("htau18",8) , rep("htau24",8), rep("nht4",12), rep("nht9",7), rep("nht18",9), rep("nht24",4), rep("nght4",12), rep("nght9",7), rep("nght18",6)))
c=as.matrix(conditionalls)
sampletable= data.frame(batchalls=batchalls, conditionalls=conditionalls)
ddsall <- DESeqDataSetFromMatrix(countData=s, colData=sampletable, design = ~batchalls+conditionalls)
ddsa=DESeq(ddsall)

##LFC shrink
##for 9 months
res.A.C <- lfcShrink(ddsa, type="ashr",  contrast=c("conditionalls","wt9","htau9"))
res.w.nlf <- lfcShrink(ddsa, type="ashr",  contrast=c("conditionalls","wt9","nlf9"))
res.w.nht <- lfcShrink(ddsa, type="ashr",  contrast=c("conditionalls","wt9","nht9"))
res.nlf.nfht <- lfcShrink(ddsa, type="ashr",  contrast=c("conditionalls","nlf9","nht9"))
res.htau.nfht <- lfcShrink(ddsa, type="ashr", contrast=c("conditionalls","htau9","nht9"))


res8 <- res.A.C[ which(res.A.C$pvalue < 2 ), ]
resdata6 <- merge(as.data.frame(res8), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata)

res9 <- res.w.nlf[ which(res.w.nlf$pvalue < 2 ), ]
resdata7 <- merge(as.data.frame(res9), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata)


res10 <- res.w.nht[ which(res.w.nht$pvalue < 2 ), ]
resdata8 <- merge(as.data.frame(res10), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata)


res11 <- res.nlf.nfht[ which(res.nlf.nfht$pvalue < 2 ), ]
resdata9 <- merge(as.data.frame(res11), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata)

res12 <- res.htau.nfht[ which(res.htau.nfht$pvalue < 2 ), ]
resdata10 <- merge(as.data.frame(res12), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata)


##for 18 months
res.A.C <- lfcShrink(ddsa, type="ashr", contrast=c("conditionalls","wt18","htau18"))
res.w.nlf <- lfcShrink(ddsa, type="ashr", contrast=c("conditionalls","wt18","nlf18"))
res.w.nht <- lfcShrink(ddsa, type="ashr", contrast=c("conditionalls","wt18","nht18"))
res.nlf.nfht <- lfcShrink(ddsa, type="ashr", contrast=c("conditionalls","nlf18","nht18"))
res.htau.nfht <- lfcShrink(ddsa, type="ashr", contrast=c("conditionalls","htau18","nht18"))

sum( res.w.nlf$padj <0.05, na.rm=TRUE)

res3 <- res.A.C[ which(res.A.C$pvalue < 2 ), ]
resdata1 <- merge(as.data.frame(res3), as.data.frame(counts(ddsa, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata1)

res4 <- res.w.nlf[ which(res.w.nlf$pvalue < 2 ), ]
resdata2 <- merge(as.data.frame(res4), as.data.frame(counts(ddsa, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata2)

res5 <- res.w.nht[ which(res.w.nht$pvalue < 2 ), ]
resdata3 <- merge(as.data.frame(res5), as.data.frame(counts(ddsa, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata3)


##doing global p adjustment
p2=as.numeric(resdata1$pvalue, resdata2$pvalue, resdata3$pvalue, resdata4$pvalue, resdata5$pvalue, resdata6$pvalue)
p3=p.adjust(p2, method= "BH")
p5=as.data.frame(split(p3, ceiling(seq_along(p3)/42771)))
write.csv(p5, "global-padjusted-wt-nlf and wt-nlfhtau-18-24.csv")


##adding the the p.adjusted columns to the corresponding table
resdata12=cbind(p5$X1, resdata1)
write.csv(resdata12, "WT-Htau-shrink-18-2.csv")

resdata22=cbind(p5$X2, resdata2)
write.csv(resdata22, "WT-nlf-shrink-18.csv")

resdata32=cbind(p5$X3, resdata3)
write.csv(resdata32, "wt-nlfhtau-shrink 18.csv")

resdata62=cbind(p5$X6, resdata6)
write.csv(resdata62, "WT-htau-shrink-9.csv")

resdata72=cbind(p5$X7, resdata7)
write.csv(resdata72, "WT-nlf-shrink-9.csv")

resdata82=cbind(p5$X8, resdata8)
write.csv(resdata82, "w-nlfht-shrink-9.csv")

resdata92=cbind(p5$X9, resdata9)
write.csv(resdata92, "nlf-nlfhtau-shrink-9.csv")

resdata102=cbind(p5$X10, resdata10)
write.csv(resdata102, "htau-nlfhtau-shrink-9.csv")

