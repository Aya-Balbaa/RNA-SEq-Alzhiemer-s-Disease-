##9month in mice corresponds to 30 in human so aging analysis is done from 9m
library (DESeq2)
library(DEGreport)
library(dplyr)
wtage=cbind(wt9m, wt14, wt18m, wt24m)
wtage2=round(wtage)
row.names(wtage2)=row.names(cts3)
conditionwt = factor(c(rep("wt9",12 ), rep("wt14", 12), rep("wt18", 12), rep("wt24", 10)))
sampleTable <- data.frame(conditionwt = conditionwt)
ddsage<- DESeqDataSetFromMatrix(countData=wtage2, colData=sampleTable, design = ~ conditionwt)
ddsw <- estimateSizeFactors(ddsage)
dat  <- counts(ddsw, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ conditionwt, colData(ddsw))
mod0 <- model.matrix(~   1, colData(ddsw))
svseq <- svaseq(dat, mod, mod0, n.sv = 3 )
svseq$sv
write.csv(dat, "wt-aging counts.csv")

ddssva <- ddsw
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

write.csv(dat, "wt-aging- aging normalised counts.csv")
design(ddssva) <- ~  SV1 + SV2 + SV3+ conditionwt
ddssva$conditionwt=relevel(conditionwt, "wt9")
dds_lrt <- DESeq(ddssva, test="LRT", reduced = ~  SV1 + SV2+ SV3)
res_LRT <- results(dds_lrt)
resultsNames(dds_lrt)
r2=lfcShrink(dds_lrt, coef=7, res=res_LRT)
sig_res_LRT <- r2 %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)
write.csv(sig_res_LRT,"wtaging lfcshrink-in mod AFTER SV=3 aging interaction.csv")



##WT-NLF age interaction.

wtnf=cbind(wtage,nlfage)
wtnf=wtnf[rowSums(wtnf)>0,]
wtnf=round(wtnf)
conditionall=factor(c(rep("wt", 46), rep("nlf", 44)))
age = factor(c(rep("9",12 ),rep("14", 12), rep("18", 12), rep("24", 10),rep("9",12 ),rep("14", 12), rep("18", 11), rep("24",9 )))
sampleTable <- data.frame( conditionall = conditionall, age=age)
ddsage <- DESeqDataSetFromMatrix(countData=wtnf, colData=sampleTable, design = ~ conditionall+age+conditionall:age)
dds3 <- estimateSizeFactors(ddsage)
dat  <- counts(dds3, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ conditionall+age+conditionall:age, colData(dds3))
mod0 <- model.matrix(~   1, colData(dds3))
svseq <- svaseq(dat, mod, mod0, n.sv = 3 )
svseq$sv

ddssva <- dds3
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

design(ddssva) <- ~  SV1 + SV2 + SV3+ conditionall+age+conditionall:age
ddssva$conditionall=relevel(conditionall, "wt")
ddssva$age=relevel(age, "9")
dds_lrt <- DESeq(ddssva, test="LRT", reduced = ~SV1 + SV2 + SV3+ conditionall+age+conditionall)
res_LRT <- results(dds_lrt)
resultsNames(dds_lrt)
r2=lfcShrink(dds_lrt, coef=11, res=res_LRT)
padj.cutoff=0.05
sig_res_LRT <- r2 %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)
write.csv(sig_res_LRT,"wt-nlf aging -lfcshrink-age interaction-in mod AFTER SV=3 with age.csv")

##wt NLGF age interaction


wtngf=cbind(wtage3,nlgfage)
wtngf=round(wtngf)
row.names(wtngf)=row.names(cts3)
conditionall=factor(c(rep("wt", 34), rep("nlgf", 29)))
age = factor(c(rep("9",12 ), rep("18", 12), rep("24", 10),rep("9",12 ), rep("18", 11), rep("24",6 )))
sampleTable <- data.frame(conditionall = conditionall, age=age)
ddsage<- DESeqDataSetFromMatrix(countData=wtngf, colData=sampleTable, design = ~ conditionall+age+ conditionall:age)
dds <- estimateSizeFactors(ddsage)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ conditionall+age+conditionall:age, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 3 )

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

design(ddssva) <- ~  SV1 + SV2+ SV3+ conditionall+ age+ conditionall:age
ddssva$conditionall=relevel(conditionall, "wt")
ddssva$age=relevel(age, "9")
dds_lrt <- DESeq(ddssva, test="LRT", reduced = ~  SV1 + SV2+ SV3+ conditionall + age)
res_LRT <- results(dds_lrt)
r2=lfcShrink(dds_lrt, coef=9, res=res_LRT)
sig_res_LRT <- r2 %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)
write.csv(dat, "wt-nlgf aging normalised counts.csv")
write.csv(sig_res_LRT,"wt-nlgf-interaction lfcshrinkage in mod-after sv=3 filter aging interaction- age.csv")


##WT-hTau
sampleTable <- data.frame(conditionall = conditionall, age=age)
ddsage<- DESeqDataSetFromMatrix(countData=wtht, colData=sampleTable, design = ~ conditionall+age+ conditionall:age)
dds <- estimateSizeFactors(ddsage)
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ conditionall+age +conditionall:age, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 3 )
svseq$sv

ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
ddssva$SV4 <- svseq$sv[,4]

write.csv(dat, "wt-htau- aging normalised counts.csv")
design(ddssva) <- ~  SV1 + SV2 + SV3+ conditionall+ age+conditionall:age
ddssva$conditionall=relevel(conditionall, "wt")
ddssva$age=relevel(age, "9")
dds_lrt <- DESeq(ddssva, test="LRT", reduced = ~  SV1 + SV2+ SV3+ conditionall+ age)
res_LRT <- results(dds_lrt)
resultsNames(dds_lrt)
r2=lfcShrink(dds_lrt, coef=9, res=res_LRT)
sig_res_LRT <- r2 %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)
write.csv(sig_res_LRT,"wt-htau- lfcshrink-interaction in mod AFTER SV=3 aging interaction.csv")


##

intnfht=cbind(nfhtage, wtage3, htauage, nlfage3)
intnfht2=round(intnfht)
app=factor(c(rep("nlf", 20), rep("wta", 60), rep("nlf",32)))
htau=factor(c(rep("ht", 20),rep("mt", 34), rep("ht", 26), rep("mt", 32)))
age=factor(c(rep("9m",7), rep("18m",9 ), rep("24m",4),rep("9m",12), rep("18m",12 ), rep("24m",10), rep("9m",10), rep("18m",8 ), rep("24m",8), rep("9m",12), rep("18m",11 ), rep("24m",9)))
sampletableintnfht= data.frame(app=app, htau=htau, age=age)
ddsintnfht <- DESeqDataSetFromMatrix(countData=intnfht2, colData=sampletableintnfht, design = ~app+htau+app:age+htau:age+age:app:htau)
ddsintnfht <- DESeqDataSetFromMatrix(countData=intnfht2, colData=sampletableintnfht, design = ~app+htau+app:age+htau:age+app:htau)
ddsintnfht$app=relevel(app, "wta")
ddsintnfht$htau=relevel(htau, "mt")
ddsintnfht$age=relevel(age, "9m")

##SVA
ddsi <- estimateSizeFactors(ddsintnfht)
dat  <- counts(ddsi, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ app+htau+app:age+htau:age+app:htau, colData(ddsi))
mod0 <- model.matrix(~   1, colData(ddsi))
svseq <- svaseq(dat, mod, mod0, n.sv = 3 )
svseq$sv

ddssva <- ddsi
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]

design(ddssva) <- ~  SV1 + SV2 + SV3+ app+htau+app:age+htau:age+app:htau

ddssva$app=relevel(app, "wta")
ddssva$htau=relevel(htau, "mt")
ddssva$age=relevel(age, "9m")
dds_lrt <- DESeq(ddssva, test="LRT", reduced = ~app+htau+app:age+htau:age)
res_LRT <- results(dds_lrt)
ddsintnfht2=DESeq(ddssva)
res=results(ddsintnfht2)
resultsNames(dds_lrt)
r2=lfcShrink(dds_lrt, coef=12, res = res_LRT)
padj.cutoff=1
sig_res_LRT <- r2 %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)
write.csv(sig_res_LRT, "nlf-htau -lfcshrinkage interaction LRT- after SV=3.csv")

resi <- lfcShrink(ddsintnfht2, type="ashr", coef = 12)
res3 <- resi[ which(resi$pvalue < 0.05 ), ]
resdatai <- merge(as.data.frame(res3), as.data.frame(counts(ddssva, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdatai)
resultsNames(ddsintnfht2)
write.csv(resdatai, "nlf-htau-interaction-lfcshrink-wald-after SV=3.csv")

##wald
ddsiw=DESeq(ddssva)
resi=results(ddsiw)
resi2 <- resi[ which(resi$padj < 0.05 ), ]
resdatai <- merge(as.data.frame(resi2), as.data.frame(counts(ddsiw, normalized=TRUE)), by="row.names", sort=FALSE)
dim(resdata1)
write.csv(resdatai, "app-htau interaction without age.csv")
