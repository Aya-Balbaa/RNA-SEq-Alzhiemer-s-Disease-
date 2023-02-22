library(ggplot2)
##using vst transformed counts with blind=TRUE "rlif"
trl=t(rlif)
pr= prcomp(trl,scale=TRUE)
pr.var<- pr$sdev^2
plot(pr$x[,1], pr$x[,2])
pr.var.per <- round(pr.var/sum(pr.var)*100, 1)
barplot(pr.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
rownames(pr$x)<- paste(r, sep="" )
pr.data <- data.frame(Sample=rownames(pr$x),
                      X=pr$x[,1],Y=pr$x[,2])

ggplot(data=pr.data, aes(x=X, y=Y, color= condition, shape= age)) +
  geom_point(size=3, alpha=0.5) + 
  xlab(paste("PC1 - ", pr.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pr.var.per[2], "%", sep="")) +
  theme_bw() 