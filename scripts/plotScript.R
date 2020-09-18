"
This script will run DESeq2 and create plots from the DESeq2 results.

SYNTAX:
   Rscript ../version_2.7/scripts/DESeq2Plot.R [SE/PE] sampleFile
   Example:
      Rscript ../version_2.7/scripts/DESeq2Plot.R SE samples.bowtie.SEsample

This script should be executed where the following files are available.
1. DESeqCountData.csv
2. DESeqCountCondition.csv
3. sampleFile (example samples.bowtie.SEsample)
4. edgeR_results.txt
5. countFile.csv
6. sampleCountNames.txt

Output:
The script would create the following plots:
1. DESeqMAplot.pdf
2. DESeqMA-LFCplot.pdf
3. DESeq-Volcano.pdf
4. pValue-histogram.pdf
5. pca.pdf
6. correlation.pdf
7. correlation_violin.pdf
8. edgeRheatMap.pdf

Author: Soumyadeep Nandi and Firoj Mahmud
Date: July 16, 2020
Date: September 08, 2020
"

library("DESeq2")
library("ggplot2")
args = commandArgs(trailingOnly=TRUE)
sample = args[2]
cd1 <- read.table(sample, sep = " ", skip = 12)
if(args[1] %in% "SE"){
   coldata1 <-  data.frame(row.names=cd1$V3,cd1$V4)
   dsn=~cd1.V4
   pcagrp="cd1.V4"
}else if(args[1] %in% "PE"){
   coldata1 <-  data.frame(row.names=cd1$V4,cd1$V5)
   dsn=~cd1.V5
   pcagrp="cd1.V5"
}else{
   print("Read type mismatch. Please specify read type, SE (single-end) or PE (paired-end).")
   quit()
}

countdataF <- read.csv("DESeqCountData.csv")
countdata <- countdataF[,-1]
dds=DESeqDataSetFromMatrix(countData=countdata,colData=coldata1,design=dsn)
dds2 <- DESeq(dds)
res <- results(dds2)
row.names(res) = countdataF$Geneid
resLFC <- lfcShrink(dds2, coef=resultsNames(dds2)[2], type="apeglm")
row.names(resLFC) = countdataF$Geneid
pdf("DESeqMAplot.pdf")
plotMA( res, ylim = c(-2, 2) )
dev.off()
pdf("DESeqMA-LFCplot.pdf")
plotMA( resLFC, ylim = c(-2, 2) )
dev.off()
pdf("DESeq-Volcano.pdf")
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(res)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()
pdf("pValue-histogram.pdf")
df <- data.frame(pn=factor(rep("P",each=dim(res)[1])), pvalue=res$pvalue)
ggplot(df, aes(x=pvalue)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.2, fill="#FF6666") + ggtitle("pValue histogram") + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))
dev.off()
rld <- rlogTransformation(dds, blind=TRUE)
data <- plotPCA(rld, intgroup=c(pcagrp),returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))
pdf("pca.pdf")
ggplot(data, aes(x=PC1, y=PC2, color=group)) + geom_point(shape=19, size=3) + ggtitle("PCA plot") + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))
dev.off()
cc <- read.csv("DESeqCountCondition.csv")
pdf("correlation.pdf")
p <- ggplot(cc, aes(x = Conditions, y = Count, fill=Conditions), at=plcd) + geom_boxplot(notch=TRUE, color="black")
p + geom_jitter(aes(color = Conditions), alpha=0.3, shape=16, position=position_jitter(0.2))
dev.off()
pdf("correlation_violin.pdf")
p <- ggplot(cc, aes(x = Conditions, y = Count, fill=Conditions), at=plcd) +  geom_violin(trim = FALSE)
p + geom_jitter(alpha=0.1, shape=16, position=position_jitter(0.2))
dev.off()

library(reshape2)
library(ggplot2)
data <- read.table("edgeR_results.txt", header=TRUE)
dataF <- data[data[, "PValue"] <=0.05,]
top <- rbind(head(dataF[order(dataF$logFC),], n = 10L),tail(dataF[order(dataF$logFC),], n = 10L))
datE <- read.csv("countFile.csv", sep=",", header=TRUE)
c=read.table("sampleCountNames.txt")
cols<-as.character(c$V1)
df <- datE[match(row.names(top),datE$Geneid),cols]
row.names(df) <- row.names(top)
mdp <- max(melt(t(df))$value)/2
pdf("edgeRheatMap.pdf")
ggplot(melt(t(df)), aes(Var1,Var2, fill=value)) + geom_raster(aes(fill = melt(t(df))$value), interpolate=TRUE) + scale_fill_gradient2(name = "Gene count", low="navy", mid="white", high="red", midpoint=mdp, limits=range(melt(t(df))$value)) + geom_tile() +  xlab(label = "Sample") + ylab(label = "Genes") + theme_classic()
dev.off()

