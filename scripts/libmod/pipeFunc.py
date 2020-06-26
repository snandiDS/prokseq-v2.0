import re
import os
import subprocess
import pandas as pd
from libmod import execCmd as ec
from libmod import errMsgFn as errMsg

def samCheck (path,reqVer):
   ver=0
   cmd = path + "/samtools --version"
   print(cmd)
   samV = subprocess.Popen([cmd], shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
   endOfPipe = samV.stdout
   for line in endOfPipe:
      if line.startswith('samtools'):
         ver = float(line.split()[1])
         continue
   if(ver>=reqVer):
      return 1
   else:
      errMsg.samtoolErrMsg()
      return 0

def checkTools (path,tool):
   if (re.match("salmon", tool)):
      cmd = path + "/salmon --version"
   if (re.match("samtools", tool)):
      cmd = path + "/samtools --version"
   if (re.match("fastqc", tool)):
      cmd = path + "/fastqc --version"
   if (re.match("bowtie", tool)):
      cmd = path + "/bowtie2 --version"
   if (re.match("afterqc", tool)):
      cmd = path + "/pypy -h"
   if (re.match("readFasta", tool)):
      cmd = path + "/readFasta test.fasta"
   if (re.match("R", tool)):
      cmd = path + "/R --version"
   if (re.match("GnBdCov", tool)):
      inFile = path + "geneBody_coverage.py"
      if not os.path.isfile(inFile):
         return 0
      else:
         return 1
   if (re.match("FeatCnt", tool)):
      inFile = path + "/featureCounts"
      if not os.path.isfile(inFile):
         return 0
      else:
         return 1
   tl = subprocess.Popen([cmd], shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
   output = tl.stdout.read()
   if output:
      return 1
   else:
      return 0

def checkRTools (tools):
   misRpkg = []
   fo = open("check.r", "w")
   fo.write("a<-installed.packages()\n")
   fo.write("packages<-a[,1]\n")
   for i in range (0,len(tools),1):
      fo.write("message(\"%s : \",is.element(\"%s\", packages))\n"%(tools[i],tools[i]))
   fo.close()
   cmd = "Rscript --vanilla check.r"
   rrun = subprocess.Popen([cmd], shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
   out = rrun.communicate()[1]
   print(out)
   for i in range(0,len(tools),1):
      rt = tools[i] + " : TRUE"
      if rt not in out: misRpkg.append(tools[i])
   cmd = "rm check.r"
   os.system(cmd)
   return misRpkg
 


def salmonCnt2csv (flNames,csv):
   #flNames = ["sal_quant1","sal_quant2"]
   transcriptLength = {}
   transcriptTPM = {}
   transcriptCount = {}
   for i in range (0,len(flNames),1):
      transcriptLength[flNames[i]] = {}
      transcriptTPM[flNames[i]] = {}
      transcriptCount[flNames[i]] = {}

   geneList = []

   for i in range (0,len(flNames),1):
      file = flNames[i] + "/quant.sf"
      fo = open(file, "r")
      genes = []
      while True:
         line = fo.readline()
         if not line:
            break
         if (re.search("^Name", line)):
            continue
         elem = line.split()
         genes.append(elem[0])
         transcriptLength[flNames[i]][elem[0]] = elem[2]
         transcriptTPM[flNames[i]][elem[0]] = elem[3]
         transcriptCount[flNames[i]][elem[0]] = round(float(elem[4]))
      if(i==0):
         geneList = genes
      else:
         geneList = list(set(geneList) & set(genes))
         #geneList = newList
      fo.close()
   
   fo = open(csv, "w")
   fo.write(",Geneid,Chr,Start,End,Strand,EffectiveLength")

   for j in range (0,len(flNames),1):
      fo.write(",%s"%flNames[j])
      # ABOVE LINE TO BE TASTED -- TO BE CHECKED
      #fo.write(",%s.count"%flNames[j])
   for j in range (0,len(flNames),1):
      fo.write(",%s.TPM"%flNames[j])
   fo.write("\n")

   for i in range (0,len(geneList),1):
      fo.write("%i,%s,NA,NA,NA,NA"%(i,geneList[i])),
      #for j in range (0,len(flNames),1):
      fo.write(",%.2f"%float(transcriptLength[flNames[0]][geneList[i]])),
      for j in range (0,len(flNames),1):
         fo.write(",%.2f"%float(transcriptCount[flNames[j]][geneList[i]])),
      for j in range (0,len(flNames),1):
         fo.write(",%.2f"%float(transcriptTPM[flNames[j]][geneList[i]])),
      fo.write("\n")


def cnt2LenNormCsv (samFiles,csvFile,method):
   #csvOut = csvFile[:len(csvFile)-3] + "LengthNorm.csv"
   csvOut = csvFile[:len(csvFile)-3] + "NucleotideAvgCount.csv"
   df = pd.read_csv(csvFile,na_filter=False)
   #dataTop = df.columns
   for i in range (0,len(samFiles),1):
     col = samFiles[i] + ".Lnorm"
     if(re.match("Bowtie",method)):
        df[col] = (df[samFiles[i]]/df['Length']).round(4)
     if(re.match("Salmon",method)):
        df[col] = (df[samFiles[i]]/df['EffectiveLength']).round(4)
        # ABOVE LINE TO BE TASTED -- TO BE CHECKED
        #df[col] = (df[samFiles[i] + ".count"]/df['EffectiveLength']).round(4)
   df.to_csv(csvOut)
   #df.to_csv(csvFile)

   
def chkFile (inFile,mod):
   if not os.path.isfile(inFile):
      print("File (%s) not found!"%inFile)
      if(re.match(mod,"bowtie")): errMsg.bowtieErrMsg()
      if(re.match(mod,"FILE")): errMsg.fileErrMsg()
      if(re.match(mod,"btBuild")): errMsg.btBuildErrMsg()
      exit()
   else:
      print("File (%s) found."%inFile)

def chkPath (inFile):
   if not os.path.exists(inFile):
      print("Directory / Path is not prepared or found")
      exit()

def samToBam (path, inSam):
   print("Running SAM to BAM")
   #inSam = "Planktonic_36h_2.sam"
   otBam = inSam[:len(inSam)-3] + "bam" 
   cmd = path + "/samtools view -bS "+inSam+" > "+otBam
   print(cmd)
   ec.executeCmds(cmd,"samtools")

   #otBamSorted = inSam[:len(inSam)-3] + "sortedBAM"
   otBamSorted = inSam[:len(inSam)-3] + "sortedBAM.bam"
   cmd = path + "/samtools sort "+ otBam + " -o " + otBamSorted
   print(cmd)
   ec.executeCmds(cmd,"samtools")

   #otBamSorted += ".bam"
   cmd = path + "/samtools index " + otBamSorted
   print(cmd)
   ec.executeCmds(cmd,"samtools")


def bmToBigWig (path,bamFile,chrSize):
   print("Running BAM to BIGWIG")
   cmd = path + "/samtools view -f 0x10 -b "+ bamFile+" | "+ path + "/samtools mpileup -S - > tbw"
   #cmd = "python2 "+paramGnBdCovPath + "bam2wig.py -s " + chr.sizes + " -i " + bamFile + " -o tbw"
   print(cmd)
   ec.executeCmds(cmd,"samtools")

   #fileName = bamFile[:len(bamFile)-4]
   #fo = open(fileName+".wig", "w")
   fo = open(bamFile+".wig", "w")
   #fo.write("track type=wiggle_0 name=%s description=%s\n"%(fileName,fileName))
   fo.write("track type=wiggle_0 name=%s description=%s\n"%(bamFile,bamFile))
   linePt=1
   lastC=""
   fi = open("tbw", "r")
   while True:
      line = fi.readline()
      if not line:
         break
      array = line.split()
      c = array[0]
      start = int(array[1])
      depth = int(array[3])
      if not (re.match(c,lastC)): fo.write("variableStep chrom=%s\n"%c)
      lastC = c
      if (((linePt % 10) ==0) and (depth>=3)): fo.write("%d\t%d\n"% (start, depth))
      linePt+=1
   fi.close()
   fo.close()

   cmd = "rm tbw"
   os.system(cmd)

   #cmd = "depend/wigToBigWig "+fileName+".wig "+chrSize+" "+fileName+".bw"
   cmd = "depend/wigToBigWig "+bamFile+".wig "+chrSize+" "+bamFile+".bw"
   print(cmd)
   ec.executeCmds(cmd,"wigtobw")

def createGeneBodyPlot(endOfPipe):

   ARRAY=[]
   rowNames=[]
   for FILE in endOfPipe:
      print("FILES %s\n"%FILE.rstrip())
      with open(FILE.rstrip(), 'r') as fo:
         AR = fo.readline().split('<- c(')
         rowNames.append(AR[0])
         print(AR[0])
         VAL = AR[1].rstrip().split(",")
         VAL[len(VAL)-1] = VAL[len(VAL)-1].split(")")[0]
         print(VAL[len(VAL)-1])
         for i in range(0,len(VAL),1):
            ARRAY.append(float(VAL[i]))

      fo.close()
   fo = open("geneBodyCoverage.curves.r","w")
   fo.write("mat <- matrix(c(")
   for i in range(0,len(ARRAY)-1,1):
      fo.write("%.2f,"%ARRAY[i])
   strRline=str(ARRAY[len(ARRAY)-1])+"), nrow="+str(len(rowNames))+", byrow=TRUE)\n"
   fo.write(strRline)
   fo.write("rownames(mat) <- c(")
   for i in range(0,len(rowNames)-1,1):
      fo.write("\"%s\","%rowNames[i])
   strRline="\""+rowNames[len(rowNames)-1]+"\")\n"
   fo.write(strRline)
   fo.write("library(\"ggplot2\")\n")
   fo.write("df2 = stack(as.data.frame(t(mat)))\n")
   fo.write("x=seq(1,100,by=1)\n")
   strRline="x = rep(x, "+str(len(rowNames))+")\n"
   fo.write(strRline)
   fo.write("pdf(\"RSeQC.geneBodyCoverage.curves.pdf\")\n")
   fo.write("qplot(x=x, y=values, color=ind, data=df2, geom=\"line\", xlab=\"Gene body percentile (5'->3')\", ylab=\"Coverage\", main=\"Gene body coverage\")+coord_fixed(80)\n")
   #fo.write("matplot(t(mat),type=\"l\")")
   fo.write("dev.off()\n")
   fo.close()

   cmd = "Rscript geneBodyCoverage.curves.r"
   os.system(cmd)

def createOp(Path):
   cmds=[]
   cmds.append("rm -rf "+Path+"/Output")
   cmds.append("mkdir "+Path+"/Output")
   cmds.append("mkdir "+Path+"/Output/plots")
   cmds.append("mkdir "+Path+"/Output/QC_preFilter")
   cmds.append("mkdir "+Path+"/Output/QC_afterFilter")
   cmds.append("mkdir "+Path+"/Output/bam")
   cmds.append("mkdir "+Path+"/Output/sam")
   cmds.append("mkdir "+Path+"/Output/AlignmentStat")
   cmds.append("mkdir "+Path+"/Output/wig")
   cmds.append("mkdir "+Path+"/Output/Results")
   for i in range(0,len(cmds),1):
      os.system(cmds[i])

def cleanUp(Path,mode):
   cmds=[]
   if(re.match("Bowtie",mode)):
      cmds.append("mv "+Path+"/*.bam "+Path+"/Output/bam/")
      cmds.append("mv "+Path+"/*.bai "+Path+"/Output/bam/")
      cmds.append("mv "+Path+"/*.sam "+Path+"/Output/sam/")
      cmds.append("mv "+Path+"/*.stat "+Path+"/Output/Results/")
      cmds.append("rm "+Path+"/*.r")
      cmds.append("mv "+Path+"/*.bw "+Path+"/Output/wig/")
      cmds.append("mv "+Path+"/*.wig "+Path+"/Output/wig/")
   cmds.append("mv "+Path+"/*.html "+Path+"/Output/QC_preFilter/")
   cmds.append("mv "+Path+"/QC "+Path+"/Output/QC_afterFilter/QCfastQ_filtered")
   cmds.append("mv "+Path+"/good "+Path+"/Output/QC_afterFilter/FilteredReads")
   cmds.append("mv "+Path+"/bad "+Path+"/Output/QC_afterFilter/RemovedReads")
   cmds.append("mv "+Path+"/*.zip "+Path+"/Output/QC_preFilter/")
   cmds.append("mv "+Path+"/*.pdf "+Path+"/Output/plots/")
   cmds.append("mv "+Path+"/*.csv "+Path+"/Output/AlignmentStat/")
   cmds.append("mv "+Path+"/*.txt "+Path+"/Output/Results/")
   for i in range(0,len(cmds),1):
      os.system(cmds[i])

#open(FH,"testbw");
#while(<FH>){
#($c, $start, undef, $depth) = split;
#if ($c ne $lastC) { print "variableStep chrom=$c\n"; };
#$lastC=$c;
#print $.," ";
#next unless $. % 10 ==0;
#print "$start\t$depth\n" unless $depth<3;
#}

def createRunR(path,fn,smpFile,parFile):
   fnFile = path + "/depend/" + fn + ".r"
   if(re.match("limma_pe",fn)): cnt = limma_pe()
   if(re.match("limma_se",fn)): cnt = limma_se()
   if(re.match("DESeq_pe",fn)): cnt = DESeq_pe()
   if(re.match("DESeq_se",fn)): cnt = DESeq_se()
   if(re.match("edgeR_pe",fn)): cnt = edgeR_pe()
   if(re.match("edgeR_se",fn)): cnt = edgeR_se()
   if(re.match("noiseq_pe",fn)): cnt = noiseq_pe()
   if(re.match("noiseq_se",fn)): cnt = noiseq_se()
   fo = open(fnFile,"w")
   fo.write(cnt)
   fo.close()
   cmd = "Rscript --vanilla "+ fnFile + " " + smpFile + " " + parFile
   print(cmd)
   ec.executeCmds(cmd,fn)
   cmd = "rm "+fnFile
   os.system(cmd)

def createRunClstPrf(path,cfP,cfN,org,gf,nf):
   fnFile = path + "/depend/clusterProf.r"
   cnt = clusterProf()
   fo = open(fnFile,"w")
   fo.write(cnt)
   fo.close()
   cmd = "Rscript --vanilla " + fnFile + " " + cfP +" "+cfN+" "+org+" "+gf+" "+nf
   print(cmd)
   ec.executeCmds(cmd,"clusterProf")
   cmd = "rm "+fnFile
   os.system(cmd)

def DESeq_se():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("DESeq2")
library("ggplot2")
#cd1 <- read.table("samples", sep = " ", skip = 1)
cd1 <- read.table(args[1], sep = " ", skip = 12)
coldata1 <-  data.frame(row.names=cd1$V3,cd1$V4)
#dat <- read.table("test.csv",row.names=1,sep=",", TRUE)
#param <- read.table("param.input", sep = " ", fill = TRUE)
param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])
dat <- read.table(fl,sep=",", TRUE)
cols = seq(8,(8+((dim(dat)[2]-7)/2))-1,1)
countdata<-dat[1:(dim(dat)[1]),cols]
lcd <- stack(log(dat[1:(dim(dat)[1]),cols]))
lcd[mapply(is.infinite, lcd)] <- NA
colnames(lcd) <- c("Count", "Conditions")
plcd=seq(1,(dim(lcd)[2] * 2),by=2)
pdf("correlation.pdf")
#p <- ggplot(lcd, aes(x = Conditions, y = Count), at=plcd) + geom_boxplot(notch=TRUE, fill='#A4A4A4', color="black")
p <- ggplot(lcd, aes(x = Conditions, y = Count, fill=Conditions), at=plcd) + geom_boxplot(notch=TRUE, color="black")
#p + geom_jitter(shape=16, position=position_jitter(0.2), col="#999999")
#p + geom_jitter(alpha=0.3, shape=16, position=position_jitter(0.2), col="#999999")
p + geom_jitter(aes(color = Conditions), alpha=0.3, shape=16, position=position_jitter(0.2))
dev.off()
pdf("correlation_violin.pdf")
p <- ggplot(lcd, aes(x = Conditions, y = Count, fill=Conditions), at=plcd) +  geom_violin(trim = FALSE)
#p + geom_jitter(color = "#999999", alpha=0.1, shape=16, position=position_jitter(0.2))
p + geom_jitter(alpha=0.1, shape=16, position=position_jitter(0.2))
dev.off()
dds=DESeqDataSetFromMatrix(countData=countdata,colData=coldata1,design=~cd1.V4)
dds2 <- DESeq(dds)
res <- results(dds2)
row.names(res) = dat$Geneid
write.table(as.matrix(res),file="DESeq2_results.txt",sep="\t",quote=F)
resultsNames(dds2)
resLFC <- lfcShrink(dds2, coef=resultsNames(dds2)[2], type="apeglm")
row.names(resLFC) = dat$Geneid
write.table(as.matrix(resLFC),file="DESeq2lfcShrink_results.txt",sep="\t",quote=F)
pdf("MA-plot.pdf")
plotMA( res, ylim = c(-1, 1) )
dev.off()

pdf("DESeq-Volcano.pdf")
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(res)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()

pdf("pValue-histogram.pdf")
df <- data.frame(pn=factor(rep("P",each=dim(res)[1])), pvalue=res$pvalue)
#hist( res$pvalue, breaks=20, col="grey" )
ggplot(df, aes(x=pvalue)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.2, fill="#FF6666") + ggtitle("pValue histogram") + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))
dev.off()
rld <- rlogTransformation(dds, blind=TRUE)
data <- plotPCA(rld, intgroup=c("cd1.V4"),returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))
pdf("pca.pdf")
#ggplot(data, aes(x=PC1, y=PC2, color=group)) + geom_point(shape=23, fill="blue", color="darkred", size=3)
ggplot(data, aes(x=PC1, y=PC2, color=group)) + geom_point(shape=19, size=3) + ggtitle("PCA plot") + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))
dev.off()'''
   return cnt

def DESeq_pe():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("DESeq2")
library("ggplot2")
cd1 <- read.table(args[1], sep = " ", skip = 12)
coldata1 <-  data.frame(row.names=cd1$V4,cd1$V5)
param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])
dat <- read.table(fl,sep=",", TRUE)
cols = seq(8,(8+((dim(dat)[2]-7)/2))-1,1)
countdata<-dat[1:(dim(dat)[1]),cols]
lcd <- stack(log(dat[1:(dim(dat)[1]),cols]))
lcd[mapply(is.infinite, lcd)] <- NA
colnames(lcd) <- c("Count", "Conditions")
plcd=seq(1,(dim(lcd)[2] * 2),by=2)
pdf("correlation.pdf")
#p <- ggplot(lcd, aes(x = Conditions, y = Count), at=plcd) + geom_boxplot(notch=TRUE, fill='#A4A4A4', color="black")
p <- ggplot(lcd, aes(x = Conditions, y = Count, fill=Conditions), at=plcd) + geom_boxplot(notch=TRUE, color="black")
#p + geom_jitter(shape=16, position=position_jitter(0.2), col="#999999")
#p + geom_jitter(alpha=0.3, shape=16, position=position_jitter(0.2), col="#999999")
p + geom_jitter(aes(color = Conditions), alpha=0.3, shape=16, position=position_jitter(0.2))
dev.off()
pdf("correlation_violin.pdf")
p <- ggplot(lcd, aes(x = Conditions, y = Count, fill=Conditions), at=plcd) +  geom_violin(trim = FALSE)
p + geom_jitter(alpha=0.1, shape=16, position=position_jitter(0.2))
#p + geom_jitter(color = "#999999", alpha=0.1, shape=16, position=position_jitter(0.2))
dev.off()
dds=DESeqDataSetFromMatrix(countData=countdata,colData=coldata1,design=~cd1.V5)
dds2 <- DESeq(dds)
res <- results(dds2)
row.names(res) = dat$Geneid
write.table(as.matrix(res),file="DESeq2_results.txt",sep="\t",quote=F)
resultsNames(dds2)
resLFC <- lfcShrink(dds2, coef=resultsNames(dds2)[2], type="apeglm")
#resLFC <- lfcShrink(dds2, coef="cd1.V5_treat_vs_control", type="apeglm")
row.names(resLFC) = dat$Geneid
write.table(as.matrix(resLFC),file="DESeq2lfcShrink_results.txt",sep="\t",quote=F)
pdf("MA-plot.pdf")
plotMA( res, ylim = c(-1, 1) )
dev.off()

pdf("DESeq-Volcano.pdf")
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(res)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()

pdf("pValue-histogram.pdf")
df <- data.frame(pn=factor(rep("P",each=dim(res)[1])), pvalue=res$pvalue)
#hist( res$pvalue, breaks=20, col="grey" )
ggplot(df, aes(x=pvalue)) + geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.2, fill="#FF6666") + ggtitle("pValue histogram") + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))
dev.off()
rld <- rlogTransformation(dds, blind=TRUE)
data <- plotPCA(rld, intgroup=c("cd1.V5"),returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))
pdf("pca.pdf")
#ggplot(data, aes(x=PC1, y=PC2, color=group)) + geom_point(shape=23, fill="blue", color="darkred", size=3)
ggplot(data, aes(x=PC1, y=PC2, color=group)) + geom_point(shape=19, size=3) + ggtitle("PCA plot") + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),axis.title = element_text(size = rel(1.25)))
dev.off()'''
   return cnt

def edgeR_se():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

samples <- read.table(args[1],skip=12)
group = factor(samples$V4)
param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])
datE <- read.csv(fl, sep=",", header=TRUE)
cols = seq(8,(8+((dim(datE)[2]-7)/2))-1,1)
x <- data.frame(datE[1:(dim(datE)[1]),cols])
row.names(x) <- datE$Geneid
idx <- which(rowSums(x)==0)
if(length(idx) > 0){ xx=x[-idx,] } else { xx = x}
library("edgeR")
y =DGEList(counts=xx,group=group)
tmm <- calcNormFactors(y,method="TMM")
tmm <- estimateCommonDisp(tmm)
tmm <- estimateTagwiseDisp(tmm)
de.com <- exactTest(tmm)
results <- topTags(de.com,n=Inf)
write.table(as.matrix(results$table),file="edgeR_results.txt",sep="\t",quote=F)

library(ggplot2)
data <- read.table("edgeR_results.txt", header=TRUE)
data$threshold = as.factor(data$FDR <= 0.05)
g <- ggplot(data=data, aes(x=logFC, y =-log10(FDR), colour=threshold)) + geom_point(alpha=0.4, size=1.75) + xlim(c(-14, 14)) + xlab("log2 fold change") + ylab("-log10 FDR") + theme_bw() + theme(legend.position="none")
pdf("edgeRvolcano.pdf")
g
dev.off()'''
   return cnt

def edgeR_pe():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

samples <- read.table(args[1],skip=12)
group = factor(samples$V5)
param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])
datE <- read.csv(fl, sep=",", header=TRUE)
cols = seq(8,(8+((dim(datE)[2]-7)/2))-1,1)
x <- data.frame(datE[1:(dim(datE)[1]),cols])
row.names(x) <- datE$Geneid
idx <- which(rowSums(x)==0)
if(length(idx) > 0){ xx=x[-idx,] } else { xx = x}
library("edgeR")
y =DGEList(counts=xx,group=group)
tmm <- calcNormFactors(y,method="TMM")
tmm <- estimateCommonDisp(tmm)
tmm <- estimateTagwiseDisp(tmm)
de.com <- exactTest(tmm)
results <- topTags(de.com,n=Inf)
write.table(as.matrix(results$table),file="edgeR_results.txt",sep="\t",quote=F)

library(ggplot2)
data <- read.table("edgeR_results.txt", header=TRUE)
data$threshold = as.factor(data$FDR <= 0.05)
g <- ggplot(data=data, aes(x=logFC, y =-log10(FDR), colour=threshold)) + geom_point(alpha=0.4, size=1.75) + xlim(c(-14, 14)) + xlab("log2 fold change") + ylab("-log10 FDR") + theme_bw() + theme(legend.position="none")
pdf("edgeRvolcano.pdf")
g
dev.off()'''
   return cnt

def noiseq_se():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("NOISeq")
samples <- read.table(args[1],skip=12)
group = factor(samples$V4)
param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])
getData <- read.table(fl,header=T, sep=",")
seqName <- getData[,2]
cols = seq(8,(8+((dim(getData)[2]-7)/2))-1,1)
getData <- getData[cols]
rownames(getData) <- seqName
mfactors <- matrix(group,nrow = length(group), ncol = 1, byrow = TRUE, dimnames = list(group,c("transcript")))
mydata <- readData(data=getData, factors=mfactors)
repl = "biological"
if(length(group) < 3) repl = "no"
getNOIseqRes <- noiseq(mydata, k = 0.1, norm = "tmm", replicates = repl, factor="transcript", pnr = 0.2, nss = 10)
#if(length(group) < 3) getNOIseqRes <- noiseq(mydata, k = 0.1, norm = "tmm", replicates = repl, factor="transcript", pnr = 0.2, nss = 10)
#if(length(group) >2) getNOIseqRes <- noiseqbio(mydata, k = 0.5, norm = "tmm", factor = "transcript", lc = 1, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345, filter = 2)
write.table(getNOIseqRes@results[[1]], file="afterNoiseq.txt", sep="\t", row.names=T, col.names=T, quote=F)'''
   return cnt

def noiseq_pe():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("NOISeq")
samples <- read.table(args[1],skip=12)
group = factor(samples$V5)

param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])

getData <- read.table(fl,header=T, sep=",")
seqName <- getData[,2]
cols = seq(8,(8+((dim(getData)[2]-7)/2))-1,1)
getData <- getData[cols]
rownames(getData) <- seqName
mfactors <- matrix(group,nrow = length(group), ncol = 1, byrow = TRUE, dimnames = list(group,c("transcript")))
mydata <- readData(data=getData, factors=mfactors)
repl = "biological"
if(length(group) < 3) repl = "no"
getNOIseqRes <- noiseq(mydata, k = 0.1, norm = "tmm", replicates = repl, factor="transcript", pnr = 0.2, nss = 10)
#if(length(group) < 3) getNOIseqRes <- noiseq(mydata, k = 0.1, norm = "tmm", replicates = repl, factor="transcript", pnr = 0.2, nss = 10)
#if(length(group) >2) getNOIseqRes <- noiseqbio(mydata, k = 0.5, norm = "tmm", factor = "transcript", lc = 1, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345, filter = 2)
write.table(getNOIseqRes@results[[1]], file="afterNoiseq.txt", sep="\t", row.names=T, col.names=T, quote=F)'''
   return cnt

def limma_se():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("limma")
cd1 <- read.table(args[1], sep = " ", skip = 12)
batch <- cd1$V4
param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])
dat <- read.table(fl,sep=",", TRUE)
cols = seq(8,(8+((dim(dat)[2]-7)/2))-1,1)
rbeDAT <- round(as.data.frame(removeBatchEffect(dat[cols], batch)))
rbeDATa <- rbeDAT + abs(min(rbeDAT))
p1 <- seq(1,cols[1]-1,1)
p2 <- seq(cols[length(cols)]+1,dim(dat)[2],1)
df <- cbind(dat[p1], rbeDATa, dat[p2])
write.csv(df,fl, row.names = FALSE, quote=FALSE)'''
   return cnt

def limma_pe():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("limma")
cd1 <- read.table(args[1], sep = " ", skip = 12)
batch <- cd1$V5
param <- read.table(args[2], sep = " ", fill = TRUE)
fl = gsub("[[:blank:]]", "", param$V2[which('COUNTFILE' == param$V1)])
dat <- read.table(fl,sep=",", TRUE)
cols = seq(8,(8+((dim(dat)[2]-7)/2))-1,1)
rbeDAT <- round(as.data.frame(removeBatchEffect(dat[cols], batch)))
rbeDATa <- rbeDAT + abs(min(rbeDAT))
p1 <- seq(1,cols[1]-1,1)
p2 <- seq(cols[length(cols)]+1,dim(dat)[2],1)
df <- cbind(dat[p1], rbeDATa, dat[p2])
write.csv(df,fl, row.names = FALSE, quote=FALSE)'''
   return cnt

def clusterProf():
   cnt = '''#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(clusterProfiler)
cutoffP=args[1]
cutoffN=args[2]
org=args[3]
TERM2GENEfile=args[4]
TERM2NAMEfile=args[5]
edger = read.table("edgeR_results.txt", sep="\t", header=TRUE)
edgGene <- data.frame(edger[which((edger$logFC <= cutoffN) | (edger$logFC >= cutoffP)),1,drop=FALSE])
geneList <- row.names(edgGene)
TERM2GENE=read.csv(TERM2GENEfile)
TERM2NAME=read.csv(TERM2NAMEfile)
ret = grep(org,search_kegg_organism(org, by='kegg_code')$kegg_code[1], value = TRUE)
ret
org
ifelse(org %in% ret, "Found organism", "Check the short names from https://www.genome.jp/kegg/catalog/org_list.html")
universe=TERM2GENE$GENE
GOenricher=enricher(geneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE=TERM2GENE,TERM2NAME= TERM2NAME)
write.table(GOenricher, file = "GOenricher.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
goPath <- data.frame(GOenricher[,1:7])
write.table(goPath, file = "GOpathways.txt", append = FALSE, quote = FALSE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
kegg <- enrichKEGG(gene= geneList,organism= 'ypy',pvalueCutoff = 0.05)
write.table(kegg, file = "KEGGenricher.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
keggPath <- data.frame(kegg[,1:7])
write.table(keggPath, file = "KEGGpathway.txt", append = FALSE, quote = FALSE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")'''
   return cnt

