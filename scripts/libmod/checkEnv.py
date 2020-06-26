#!/usr/bin/python

import re
from libmod import pipeFunc as fn
from libmod import execCmd as ec

def check(paramFile):
   paramFastqcPath=""
   paramBowtiePath=""
   paramRootPath=""
   paramPYPYPath=""
   paramSAMTOOLSPath=""
   paramSalmonPath=""
   paramGnBdCovPath=""
   paramFeatCntPath=""
   paramReadFastaPath=""
   default=0

   """    READING PARAMETER FILE
          ----------------------"""
   fo = open(paramFile, "r")
   while True:
      line = fo.readline()
      if not line:
         break
      if line.startswith('#'):
         continue
      if (re.search("^PATH FASTQC", line)):
         param = line.split()
         paramFastqcPath = param[2]
      if (re.search("^PATH BOWTIE", line)):
         param = line.split()
         paramBowtiePath = param[2]
      if (re.search("^PATH ROOT", line)):
         param = line.split()
         paramRootPath = param[2]
      if (re.search("^PATH PYPY", line)):
         param = line.split()
         paramPYPYPath = param[2]
      if (re.search("^PATH SAMTOOLS", line)):
         param = line.split()
         paramSAMTOOLSPath = param[2]
      if (re.search("^PATH SALMON", line)):
         param = line.split()
         paramSalmonPath = param[2]
      if (re.search("^PATH geneBody_coverage", line)):
         param = line.split()
         paramGnBdCovPath = param[2]
      if (re.search("^PATH FEATURECOUNTS", line)):
         param = line.split()
         paramFeatCntPath = param[2]
      if (re.search("^PATH READFASTA", line)):
         param = line.split()
         paramReadFastaPath = param[2]
      if (re.search("^PATH DEFAULT", line)):
         default=1
         break
   fo.close()
   """-------END OF READING PARAMETER FILE--------"""

   if(default):
      #paramFastqcPath = paramBowtiePath = paramPYPYPath = paramRootPath + "/depend"
      paramFastqcPath = paramRootPath + "/depend/FastQC"
      paramBowtiePath = paramRootPath + "/depend/bowtie2/bowtie2-2.3.5.1-linux-x86_64"
      paramPYPYPath = paramRootPath + "/depend/pypy2.7-v7.2.0-linux64/bin"
      paramReadFastaPath = paramRootPath + "/depend"
      paramSAMTOOLSPath = paramRootPath + "/depend/samtools/bin"
      paramSalmonPath = paramRootPath + "/depend/salmon-latest_linux_x86_64/bin"
      paramGnBdCovPath = paramRootPath + "/depend/RSeQC-2.6.2/scripts/"
      paramFeatCntPath = paramRootPath + "/depend/subread-1.4.6-p5-Linux-i386/bin/"

   #root=['/usr/local/bin','/usr/bin']
   root=['/usr/bin','/usr/local/bin']

   tools = {}
   tools['samtools']=1
   tools['salmon']=0
   tools['fastqc']=1
   tools['bowtie']=0
   tools['afterqc']=1
   tools['readFasta']=1
   tools['R']=1
   tools['GnBdCov']=1
   tools['FeatCnt']=1

   print("\nChecking required tools or packages:\n")
   if(fn.checkTools("/usr/bin","samtools")==0):
      if(fn.checkTools("/usr/local/bin","samtools")==0):
         tools['samtools'] = fn.checkTools(paramSAMTOOLSPath,"samtools")
   tools['salmon'] = fn.checkTools(paramSalmonPath,"salmon")
   tools['fastqc'] = fn.checkTools(paramFastqcPath,"fastqc")
   tools['bowtie'] = fn.checkTools(paramBowtiePath,"bowtie")
   tools['afterqc'] = fn.checkTools(paramPYPYPath,"afterqc")
   tools['GnBdCov'] = fn.checkTools(paramGnBdCovPath,"GnBdCov")
   tools['FeatCnt'] = fn.checkTools(paramFeatCntPath,"FeatCnt")
   fo = open("test.fasta", "w")
   fo.write(">test\n")
   fo.write("ATGCGATCGTA\n")
   fo.close()
   tools['readFasta'] = fn.checkTools(paramReadFastaPath,"readFasta")
   for i in range (0,len(root),1):
      tools['R'] = fn.checkTools(root[i],"R")
      if tools['R']:
         break

   pkgs = []
   bwtsal = 0
   misPkg = 0
   misRpkg = []
   for x in tools:
      if(tools[x]):
         print("%15s : Found"%x)
      else:
         print("%15s : Not found"%x)
         pkgs.append(x)
         if(re.match("bowtie",x)):
            bwtsal += 1
            if(tools['salmon']==0):
               misPkg += 1
         elif(re.match("salmon",x)):
            bwtsal += 1
            if(tools['bowtie']==0):
               misPkg += 1
         else:
            misPkg += 1

   print("\nChecking R packages:\n")
   rPkgs = ['DESeq2','ggplot2','edgeR','NOISeq','limma','clusterProfiler','apeglm']
   misRpkg = fn.checkRTools(rPkgs)


   for i in range(0,len(pkgs),1):
      pkNotFound(pkgs[i])

   if(bwtsal==2):
      print("""
         BOTH bowtiw2 AND salmon are missing.
         ATLEAST YOU SHOULD HAVE ONE.""")
      return 0
   if(misPkg>=1):
      print("""
         Some of the packages are missing.
         Try to install these packages.
         """)
      return 0
   if(len(misRpkg)>=1):
      print("Missing R packages")
      for i in range(0,len(misRpkg),1):
         print("%15s : Not found"%misRpkg[i])
      print("\n")
      return 0

   return 1



def pkNotFound(pkg):
   if(re.match("samtools",pkg)):
      print('''

SAMTOOLS:
        samtools IS NOT FOUND!

        Some features of the quality check will not be done.
        Features such as, ... TO BE FILLED BY FIROJ
        This will not hamper the overall pipeline.

        User can also install the samtools and can specify the path
        in the parameter file as
        PATH SAMTOOLS path_to_sam_folder/samtools/bin

        However, you can install samtool as follows:
        sh scripts/samtools.sh
        After successful execution of the script you will get
        the following output on your screen.

        ...
        ...
        Program: samtools (Tools for alignments in the SAM format)
        Version: 1.10 (using htslib 1.10)

        Usage:   samtools <command> [options]

        Commands:
           -- Indexing
        ...
        ...
           -- Viewing
              flags          explain BAM flags
              tview          text alignment viewer
              view           SAM<->BAM<->CRAM conversion
              depad          convert padded BAM to unpadded BAM



     ''')

   if(re.match("bowtie",pkg)):
      print('''

BOWTIE2:
        bowtie2 IS NOT FOUND!


     ''')

   if(re.match("salmon",pkg)):
      print('''

SALMON:
        salmon IS NOT FOUND!


     ''')

   if(re.match("afterqc",pkg)):
      print('''

AFTERQC:
        afterqc IS NOT FOUND!


     ''')

   if(re.match("fastqc",pkg)):
      print('''

FASTQC:
        fastqc IS NOT FOUND!


     ''')

   if(re.match("FeatCnt",pkg)):
      print('''

featureCount:
        featureCount IS NOT FOUND!


     ''')


   if(re.match("readFasta",pkg)):
      print('''

READFASTA:
        readFasta IS NOT FOUND!


     ''')

   if(re.match("GnBdCov",pkg)):
      print('''

geneBody_coverage:
        geneBody_coverage.pl IS NOT FOUND!


     ''')


