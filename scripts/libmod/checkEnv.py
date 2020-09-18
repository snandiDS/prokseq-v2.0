#!/usr/bin/python

import re
from libmod import pipeFunc as fn
from libmod import execCmd as ec

def check(Inp):
   '''
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
   '''

   if(Inp['default']):
      #paramFastqcPath = paramBowtiePath = paramPYPYPath = paramRootPath + "/depend"
      Inp['paramFastqcPath'] = Inp['paramRootPath'] + "/depend/FastQC"
      Inp['paramBowtiePath'] = Inp['paramRootPath'] + "/depend/bowtie2/bowtie2-2.3.5.1-linux-x86_64"
      Inp['paramPYPYPath'] = Inp['paramRootPath'] + "/depend/pypy2.7-v7.2.0-linux64/bin"
      Inp['paramReadFastaPath'] = Inp['paramRootPath'] + "/depend"
      Inp['paramSAMTOOLSPath'] = Inp['paramRootPath'] + "/depend/samtools/bin"
      Inp['paramSalmonPath'] = Inp['paramRootPath'] + "/depend/salmon-latest_linux_x86_64/bin"
      Inp['paramGnBdCovPath'] = Inp['paramRootPath'] + "/depend/RSeQC-2.6.2/scripts/"
      Inp['paramFeatCntPath'] = Inp['paramRootPath'] + "/depend/subread-1.4.6-p5-Linux-i386/bin/"

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
         tools['samtools'] = fn.checkTools(Inp['paramSAMTOOLSPath'],"samtools")
   tools['salmon'] = fn.checkTools(Inp['paramSalmonPath'],"salmon")
   tools['fastqc'] = fn.checkTools(Inp['paramFastqcPath'],"fastqc")
   tools['bowtie'] = fn.checkTools(Inp['paramBowtiePath'],"bowtie")
   tools['afterqc'] = fn.checkTools(Inp['paramPYPYPath'],"afterqc")
   tools['GnBdCov'] = fn.checkTools(Inp['paramGnBdCovPath'],"GnBdCov")
   tools['FeatCnt'] = fn.checkTools(Inp['paramFeatCntPath'],"FeatCnt")
   fo = open("test.fasta", "w")
   fo.write(">test\n")
   fo.write("ATGCGATCGTA\n")
   fo.close()
   tools['readFasta'] = fn.checkTools(Inp['paramReadFastaPath'],"readFasta")
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
   rPkgs = ['DESeq2','ggplot2','edgeR','NOISeq','limma','clusterProfiler','apeglm','RUVSeq','RColorBrewer']
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


