#!/usr/bin/python

'''-------------------------------------------------------------------------------------------------
Program for automated RNA-seq data analysis package for Prokaryotic. The program takes care of the
RNA-seq data analysis from quality control to pathway enrichment analysis.

    Copyright (C) .....

    Author: Soumyadeep Nandi <snandi@ggn.amity.edu>
            A K M Firoj Mahmud <firoj.mahmud@umu.se>
--------------------------------------------------------------------------------------------------'''


import os
import subprocess
import re
import pandas as pd
from libmod import pipeFunc as fn
from libmod import execCmd as ec
from libmod import checkEnv as ck
from libmod import parseInput as ps
import optparse
import threading
import math

def extract():
   file="exampleFiles.tar.gz"
   print("\n Extracting files...\n")
   if not os.path.isfile(file):
      print("\n exampleFiles.tar.gz not found!\n Please copy exampleFiles.tar.gz to this current location.\n")
      exit()
   cmd = "tar -xzf "+file
   os.system(cmd)
   print(cmd)

def runCheck(sampleFile,paramFile,np,skp):

   fn.chkFile(paramFile,"FILE")
   fn.chkFile(sampleFile,"FILE")

   reportfo = open("report","a")

   inp = {}
   inp = ps.parseIn(paramFile)
   #inp={'paramBowtie': {'-I': 0, '-X': 500, '-k': 1, '-p': 1}, 'paramSalmonIDX': {'-k': 29}, 'paramSalmonQNT': {'-l': 'A', '-p': 2, '--validateMappings': 'TRUE'}, 'paramFeatCnt': {'a': 'oldAnnotationGFF.gtf', 'o': 'FeatCount'}, 'paramCntFile': 'countFile.csv', 'paramGeneBdyCov': {'r': 'oldAnnotationGFF.bed'}, 'paramFastqcPath': pwd+'/depend/FastQC', 'paramBowtiePath': pwd+'/depend/bowtie2/bowtie2-2.3.5.1-linux-x86_64', 'paramRootPath': pwd, 'paramPYPYPath': pwd+'/depend/pypy2.7-v7.2.0-linux64/bin', 'paramSAMTOOLSPath': pwd+'/depend/samtools/bin', 'paramSalmonPath': pwd+'/depend/salmon-latest_linux_x86_64/bin', 'paramGnBdCovPath': pwd+'/depend/RSeQC-2.6.2/scripts/', 'paramFeatCntPath': pwd+'/depend/subread-1.4.6-p5-Linux-i386/bin/', 'paramReadFastaPath': pwd+'/depend', 'default': 0, 'geneBdyCov': 1, 'transcriptAl': 1, 'batchEfRm': 0, 'pathway': {'cutoffP': 2.0, 'cutoffN': -2.0, 'org': 'ypy', 'term2geneFile': 'data/TERM2GENE.csv', 'term2nameFile': 'data/TERM2NAME.csv'}, 'pathwayAna': 1, 'paramAfterqc': {}}

   pwd = return_value = os.popen('pwd').read().rstrip()
   inp['paramFastqcPath']=pwd+'/depend/FastQC'
   inp['paramBowtiePath']=pwd+'/depend/bowtie2/bowtie2-2.3.5.1-linux-x86_64'
   inp['paramRootPath']=pwd
   inp['paramPYPYPath']=pwd+'/depend/pypy2.7-v7.2.0-linux64/bin'
   inp['paramSAMTOOLSPath']=pwd+'/depend/samtools/bin'
   inp['paramSalmonPath']=pwd+'/depend/salmon-latest_linux_x86_64/bin'
   inp['paramGnBdCovPath']=pwd+'/depend/RSeQC-2.6.2/scripts/'
   inp['paramFeatCntPath']=pwd+'/depend/subread-1.4.6-p5-Linux-i386/bin/'
   inp['paramReadFastaPath']=pwd+'/depend'
   if(ck.check(inp)==0): exit()
   
   reportfo.write("Required dependencies : OK\n")

   paramBowtie={}
   paramBowtie=inp['paramBowtie']
   paramSalmonIDX={}
   paramSalmonIDX=inp['paramSalmonIDX']
   paramSalmonQNT={}
   paramSalmonQNT=inp['paramSalmonQNT']
   paramFeatCnt={}
   paramFeatCnt=inp['paramFeatCnt']
   paramCntFile=inp['paramCntFile']
   paramGeneBdyCov={}
   paramGeneBdyCov=inp['paramGeneBdyCov']
   pathway={}
   pathway=inp['pathway']
   paramAfterqc={}
   paramAfterqc=inp['paramAfterqc']
   if(inp['default']):
      #paramFastqcPath = paramBowtiePath = paramPYPYPath = paramRootPath + "/depend"
      inp['paramFastqcPath'] = inp['paramRootPath'] + "/depend/FastQC"
      inp['paramBowtiePath'] = inp['paramRootPath'] + "/depend/bowtie2/bowtie2-2.3.5.1-linux-x86_64"
      inp['paramPYPYPath'] = inp['paramRootPath'] + "/depend/pypy2.7-v7.2.0-linux64/bin"
      inp['paramReadFastaPath'] = inp['paramRootPath'] + "/depend"
      inp['paramSAMTOOLSPath'] = inp['paramRootPath'] + "/depend/samtools/bin"
      inp['paramSalmonPath'] = inp['paramRootPath'] + "/depend/salmon-latest_linux_x86_64/bin"
      inp['paramGnBdCovPath'] = inp['paramRootPath'] + "/depend/RSeQC-2.6.2/scripts/"
      inp['paramFeatCntPath'] = inp['paramRootPath'] + "/depend/subread-1.4.6-p5-Linux-i386/bin/"


   samFiles=[]
   samFilesT=[]
   samFilesC=[]
   indexedGenome=""
   #indexed = 0
   fn.createOp(inp['paramRootPath'])


   sam = {}
   sam = ps.parseSample(sampleFile)
   SAMPLE=sam['SAMPLE']
   conditions=sam['conditions']
   pe=sam['pe']

   if(len(conditions.keys()) > 2):
      print("Can't handle more than two conditions")
      exit(1)

   samples=[]
   for i in SAMPLE:
      samples.append(i)

   """ ---------END OF READING SAMPLE FILE----------"""


   """      INDEXING GENOME
          -------------------"""


   """   INDEXING GENOME FOR SALMON
        ----------------------------"""
   if(inp['transcriptAl']):
      fn.chkFile(SAMPLE['GENOME'][0],"FILE")
      cmd = inp['paramSalmonPath'] + "/" + "salmon index -t " + SAMPLE['GENOME'][0] + " -i " + SAMPLE['GENOME'][1]
      for x, y in paramSalmonIDX.items():
         cmd = cmd + " " + str(x) + " " + str(y)
      cmd = cmd +" /dev/null 2>&1"
      print(cmd)
      ec.executeCmds(cmd,"salmon")
      indexedGenome=SAMPLE['GENOME'][1]

   reportfo.write("Indexing genome (salmon): OK\n")

   """   INDEXING GENOME FOR BOWTIE2
       ------------------------------"""
   if(inp['transcriptAl']==0):
      print("Checking the annotation file for feature count provided in %s as [FEATURECOUNTS a]"%paramFile)
      fn.chkFile(paramFeatCnt['a'],"FILE")
      if(inp['geneBdyCov']):
         print("Checking the annotation file for the program geneBody_coverage from RSeQC package provided in %s as [geneBody_coverage r]"%paramFile)
         fn.chkFile(paramGeneBdyCov['r'],"FILE")

      chrSz = open("chr.sizes", "w")
      fn.chkFile(SAMPLE['GENOME'][0],"FILE")
      cmd = inp['paramReadFastaPath'] + "/readFasta " + SAMPLE['GENOME'][0]
      cmd = cmd +" /dev/null 2>&1"
      print(cmd)
      chrSeq = subprocess.Popen([cmd], shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
      endOfPipe = chrSeq.stdout
      chrSize = {}
      for line in endOfPipe:
         if line.startswith('>'):
            chr = (line.rstrip())[1:len(line.rstrip())]
            continue
         else:
            chrSize[chr] = len(line.rstrip())
            chr=''
      for chr in chrSize:
         chrSz.write("%s\t%d\n"% (chr,chrSize[chr]))
      chrSz.close()

      cmd = "mkdir "+inp['paramRootPath']+"/bowtie2_genome"
      os.system(cmd)
      fn.chkPath(inp['paramRootPath']+"/bowtie2_genome")

      print("\n\nRunning genome indexing on %s; and the output is stored as %s\n"%(SAMPLE['GENOME'][0],SAMPLE['GENOME'][1]))
      cmd = inp['paramBowtiePath'] + "/" + "bowtie2-build "+inp['paramRootPath'] + "/" +SAMPLE['GENOME'][0]+" "+inp['paramRootPath']+"/"+SAMPLE['GENOME'][1]
      cmd = cmd +" /dev/null 2>&1"
      print(cmd)
      ec.executeCmds(cmd,"bowtie")
      indexedGenome=SAMPLE['GENOME'][1]
      #indexed=1
      fn.chkFile(SAMPLE['GENOME'][1]+".1.bt2","btBuild")

   reportfo.write("Indexing genome (bowtie) : OK\n")

   """ -------END OF INDEXING GENOME --------"""


   """ RUNNING FASTQC
       --------------"""

   if(int(skp)==0):
      cmds=[]
      for fastqIDX in range(1,len(samples),1):
         if(pe):
            for j in range(0,2,1):
               fn.chkFile(SAMPLE[samples[fastqIDX]][j],"FILE")
               #cmd = paramFastqcPath +"/fastqc  -o results/qc "+ SAMPLE[samples[fastqIDX]][j]
               cmd = inp['paramFastqcPath'] +"/fastqc "+ SAMPLE[samples[fastqIDX]][j]
               #cmd = cmd +" /dev/null 2>&1"
               cmds.append(cmd)
         if(pe==0):
               fn.chkFile(SAMPLE[samples[fastqIDX]][0],"FILE")
               #cmd = paramFastqcPath +"/fastqc  -o results/qc "+ SAMPLE[samples[fastqIDX]][0]
               cmd = inp['paramFastqcPath'] +"/fastqc "+ SAMPLE[samples[fastqIDX]][0]
               #cmd = cmd +" /dev/null 2>&1"
               cmds.append(cmd)


      threads = list()
      for i in range(0,len(cmds),int(np)):
         for k in range(i,i+int(np),1):
            try:
               cmds[k]
            except IndexError:
               print("")
            else:
               print("\n\nRunning FastQC as %s\n"%cmds[k])
               pt = threading.Thread(target=ec.executeCmds,args=(cmds[k],"fastqc"))
               threads.append(pt)
               print("starting thread %s"%pt)
               pt.start()
         for j, thread in enumerate(threads):
            thread.join()

   reportfo.write("Running FastQC : OK\n")

   """ ----------END RUNNING FASTQC---------------"""

   # Creating a list of disabled parameters for AfterQC
   afqcDisPar = ["-1","READ1_FILE","-2","READ2_FILE","-7","INDEX1_FILE","-5","INDEX2_FILE","-d","INPUT_DIR","-g","GOOD_OUTPUT_FOLDER","-b","BAD_OUTPUT_FOLDER","--read1_flag","--read2_flag","--index1_flag","--index2_flag"]
   # RUNNING AFTERQC and ALIGNMENT FOR SINGLE END SAMPLES
   # ----------------------------------------------------

   if(pe==0):
      print("\n\n\nWorking with single end files\n")

   if(pe==1):
      print("\n\n\nWorking with paired-end files\n")

      for fastqIDX in range(1,len(samples),int(np)):
         sampleEXT=[]
         #CHECK IF WE HAVE ENOUGH SAMPLES i.e. LINES IN THE SAMPLE FILE
         for i in range(0,int(np),1):
            try:
               SAMPLE[samples[fastqIDX+i]][0]
            except IndexError:
               sampleEXT.append(0)
            else:
               sampleEXT.append(1)
         """ RUNNING AFTERQC
             ---------------"""
         if(int(skp)==0):
            threads = list()
            for i in range(0,len(sampleEXT),1):
               if(sampleEXT[i]):
                  print("\nRunning AfterQC with %s and %s\n"%(SAMPLE[samples[fastqIDX+i]][0],SAMPLE[samples[fastqIDX+i]][1]))
                  cmd = inp['paramPYPYPath'] + "/pypy " + inp['paramRootPath'] + "/depend/afterqc/AfterQC-master/after.py -1 "+SAMPLE[samples[fastqIDX+i]][0]+" -2 "+SAMPLE[samples[fastqIDX+i]][1]
                  for x, y in paramAfterqc.items():
                     if str(x) in afqcDisPar: continue
                     elif(y == "FALSE"): continue
                     elif(y == "TRUE"):
                        cmd = cmd + " " + str(x) + " "
                     else:
                        cmd = cmd + " " + str(x) + " " + str(y)
                  cmd = cmd +" /dev/null 2>&1"
                  print(cmd)
                  pt = threading.Thread(target=ec.executeCmds,args=(cmd,"afterqc"))
                  threads.append(pt)
                  print("starting thread %s"%pt)
                  pt.start()
            for j, thread in enumerate(threads):
               thread.join()

         reportfo.write("Running AfterQc : OK\n")

         """ END OF RUNNING AFTERQC """

         if(inp['transcriptAl']):
            print("TRANSCRIPT")
            threads = list()
            for i in range(0,len(sampleEXT),1):
               nFiles = []
               if(sampleEXT[i]):
                  samFiles.append(samples[fastqIDX+i])
                  #../bin/salmon quant -i athal_index -l A -1 reads_1.fastq -2 reads_2.fastq -p 2 --validateMappings -o sal_quant
                  nFiles.append(SAMPLE[samples[fastqIDX+i]][0])
                  if(re.search("fq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fq','good.fq', nFiles[len(nFiles)-1])
                  if(re.search("fastq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fastq','good.fq', nFiles[len(nFiles)-1])
                  nFiles.append(SAMPLE[samples[fastqIDX+i]][1])
                  if(re.search("fq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fq','good.fq', nFiles[len(nFiles)-1])
                  if(re.search("fastq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fastq','good.fq', nFiles[len(nFiles)-1])
                  #nFiles.append(re.sub('fastq','good.fq', SAMPLE[samples[fastqIDX+i]][0]))
                  #nFiles.append(re.sub('fastq','good.fq', SAMPLE[samples[fastqIDX+i]][1]))
                  cmd = inp['paramSalmonPath'] + "/" + "salmon quant -i " + SAMPLE['GENOME'][1] + " -l " + paramSalmonQNT['-l']  + " -1 good/"+ nFiles[0] + " -2 good/"+nFiles[1]
                  for x, y in paramSalmonQNT.items():
                     if(x == "-l"): continue
                     elif (y == "FALSE"): continue
                     elif(y == "TRUE"):
                        cmd = cmd + " " + str(x) + " "
                     else:
                        cmd = cmd + " " + str(x) + " " + str(y)
                  cmd = cmd + " -o " + samples[fastqIDX+i]
                  cmd = cmd +" /dev/null 2>&1"
                  print(cmd)
                  pt = threading.Thread(target=ec.executeCmds,args=(cmd,"salmon"))
                  threads.append(pt)
                  print("starting thread %s"%pt)
                  pt.start()
            for j, thread in enumerate(threads):
               thread.join()

         reportfo.write("Running salmon : OK\n")


         if(inp['transcriptAl']==0):
            threads = list()
            for i in range(0,len(sampleEXT),1):
               nFiles = []
               if(sampleEXT[i]):
                  nFiles.append(SAMPLE[samples[fastqIDX+i]][0])
                  if(re.search("fq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fq','good.fq', nFiles[len(nFiles)-1])
                  if(re.search("fastq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fastq','good.fq', nFiles[len(nFiles)-1])
                  nFiles.append(SAMPLE[samples[fastqIDX+i]][1])
                  if(re.search("fq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fq','good.fq', nFiles[len(nFiles)-1])
                  if(re.search("fastq$",nFiles[len(nFiles)-1],re.IGNORECASE)):
                     nFiles[len(nFiles)-1] = re.sub('fastq','good.fq', nFiles[len(nFiles)-1])
                  #nFiles.append(re.sub('fastq','good.fq', SAMPLE[samples[fastqIDX+i]][0]))
                  #nFiles.append(re.sub('fastq','good.fq', SAMPLE[samples[fastqIDX+i]][1]))
                  statFile = samples[fastqIDX+i]+".stat"
                  print("\n\nRunning bowtie with %s and %s\n"%(nFiles[0],nFiles[1]))
                  cmd = "rm -rf "+samples[fastqIDX+i]
                  os.system(cmd)
                  #cmd = inp['paramBowtiePath'] + "/bowtie2 -I "+str(paramBowtie['-I'])+" -X "+str(paramBowtie['-X'])+" -x "+indexedGenome+" -1 good/"+nFiles[0]+" -2 good/"+nFiles[1]+" -k "+str(paramBowtie['-k'])+" -p "+str(paramBowtie['-p'])+" -S "+samples[fastqIDX+i]+" 2> "+statFile
                  cmd = inp['paramBowtiePath'] + "/bowtie2 -x "+indexedGenome+" -1 good/"+nFiles[0]+" -2 good/"+nFiles[1]
                  for x, y in paramBowtie.items():
                     if(y == "TRUE"):
                        cmd = cmd + " " + str(x) + " "
                     elif(y == "FALSE"): continue
                     else:
                        cmd = cmd + " " + str(x) + " " + str(y)
                  cmd = cmd + " -S "+samples[fastqIDX+i]+" 2> "+statFile
                  print(cmd)
                  pt = threading.Thread(target=ec.executeCmds,args=(cmd,"bowtie"))
                  threads.append(pt)
                  print("starting thread %s"%pt)
                  pt.start()
            for j, thread in enumerate(threads):
               thread.join()

            reportfo.write("Running bowtie : OK\n")


            for i in range(0,len(sampleEXT),1):
               if(sampleEXT[i]):
                  samFiles.append(samples[fastqIDX+i])
                  fn.chkFile(samples[fastqIDX+i],"bowtie")

            threads = list()
            for i in range(0,len(sampleEXT),1):
               if(sampleEXT[i]):
                  pt = threading.Thread(target=fn.samToBam, args=(inp['paramSAMTOOLSPath'],samples[fastqIDX+i]))
                  threads.append(pt)
                  print("starting thread %s"%pt)
                  pt.start()
            for j, thread in enumerate(threads):
               thread.join()

            threads = list()
            for i in range(0,len(sampleEXT),1):
               if(sampleEXT[i]):
                  otBam = samples[fastqIDX+i][:len(samples[fastqIDX+i])-3] + "sortedBAM.bam"
                  fn.bmToBigWig(inp['paramSAMTOOLSPath'],otBam,"chr.sizes")
                  if(inp['geneBdyCov']):
                     fn.chkFile(paramGeneBdyCov['r'],"FILE")
                     cmd = "rm "+otBam+"_RSeQC.geneBodyCoverage.*"
                     os.system(cmd)
                     #cmd = "python2 "+paramRootPath + "/depend/RSeQC-2.6.2/scripts/geneBody_coverage.py -r "+paramGeneBdyCov['r']+" -i " + otBam + " -o "+otBam+"_RSeQC"
                     cmd = "python2 "+ inp['paramGnBdCovPath'] + "geneBody_coverage.py -r "+paramGeneBdyCov['r']+" -i " + otBam + " -o "+otBam+"_RSeQC"
                     cmd = cmd +" /dev/null 2>&1"
                     print(cmd)
                     pt = threading.Thread(target=ec.executeCmds,args=(cmd,"rseqc"))
                     threads.append(pt)
                     print("starting thread %s"%pt)
                     pt.start()
            for j, thread in enumerate(threads):
               thread.join()
   if(inp['transcriptAl']==0):
      rList = []
      FILES = subprocess.Popen(["ls *_RSeQC.geneBodyCoverage.r"], shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
      endOfPipe = FILES.stdout
      for line in endOfPipe:
         rList.append(line.rstrip())
      fn.createGeneBodyPlot(rList)

      print("\n\nRunning featureCounts\n")
      #cmd = paramRootPath + "/depend/subread-1.4.6-p5-Linux-i386/bin/featureCounts -a "+paramFeatCnt['a']+" -o "+paramFeatCnt['o']
      cmd = inp['paramFeatCntPath'] + "featureCounts -a "+paramFeatCnt['a']+" -o "+paramFeatCnt['o']
      for file in samFiles:
         cmd = cmd+" "+file
      cmd = cmd+" 2> /dev/null"
      cmd = cmd +" /dev/null 2>&1"
      print(cmd)
      ec.executeCmds(cmd,"featcnt")
      fn.chkFile(paramFeatCnt['o'],"FILE")

      reportfo.write("Running gene body cooverage : OK\n")

      df = pd.read_csv(paramFeatCnt['o'],sep='\t', skiprows=1)
      dataTop = df.columns
      #print(dataTop)
      for i in range (0,len(samFiles),1):
         colNam = samFiles[i]+".tpm"
         df[colNam] = df[samFiles[i]] / (df['Length']/1000)
         pmSF=(df[colNam].sum())/10**6
         df[colNam] = (df[colNam] / pmSF).round(4)
      #print(df)
      #paramFeatCnt['o'] = paramFeatCnt['o']+".csv"
      #df.to_csv(paramFeatCnt['o'])
      df.to_csv(paramCntFile)
      print("Wrote TPM values in %s"%paramCntFile)
      for i in range (0,len(samFiles),1):
         Total=df[samFiles[i]].sum()
         colNam = samFiles[i]+".cpm"
         df[colNam] = ((df[samFiles[i]] / Total ) * 1000000).round(4)
      otCPM = paramCntFile[:len(paramCntFile)-4] + "_CPM.csv"
      df.to_csv(otCPM)
      print("Wrote CPM values in %s"%otCPM)
      cmd = "mv "+paramFeatCnt['o']+" "+inp['paramRootPath']+"/Output/AlignmentStat"
      os.system(cmd)
      cmd = "mv "+paramFeatCnt['o']+".summary "+inp['paramRootPath']+"/Output/AlignmentStat"
      os.system(cmd)

      fn.cnt2LenNormCsv(samFiles,paramCntFile,"Bowtie")


   if(inp['transcriptAl']):
      fn.salmonCnt2csv(samFiles,paramCntFile)
      fn.cnt2LenNormCsv(samFiles,paramCntFile,"Salmon")

   if(inp['batchEfRm']):
      print("\n\nRemoving batch effect\n")
      if(pe==1): fn.createRunR(inp['paramRootPath'],"limma_pe",sampleFile,paramFile)
      if(pe==0): fn.createRunR(inp['paramRootPath'],"limma_se",sampleFile,paramFile)


   print("\n\nDESeq2\n")
   if(pe==1): fn.createRunR(inp['paramRootPath'],"DESeq_pe",sampleFile,paramFile)
   if(pe==0): fn.createRunR(inp['paramRootPath'],"DESeq_se",sampleFile,paramFile)
   reportfo.write("Running DESeq2 : OK\n")


   if(len(samFiles)>=4):
      print("\n\nEdgeR\n")
      if(pe==1): fn.createRunR(inp['paramRootPath'],"edgeR_pe",sampleFile,paramFile)
      if(pe==0): fn.createRunR(inp['paramRootPath'],"edgeR_se",sampleFile,paramFile)
   else:
      print("There is no replication. EdgeR will not run.\n")
   reportfo.write("Running EdgeR : OK\n")

   print("\n\nNOISeq\n")
   if(pe==1): fn.createRunR(inp['paramRootPath'],"noiseq_pe",sampleFile,paramFile)
   if(pe==0): fn.createRunR(inp['paramRootPath'],"noiseq_se",sampleFile,paramFile)
   reportfo.write("Running NOISeq : OK\n")

   if(inp['pathwayAna']):
      print("\n\nPATHWAY ANALYSIS\n")
      fn.createRunClstPrf(inp['paramRootPath'],str(pathway['cutoffP']),str(pathway['cutoffN']),pathway['org'],pathway['term2geneFile'],pathway['term2nameFile'])
   reportfo.write("Running pathway analysis : OK\n")

   reportfo.close()

def main():

   input=['param.salmon.yaml','samples.salmon.SEsample','samples.salmon.PEsample','param.bowtie.yaml','samples.bowtie.SEsample','samples.bowtie.PEsample']

   for ifile in input:
      if not os.path.isfile(ifile):
         extract()

   for i in range (0,5,3):
      sam = {}
      sam = ps.parseSample(input[i])
      SAMPLE=sam['SAMPLE']
      for x, y in SAMPLE.items():
         if(re.match("GENOME",x)):
            if not os.path.isfile(SAMPLE['GENOME'][0]):
               extract()
         else:
            for ifile in y:
               if not os.path.isfile(ifile):
                  extract()
                  #fn.chkFile(i,"FILE")


   reportfo = open("report","w")
   reportfo.write("\nRequired files : OK\n")
   reportfo.close()

   pwd = return_value = os.popen('pwd').read().rstrip()
   print(pwd)


   print("\nPlease wait. Checking SALMON. This may take several minuits. For details please visit test.log\n")
   runCheck("samples.salmon.PEsample","param.salmon.yaml",2,0)
   print("\nPlease wait. Checking BOWTIE. This may take several minuits. For details please visit test.log\n")
   runCheck("samples.bowtie.PEsample","param.bowtie.yaml",2,1)
   fn.cleanUp(pwd,"Bowtie")
   fn.cleanUp(pwd,"Salmon")

   rep={}
   reportfo = open("report","r")
   while True:
      line = reportfo.readline()
      if not line:
         break
      lne = line.rstrip()
      rep[lne]=1
   reportfo.close()
   os.system("rm report")
   for x in rep:
      print("%35s"%x)
   print("\nLooks good! \n Found all the dependencies working.\n")


if __name__ == '__main__':
        main()


