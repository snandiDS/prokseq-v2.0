import re
import yaml
def parseIn(file):

   param={}
   param['paramBowtie']={} 
   param['paramBowtie']['I']=0
   param['paramBowtie']['X']=500
   param['paramBowtie']['k']=1
   param['paramBowtie']['p']=1
   param['paramSalmonIDX']={}
   param['paramSalmonIDX']['-k']=31
   param['paramSalmonQNT']={}
   param['paramSalmonQNT']['-l']="A"
   param['paramSalmonQNT']['-p']=2
   param['paramSalmonQNT']['--validateMappings']=""
   #paramSalmonQNT['-o']="sal_quant"
   param['paramFeatCnt']={}
   param['paramCntFile']="test"
   param['paramGeneBdyCov']={}
   param['paramFastqcPath']=""
   param['paramBowtiePath']=""
   param['paramRootPath']=""
   param['paramPYPYPath']=""
   param['paramRUVseq']={}
   param['paramRUVseq']['nc']=0
   param['paramRUVseq']['GeneListFile'] = "F"
   param['paramSAMTOOLSPath']=""
   param['paramSalmonPath']=""
   param['paramGnBdCovPath']=""
   param['paramFeatCntPath']=""
   param['paramReadFastaPath']=""
   param['default']=0  
   param['geneBdyCov']=0 
   param['transcriptAl']=0
   param['batchEfRm']=0
   param['pathway']={} 
   param['pathwayAna']=0
   param['transcriptAl']=0
   param['paramAfterqc']={}

   with open(file) as f:
      data = yaml.load(f, Loader=yaml.FullLoader)

   #for key in data:
   if data.get('BOWTIE') is not None:
      param['paramBowtie'] = data['BOWTIE']
   if data.get('SALMON') is not None:
      param['paramSalmonIDX'] = data['SALMON']['SALMONINDEX']
      param['paramSalmonQNT'] = data['SALMON']['SALMONQUANT']
      param['transcriptAl']=1
   if data.get('AFTERQC') is not None:
      param['paramAfterqc'] = data['AFTERQC']
   if data.get('FEATURECOUNTS') is not None:
      param['paramFeatCnt'] = data['FEATURECOUNTS']
   if data.get('COUNTFILE') is not None:
      param['paramCntFile'] = data['COUNTFILE']
   if data.get('geneBody_coverage') is not None:
      param['paramGeneBdyCov'] = data['geneBody_coverage']
      param['geneBdyCov'] = 1
   if data.get('BATCH_EFFECT_REMOVE') is not None:
      if(re.match("TRUE",data['BATCH_EFFECT_REMOVE'])):
         param['batchEfRm'] = 1
   if data.get('RUVSeq') is not None:
      param['paramRUVseq'] = data['RUVSeq']
   if data.get('PATHWAY') is not None:
      param['pathwayAna']=1
      param['pathway']["cutoffP"] = data['PATHWAY']['cutoffPositive']
      param['pathway']["cutoffN"] = data['PATHWAY']['cutoffNegative']
      param['pathway']["org"] = data['PATHWAY']['Organism']
      param['pathway']["term2geneFile"] = data['PATHWAY']['TERM2GENE']
      param['pathway']["term2nameFile"] = data['PATHWAY']['TERM2NAME']
   if data.get('PATH FASTQC') is not None:
      param['paramFastqcPath'] = data['PATH FASTQC']
   if data.get('PATH BOWTIE') is not None:
      param['paramBowtiePath'] = data['PATH BOWTIE']
   if data.get('PATH ROOT') is not None:
      param['paramRootPath'] = data['PATH ROOT']
   if data.get('PATH PYPY') is not None:
      param['paramPYPYPath'] = data['PATH PYPY']
   if data.get('PATH SAMTOOLS') is not None:
      param['paramSAMTOOLSPath'] = data['PATH SAMTOOLS']
   if data.get('PATH SALMON') is not None:
      param['paramSalmonPath'] = data['PATH SALMON']
   if data.get('PATH geneBody_coverage') is not None:
      param['paramGnBdCovPath'] = data['PATH geneBody_coverage']
   if data.get('PATH FEATURECOUNTS') is not None:
      param['paramFeatCntPath'] = data['PATH FEATURECOUNTS']
   if data.get('PATH READFASTA') is not None:
      param['paramReadFastaPath'] = data['PATH READFASTA']
   if data.get('PATH DEFAULT') is not None:
      param['default']=1

   return(param)

def parseSample(fl):

   #    READING SAMPLE FILE
   #    -------------------
   fo = open(fl, "r")

   param={}
   param['SAMPLE']={}
   param['conditions']={}
   param['pe']=0

   while True:
      line = fo.readline()
      # check if line is not empty
      if not line:
         break
      if re.search("^GENOME", line):
         files = line.split()
         param['SAMPLE'][files[0]]=[]
         param['SAMPLE'][files[0]].append(files[1])
         param['SAMPLE'][files[0]].append(files[2])
         #print("%s %s %s"% (files[0],files[1],files[2]))
      if re.search("^FASTQ", line):
         files = line.split()
         if(len(files) == 4):
            param['SAMPLE'][files[2]] = []
            param['SAMPLE'][files[2]].append(files[1])
            param['conditions'][files[3]]=1
         if(len(files) == 5):
            param['pe']=1
            param['SAMPLE'][files[3]] = []
            param['SAMPLE'][files[3]].append(files[1])
            param['SAMPLE'][files[3]].append(files[2])
            param['conditions'][files[4]]=1
   fo.close()

   return(param)
