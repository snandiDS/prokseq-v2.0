def fileErrMsg():
   fileErrorMsg = """
                    File not created or found!
                    Please check if the file exists.
                    Please check the permission of the directory and give permission.
                 """

   print(fileErrorMsg)
   return 1 

def btBuildErrMsg():
   btBuildErrorMsg = """
                    IN CASE OF ERROR
                    Indexing of the genome is not done! 
                       Suggestions for most common errors:
                          1. Please check the permission of the directory bowtie2_genome.
                          2. Pleae check the genome file is correct.
                  
                  """

   print(btBuildErrorMsg)
   return 1

def pypyErrMsg():
   pypyErrorMsg = """
                    IN CASE OF ERROR
                    Suggestions for most common errors:
                    Error: Status : FAIL 127 /bin/sh: ... ... ...pypy: No such file or directory
                    Fix: Please check the path of the pypy executables.
                    pypy: error while loading shared libraries: libbz2.so.1.0: cannot open shared object file: No such file or directory
                    Fix:         ln -s /usr/lib64/libbz2.so.1.0.6 /usr/lib64/libbz2.so.1.0
                    pypy: error while loading shared libraries: libtinfo.so.5: cannot open shared object file: No such file or directory
                    Fix:         ln -s libtinfo.so.6.1 libtinfo.so.5
                  
                  """

   print(pypyErrorMsg)
   return 1

def bowtieErrMsg():
   bowtieErrorMsg = """
                    IN CASE OF ERROR
                    Suggestions for most common errors:
                    Error: Status : FAIL 127 /bin/sh: ... ... ...bowtie2 / bowtiebuild: No such file or directory
                    Fix: Please check the path of the bowtie executables.
                    Error: Permission denied.
                    Fix: Check the permission of the directory and give permission.
                    """

   print(bowtieErrorMsg)
   return 1

def salmonErrMsg():
   salmonErrorMsg = """
                    IN CASE OF ERROR
                    Suggestions for most common errors:
                    Error: Status : FAIL 127 /bin/sh: ... ... ...salmon / build / quant: No such file or directory
                    Fix: Please check the path of the salmon executables.
                    Error: Permission denied.
                    Fix: Check the permission of the directory and give permission.
                    """

   print(salmonErrorMsg)
   return 1

def featCntErrMsg():
   featCntErrorMsg = """
                    IN CASE OF ERROR
                    Suggestions for most common errors:
                    Error: Status : FAIL 127 /bin/sh: ... ... ... featureCounts: No such file or directory
                    Fix: Please check the path of the featureCounts executables. 

                    """

   print(featCntErrorMsg)
   return 1

def fastqcErrMsg():
   fastqcErrorMsg = """
                    IN CASE OF ERROR
                    Suggestions for most common errors:
                    Error: Status : FAIL 127 /bin/sh: ... ... ... fastqc: No such file or directory
                    Fix: Please check the path of the fastqc executables.

                    """

   print(fastqcErrorMsg)
   return 1

def DESeqErrMsg():
   deseqErrorMsg = """
                   COMMON ERRORS:
                   Please chack library("DESeq2")  and library("ggplot2") are installed.

                   """
   print(deseqErrorMsg)
   return 1

def EdgeRErrMsg():
   edgerErrorMsg = """
                   COMMON ERRORS:
                   Please chack library("egdeR") is installed.

                   """
   print(edgerErrorMsg)
   return 1

def NOISeqErrMsg():
   noiseqErrorMsg = """
                   COMMON ERRORS:
                   Please chack library("egdeR") is installed.

                   """
   print(noiseqErrorMsg)
   return 1

def limmaErrMsg():
   limmaErrorMsg = """
                   COMMON ERRORS:
                   Please chack library("limma") is installed.

                   """
   print(limmaErrorMsg)
   return 1

def samtoolErrMsg():
   samtoolErrorMsg = '''

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

   '''
   print(samtoolErrorMsg)
   return 1

