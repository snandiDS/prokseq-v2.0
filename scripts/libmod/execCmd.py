import re
import os
import subprocess
from libmod import errMsgFn as errMsg
import time

def print_cube(num):
    """
    function to print cube of given num
    """
    time.sleep(2)
    print("Cube: {}".format(num * num * num))

def print_square(num):
    """
    function to print square of given num
    """
    time.sleep(6)
    print("Square: {}".format(num * num))

def executeCmds (cmd,mod):
   try:
      output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, universal_newlines=True)
   except subprocess.CalledProcessError as exc:
      print("Status : FAIL", exc.returncode, exc.output)
      if(re.search("bowtie",mod)): errMsg.bowtieErrMsg()
      if(re.search("fastqc",mod)): errMsg.fastqcErrMsg()
      if(re.search("afterqc",mod)): errMsg.pypyErrMsg()
      if(re.search("featcnt",mod)): errMsg.featCntErrMsg()
      if(re.search("DESeq",mod)): errMsg.DESeqErrMsg()
      if(re.search("edgeR",mod)): errMsg.EdgeRErrMsg()
      if(re.search("noiseq",mod)): errMsg.DESeqErrMsg()
      if(re.search("salmon",mod)): errMsg.salmonErrMsg()
      if(re.search("limma",mod)): errMsg.limmaErrMsg()
      exit()
   else:
      print("Output: \n{}\n Success.\n".format(output))
      #os.system(cmd)
