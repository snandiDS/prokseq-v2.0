#!/usr/bin/python
#
#Copyright (c) Soumyadeep Nandi and Firoj Mahmud
#
#Permission is hereby granted, free of charge, to any person obtaining
#a copy of this software and associated documentation files (the
#"Software"), to deal in the Software without restriction, including
#without limitation the rights to use, copy, modify, merge, publish,
#distribute, sublicense, and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
#    Author: Soumyadeep Nandi <snandi@ggn.amity.edu>
#            Firoj Mahmud <firoj.mahmud@umu.se>
#
#
# This script transforms the gene bank file from embl from
# ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/embl/ to a comma seperated file.
# The output file will have the following format.
# GeneOntology,GeneName
#
# Options:
#     -h  help message
#     -i  inputfile name
#     -o  outputfile name
# 
# Syntax:
# python parseEMBLgeneFile.py -i inputFile_name -o outputFile_name
# 

import re, sys, getopt

def main(argv):
   inile = ''
   outfile = ''
   gene = '' 
   go = '' 
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifl=","ofl="])
   except getopt.GetoptError:
      print('parseEMBLgeneFile.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   if(len(opts)<2):
      print('parseEMBLgeneFile.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('parseEMBLgeneFile.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifl"):
         infile = arg
      elif opt in ("-o", "--ofl"):
         outfile = arg
      else:
         print('parseEMBLgeneFile.py -i <inputfile> -o <outputfile>')
         sys.exit()
   fl = open(infile, 'r') 
   fo = open(outfile, 'w') 
   while True: 
      line = fl.readline() 
      if(re.match("^FT\s+/gene=\"\S+",line)):
         LINE = re.findall("^FT\s+/gene=\"(\S+)\"",line)
         gene=LINE[0]
      if(re.match("^FT\s+/db_xref=\"GO:\S+",line)):
         LINE = re.findall("^FT\s+/db_xref=\"(GO:\S+)\"",line)
         go=LINE[0]
         fo.write(str(go)+","+str(gene)+"\n")
      if not line: 
           break
   fl.close()
   fo.close()

if __name__ == "__main__":
   main(sys.argv[1:])
