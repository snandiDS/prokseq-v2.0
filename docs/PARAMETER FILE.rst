PARAMETER FILE:
===============
**There should be one parameter file. The entries of the file should be as follows.**

        
**File "param.yaml"** - definition file in YAML format. 

Define the paths and parameters to run the external packages. All the default parameters are considered. However, if users wish to modify any default parameters, please specify here.

    NOTE: Any flag/parameter which do not require any value, specify **"TRUE"/"FALSE."** Quotes are necessary. "**TRUE**" to invoke the flag and "**FALSE**" to suppress.
        

**Describe the Bowtie options below as:**
        
BOWTIE:
           	``-I : 0,`` 

           	``-X : 500,``

           	``-k : 1,``

           	``-p : 1,``

**In case the package salmon is to be run, uncomment(remove #) the following options. These are the default values. These options will overwrite the BOWTIE parameters.**
        
**SALMON:**

        	``SALMONINDEX:-k : 29``

        	``SALMONQUANT:``

        	``-l : A-p : 2``

        	``--validateMappings : "TRUE"``

        

**Define the parameters for AfterQC. Example default parameters are shown as bellow**. All the default parameters are taken. Specify the parameter and value in case of any changes. One cana add the other parameters with value.

**AFTERQC:**
        	``-s : 35``

        	``-n : 5``

        	``-p : 35``

        	``--debubble : "TRUE"``

        	``-a : 5``



        
**Define the Featurecounts options as bellow:**

**FEATURECOUNTS:**

           ``a : oldAnnotationGFF.gtf`` #Define the Featurecounts input GTF file.

           ``o : FeatCount`` #Define the Featurecounts output file name.




**Specify a count file name**

COUNTFILE : **countFile.csv**

        

**geneBody coverage.r require a bed file. Specify the name of bed file as bellow:**

geneBody_coverage:
           ``r : oldAnnotationGFF.bed``
        


**Specify if batch effect removal is required. FALSE if not required, else TRUE.**

BATCH_EFFECT_REMOVE : "**FALSE**"



        
Remove Unwanted Variation from RNA-Seq Data
-------------------------------------------   
**To run RUVSeq, uncomment the following two lines.**

RUVSeq:
           ``nc : 100`` #Integer - Top 100 genes as ranked by edgeR p-values. Negative control genes to estimate the factors of unwanted variation.
        



Parameters for pathway analysis.
-------------------------------

**For pathway analysis, define the organism in three alphabets as bellow.**

ypy = *Yersinia pseudotuberculosis*

-**User need to change the keg abbreviation of their genome which can be found in** https://www.genome.jp/kegg/catalog/org_list.html. Here ypy is the *Yersinia pseudotuberculosis* YPIII



**If PATHWAY analysis is not required comment out (add # infront) the PATHWAY section.**

PATHWAY:
           ``cutoffPositive : 2.0``  #logfold upper limit
           ``cutoffNegative : -2.0`` #logfold lower limit
           ``Organism : ypy``


**For Gene Ontology of the pathway analysis, define GO term and gene name file.**

TERM2GENE : data/TERM2GENE.csv

TERM2NAME : data/TERM2NAME.csv



Select the path of tools
------------------------

**Specify the root path. That means where the ProkSeq bundle is unpacked.**

The PATH ROOT : should indicate the current working directory.

**The location should have the following folders**

      1. depend - contains all the binaries of the external packages

      2. scripts - contains all the modules required for ProkSeq, and prokseq-vx.x.sh

      3. data - Contains Gene ontology files for pathway analysis.

      4. Fastq files, and other genome/transcript files.


**If the package is stored in the folder/path /home/user/PROKSEK**

PATH ROOT : /home/user/PROKSEK


**If the above environment (depend, scripts, data) is true, the following line may be uncommented.**
        
PATH DEFAULT : "TRUE"

PATH SAMTOOLS : /home/user/PROKSEK/depend/samtools/bin #Specify the path to samtools
       
PATH geneBody_coverage : /home/user/PROKSEK/depend/RSeQC-2.6.2/scripts/  # Specify the path to geneBody_coverage

PATH FEATURECOUNTS : /home/user/PROKSEK/depend/subread-1.4.6-p5-Linux-i386/bin/  #Specify the path to FEATURECOUNTS
        
PATH FASTQC : /home/user/PROKSEK/depend/FastQC # Specify the path to fastqc 
        
PATH BOWTIE : /home/user/PROKSEK/depend/bowtie2/bowtie2-2.3.5.1-linux-x86_64 #Specify the path to bowtie
        
PATH SALMON : /home/user/PROKSEK/depend/salmon-latest_linux_x86_64/bin  #Specify the path to salmon if salmon is required

PATH PYPY : /home/user/PROKSEK/depend/pypy2.7-v7.2.0-linux64/bin #Specify the path to pypy required for running afterqc
        
PATH READFASTA : /home/user/PROKSEK/depend  #Specify the path to readfasta
        
**End of file "param.input"**
        


In general, the entries starting with BOWTIE instructs the program to run the tool with the additional parameter. Similarly SALMON, AFTERQC, FEATURECOUNTS, etc. Entries beginning with PATH indicates the path to the executables of the external tools.  If one uses the bundled packages in the depend folder, `PATH DEFAULT : TRUE` line should be uncommented. 

**Note: PATH ROOT : path_to_the_current_working_directory should be mentioned.**


**NOTE: The program will override the Bowtie options, and the package salmon will be run if both bowtie and salmon's options are provided. Therefore, if salmon is required, please comment on the Bowtie option lines, and vice versa.**


