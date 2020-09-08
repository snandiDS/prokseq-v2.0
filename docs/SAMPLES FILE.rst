SAMPLE FILE:
============
Sample file is required for the program. This file should have the following format. Please don't change the format of the file. Simply replace the fastq, sam and the conditions of the sample.

        
**File "sample" - sample description file.Specify the names of the sample files and tag them as "treat" and "control".**
     

# Specify the genome file, and specify the path where the indexed file will be stored, and the prifex of the indexed genome. Default is 'bowtie2_genome'.
        
GENOME SequenceChromosome.fasta bowtie2_genome/sequenceChr

-Specify the fastq files

-Specify the output name of the sam files.

-Followed by the tag/class/condition of the sample (treated or control)

**List all the fastq files as below.**

FASTQ sampleTreat_1.R1.fq sampleTreat_1.R2.fq sampleTreat_1.sam treat

FASTQ sampleTreat_2.R1.fq sampleTreat_2.R2.fq sampleTreat_2.sam treat

FASTQ sampleTreat_3.R1.fq sampleTreat_3.R2.fq sampleTreat_3.sam treat

FASTQ sampleCtrl_1.R1.fq sampleCtrl_1.R2.fq sampleCtrl_1.sam control

FASTQ sampleCtrl_2.R1.fq sampleCtrl_2.R2.fq sampleCtrl_2.sam control

FASTQ sampleCtrl_3.R1.fq sampleCtrl_3.R2.fq sampleCtrl_3.sam control

**End of file "sample"**
      

In general, this file starts with GENOME entry. As in the example, the genome file to be used in the analysis is **SequenceChromosome.fasta**. The next argument indicates what would be the prefix of the **indexed genome**, and where to store.<br/>

The following lines are FASTQ. The first argument in the entry is **forward fastq file**, and second is the **reverse fastq** file. **However, if one has only single-end reads, only one entry may be there**. 

The subsequent argument is the **SAM** file name. This argument specifies the **output** name of the sam files after running the bowtie. In this example, 
**sampleTreat_1.sam**, **sampleTreat_2.sam**, etc. are the sam files for **sampleTreat_1.R1/R2.fq**, **sampleTreat_2.R1/R2.fq**, etc. 

The last argument is the condition/class of the sample file (example: "**treat**" and "**control**").


In case of SALMON:
-----------------

        
**File "sample" - sample description file**

-Specify the names of the sample files and tag them as "**treat**" and "**control**".

-Specify the genome file and specify the path where the indexed file will be stored, and the prifex of the indexed genome. Default is 'transcripts_index'.
        
-**GENOME orf_coding_all.fasta transcripts_index**

-Specify the fastq files

-Specify the output name of the quant files.Followed by the tag/class/condition of the sample (**treated** or **control**)


-**List all the fastq files as bellow.**
        FASTQ sampleTreat_1.R1.fq sampleTreat_1.R2.fq sal_quant1 treat

        FASTQ sampleTreat_2.R1.fq sampleTreat_2.R2.fq sal_quant2 treat

        FASTQ sampleTreat_3.R1.fq sampleTreat_3.R2.fq sal_quant3 treat

        FASTQ sampleCtrl_1.R1.fq sampleCtrl_1.R2.fq sal_quant4 control

        FASTQ sampleCtrl_2.R1.fq sampleCtrl_2.R2.fq sal_quant5 control

        FASTQ sampleCtrl_3.R1.fq sampleCtrl_3.R2.fq sal_quant6 control

**End of file "sample"**
       

In general, this file starts with GENOME entry. As in the example, the transcript file to be used in the analysis is **transcripts.fasta**. The next argument indicates what would be the prefix of the **indexed** file, and where to store.<br/>'

The following lines are FASTQ. The first argument in the entry is **forward fastq** file and second is the **reverse fastq** file. However, if one has only single-end reads, only one entry may be there. 

The subsequent argument is the directory where the **alignment file has to be stored**. The last argument is the condition/class of the sample file (example: "**treat**" and "**control**").


DATA FILES:
===========
1. Samples files in fastq format.

2. Pathway analysis (Optional):

    1. TERN2NAME.csv :Gene Ontology to terms mapping csv file (Eg: GO:0000001,mitochondrion inheritance).This is the GO term classification which is common for all organisms.
  
    2. TERM2GENE.csv :Gene Ontology to gene mapping csv file (Eg: GO:0003688,YPK_0001). This is genome specific Gene ontology file. **TERM2GENE.csv** is a comma delimited 2 column file. First column is the **GO term** and second column is the **gene name**. 
          :User can download the GO file or GFF annotation file from Genome2D webserver
           (http://genome2d.molgenrug.nl/g2d_core_select_genbank.php)

    3. For Bowtie implementation:***Genome file in fasta format**

       For Salmon implementation: **Transcript file in fasta format**


    4. GTF file: Genome annotation in .gtf format. GTF file for bacteria can be downloaded from this link https://bacteria.ensembl.org/info/data/ftp/index.html. Or GFF3 annotation file can be converted to .gtf file by the gff2gtf.sh script provided in ProkSeq folder  

    5. BED file: BED file can be created by using gtf2bed.sh script which is also included in ProkSeq folder

**All these files should be declared in SAMPLE FILE and PARAMETER file.**
