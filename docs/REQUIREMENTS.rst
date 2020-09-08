REQUIRMENTS:
============
Users can run the ProkSeq program to see which depending packages are missing. By default, ProkSeq will search for required programs and 
print the availability of the program.

ProkSeq requires the following packages:

EXTERNAL TOOLS:
---------------
**Package/program   : Purpose**

**FastQC**          : Quality check

**Bowtie**          : Aligning the reads

**Pypy**            : For speed and memory usage we sometimes uses pypy an alternative implementation of python 3.6

**featureCounts**   : Counting reads to genomic features such as genes, exons, promoters, and genomic bins.


**AfterQc**         : Automatic filtering trimming of the fastq sequences.

**Samtools**        : For post-processing of the SAM and BAM reads alignment files.

**Salmon**          : A tool for wicked-fast transcript quantification from RNA-seq data.

The dependencies mentioned above are essential. However, the executable binaries
are bundled in the folder "depend". If the user is fetching the package from

github [https://github.com/snandiDS/prokseq-v2.1], then the user will get a

script [setup.sh] inside the depend folder.

**Please run this script as**

        ``sh setup.sh``

This script will fetch the required dependencies from http://www.fallmanlab.org.
The script will also ask if the user wants to compile samtools.
After running the script and accepting Y for compiling the samtools, the output
on the screen would be as follows:

        **-- Viewing**

           flags         : explain BAM flags

           tview         : text alignment viewer

           view          : SAM<->BAM<->CRAM conversion

           depad         : convert padded BAM to unpadded BAM

**This means the program ran successfully.**

R packages:
-----------

**R packages**  :**Purpose**

    
ggplot2         : A system for declaratively creating graphics, based on The Grammar of Graphics.


**Bioconductor Packages:  Purpose**
        
DESeq2          :       Differential gene expression analysis based on the negative binomial distribution.

edgeR           :       Package for examining differential expression of replicated count data.

NOISeq          :       A non-parametric approach for the differential expression analysis of RNseq-data.

limma           :       A package for the analysis of gene expression data arising from microarray or RNA-seq technologies.

clusterProfiler :       To analyze functional profiles of genomic coordinates (supported by ChIPseeker), gene and gene clusters.

apeglm          :       The adaptive t prior shrinkage estimator used to Shrink log2 fold changes.

RUVSeq          :       Remove Unwanted Variation from RNA-Seq Data

RColorBrewer    :       Required to create nice looking color palettes especially for thematic maps, used with RUVSeq.


PYTHON LIBRARIES:
----------------
**This program is written in python 3.6, and uses the following python libraries. os, subprocess, re, pandas, optparse, math, threading. These libraries may be installed using pip.**

