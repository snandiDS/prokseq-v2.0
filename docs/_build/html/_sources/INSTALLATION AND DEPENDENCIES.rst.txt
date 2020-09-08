INSTALLATION AND DEPENDENCIES:
=============================================
**If one or any of the above dependencies are missing user can install it  by following the instructions below.**

Pyhton3:
--------
**Ubuntu**
Ubuntu 17.10, Ubuntu 18.04 (and above) come with Python 3.6 by default.
Ubuntu 16.10 and 17.04 do not come with Python 3.6 by default, but it is in the Universe repository.
You should be able to install it with the following commands:

        ``- sudo apt-get update``
        ``- sudo apt-get install python3.6``

For Ubuntu 14.04 or 16.04, Python 3.6 is not in the Universe repository, and user do not need to get
it from a Personal Package Archive (PPA). For example, to install Python from the “deadsnakes” PPA,
do the following:
       ``- sudo apt-get update``
       ``- sudo apt-get install python3.6``
**CentOS**
User should first update the system with the yum package manager:
        ``- sudo yum update``
        ``- sudo yum install yum-utils``

Then install the CentOS IUS package
        ``- sudo yum install https://centos7.iuscommunity.org/ius-release.rpm``

Then install Python and Pip:
        ``- sudo yum install python36u``
        ``- sudo yum install python36u-pip``

Installation of R:
------------------
The pipeline is tested on R version 3.6.0.
**Installing R on Ubuntu 19.04/18.04/16.04**
        
Before installing R, user need to update the system package index and upgrade all installed packages
using the following two commands:
        ``-sudo apt update``
        ``-sudo apt -y upgrade``
        
After that, run the following in the command line to install base R.
        ``-sudo apt -y install r-base``
        
**Install R on CentOS 7**
R packages are available in the EPEL repositories. It can be installed by typing:
        ``-sudo yum install epel-release``

Once the repository is added, install R by typing:
        ``-sudo yum install R``

Installation of R Bioconductor packages:
----------------------------------------
        if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")

        **edgeR:**
        ``BiocManager::install("edgeR")``

        **DESeq2:**
        ``BiocManager::install("DESeq2")``

        **NOISeq:**
        ``BiocManager::install("NOISeq")``

        **limma:**
        ``BiocManager::install("limma")``

        **clusterProfiler:**
        ``BiocManager::install("clusterProfiler")``

        **apeglm:**
        ``BiocManager::install("apeglm")``

        **RUVSeq:**
        ``BiocManager::install("RUVSeq")``

Samtools:
---------
**Ubuntu 18.04 or higher**
Install samtools by entering the following commands in the terminal:
        ``-sudo apt update``
        ``-sudo apt install samtools``
For the other version of Ubuntu or centose use the samtools.sh script in the package folder. User can
go to the PorkSeq folder and open a terminal and write sh samtools.sh. The program will install samtools
in the samtools directory.

EXTERNAL TOOLS:
---------------
**This program uses the following tools.**

        1.FastQC : This package runs the quality check

        2.Bowtie : Needed for aligning the reads

        3.Pypy : For speed and memory usage we sometime uses pypy an alternative implementation of python 3.6

        4.featureCounts from subread: a software program developed for counting reads to genomic features such
        as genes, exons, promoters and genomic bins.

        5.AfterQc: Tools for automatic filtering trimming of the fastq sequences.

To run the program the above mentions dependencies are essential. However, the executable binaries are bundled in the folder depend.
For ubuntu 18.04 or higher version use

``sudo apt-get update | apt-get install python3-pandas to install pandas``

User can download the depend folder along with all the dependencies from the following the Link: 
https://umeauniversity-my.sharepoint.com/:u:/g/personal/aakk0004_ad_umu_se/EZ6UF28lCcJGiuPOWQ8oVr0BtQAK1caGUEdVHuP29_I01g?e=o1K0mh

OR, follow the following:
        1. Create a folder named depend
        2. cd depend
        3. cp ../scripts/setup.sh .
        4. sh setup.sh
The depend foldr will be populated.
