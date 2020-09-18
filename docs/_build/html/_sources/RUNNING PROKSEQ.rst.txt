RUNNING PROKSEQ:
===============

SYNTAX:
-------
        Usage: prokseq.py [options] arg

        Options:
          ``-h, --help``        
                                **show this help message and exit**
          ``-s SAMPLE_FILE_NAME, --sample=SAMPLE_FILE_NAME``
                                **provide the sample file**
          ``-p PARAMETER_FILE_NAME, --param=PARAMETER_FILE_NAME``
                                **provide the parameter file**
          ``-n NUMBER OF PROCESSORS, --numproc=NUMBER OF PROCESSORS``
                                **provide the number of processors**
          ``-v, --version``         
                                **Version of the package**


EXAMPLE:
--------
        ``python scripts/prokseq.py -s samples.bowtie.PEsample -p param.bowtie.yaml -n 4``


The program will run with sample file **"samples.bowtie.PEsample"**, and parameter file **"param.bowtie.yaml"**. The program will also utilize **4 processors**.

To run the program, the dependencies mentioned above are essential. 

However, the executable binaries are bundled in the folder **"depend"**. 


**An example parameter (param.bowtie.yaml) and sample file (samples.bowtie.PEsample) are bundled together with the package in the exampleFiles.tar.gz.**



Default directory layout should look like below:
------------------------------------------------
        
        ./README.md

        ./depend/afterqc

        ./depend/bowtie2

        ./depend/FastQC

        ./depend/pypy2.7-v7.2.0-linux64

        ./depend/readFasta

        ./depend/RSeQC-2.6.2

        ./depend/salmon-latest_linux_x86_64

        ./depend/samtools

        ./depend/samtools-1.10

        ./depend/subread-1.4.6-p5-Linux-i386

        ./depend/wigToBigWig

        ./scripts/prokseq.py

        ./scripts/runCheck.py

        ./scripts/prokseq-vT1.py

        ./scripts/plotScript.R

        ./scripts/gff3_2_gtf.sh

        ./scripts/gtf2bed.sh

        ./scripts/setup.sh

        ./scripts/samtools.sh

        ./scripts/libmod

        ./scripts/libmod/checkEnv.py

        ./scripts/libmod/errMsgFn.py

        ./scripts/libmod/execCmd.py

        ./scripts/libmod/__init__.py

        ./scripts/libmod/pipeFunc.py

        ./scripts/libmod/__pycache__

Example files layout:
---------------------
        ./sampleCtrl_1.R1.fq

        ./sampleCtrl_1.R2.fq

        ./sampleCtrl_2.R1.fq

        ./sampleCtrl_2.R2.fq

        ./sampleCtrl_3.R1.fq

        ./sampleCtrl_3.R2.fq

        ./sampleTreat_1.R1.fq

        ./sampleTreat_1.R2.fq

        ./sampleTreat_2.R1.fq

        ./sampleTreat_2.R2.fq

        ./sampleTreat_3.R1.fq

        ./sampleTreat_3.R2.fq

        ./samples.bowtie.PEsample

        ./samples.bowtie.SEsample

        ./samples.salmon.PEsample

        ./samples.salmon.SEsample

        ./param.bowtie.yaml

        ./param.salmon.yaml

        ./oldAnnotationGFF.bed

        ./oldAnnotationGFF.gtf

        ./orf_coding_all.fasta

        ./SequenceChromosome.fasta


Check test run:
---------------

After setting up of the depend directory, one can check if the environment is all setup. The required fastq and other required files for the check run is bundled in **exampleFiles.tar.gz.** Therefore, this file should be untared or at least should be there in the current working directory. The check script can be run as follows.

``> python scripts/runCheck.py``

**Note: This will require the files from exampleFiles.tar.gz.**




