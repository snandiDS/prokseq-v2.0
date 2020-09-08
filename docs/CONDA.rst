CONDA:
------
**Step 1:** Fetch the package.

``> conda install -p <PATH_TO_DOWNLOAD> -c <CHANNEL> prokseq``

Example:

``> mkdir testPrseq``

``> conda install -p /home/user/testPrseq -c snandids prokseq``

**Step 2:** Once the package is obtained, run the following commands.

``> tar -xvzf exampleFiles.tar.gz``

**Step 3:** Install dependencies.

``>tar -xvzf depend.tar.gz``

A depend folder will be created with all the required dependencies. Most of these packages were compiled on architecture x86_64 Fedora 31. Users may have to recompile some of them.

**Install the R and the R bioconductor packages.**
	   
Though the pipeline is written in Python3.6, but some packages included in the pipeline require Python2.7. Therefore, it is advised to install Python2. The program should find python2 and python (python3) in the env PATH. To make life easier, we recommend create environment (Step 4).

**Step 4:** Create virtual environment

``> conda create -n yourenvname python=3.6``

``> conda activate yourenvname``

``> conda install pandas``

``> pip2.7 install numpy``

``> pip2 install qc bitsets RSeQC``

``> pip2 install --upgrade cython bx-python pysam RSeQC numpy``

**Step 5:**

Once all the dependencies and R packages are installed, and the example files are untared, change the PATH in parameter file (Eg. param.input.bowtie). The PATH should point to the packages.

**For example:** If you are using the above-mentioned path `[/home/user/testPrseq]` from **Step 1**, then specify the path as below for all the packages in the parameter file.


Specify the path to pypy required for running afterqc

``PATH PYPY /home/path/testPrseq/depend/pypy2.7-v7.2.0-linux64/bin``

Then run the following command to test run the pipeline.

``> python scripts/pipeline-v2.8.py -s samples.bowtie.PEsample -p param.bowtie.yaml -n 4``

**Description:** The script is running with PE (paired-end) samples described in **samples.bowtie.PEsample**, and with the parameters defined in **param.bowtie.yaml**. The program is submitted with four processors.

**To remove:**

``> conda remove -p /home/path/testPrseq prokseq``


