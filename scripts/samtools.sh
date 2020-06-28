#!/bin/sh

# Author : Soumyadeep Nandi
# The script will compile the samtools version 1.10 in the depend
# folder. The script will create a folder "samtools" to store the
# compiled programs.

current=$(pwd)
echo $current
dest="${current}/depend/samtools"
echo $dest
mkdir $dest
cd depend/samtools-1.10
make clean
./configure --prefix=$dest
make
make install
cd $current
./depend/samtools/bin/samtools
