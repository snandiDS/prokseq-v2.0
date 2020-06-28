#!/bin/sh

# Author : Soumyadeep Nandi
# The script will download the packages from http://www.fallmanlab.org 
# to the local system and store them in the folder "depend".
# The executable binaries are bundled in this folder.
# By default this folder should have the setup.sh script.

echo "Running depend setup"
samtools(){
  current=$(pwd)
  echo $current
  dest="${current}/samtools"
  echo $dest
  mkdir $dest
  cd samtools-1.10
  chmod 755 config.status
  chmod 755 configure
  chmod 755 install-sh
  chmod 755 version.sh
  make clean
  ./configure --prefix=$dest
  make
  make install
  cd $current
  ./samtools/bin/samtools
}
echo "Fetching bowtie2"
mkdir bowtie2
cd bowtie2
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/bowtie2-2.3.5.1-linux-x86_64.zip
unzip bowtie2-2.3.5.1-linux-x86_64.zip
cd ../
echo "Fetching SALMON"
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/salmon-1.1.0_linux_x86_64.tar.gz
tar -xvzf salmon-1.1.0_linux_x86_64.tar.gz
echo "Fetching pypy"
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/python2.7.zip
unzip python2.7.zip
mv python2.7/pypy2.7-v7.2.0-linux64.tar .
tar -xvf pypy2.7-v7.2.0-linux64.tar
rm -rf python2.7
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/pyhton3.6.zip
unzip pyhton3.6.zip
mv pyhton3.6/pypy3.6-v7.2.0-linux64.tar .
tar -xvf pypy3.6-v7.2.0-linux64.tar
rm -rf pyhton3.6
echo "Fetching RSeQC"
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/RSeQC-2.6.2.zip
unzip RSeQC-2.6.2.zip 
echo "Fetching subread"
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/subread-1.4.6-p5-Linux-i386.tar.gz
tar -xvzf subread-1.4.6-p5-Linux-i386.tar.gz
echo "Fetching samtools-1.10"
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/samtools-1.10.zip
unzip samtools-1.10.zip
echo "Fetching afterQC"
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/afterqc.zip
unzip afterqc.zip
echo "Fetching FastQC"
wget http://www.fallmanlab.org/wp-content/uploads/2020/06/FastQC.zip
unzip FastQC.zip
while true; do
   read -p "Do you wish to compile samtools [Y/N]?" yn
   case $yn in
      [Yy]* ) samtools; break;;
      [Nn]* ) exit;;
      * ) echo "Please answer [Y] yes or [N] no.";;
   esac
done

