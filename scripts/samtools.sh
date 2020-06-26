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
