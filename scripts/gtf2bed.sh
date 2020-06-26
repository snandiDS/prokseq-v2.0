#!/bin/sh
#
#    Copyright .
#
#    Author: Soumyadeep Nandi <snandi@ggn.amity.edu>
#	     Firoj ...
#
#
# This script transform the GTF file to a bed file.
# The GTF file shoule have following format.
# NC_0XXXXX       TEST    exon    37      1389    .       +       .       gene_id "XXX_0001"; transcript_id "TR:XXX_0001"
#
# Each entries should be tab delimited.
#
usage="
$(basename "$0") [-h] [-f fileName.gtf] -- script to convert GTF format to BED format

Options:
    -h  help message
    -f	GTF file name

Syntax:
sh gtf2bed.sh -f gtfFile_name > bedFile_name

Example:
sh gtf2bed.sh -f annotation.gtf > annotation.bed

Format:
GTF file:
NC_0XXXXX       TEST    exon    37      1389    .       +       .       gene_id "XXX_0001"; transcript_id "TR:XXX_0001"
BED file:
NC_0XXXXX       37      1389    XXX_0001        0       +       37      1389    0       1       1352,   0,
"

if [ $# -eq 0 ]; then
    echo "$usage" >&2
    echo ""
    echo "!!!!Please provide argument with option!!!!"
    echo ""
    exit 1
fi

while getopts ':hf:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) gtfF=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

#gtfF=$1
awk '{print $1"\t"$4"\t"$5"\t"$10"\t0\t"$7"\t"$4"\t"$5"\t0\t1\t"$5-$4",\t0,"}' $gtfF > _tempBED
sed -i -e 's/"//g' -e 's/;//g' _tempBED
cat _tempBED
rm _tempBED
