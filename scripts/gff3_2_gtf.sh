#!/bin/sh
#
#Copyright (c) Soumyadeep Nandi and Firoj Mahmud
#
#Permission is hereby granted, free of charge, to any person obtaining
#a copy of this software and associated documentation files (the
#"Software"), to deal in the Software without restriction, including
#without limitation the rights to use, copy, modify, merge, publish,
#distribute, sublicense, and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
#    Author: Soumyadeep Nandi <snandi@ggn.amity.edu>
#
#
# This script transform the GFF3 file to a GTF file.
# The GTF file shoule have following format.
# NC_0XXXXX       TEST    exon    37      1389    .       +       .       gene_id "XXX_0001"; transcript_id "TR:XXX_0001"
#
# Each entries should be tab delimited.
#
usage="
$(basename "$0") [-h] [-f fileName.gff] -- script to convert GFF format to GTF format

Options:
    -h  help message
    -f  GFF file name

Syntax:
sh gff3_2_gtf.sh -f gffFile_name > gtfFile_name

Example:
sh gff3_2_gtf.sh -f annotation.gff > annotation.gtf

Format:
GFF file:
##sequence-region NC_010465.1 1 4689441
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=502800
NC_010465.1     RefSeq  region  1       4689441 .       +       .       ID=NC_010465.1:1..4689441;Dbxref=taxon:502800;Is_circular=true;Name=ANONYMOUS;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=YPIII
NC_010465.1     RefSeq  gene    1       1389    .       +       .       ID=gene-YPK_RS00005;Name=dnaA;gbkey=Gene;gene=dnaA;gene_biotype=protein_coding;locus_tag=YPK_RS00005;old_locus_tag=YPK_0001
NC_010465.1     Protein Homology        CDS     1       1389    .       +       0       ID=cds-WP_002220732.1;Parent=gene-YPK_RS00005;Dbxref=Genbank:WP_002220732.1;Name=WP_002220732.1;gbkey=CDS;gene=dnaA;inference=COORDINATES: similar to AA sequence:RefSeq:YP_002348949.1;locus_tag=YPK_RS00005;product=chromosomal replication initiator protein DnaA;protein_id=WP_002220732.1;transl_table=11
NC_010465.1     RefSeq  gene    1394    2494    .       +       .       ID=gene-YPK_RS00010;Name=dnaN;gbkey=Gene;gene=dnaN;gene_biotype=protein_coding;locus_tag=YPK_RS00010;old_locus_tag=YPK_0002
NC_010465.1     Protein Homology        CDS     1394    2494    .       +       0       ID=cds-WP_002209645.1;Parent=gene-YPK_RS00010;Dbxref=Genbank:WP_002209645.1;Name=WP_002209645.1;gbkey=CDS;gene=dnaN;inference=COORDINATES: similar to AA sequence:RefSeq:WP_005976672.1;locus_tag=YPK_RS00010;product=DNA polymerase III subunit beta;protein_id=WP_002209645.1;transl_table=11

GTF file:
NC_010465.1     RefSeq  gene    1       1389    .       +       .       gene_id YPK_RS00005 Name        dnaA
NC_010465.1     RefSeq  gene    1394    2494    .       +       .       gene_id YPK_RS00010 Name        dnaN
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
    f) gffF=$OPTARG
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

#awk '{if($3=="gene") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' $gffF |sed -e 's/ID=gene-/gene_id\t/' -e 's/;Name=/ Name\t/' -e 's/;.*//' > _tempGTF
#awk '{if($3=="gene") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' $gffF |sed -e 's/ID=/gene_id\t/' -e 's/;Name=/ Name\t/' -e 's/;.*//' > _tempGTF
awk '{if($3=="gene") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' $gffF |sed -e 's/ID=/gene_id\t/' -e 's/;Name=/ Name\t/' -e 's/;.*//' |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$12"\t"$11"\t"$10}'> _tempGTF
#sed -i -e 's/"//g' -e 's/;//g' _tempGTF
cat _tempGTF
rm _tempGTF

