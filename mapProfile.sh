#!/bin/bash
###
# Description:
#   This script generates mappability profiles to be used in PICS and PING packages
#
# 16 July 2012
###


usage="
Usage: mapProfile.sh -l readLength -g genomeName [options]

Options:
-c : Specify a chromosome to download. If the entire genome is not needed.
-h : Print this help message and exit.

Examples:
./mapProfile.sh -l 50 -g hg19          #For the entire human genome with 50bp reads
./mapProfile.sh -l 36 -g sacCer2 -c I  #For the chromosome I of saccharomyces cerevisiae with 36bp reads
"

scriptdir=`dirname $0`
genome=
length=
chr="*" #default is all chr selected

function printHelp()
{
  echo "$usage"
}

while [ "$1" != "" ]
do
  case $1 in
    -l|--length) length=$2;shift;;
    -g|--genome) genome=$2;shift;;
    -c|--chr) chr=$2;shift;;
    -h|--h|-help|--help) printHelp; exit 1;;
    -*|--*) echo "Error! Wrong argument $1"; printHelp; exit 1;;
  esac
  shift
done

if [ "$genome" == "" ]
  then
  echo "Missing genome argument. Use -g <genomeName> to specify the genome."
  exit 1
fi
if [ "$length" == "" ]
  then
  echo "Missing length argument. Use -l <length> to specify the length of the reads."
  exit 1
fi

# Get the fasta for each chromosome
wget --timestamping "ftp://hgdownload.cse.ucsc.edu/goldenPath/${genome}/chromosomes/chr${chr}.fa.gz" #-P to specify DL directory

# Unzip
gunzip *fa.gz
# Merge
cat *.fa > $genome.fa

# bwa
bwa index -a bwtsw $genome.fa
samtools faidx $genome.fa

## pileup and bedFile creation
#if [ "$chr" == "*" ]
#then
#  for i in chr*.fa
#  do
#    $scriptdir/proMap.sh -l $length -g $genome.fa -c $i #& #run on every cores simultaneously
#  done
#  cat chr*.fa-$length/chr*.fa-$length.bed > $genome-$length.bed
#  #clean dir
#  rm -r chr*.fa-$length
#  rm -r $genome.fa.*
#  echo "Mappability profile for the genome $genome with reads of length $length: $genome-$length.bed"
#else
#  $scriptdir/proMap.sh -l $length -g $genome.fa -c chr$chr.fa
#  echo "Mappability profile for chr$chr of the genome $genome with reads of length $length: chr$chr.fa-$length/chr$chr.fa-$length.bed"
#fi
