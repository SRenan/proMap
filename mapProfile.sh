#!/bin/bash
###
# Description:
#   This script generates mappability profiles to be used in PICS and PING packages
#   ref: 
#
# 06 July 2012
###


usage="
Usage: mapProfile.sh -l readLength -g genomeName [options]

Options:
-d : Directory where the files will be downloaded and extracted.
-c : Specify a chromosome to download. If the entire genome is not needed.
-h : Print this help message and exit.
"

scriptdir=`dirname $0`
workdir=./
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
    -d|--dir) workdir=$2;shift;;
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
wget -P $workdir  --timestamping "ftp://hgdownload.cse.ucsc.edu/goldenPath/${genome}/chromosomes/chr${chr}.fa.gz" #-P to specify DL directory

# Unzip
gunzip $workdir/*fa.gz
# Merge
cat *.fa > $genome.fa

# bwa
bwa index -a bwtsw $genome.fa
samtools faidx $genome.fa

# pileup and bedFile creation
#python pilup.py 
$scriptdir/proMapAll.sh -l $length -g $genome.fa
if [ "$chr" == "*" ]
then
  echo "Mappability profile for the genome $genome with reads of length $length: $genome-$length.bed"
else
  echo "Mappability profile for chr$chr of the genome $genome with reads of length $length: $genome-$length.bed"
fi
