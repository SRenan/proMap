#!/bin/bash

##################################
## Command line
# proMapAll -l 36 -g "hg18.fa"
# proMapAll --length 36 --genome "hg18.fasta"
# -l or --length: sequence length
# -g or --genome: genome name
###################################
## Before running this command
## You need to index your genome with the following two commands
# bwa index -a bwtsw genome
# samtools faidx genome
###################################

scriptdir=`dirname $0`
args=("$@") # Input arguments
nargs="${#args[@]}" # Number of arguments

index=0          # Initialize count.
for (( i=1; i<=$nargs; i++ ))
do
  if [ ${args[$i-1]} = '-l' ] || [ ${args[$i-1]} = '--length' ]; then #  || "$arg[i]"=="--length"
    length=${args[$i]}
    let index++
    let i++
  else
    if [ ${args[$i-1]} = '-g' ] || [ ${args[$i-1]} = '--genome' ]; then 
      genome=${args[$i]}
      let index++
      let i++
    fi
  fi
done

if [ $index -ne 2 ]
  then
  echo "Error, usage should be of the form: `basename $0` -l length -g genome"
  exit 1
else
  echo "Starting mapping for all chromosomes on genome $genome with mers of length $length."
fi

for i in chr*.fa
do
   $scriptdir/proMap.sh --length $length --genome $genome --chr $i &
done

# For large genomes (e.g. human), you may want to split the chromosomes, see below
# for i in chr[1]*.fa
# do
#    proMap --length $length --genome $genome --chr $i &
# done
# 
# wait
# 
# for i in chr[!1]*.fa
# do
#    proMap --length $length --genome $genome --chr $i &
# done

wait

for chr in chr*.fa
do
  mv "$chr-$length/$chr-$length.bed" ./
  # mv "$chr-$length/$chr-$length.map2" ./
  rm -r "$chr-$length"
done

cat chr*-$length.bed > $genome-$length.bed
# cat chr*-$length.map2 > $genome-$length.bed2
rm chr*-$length.bed
# rm chr*-$length.map2
