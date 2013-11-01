#!/bin/bash

script_dir=`dirname $0`
perl_dir=$script_dir/perl

##################################
## Command line
# proMap -l 36 -g "hg18.fa" -c "chr1"
# proMap --length 36 --genome "hg18.fasta" --chr "chr1"
# -l or --length: sequence length
# -g or --genome: genome name
# -c or --chr: chromosome to be mapped
### Hard Coded parameters

# freq is hard coded at 0.05 in awk command

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
    else 
      if [ ${args[$i-1]} = '-c' ] || [ ${args[$i-1]} = '--chr' ]; then
        chr=${args[$i]}
        let index++
        let i++
      fi
    fi
  fi
done

if [ $index -ne 3 ]
  then
  echo "Error, usage should be of the form: `basename $0` -l length -g genome -c chr"
  exit 1
else
  echo "Starting mapping for $chr on genome $genome with mers of length $length ."
fi

mkdir $chr-$length
perl "$perl_dir/enumeratenmer" -nmer $length -skip_gappy -fasta=$chr -raw -fastq > "$chr-$length/$chr-$length+.txt" 2> "$chr-$length/$chr-$length-.txt"
bwa aln $genome "$chr-$length/$chr-$length+.txt" > "$chr-$length/$chr-$length.sai"
bwa samse $genome "$chr-$length/$chr-$length.sai" "$chr-$length/$chr-$length+.txt" > "$chr-$length/$chr-$length.sam"
##cleanup
#rm $chr-$length/$chr-$length*.txt
#rm $chr-$length/$chr-$length*.sai
grep "XT:A:U" "$chr-$length/$chr-$length.sam" > "$chr-$length/$chr-$length-U.sam"
##cleanup
#rm $chr-$length/$chr-$length.sam
samtools view -b -S -t "$genome.fai" -o "$chr-$length/$chr-$length.bam" "$chr-$length/$chr-$length-U.sam"
##cleanup
#rm $chr-$length/$chr-$length-U.sam
samtools sort "$chr-$length/$chr-$length.bam" "$chr-$length/$chr-$length-sorted"
##rm $chr-$length/$chr-$length.bam
samtools mpileup -f $genome "$chr-$length/$chr-$length-sorted.bam" > "$chr-$length/$chr-$length.pileup"
##cleanup
#rm $chr-$length/$chr-$length*.b*
echo | awk -v var=${length} '{if($4>=0.05*var+1){print$1"\t"$2"\t"$4}}' "$chr-$length/$chr-$length.pileup" > "$chr-$length/$chr-$length.map"
## cleanup
#rm $chr-$length/$chr-$length.pileup
perl "$perl_dir/bedcompact.pl" "$chr-$length/$chr-$length.map" "$chr-$length/$chr-$length.map2"
## cleanup
#rm $chr-$length/$chr-$length.map
awk '{if($4==0){print$1"\t"$2"\t"$3"\t"$4}}' "$chr-$length/$chr-$length.map2" > "$chr-$length/$chr-$length.bed"
## cleanup
#rm $chr-$length/$chr-$length.map2

echo "Mappability file: $chr-$length/$chr-$length.bed"
echo "Done!"

