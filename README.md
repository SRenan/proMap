Requirements:
-------------
You will need the following software installed: 
- samtools (http://samtools.sourceforge.net/)
- bwa (http://bio-bwa.sourceforge.net/bwa.shtml)

You will also need bash and perl, and the following perl packages installed on your machine:  
Config-General: http://search.cpan.org/~tlinden/Config-General-2.49/General.pm  
Math-VecStat: http://search.cpan.org/~aspinelli/Math-VecStat-0.08/VecStat.pm  
Set-IntSpan: http://search.cpan.org/~swmcd/Set-IntSpan-1.14/IntSpan.pm  
Statistics-Descriptive: http://search.cpan.org/~shlomif/Statistics-Descriptive-3.0200/lib/Statistics/Descriptive.pm

    perl -MCPAN -e shell
    install Config::General
    install Math::VecStat
    install Set::IntSpan
    install Statistics::Descriptive


Usage:
------
To make new mappability profiles, use mapProfile.sh which is a simple script that will launch the appropriate commands.

Here are a few steps for constructing your own mappability file:
We use the sacCer2 genome which is small enough for running a test. We select a reads length of 50bp.

    mkdir sacCer2
    cd sacCer2

For the entire genome:

    ../mapProfile.sh -l 50 -g sacCer2

For the chromosome I only:

    ../mapProfile.sh -l 50 -g sacCer2 -c I

At the end, you should have a file named sacCer.fa-50.bed containing all the non-mappable intervals.
You can then simply load it up in R using read.table or the import function of the rtracklayer package and use it in PICS or PING.


Notes:
------
- You will have to specify a genome name available on the UCSC ftp. ("ftp://hgdownload.cse.ucsc.edu/goldenPath/")
- The working directory is the current directory 
- By default, the script will launch chromosomes in parallel which requires multiple cpus/cores and a large amount of memory. If you do not have a large cpu/memory machine you could simply do it one chromosome at a time. 
- For a list of available options and examples, run ./mapProfile.sh -h
