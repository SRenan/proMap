#!/usr/bin/perl


use strict;

my $ligne;
my $ligne2;
my $flagScore=0;
my $increment=0;
my $scoreTemp;
my $startTemp;
my $flagregroupe=0;
my $endTemp;
my $chrTemp;
my $flag=0;
my $starPosition;
my $flagflag=0;
my $nombreLigne=0;

my $file=$ARGV[0];
my $out_file=$ARGV[1];

open (IN,$file);
open (out,">$out_file");


open (INN,$file);

while($ligne2=<INN>)
{
  $nombreLigne++;
}

while($ligne=<IN>)
{
  $ligne=~ /^(\S+)\s+(\S+)\s+(\S+)/;

  if($flagScore==0)
  {
    $scoreTemp=$3;
    $startTemp=$2;
    $flagScore=1;
    $endTemp=$2;
    $chrTemp=$1;
    $starPosition=$2
  }


  #print "****************************\n";
  #print "valeur temporaire : $scoreTemp\n";
  #print "valeur  du score : $3\n";
  #print "valeur  du start :$startTemp\n";
  #print "valeur  du end : $endTemp\n";
  #print "****************************\n";

  if(($scoreTemp!=$3)&&($increment>0))
  {	
    print out "$chrTemp\t$startTemp\t$endTemp\t$scoreTemp\n";
    $scoreTemp=$3;
    $startTemp=$2;
    $flagregroupe=0;
    #sleep(1);
  }

  #print "$2\t $starPosition\n";
  if($2!=($starPosition+1)&&($increment>0))
  {
    if(($endTemp-$startTemp)==0)
    {
      print out "$chrTemp\t$startTemp\t$endTemp\t$scoreTemp\n";
    }

    my $valeurstart=$starPosition+1;
    my $valeurend=$2-1;
    #print "Probleme**********************\n\n\n\n";	
    print out "$chrTemp\t$valeurstart\t$valeurend\t0\n";
    $startTemp=$2;
    #sleep(2);
  }

  if(($startTemp!=1)&&($increment==0))
  {
    my $temppo=($2-1);
    print out "$chrTemp\t1\t$temppo\t0\n";
  }

  if (($nombreLigne==($increment+1)))
  {
    # print "rentre";
    print out "$chrTemp\t$startTemp\t$2\t$scoreTemp\n";
  }

  $increment++;
  $endTemp=$2;
  $chrTemp=$1;
  $starPosition=$2;
}


close(out);
close(IN);