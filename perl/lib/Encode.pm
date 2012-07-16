
=pod

=head1 NAME

Encode.pm - encodes nmer sequences

=head1 SYNOPSIS
 
  use lib "/home/martink/cvs/encode/current";
  use Encode;

  $s            = Encode::encode_nmer($nmer)
  $nmer_decoded = Encode::decode_nmer($s)

=head1 DESCRIPTION

Encodes sequence strings into a denser form using printable filesystem-compatible characters.

=head2 Encoding Scheme

The encoding scheme is not sensitive to the case of the input
sequence. The input will be convered to uppercase before encoding.

The encoding scheme associates a printable ASCII character with each 3mer.

 AAA !   AAC +   AAG ,   AAT 0   ACA 1   ACC 2   ACG 3   ACT 4   AGA 5   AGC 6 
 AGG 7   AGT 8   ATA 9   ATC :   ATG =   ATT @   CAA B   CAC D   CAG E   CAT F
 CCA H   CCC I   CCG J   CCT K   CGA L   CGC M   CGG O   CGT P   CTA Q   CTC R
 CTG S   CTT U   GAA V   GAC W   GAG X   GAT Y   GCA Z   GCC _   GCG a   GCT b
 GGA c   GGC d   GGG e   GGT f   GTA g   GTC h   GTG i   GTT j   TAA k   TAC l
 TAG m   TAT n   TCA o   TCC p   TCG q   TCT r   TGA s   TGC t   TGG u   TGT v
 TTA w   TTC x   TTG y   TTT z   NNN -

The characters A, C, G, T, N are excluded from the code vocabulary so
that sequences that are not multiples of 3mers can be encoded. Given a
sequence 

  s = a1+a2+a3+...+an+b

where (a1..an) are 3mers and b is the remaining 1mer or 2mer, the
encoded sequence is

  s' = code(a1)+code(a2)+code(a3)+...+code(an)+b

In other words, the trailing 1mer or 2mer is not encoded.

Some 9mer encodings

 TCTTGAAAT rs0
 TAAGGAGGG kce
 TATCTACTC nQR

Some 10mer encodings

 CATCGGCGTA FOPA
 TTCGTACCCT xgIT
 TACGGGCTAT leQT

Some 11mer encodings

 CGAAAGCTTGT L,UGT
 TTCAATAATAC x00AC
 CAGAACTGGCC E+uCC

You can immediately tell that the encoded string does not correspond
to a 3mer multiple because it contains one of ATGC.

=head2 Code Vocabulary

All code characters are printable and allowed in file and directory
names. Each code can for a single character file or directory. Certain
characters were excluded (e.g. . < > &) to allow for this and avoid
complications with special shell characters.

=head2 Encoding Gaps

Gaps in the input sequence (Ns) are not encoded. If the sequence
contains gaps then it is treated as a concatenate of encodable
sequences as shown above

  s_gappy = s1 + gap + s2 + gap +

  s_gappy_encoded = encode(s1) + gap + encode(s2) + gap + ...

If the gap is a 3mer, or larger, any NNN will be encoded by the
missing code (-), with remaining N or NN not encoded.

The presence of 1- or 2-base gaps makes encoding inefficient.

Below are some examples of encoded gappy nmers.

 CGANNNNNNAGTACNNAACNAGGCNG L--8ACNN+N7CNG 
 CNAANCGTNNNNNNTCTCTCNCNNCN CNAANP--rRNCNNCN 
 AGCGANNNNNNCNAATGTTGGCCCCC 6GA--CN0jdIC 

=head2 Missing Code

Any 3mers that contain characters other than ATGC and contain no gap
Ns, will be mapped to the missing code "^". The missing code decodes
into "???". This is not an invertable operation!

 GNTGNNNAQGNNNNNNNTCACTCGCN GNTG-^--NoRGCN 

Please make sure that the sequence data fed into the encoder doesn't
contain characters other than ATGCN, and if it does, make sure that
you understand what the encoding is doing.

=head2 Encoding Limitations

The encoding scheme is limited to the domain of printable, common
characters. Thus, while all 3mers of ATGC can be encoded (64
combinations), not all 3mers of ATGCN (125 combinations) can be
without resorting to the extended ASCII set.

=head1 HISTORY

=over

=item * 2 Feb 2009 v0.3

Adjusted encode and decode routines to use substr() for faster performance.

Removed temporary variables for faster performance.

=item * 6 Nov 2007 v0.2

Added support for sequences that are not multiples of 3mer. Altered
encoding to suit.

=item * 5 Nov 2007 v0.1

Started.

=back 

=head1 BUGS

=head1 AUTHOR

Martin Krzywinski

=head1 CONTACT

  Martin Krzywinski
  Genome Sciences Centre
  Vancouver BC Canada
  www.bcgsc.ca
  martink@bcgsc.ca

=cut

#######

package Encode;

use strict;
use vars qw(@seqchars $gap_code $missing_code $debug $huffman_lookup $missing_inverse @excludedchars $excludedrx);

$debug           = 0;
$huffman_lookup  = undef;
@seqchars        = qw(A C G T);
@excludedchars   = (@seqchars,"N");
$excludedrx      = join("",@excludedchars);
$gap_code        = "-";
$missing_code    = "^";
$missing_inverse = "???";

sub init_nmer {
  my @huffman_codes = grep($_ !~ /[$excludedrx]/ &&
			   $_ =~ /[\w,=+:@!]/i, map { chr($_) } (33..122));
  my $d = 3;
  for my $i (0..@seqchars**$d-1) {
    my $s;
    
    for my $j (0..$d-1) {
      my $ii = ( $i >> 2*$j ) & 3;
      $s = $seqchars[$ii] . $s;
    }
    my $code = defined $huffman_codes[$i] ? $huffman_codes[$i] : die "ran out of codes";
    print join(" ","lookuptable",$s,$code),"\n" if $debug;
    $huffman_lookup->{nmer}{$s} = $code;
    $huffman_lookup->{nmer}{$code} = $s;
  }
  $huffman_lookup->{nmer}{NNN} = $gap_code;
  $huffman_lookup->{nmer}{$gap_code} = "NNN";
  $huffman_lookup->{nmer}{$missing_code} = $missing_inverse;
  $huffman_lookup->{nmer}{$missing_inverse} = $missing_code;
}

sub encode_nmer {
  my $s = shift;
  my $s = uc $s;
  my $s_encoded;
  if(! defined $huffman_lookup) {
    init_nmer();
  }
  if($s =~ /n/i) {
    while($s =~ /(([^N]){1,3})|(N{1,3})/g) {
      my $string = $1 || $3;
      my $code;
      if(length($string) == 3) {
	$s_encoded .= defined $huffman_lookup->{nmer}{$string} ? $Encode::huffman_lookup->{nmer}{$string} : $missing_code;
      } else {
	$s_encoded .= $string;
      }
    }
  } else {
    for my $i ( 0 .. length($s)/3 ) {
      my $string = substr($s,$i*3,3);
      my $code;
      if(length($string) == 3) {
	$s_encoded .= defined $huffman_lookup->{nmer}{$string} ? $Encode::huffman_lookup->{nmer}{$string} : $missing_code;
      } else {
	$s_encoded .= $string;
      }
    }
  }
  return $s_encoded;
}

sub decode_nmer {
  my $s_encoded = shift;
  my $s;
  if(! defined $huffman_lookup) {
    init_nmer();
  }
  for my $i (0..length($s_encoded)-1) {
    my $code = substr($s_encoded,$i,1);
    $s .= defined $huffman_lookup->{nmer}{$code} ? $huffman_lookup->{nmer}{$code} : $code;
  }
  # old from v0.2
  #while($s_encoded =~ /(.)/g) {
  #  my $code = $1;
  #  my $inverse = defined $huffman_lookup->{nmer}{$1} ? $huffman_lookup->{nmer}{$1} : $code;
  #  $s .= $inverse;
  #}
  return $s;
}
