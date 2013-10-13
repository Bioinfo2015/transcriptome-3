use warnings;
use strict;
use List::Util qw[min max];
use Getopt::Std;
use File::Basename;


die(qq/
Mining_Transcripts.pl [options] 

external dependencies: exonerate

options:
-f  FILE  Annotated transcripts
-b  FILE  combined reference protein database
-c  INT   number of bp to keep outside the orf [500]
-o  CHAR  prefix for outfile

\n\n/) unless (@ARGV);


my %opts = (f=>undef, b=>undef, c=>500, o=>undef);
getopts('f:b:c:o:', \%opts);


my $trans = $opts{f};
my $ref = $opts{b};
my $cut = $opts{c};
my $out = dirname($trans) . "/" . $opts{o} . ".txt";

open (PRO, "<", $ref);
my $id;
my %refpro;
my ($pos,$pos2);

while (<PRO>) {
  chomp (my $line = $_);
  if ($line =~ m/^>(\S+)/) {    
    $id = $1; 
    $pos = tell (PRO);
    
  }
  unless ($line =~ m/^>(\S+)/) {
    push @{$refpro{$id}}, $pos;
    $pos2 = tell (PRO);
    $pos = $pos2;    
  }    
}	

my $ID;
my %seq;
open (SEQ, "<", $trans);
while (<SEQ>) {
  chomp (my $line =$_);
  if ($line =~ m/^>(\S+)/) {
    $ID = $1;
    my @d = split /\s+/, $line;
    my $prot = $d[2];
    my $pos = tell (SEQ);
    chomp (my $seq = <SEQ>);
    $seq{$prot} = {'pos' => $pos, 'contig' => $ID, 'rest' => join ("\t", @d[2..$#d])};
  }  
} 

my %new_seq;
foreach my $contig (sort {$a cmp $b} keys %seq) {
  if ($refpro{$contig}) {
    my $e;
    my $p;
    my $s;
    for (my $i = 0 ; $i< scalar @{$refpro{$contig}}; $i++) {
      seek PRO, $refpro{$contig}[$i], 0;
      chomp ($e = <PRO>);
      $p .= $e;
    }
    seek SEQ, $seq{$contig}{'pos'}, 0; 
    chomp ($s = <SEQ>);

    
    my $target = "target.fa";
    my $query = "query.fa";
  
    open (T, ">", $target);
    open (Q, ">", $query);
    
    print Q ">seq","\n", $s, "\n";
    print T ">protein", "\n", $p, "\n";  
    my $contiglength = length($s);
    
    my @exon = `exonerate $target $query -m protein2genome --showalignment no --showcigar 0`;
    unlink($query); unlink($target);
    my @match;
    close T; close Q;
  my $start;
  my $end;
  foreach (@exon) {
    chomp (my @line = split /\s+/, $_);
    if ($line[0] =~m /vulgar/) {
      if ($line [8] eq '-') {
	print "WTF?? still has errors at this point?? I am gonna quit science and I want to be a super hero like spiderman instead!", "\n";
	print "the error was ", $line[5], "\n";
      }
      push @match, $line[6]+1;
      push @match, $line[7]+1;
    } # if ($line[0] =~m /vulgar/) {
  }
    $start = min (@match);
    $end = max (@match);
    
  if ($contiglength - $end <= $cut ) {
    if ($start <= $cut) {
      $new_seq{$contig} = {'seq' => $s, 'contig' => $seq{$contig}{'contig'}, 'rest' =>$seq{$contig}{'rest'}};
    }
    if ($start > $cut) {
      my $seq2 = substr ($s, $start-$cut-1);
      $new_seq{$contig} = {'seq' => $seq2, 'contig' => $seq{$contig}{'contig'}, 'rest' =>$seq{$contig}{'rest'}};
    }
  }
  if ($contiglength - $end > $cut) {
    if ($start <= $cut) {
      my $seq2 = substr ($s, 0, $end+$cut);
      $new_seq{$contig} = {'seq' => $seq2, 'contig' => $seq{$contig}{'contig'}, 'rest' =>$seq{$contig}{'rest'}};
    }
    if ($start > $cut) {
      my $seq2 = substr ($s, $start-$cut-1, $end-$start+$cut+$cut+1);
      $new_seq{$contig} = {'seq' => $seq2, 'contig' => $seq{$contig}{'contig'}, 'rest' =>$seq{$contig}{'rest'}};
    }
  }
    
    
  }  
  close SEQ; close PRO;
}
open (OUT, ">", $out);
foreach my $name ( sort{$a cmp $b} keys %new_seq) {
  print OUT ">", $new_seq{$name}{'contig'}, "\t", $new_seq{$name}{'rest'},"\n",  $new_seq{$name}{'seq'},"\n";
}
close OUT;
