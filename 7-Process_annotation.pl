#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use Getopt::Std;
use Getopt::Long;
#use File::Temp;
use List::Util qw( min max );

#Ke Bi (kebi@berkeley.edu) Sept. 3 2013



&main;
exit;

sub main {
        &usage if (@ARGV<1);
        my $command = shift(@ARGV);
        my %fun = (merge=>\&merge, basic=>\&basic, masking=>\&masking, rmgene=>\&rmgene);
        die("Unknown command \"$command\"\n") if (!defined($fun{$command}));
        &{$fun{$command}};
      }


sub usage {
  die(qq/
Usage:  Process_annotation.pl <command> [<arguments>]\n

Command: 

merge:          merge annotated results by blasting against various references. It also filters out redundancies by self-blasting

basic:          basic filters on length and GC content

masking:        repeat masking 

rmgene:         remove genes by names

\n/);
}



####################################################################################################################################################


sub merge {
  
  die(qq/
Process_annotation.pl merge [options] 

options:
-a FILE   master annotation file (e.g. anolis.fasta).
-b CHAR   directory with all other annotation files (e.g. a folder containing: chicken.fasta, human.fasta, fish.fasta, etc.)
-c FILE   a list of fasta sequences that are needed to be filtered out [null]. Note if you don't have any then just don't use -c
-d CHAR   output prefix

\n\n/) unless (@ARGV);
  
  my %opts = (a=>undef, b=>undef, c=>undef, d=>undef);
  getopts('a:b:c:d:', \%opts);
  
  my $dir;
  
  if ($opts{b} =~ m/\/$/ ){
    $dir = $opts{b}; 
  }
  else {
    $dir = $opts{b} . "/";
  }
  
  my @all_data = <$dir*fasta>;
  my $anolis = $opts{a};
  
  my $name = $opts{d};
  
  my $new_anolis = dirname($anolis). "/" . $name . "_final_original_names.txt";   
  readAnnotation ($anolis, $new_anolis);
  my $out = $anolis . "_only2.txt";  
  my $new_anole = readSeq ($new_anolis, $out);
  
  
  foreach (@all_data) {
    my $annotation = $_ . "_only.txt"; 
    my $out1 = $_ . "_only2.txt";
    my $no_match =  $_ . "_no_match.txt";
    
    readAnnotation($_, $annotation);
    my $new = readSeq($annotation, $out1);
    my %new = %{$new};
    my $out2 = "anolis_annotation_combined.txt"; 
    print "\n\n";
    print "now comparing ", basename($out), " with ", basename($out1), "!\n";
    PWblast($out, $out1, $out2,\%new, $no_match);
    system ("mv $out2 $out");
    my $final1 = 'temp';
    system ("cat $new_anolis $no_match > $final1");
    system ("rm $new_anolis");
    system ("mv $final1 $new_anolis");
    system ("rm $annotation");
  }
  
  unlink ($out);
  my $self1 = $new_anolis . "1";
  my $hash = renameSeq($new_anolis, $self1);
  my %hash = %{$hash};
  my $self2 = $new_anolis . "2";
  system ("cp $self1 $self2");
  my $final_nr = dirname($anolis). "/" . $name. '_annotation.final';
  
  print "\n\n";
  print "now doing self-blasting!","\n";
  self($self1, $self2, \%hash, $final_nr);
  
  if ($opts{c}) {
  my $clean = $final_nr . "1";
  my $hash1 = readSeq($final_nr, $clean);
  my %hash1 = %{$hash1};
  my $final_out = dirname($final_nr). "/1-" . $name. '_annotation_filtered.final';
  
  print "\n\n";
  print "now filtering out stuff that you do not want!","\n";
  filter($clean, $opts{c}, \%hash1, $final_out);
  system ("rm $final_nr");
  my $s = size($final_out);
  print "\n\n";
  print "the size of final dataset is ", $s, "bp!", "\n"; 
  print "\n\n";
}
  
  else {
    my $s = size($final_nr);
    print "\n\n";
    print "the size of final dataset is ", $s, "bp!", "\n";
    print "\n\n";
  }

  
  
  #############################################################
  
  sub filter {
    my ($db, $query, $hash, $out) = @_;
    my %clean = %{$hash};
    my $blastout = 'filter_blast.out';
    my $call1 = system("makeblastdb -in $db -dbtype nucl");
    my $call2 = system("blastn -db $db -query $query -evalue 1e-20 -outfmt 6 -out $blastout");
    system("rm $db.n*");
    
    my %tmp;
    open(IN, "<$blastout");
    while (<IN>) {
      chomp(my $line = $_);
      my @d = split(/\s+/,$line);
      $tmp{$d[1]}++;
    }
    close(IN);
    system ("rm $blastout");
    open (OUT, ">", $out);
    foreach my $anno (sort {$a cmp $b} keys %clean) {
      unless ($tmp{$anno}) {
	print OUT  ">", $anno, "\t", $clean{$anno}{'rest'}, "\n", $clean{$anno}{'seq'},"\n"; 
      }   
    }
    close OUT;
    unlink ($db);
    
  } 
  
  
  sub self {
    my ($query, $db, $final, $out) = @_;
    my %final = %{$final};
    my $blastout = 'blast.out';
    my $call1 = system("makeblastdb -in $db -dbtype nucl");
    my $call2 = system("blastn -db $db -query $query -evalue 1e-50 -outfmt 6 -out $blastout");
    system("rm $db.n*");
    open(IN, "<$query");
    my %seq;
    while (<IN>) {
      chomp(my @line = split /\s+/, $_);
      if ($line[0] =~ m/>(\S+)/) {
	my $id = $1;
	chomp(my $seq = <IN>);
	$seq{$id} = {'seq'=>$seq, 'len'=> length ($seq)};
      }
    }
    close(IN);
    
    my %tmp;
    open(IN, "<$blastout");
    while (<IN>) {
      chomp(my $line = $_);
      my @d = split(/\s+/,$line);
      push(@{$tmp{$d[0]}},\@d);
    }
    close(IN);
    system ("rm $blastout");
    
    foreach my $id (sort {$a cmp $b} keys %tmp) {
      if (scalar(@{$tmp{$id}}) > 1 ) {
	my %match;
	for (my $i = 0; $i < scalar(@{$tmp{$id}}); $i++) {   
	  $match{$tmp{$id}[$i][1]} ={'leng'=> $seq{$tmp{$id}[$i][1]}{'len'}};      
	  for (my $j = $i+1; $j < scalar(@{$tmp{$id}}); $j++) {
	    $match{$tmp{$id}[$j][1]} = {'leng'=> $seq{$tmp{$id}[$j][1]}{'len'}};
	    
	    if ($match{$tmp{$id}[$i][1]}{'leng'} > $match{$tmp{$id}[$j][1]}{'leng'}) {
	      delete $final{$tmp{$id}[$j][1]}  if ($final{$tmp{$id}[$j][1]});	   
	    }
	    
	    if ($match{$tmp{$id}[$i][1]}{'leng'} < $match{$tmp{$id}[$j][1]}{'leng'}) {
	      delete $final{$tmp{$id}[$i][1]} if ($final{$tmp{$id}[$i][1]});	 	   
	    }	 
	    
	    if ($match{$tmp{$id}[$i][1]}{'leng'} == $match{$tmp{$id}[$j][1]}{'leng'}) {
	      if ($tmp{$id}[$i][1] ne $tmp{$id}[$j][1]) {
		delete $final{$tmp{$id}[$i][1]} if ($final{$tmp{$id}[$i][1]});	     
	      }
	    else {
	      next;
	    }
	    }
	  } 
	} 
      }
    }   
  
    open (OUT, ">", $out);
    foreach my $anno (sort {$a cmp $b} keys %final) {
      print OUT  ">", $anno, "\t", $final{$anno}{'rest'}, "\n", $final{$anno}{'seq'},"\n"; 
    }   
    
    close OUT;  
    unlink ($query);
    unlink ($db);
  }
  
  
  sub PWblast {
    my ($anolis, $next, $out, $new, $no_match) = @_;
    my %new =%{$new};
    my $blastout = dirname($anolis) . "/". basename($anolis) . "_" . basename($next). '.blast.out';
    my $call1 = system("makeblastdb -in $next -dbtype nucl");
    my $call2 = system("blastn -db $next -query $anolis -evalue 1e-50 -outfmt 6 -out $blastout");
    
    my %r;
    open(IN, "<$blastout");
    while(<IN>) {
      chomp(my @line = split /\s+/,$_);
      $r{$line[1]}++;
    }
    close(IN);
    unlink($blastout);
    system("rm $next.n*");
    
    open (NEW, "<", $next);
    my $left = dirname($anolis) . "/". "no_match_to_anole_in_" . basename($next). ".txt";
    open (OUT, ">", $left);
    open (OUT2, ">", $no_match);
    while (<NEW>) {
      chomp(my @line = split /\s+/,$_);
      my $id = $1 if $line[0] =~ m/>(\S+)/;
      chomp(my $seq = <NEW>);
      unless ($r{$id}) {
	print OUT join ("\t", @line),  "\n", $seq, "\n";
	print OUT2 ">", $id, "\t", $new{$id}{'rest'}, "\n", $new{$id}{'seq'},"\n";
      }
    }
    
    close NEW; close OUT; close OUT2;   
    system ("cat $anolis $left > $out");
    unlink ($left);
    unlink ($next);
    
  }
  
  
  sub renameSeq {
    my ($file, $out1) = @_;
    my %hash;
    my $out = $file . '2';
    open(IN, "<$file");
    open(OUT, ">$out");
    
    my $number = 1;
    while (<IN>) {
      chomp(my @line = split /\s+/, $_);
      
      if ($line[0] =~ m/>/) {
	chomp(my $seq = <IN>);
	my $id = 'Contig' . $number;
	$hash{$id} = {'seq'=>$seq, 'rest'=>join ("\t",@line[1..$#line])};
	print OUT ">", $id, "\n", $seq,"\n";
	$number++;
      }
    }
    close(IN); close(OUT);
    rename($out,$out1);
    return (\%hash);
  }       
  
  sub readSeq {
    my ($seqfile, $out) = @_;
    my %seq;
    open(IN, "<$seqfile");
    
    open (OUT, ">", $out);
    while (<IN>) {
      chomp(my @line = split /\s+/, $_);
      if ($line[0] =~ m/>(\S+)/) {
	my $id = $1;
	chomp(my $seq = <IN>);
	$seq{$id} = {'seq'=>$seq, 'rest'=>join ("\t",@line[1..$#line])};
	print OUT ">",$id, "\n", $seq,"\n";
      }
    }
    close IN;
    close OUT;
    return(\%seq); 
  }
  
  sub readAnnotation {
    my ($seqfile, $annotation) = @_; 
    open(IN, "<$seqfile");
    open (OUT, ">", $annotation);
    while (<IN>) {
      chomp(my $line =  $_);
      if ($line =~ m/ENS/) {
	chomp(my $seq = <IN>);
	print OUT $line, "\n", $seq,"\n";
      }
      
    }
    close IN; close OUT;
  }
  
  
  sub size {
    my ($seq) = @_;
    open (IN, "<", $seq); 
    my $len;
    while (<IN>) {
      chomp(my $line = $_);
      if ($line=~ m/>(\S+)/) {
	chomp(my $seq = <IN>);
	$len += length($seq);    
      }   
    }
    close IN;
    return ($len); 
  }
  
}

###########################################################################

sub basic {
  
  die(qq/
Process_annotation.pl basic [options] 

options: 
-a INT    min length of contig [201]
-b INT    max length of contg [9999999]
-c INT    min GC content (percent) [35]
-d INT    max GC content (percent) [70]
-f file   annotated fasta sequences (1-***) outputted by "Process_annotation.pl merge" 
-g CHAR   outfile prefix
\n\n/) unless (@ARGV);
  
  my %opts = (a=>201, b=>9999999, c=>35, d=>70, f=>undef, g=>undef);
  getopts('a:b:c:d:f:g:', \%opts);
  
  my $annotation = $opts{f};
  
  my ($gc_out, $gc) = GC($annotation, $opts{c}, $opts{d}, $opts{g});
  my $s = size($gc_out);
  print "\n\n";
  print "the size of transcripts after GC filtering is ", $s, "bp!", "\n"; 
  print "\n\n";
  
  my ($len_out, $len) = cut_len($gc_out,$opts{a}, $opts{b}, $opts{g});
  my $le = size($len_out);
  print "\n\n";
  print "the size of transcripts after GC and length filtering is ", $le, "bp!", "\n"; 
  print "\n\n";
  
  
  
  sub GC {
    my ($file, $min, $max, $prefix) = @_;
    open (IN, "<", $file);
    my $out_final = dirname ($file) . "/". $prefix. "_GC_filtered.txt";
    my $out_GC =  dirname ($file) . "/". $prefix. "_contigs_with_extreme_GC.txt";
    open (FINAL, ">",  $out_final);
    open (GC, ">", $out_GC);
    
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ m/^>/) {
	chomp (my $seq = <IN>);
	my $SEQ = uc($seq);
	my $l = ($SEQ =~ tr/ATGCN//);
	my $GC = ($SEQ =~ tr/GC/gc/);
	my $percentGC = ($GC/$l)*100;
	my $rounded = sprintf("%.2f", $percentGC);	  	     
	if ($percentGC >= $min &&  $percentGC <= $max ) {	  
	  print FINAL $line, "\n", $seq,"\n";
	}
	else {
	  print GC $line, "\n", $seq, "\n", $rounded, "\n";
	}
      }
      
    }
    close IN; close FINAL; close GC;
    return ($out_final, $out_GC);	
  }
  
  sub cut_len {
    my ($file, $min, $max, $prefix) = @_;
    open (IN, "<", $file);
    my $len_final = dirname ($file) . "/2-". $prefix. "_GC_length_filtered.txt";
    my $len_out =  dirname ($file) . "/". $prefix. "_contigs_with_extreme_lengths.txt";
    open (FINAL, ">",  $len_final);
    open (LEN, ">", $len_out);
    
    while (<IN>) {
      chomp (my $line = $_);
      if ($line =~ m/^>/) {
	chomp (my $seq = <IN>);
	my $size = length ($seq);	  	     
	if ($size >= $min && $size <= $max ) {	  
	  print FINAL $line, "\n", $seq,"\n";
	}
	else {
	  print LEN $line, "\n", $seq, "\n", $size, "\n";
	}
      }
      
    }
    close IN; close FINAL; close LEN;
    return ($len_final, $len_out);	
  } 
}

####################################################################################


sub masking {
  
  die(qq/
Process_annotation.pl masking [options] 

options: 
-a FILE   repeat_masking file (repeat masker format)
-b FILE   annotated fasta sequences (2-***) outputted by "Process_annotation.pl basic" 
-g CHAR   outfile prefix
\n\n/) unless (@ARGV);
  
  my %opts = (a=>undef, b=>undef, g=>undef);
  getopts('a:b:g:', \%opts);
  
  my $repeat = $opts{a};
  my $annotation = $opts{b};
  my $name = $opts{g};
    
  my $clean = $annotation . "1";
  
  my $hash1 = readSeq($annotation, $clean);
  my %hash1 = %{$hash1};

  my $outfile = repeatmasking($repeat, $clean, $name, \%hash1);
  
  system ("rm $clean"); 
  
  my $le = re_size($outfile);
  print "\n\n";
  print "the size of transcripts after repeat masking is ", $le, "bp!", "\n"; 
  print "\n\n";

  
  sub re_size {
    my ($seq) = @_;
    open (IN, "<", $seq); 
    my $len;
    while (<IN>) {
      chomp(my $line = $_);
      if ($line=~ m/>(\S+)/) {
	chomp(my $s = <IN>);
	my $sz = ($s =~ tr/ATCGN//);
	
	$len += $sz;
      }
      
    }
    close IN;
    return ($len); 
  }

  sub repeatmasking {
    my ($repeat, $anno, $prefix, $hash) = @_; 
    my %h = %{$hash};
    my %mask;
    open (RE, "<", $repeat); 
    while (<RE>) {
      chomp(my $line = $_);
      next if ($line !~ /^\s*\d+/);
      next if ($line =~ m/Low_complexity/);
      my @d = split;   
      push(@{$mask{$d[4]}}, [$d[5], $d[6]]);
    }
    close RE;
    
    my %seq;
    open (IN, "<", $anno); 
    while (<IN>) {
      chomp (my $line = $_);
      my $ID = $1 if ($line =~ m/^>(\S+)/);
      chomp (my $seq = <IN>);
      $seq{$ID} = $seq;
    }
    close IN;
    
    my $out = dirname($anno) . "/3-". $prefix ."_repeats_masked.txt";
    open (OUT, ">", $out);
 
    foreach my $id (sort {$a cmp $b} keys %seq) {
      if ($mask{$id}) {	
	my $seq = $seq{$id};
	for (my $i = 0; $i < scalar(@{$mask{$id}}); $i++) {   
	  my $start = $mask{$id}[$i][0];
	  my $end = $mask{$id}[$i][1];
	  my $s = min ($start, $end);
	  my $e = max ($start, $end);	  	  
	  my $masked = 'n' x ($e-$s+1);	  
	  my $length = length ($masked);
	  substr ($seq, $s-1, $length) = $masked;
	}
	print OUT ">", $id, "\t", $h{$id}{'rest'}, "\n", $seq, "\n"; 
      }     
      else {
	my $seq = $seq{$id};
	print OUT ">", $id, "\t", $h{$id}{'rest'}, "\n", $seq, "\n";	
      }
      
    }
    close OUT;
    return ($out);

  }

}

##############################################################################################

sub rmgene {
  die (qq/
Usage: Process_annotation.pl rmgene [options]

Options:
        -genes       CHAR.....         genes to be removed (e.g. rRNA snRNA 18S ribosomal Histone). You can give as many as possible (case insensitive)
        -annot       FILE              annotated fasta sequences (3-***) outputted by "Process_annotation.pl masking"
        -out         CHAR              outfile prefix

\n/) if !@ARGV;
  
  
  my ($annot, $genes, $out) = (undef, undef, undef);
  GetOptions('genes=s@{1,}' => \$genes,'annot=s@{1,1}' => \$annot,'out=s@{1,1}' => \$out);
 
  my $dir = dirname (@{$annot}[0]) . "/";
  my $file =  @{$annot}[0]; 
  my $prefix = @{$out}[0]; 
  

  my %seq;
  open(IN, "<", $file);
  while (<IN>) {
    chomp(my @line = split /\s+/, $_);
    if ($line[0] =~ m/>(\S+)/) {
      my $id = $1;
      chomp(my $seq = <IN>);
      $seq{$id} = {'seq'=>$seq, 'rest'=>join ("\t",@line[1..$#line])};
    }
  }
  close IN;

  my $out_seq = $dir . "4-". $prefix. "_with_unwanted_genes_being_removed.txt";
  my $genes_removed = $dir . $prefix. "_removed_genes.txt";
  open (FINAL, ">", $out_seq);
  open (REMOVE,">", $genes_removed);

  foreach (@{$genes}) {
    my $gene = $_;
    foreach my $id (sort {$a cmp $b} keys %seq) { 
      if ($seq{$id}{'rest'} =~ m/$_/i) {
      print REMOVE ">", $id, "\t", $seq{$id}{'rest'}, "\n", $seq{$id}{'seq'}, "\n";
      delete $seq{$id};  
      }    
    }
  }
  close REMOVE;
  
  foreach my $id (sort {$a cmp $b} keys %seq) { 
    print FINAL ">", $id, "\t", $seq{$id}{'rest'}, "\n", $seq{$id}{'seq'}, "\n";
  }
  close FINAL;

}
  
    











