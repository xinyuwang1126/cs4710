#!/usr/bin/perl -w

use strict;
use warnings;
use List::Util qw/sum/;
use List::Util qw/max min/;
# Reading protein sequence data from fake_protein file
my @PDB_data = get_file_data($ARGV[0]);
my %recordtypes = parsePDBrecordtypes(@PDB_data);
my @atoms = parseATOM ($recordtypes{'ATOM'});
#PART1
my @x; my @y; my @z;
foreach my $atom (@atoms){
  my $x = substr($atom, 0, 2);
  my $y = substr($atom, 3, 2);
  my $z = substr($atom, 6, 2);
  push @x, $x;
  push @y, $y;
  push @z, $z;
}
my $n = @x;
my @d;
for (my $i=0; $i+1<$n; $i++){
  for (my $j=$i+1; $j<$n; $j++){
    my $x2 = ($x[$j]-$x[$i])*($x[$j]-$x[$i]);
    my $y2 = ($y[$j]-$y[$i])*($y[$j]-$y[$i]);
    my $z2 = ($z[$j]-$z[$i])*($z[$j]-$z[$i]);
    my $d2 = $x2 + $y2 + $z2;
    my $d = sqrt($d2);
    push @d, $d;
  }
}
my $diameter = max @d;
print "\n","PART1","\n";
print "The distances of all pairs of atoms are","\n";
foreach my $d (@d){
  print $d, ",";
}
print "\n\n";
print "The diameter of the protein is ";
print $diameter,"\n";

#PART2
print "\n","PART2","\n";
my $min = min @d; my $max = max @d;
print "The minimum distance is ",$min, "\n";
print "The maximum distance is ",$max, "\n";
my $f1=0; my $f2=0; my $f3=0; my $f4=0;
foreach my $d (@d){
  if (1<=$d && $d<1.2) {$f1++;}
  elsif (1.2<=$d && $d<1.4) {$f2++;}
  elsif (1.4<=$d && $d<1.6) {$f3++;}
  elsif (1.6<=$d && $d<1.8) {$f4++;}
}
print "Bin size is 0.2, starting from 1.0 to 1.8","\n";
print "The coresponding frequency for each bin is","\n";
print "1.0<=f[i]<1.2: ",$f1,"\n";
print "1.2<=f[i]<1.4: ",$f2,"\n";
print "1.4<=f[i]<1.6: ",$f3,"\n";
print "1.6<=f[i]<1.8: ",$f4,"\n";
#get_file_data
sub get_file_data {
  my($filename) = @_;
  use strict;
  use warnings;
# Initialize variables
  my @filedata = ( );
  unless( open(GET_FILE_DATA, $filename) ) {
  print STDERR "Cannot open file \"$filename\"\n\n";
  exit;
  }
  @filedata = <GET_FILE_DATA>;
  close GET_FILE_DATA;
  return @filedata;
  }
#parsePDBrecordtypes
sub parsePDBrecordtypes {
  my @file = @_;
  use strict;
  use warnings;
  my %recordtypes = ( );
  foreach my $line (@file) {
    my($recordtype) = ($line =~ /^(\S+)/);
    if (defined $recordtypes{$recordtype} ) {
    $recordtypes{$recordtype} .= $line;}
    else {$recordtypes{$recordtype} = $line;}
   }
  return %recordtypes;
}
#parseATOM
sub parseATOM {
  my($atomrecord) = @_;
  use strict;
  use warnings;
  my @results;
  my(@atomrecord) = split(/\n/, $atomrecord);
  foreach my $record (@atomrecord) {
    my $number = substr($record, 6, 5);
    my $x = substr($record, 31, 2);
    my $y = substr($record, 33, 2);
    my $z = substr($record, 35, 2);
    push @results, "$x $y $z";
  }
  return @results;
}

exit;
