#!/usr/bin/perl -w
# Reading protein sequence data from fake_protein file
use strict;
use warnings;
use List::Util qw/sum/;
#use BeginPerlBioinfo;
my @PDB_data = get_file_data($ARGV[0]);
my %recordtypes = parsePDBrecordtypes(@PDB_data);
my @atoms = parseATOM ($recordtypes{'ATOM'});
my @x; my @y; my @z;
foreach my $atom (@atoms){
  my $x = substr($atom, 0, 2);
  my $y = substr($atom, 3, 2);
  my $z = substr($atom, 6, 2);
  push @x, $x;
  push @y, $y;
  push @z, $z;
}
my $x = sum @x;
my $y = sum @y;
my $z = sum @z;
my $n = @x;
my $xc = $x/$n;
my $yc = $y/$n;
my $zc = $z/$n;
print "The centroid of the atoms is ";
print "($xc,$yc,$zc)";
#print @x, "\n";
#print @y, "\n";



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
