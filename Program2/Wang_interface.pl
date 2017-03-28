#!/usr/bin/perl

use strict;

# Extract different chains
my @chain1 = extractAtom($ARGV[1],$ARGV[0]);
my @chain2 = extractAtom($ARGV[2],$ARGV[0]);
my $name1 = $ARGV[1];
my $name2 = $ARGV[2];
my @dist = ();
my $threshold = $ARGV[3];
# Compute distance and output interface amino acids.
for(my $i=0; $i < scalar @chain2; $i++)
{
  for(my $j=0; $j < scalar @chain1; $j++)
  {
    my $distance = computedist($chain2[$i][0],$chain2[$i][1],$chain2[$i][2],$chain1[$j][0],$chain1[$j][1],$chain1[$j][2]);
    #store the distance
    $dist[scalar @dist] = $distance;
    if($distance<$threshold)
    {
      print $name1,":",$chain1[$j][3],"(",$chain1[$j][4],")"," interacts with ",$name2,":",$chain2[$i][3],"(",$chain2[$i][4],")","\n";
    }
  }
}


# Subroutine computing distance
sub computedist
{
  my @co = @_;
  my $distance = sqrt(($co[0]-$co[3])**2+($co[1]-$co[4])**2+($co[2]-$co[5])**2);
  return $distance;
}

# Subroutine extrating C-alpha atom in chain E
sub extractAtom{
my $identifier = shift(@_);
my ($filename) = @_;
unless (open(File, $filename))
{
  print "Could not open file $filename!\n";
  exit;
}
my @E=();
while(my $line = <File>){
	chomp($line);
	if($line=~m"^ATOM"){
		# get substring
		my $x = substr($line, 30, 8);  # skip the first 30 characters, read next 8 characters
		my $y = substr($line, 38, 8);
		my $z = substr($line, 46, 8);
    my $amino = substr($line, 17, 3);
    my $aanum = substr($line, 22, 4);
		# convert string to float
		$x=$x*1;
		$y=$y*1;
		$z=$z*1;
    $aanum=$aanum*1;
    my $name = substr($line, 12, 4);
    my $chain = substr($line, 21, 1);
    if($name=~"CA")
    {
		# store data
    if($chain=~"$identifier"){
		my $length = scalar @E;
		$E[$length][0]=$x;
		$E[$length][1]=$y;
		$E[$length][2]=$z;
    $E[$length][3]=$amino;
    $E[$length][4]=$aanum;
    }
	  }
  }
 }
close(File);
return @E;
}
