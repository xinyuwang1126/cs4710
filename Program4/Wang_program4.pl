#!/usr/bin/perl
use strict;
use List::MoreUtils qw/uniq/;

my ($protein_ref,$n,$matrix_ref) = Create_Adjacency_Matrix('kshv.sif');
my @protein = @$protein_ref;
my %matrix = %$matrix_ref;

foreach(@protein)
{
  my $i = 0;
  while($i<$n){
#  print $matrix{$_}{$protein[$i]};
  $i++;
}
#print "\n";
}

my %degree; #store the degree of each protein
foreach(@protein)
{
  $degree{$_} = 0;
}
foreach(@protein)
{
  my $i = 0;
  while($i<$n){
  if ($matrix{$_}{$protein[$i]}==1)
  {$degree{$_}++;}
  $i++;
}
}
#print $degree{$protein[4]};
#print "\n";

my @degree;
my $i=0;
while ($i<$n)
{
  push @degree,$degree{$protein[$i]};
  $i++;
}
my @uniquedegree = uniq @degree;
#print @degree;
my %stat;
foreach(@uniquedegree)
{$stat{$_}=0;}
foreach(@degree)
{$stat{$_}++;}

#foreach(%stat) {print "$_\n";}
#print "k=","$_\n" for keys %stat;
my $key;
foreach $key (keys %stat)
{
  print "k=$key, pk=$stat{$key}\n";
}
print "\n";
my %neighbors;

foreach(@protein){
  $neighbors{$_}=[];
}

foreach(@protein) {
  my $i = 0;
  while($i<$n){
  if ($matrix{$_}{$protein[$i]}==1)
  {push $neighbors{$_},  $protein[$i]; }
  $i++;
 }
}
my %n;
foreach(@protein){$n{$_}=0;}
foreach(@protein){
my $neighbor_ref = $neighbors{$_};
my @array = @$neighbor_ref;
#print @array,"\n";
my $size = scalar @array;
  my $count=0;
  my $i=0;
  while($i<$size-1)
  {
    my $j = $i+1;
    while($j<$size){
      if($matrix{$array[$i]}{$array[$j]}==1){
        $count++;
      }
      $j++;
    }
    $i++;
  }
$n{$_}=$count;
}

my %coefficient;
foreach(@protein){$coefficient{$_}=0;}
foreach(@protein)
{
  if($degree{$_}==0||$degree{$_}==1){$coefficient{$_}=0;}# if it has degree 0 or 1, the coefficient is 0.
  else {$coefficient{$_}=2*$n{$_}/($degree{$_}*($degree{$_}-1));}
}
my $sum=0;
my $key;
foreach $key (keys %coefficient)
{
  $sum += $coefficient{$key};
  #print "protein=$key, coefficient is $coefficient{$key}\n";
}
my $AVG_C= $sum/$n*1.0;
print "AVERAGE_C is $AVG_C\n\n";

my %avgc_k;
foreach(@uniquedegree)
{$avgc_k{$_}=0;}
my $key;
foreach $key (keys %avgc_k)
{
  my $sum = 0;
  foreach(@protein){
    if($key==$degree{$_}){
     $sum += $coefficient{$_};
    }
  }
  $avgc_k{$key} = $sum/$stat{$key}*1.0;
  print "k=$key, AVERAGE_C(k) is $avgc_k{$key}\n";

}
print "\n";
my %compute_max_nodes = %coefficient;
my $i=0;
print "Top 5 nodes with highest cluster coefficient:\n";
while($i<5){
my $max_coefficient = 0;
my $key; my $node;
foreach $key (keys %compute_max_nodes)
{
  if($compute_max_nodes{$key}>$max_coefficient)
  {
    $node = $key;
    $max_coefficient = $compute_max_nodes{$key};
  }
}
print "$node, coefficient is $max_coefficient\n";
$compute_max_nodes{$node}=0;
$i++;
}
print "\n";
my %A = {%matrix};
foreach(@protein){
  my $i = 0;
  while($i<$n){
  if ($matrix{$_}{$protein[$i]}!=1)
  {$A{$_}{$protein[$i]}=500;}
  else {$A{$_}{$protein[$i]}=1;}
  $i++;
}
}

my $k=0;
while($k<$n){
  my $i=0;
  while($i<$n){
    my $j=0;
    while($j<$n){
      if($A{$protein[$i]}{$protein[$k]}+$A{$protein[$k]}{$protein[$j]}<$A{$protein[$i]}{$protein[$j]})
      {$A{$protein[$i]}{$protein[$j]}=$A{$protein[$i]}{$protein[$k]}+$A{$protein[$k]}{$protein[$j]};}
      $j++;
    }
    $i++;
  }

  $k++;
}

my %closeness;
foreach(@protein)
{
  my $i = 0;
  my $sumlength=0;
  while($i<$n){
  if($protein[$i]ne$_){$sumlength+=$A{$_}{$protein[$i]};}
  $i++;
}
  $closeness{$_}=1.0/$sumlength;
}

my $key;
foreach $key (keys %closeness)
{
#  print "protein=$key, closeness is $closeness{$key}\n";
}
my $i=0;
print "Top 5 nodes with highest closeness centrality:\n";
while($i<5){
my $max_closeness = 0;
my $key; my $node;
foreach $key (keys %closeness)
{
  if($closeness{$key}>$max_closeness)
  {
    $node = $key;
    $max_closeness = $closeness{$key};
  }
}
print "$node, closeness centrality is $max_closeness\n";
$closeness{$node}=0;
$i++;
}






sub Create_Adjacency_Matrix {
my ($filename) = @_;
open(my $fh, $filename);
my @proteinlist = ();
while (my $line = <$fh>)
  {
    chomp($line);
    my ($protein1,$protein2) = split /1.0/, $line;
    $protein1 =~ s/^\s+|\s+$//g;
    $protein2 =~ s/^\s+|\s+$//g;
    push @proteinlist, $protein1;
    push @proteinlist, $protein2;
  }
my @uniqueproteinlist = uniq @proteinlist;
#foreach(@uniqueproteinlist)
#{
#  print "$_\n";
#}
my $n = @uniqueproteinlist;
my %adjacencyMatrix;
foreach(@uniqueproteinlist)
{
  my $i = 0;
  while($i<$n){
  $adjacencyMatrix{$_}{$uniqueproteinlist[$i]}=0;
  $i++;
}
}
open(my $fh, $filename);
while (my $line =<$fh>)
{
  chomp($line);
  my ($protein1,$protein2) = split /1.0/, $line;
  $protein1 =~ s/^\s+|\s+$//g;
  $protein2 =~ s/^\s+|\s+$//g;
  $adjacencyMatrix{$protein1}{$protein2} = 1;
  $adjacencyMatrix{$protein2}{$protein1} = 1;
}

return (\@uniqueproteinlist,$n,\%adjacencyMatrix);
}
