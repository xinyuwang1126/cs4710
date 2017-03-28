use strict;
use warnings;
use File::Find;

my $directory = "proteins";
my @files = read_dir($directory);
my $n = @files;
my $sumRh = 0.0; my $sumRe = 0.0; my $sumQ3 = 0.0;
foreach my $e(@files)
{
	my $name = substr($e, 0, -4);
	my @result = call_stride("$directory/$e");
	my ($str,$both) = parse_stride(@result);
  print "Name: ", $name; print "\n","\n";
# calculate Nh, Ne
  open(my $FH, "proteins/$e");
	my $Nh = 0; my $Ne = 0;
	while(my $line = <$FH>)
	{
		chomp($line);
		if($line=~m"^HELIX")
		{
			# get substring
			my $h = substr($line, 71, 5);
			# convert string to float
			$h = $h*1;
			# sum data
			$Nh += $h;
		}
		if($line=~m"^SHEET")
		{
			# get substring
			my $start = substr($line, 22, 4);
			my $terminal = substr($line, 33, 4);
			my $e = $terminal*1 - $start*1 +1;
			# sum data
			$Ne += $e;
		}
	}

# get the array of predicted structure
	my @str = split //, $str;

# use hash
  my %stat = ('Rh'=> 0, 'Re'=> 0, 'Q3'=> 0);
# count Rh, Re, Q3
  open($FH, "$directory/$e");
  while(my $line = <$FH>)
  {
	  chomp($line);
	  if($line=~m"^HELIX")
	  {
	  	# get substring
		  my $start = substr($line, 21, 4)*1;
			my $terminal = substr($line, 33, 4)*1;
      while($start <= $terminal)
			{
				if($str[$start -1] eq 'H')
				 {$stat{'Rh'} += 1.0/$Nh;
			    $stat{'Q3'} += 1.0/($Nh+$Ne);}
				elsif($str[$start -1] eq 'G')
				 {$stat{'Rh'} += 1.0/$Nh;
				  $stat{'Q3'} += 1.0/($Nh+$Ne);}
				elsif($str[$start -1] eq 'I')
				 {$stat{'Rh'} += 1.0/$Nh;
					$stat{'Q3'} += 1.0/($Nh+$Ne);}

				$start++;
			}

	  }
	  if($line=~m"^SHEET")
	  {
		  # get substring
		  my $start = substr($line, 22, 4)*1;
	  	my $terminal = substr($line, 33, 4)*1;
			while($start <= $terminal)
			{
				if($str[$start -1] eq 'E')
				{$stat{'Re'} += 1.0/$Ne;
			   $stat{'Q3'} += 1.0/($Nh+$Ne);}
				$start++;
			}
	  }
  }

  $sumRh += $stat{'Rh'}; $sumRe += $stat{'Re'}; $sumQ3 += $stat{'Q3'};
	print $both; print "\n";
	print "Rh = "; printf("%6.4f", $stat{'Rh'}); print "\n";
	print "Re = "; printf("%6.4f", $stat{'Re'}); print "\n";
	print "Q3 = "; printf("%6.4f", $stat{'Q3'}); print "\n";
	print "\n";
};
my $avgRh = $sumRh/$n*1.0; print "avgRh = "; printf("%6.4f",$avgRh); print "\n";
my $avgRe = $sumRe/$n*1.0; print "avgRe = "; printf("%6.4f",$avgRe); print "\n";
my $avgQ3 = $sumQ3/$n*1.0; print "avgQ3 = "; printf("%6.4f",$avgQ3); print "\n";

sub read_dir
{
	my ($directory) = @_;
	opendir(FOLDER, $directory) || die "Error in opening directory.";
	my @files = grep(/\.pdb/, readdir(FOLDER));
	closedir(FOLDER);
	return @files;
}

sub call_stride
{
  my($file) = @_;
  my $stride = "/usr/local/bin/stride";
  my $options = "";
  my @result = `$stride $options $file`;
  return @result;
}
sub parse_stride
{
	my @result = @_;
  # Extract the lines of interest
	my @seq = grep(/^SEQ/,@result);
	my @str = grep(/^STR/,@result);
	my @both = grep(/^SEQ|^STR/,@result);

	my $both = join("", @both);
  # Process those lines to discard all but the sequence or structure information
	for(@seq){$_ = substr($_,10,50);}
	for(@str){$_ = substr($_,10,50);}

  # Return the information as an array of two strings
	my $seq = join("", @seq);
	my $str = join("", @str);
	#print "$str\n";
	return ($str, $both);
}
