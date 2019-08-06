#/usr/bin/perl
use strict;
use warnings;
use Cwd;

#Usage perl NexustoPhylip.pl <file_list>#
#It takes a list of nexus alignments and convert each file into phylip alignment file#

my $dir = getcwd();
my $file_list = shift;

my @loci = read_file($file_list);

foreach my $locus (@loci)
{
my $input_file = $locus . ".nex";
my $out_file = $locus . ".phy";
my @taxon =();
my @seq = ();

open(my $FH,"<","$dir/$input_file") or die "couldn't read the file $!\n";
{
open( my $FH1,">>","$dir/$out_file") or die "Couldn't create $!\n";
	{
		while(my $line = <$FH>)
		{
		chomp($line);
		if ($line =~ m/\bdimensions|Dimensions\b/)
		{
		my ($word11,$word22,$word33) = split(/ /,$line);
		#print"$word11\t,$word22\t,$word33\n";
		my ($chr11,$chr22) = split(/=/,$word22);
		#print "$chr22$chr11\n";
		my ($chr1,$chr2) = split(/=/,$word33);
		#chop($chr22);
		chop($chr2);
		#print "$chr22";
		print $FH1 "$chr22\t$chr2 \n";
		}
		if ($line=~/Nexus/i || $line=~/Matrix/i || $line=~/\;/ || $line=~/\,/ || $line=~/format/i || $line eq /\s+/|| $line =~/\bTREE\b/i )  #|| $line =~/TREE/i
				{ 
				next;
				}
			
		else
				{
					my ($taxon,$seq) = split(/\s+/,$line);
					print $FH1 "$taxon\t$seq\n";
					
				}
		}
	}close($FH);

print "$input_file\n";
}
}

sub read_file  #This function is to read the text file and return a list of contents
{
	my $file = shift;
	chomp($file);
	my @content = ();
	open(my $fh1,"<","$dir/$file") or die "Can't open $file\n"; #opens the input file 
	{
	while(my $line = <$fh1>)
	{
	chomp($line);
	
		push(@content,$line);
	
	}
	}
		return @content;
}