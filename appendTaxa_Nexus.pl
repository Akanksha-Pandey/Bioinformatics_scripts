#/usr/bin/perl
use strict;
use warnings;
use Cwd;

my $working_dir = getcwd();

if ($#ARGV != 0 ) { # This part of the code is to check whether you entered all the arguments at the command line
	print "usage: perl appendTaxa_Nexus.pl <File_list> \n"; 
	exit;
}
my $file_list = shift;
print "You entered the following arguments at the command line\n";
foreach my $temp(@ARGV) #prints arguments entered at the command line 
{
print"$temp\n";
}
my @files = ();
print "Checking query files in $working_dir\n";
my @input_files = read_file($file_list);
	foreach my $temp(@input_files)
	{
		my $file  = $temp . ".nex";
		if (-s "$working_dir/$file")
		{
			print "File $file is present in $working_dir\n";
			chomp($temp);
			push(@files,$temp); 
			next;
		}
		else
		{
			die ("File $file is not found! Exiting the pipeline\n");
		}
	}

foreach my $test(@files)
{
read_Nexus_out_fasta($test);
}

########################################Subroutines################################################
sub read_file  #This function is to read the text file and return a list of contents
{
	my $file = shift;
	chomp($file);
	my @content = ();
	open(my $fh1,"<","$working_dir/$file") or die "Can't open $file\n"; #opens the input file 
	{
	while(my $line = <$fh1>)
	{
	chomp($line);
	
		push(@content,$line);
	
	}
	}
		return @content;
}

sub read_Nexus_out_fasta  #This function is to read the text file and return a list of contents
{
	my $temp = shift;
	my $file = $temp . ".nex";
	chomp($file);
	#my @content = ();
	my @taxa = ();
	my @sequence = ();
	my $num_taxa = 0;
	my $num_sites = 0;
	my $out_file = $temp . ".fasta";
	open( my $FH1,"<","$working_dir/$file") or die "Couldn't read $file $!\n";
	{
	open(my $FH2,">>","$working_dir/$out_file") or die "Couldn't write $out_file $!\n";
	{
		while(my $line = <$FH1>)
		{
		chomp($line);
		if ($line =~ m/\bDimensions\b/)
		{
		my ($word11,$word22,$word33) = split(/ /,$line);
		#print"$word11\t,$word22\t,$word33\n";
		my ($chr11,$chr22) = split(/=/,$word22);
		#print "$chr22$chr11\n";
		my ($chr1,$chr2) = split(/=/,$word33);
		#chop($chr22);
		chop($chr2);
		#print "$chr22";
		$num_sites = $chr2; 
		}
		if ($line=~/Nexus/i || $line=~/Matrix/i || $line=~/\;/ || $line=~/\,/ || $line=~/format/i || $line eq /\s/|| $line =~/\bTREE\b/i )  #|| $line =~/TREE/i
				{ 
				next;
				}
			
		else
				{
					my ($taxon,$seq) = split(/\s+/,$line);
					#print $FH1 ">$temp$taxon\n$seq\n";
					my $label = $temp . "_" .$taxon;
					#print ">$label\n$seq\n";
					print $FH2 ">$label\n$seq\n\n";
					#push(@taxa,$taxon);
					#push(@sequence,$seq); 
					
				}
		}
		$num_taxa = scalar(@taxa);
	}
	}close($FH1);
}
