#/usr/bin/perl
use strict;
use warnings;
use Cwd;
# Usage remove_gaps.pl <filename>#
#It removes gaps from each sequence in a MSA file#

my $dir = getcwd();

my $input_file = shift; 
my $out_file = $input_file . ".degapped.fasta";

open(my $FH1, "<", "$dir/$input_file");
{
open(my $FH2, ">>", "$dir/$out_file");
{
while(my $line = <$FH1>)
	{
	chomp($line);
	if($line =~ /^>/)
	{
	#print "$line\n";
	print $FH2 "$line\n";
	}
	else
	{
	#$line =~ s/[a-z]+?/N/g;
	#$line =~ s/N{3,10000}/NNNNN/g;
	$line =~ s/-+?//g;
	#print "$line\n";
	print $FH2 "$line\n";
	}
	}
}
}

