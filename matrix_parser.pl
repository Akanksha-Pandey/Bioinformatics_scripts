#/usr/bin/perl
use strict;
use warnings;
use Cwd;

# Usage perl matrix_parser.pl <filename>#

#this script takes an input half matrix file (used in programs like Iqtree, raxml etc. ) and creates a column wise output

my $dir = getcwd();

my $file = shift;

my %matrix_hash =();
my @lin_arr =();
open(my $FH,"<","$file");
{
	my @content = <$FH>;
	my $len = scalar(@content);
	for(my $i=0;$i<$len;$i++)
	{
	chomp($content[$i]);
	my $line = $content[$i];
	chomp($line);
	@lin_arr = split(/\s+/,$line);
	my $len_sc = scalar(@lin_arr);
	for(my $j=0;$j<=$len_sc;$j++)
	{
	$matrix_hash{$i}[$j] = $lin_arr[$j];
	}
	}
		for (my $i=0;$i<$len;$i++)
	{
	for my $key ( sort {$a<=>$b} keys %matrix_hash)
	{
	if ((exists ${$matrix_hash{$key}}[$i]) and (${$matrix_hash{$key}}[$i] =~ /\w/))
	{
	print"${$matrix_hash{$key}}[$i]\n";
	}
	else
	{
	next;
	}
	}
	
	}
}