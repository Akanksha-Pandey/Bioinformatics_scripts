#/usr/bin/perl
use strict;
use warnings;
#use Math::Random::MT::Perl;
use Cwd;

my $dir = getcwd();
#This code creates random numbers from a given range for different size of subsets and multiple also generates multiple replicates of those subsets

my $range =401632;
my @subsets = (47986,54014,161897,194117);
my $run =1;
while($run <=10)
{
srand();
my @maskarray = ();
foreach my $total_sample (@subsets)
{
my $catagory = "test" . "_" . $run ."_$total_sample";
	for(my $i=0;$i<= $range;$i++)
	{
	$maskarray[$i]=0;
	}
	
my $counter = 0;
	while($counter < $total_sample)
	{
	my $r_num = int(rand($range));
	if( $maskarray[$r_num] ==0)
	{
	$maskarray[$r_num] =1;
	$counter++;
	}
	}
open(my $FH1,">>","$dir/Simion_consensus.nex") or die "couldn't read the file $!\n"; # This part ids used to print those random numbers as charsets to a nexus alignment#
{
	print $FH1 "CHARSET $catagory =";
	#print "CHARSET $catagory =";
	for(my $j=0;$j<= $range;$j++)
	{
		if( $maskarray[$j] ==1)
		{
		my $k =$j+1;
		print $FH1 "$k ";
		#print "$k ";
		}
		else
		{
		next;
		}
	}
print $FH1 ";\n";
#print ";\n";
}
}
$run++;
}
	
