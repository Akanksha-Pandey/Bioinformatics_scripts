#/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Cwd;

###Usage perl create_consensus.pl <file_list>###
#It takes a list of multiple sequence alignments in nexus format and creates a consensus sequence based on henikoff and henikoff weights in fasta format per file
# <file_list> is a a text file with names of nexus alignment file without any suffix##
my $dir = getcwd();

my ($Query_list) = shift;

my @loci = ();
my @new =();
my $count_gene=1;

if (-s "$dir/$Query_list") # This part of the code checks whether query list files are present in current working directory
	{
		print "$Query_list is present in $dir\n";
		my @content = read_file($Query_list);
		foreach my $temp(@content)
		{
			my $file  = $temp . ".upd.nex";
			if (-s "$dir/$file")
			{
				print "File $file is present in $dir\n";
				chomp($temp);
				push(@loci,$temp); 
			}
		}
	}
foreach my $locus (@loci)  # Get a file with all loci and add them here in for loop(it initializes names for all other files
{
my $file = $locus . ".upd.nex";
my $file_fasta = "$locus.fasta";
my $updated_nex = $locus . ".updated" . ".txt";
###############Parsing only sequences from nexus file into a scratch file###############
	open (my $FH_input,"$dir/$file") or die "Couldn't read file $!";
	{
		while (my $t = <$FH_input>)
		{
			chomp($t);
			if ($t=~/Nexus/i || $t=~/Matrix/i || $t=~/\;/ || $t=~/\,/ || $t=~/format/i || $t eq /\s/) 
				{ 
				next;
				}
			else
				{
					push (@new,$t);
					open(my $FH1,">>","$dir/$updated_nex")or die "Couldn't read file $!";
					{
					print $FH1 "$t\n";
					}close($FH1);
				}
		}
	}close($FH_input);

create_consensus($updated_nex,$file_fasta,$locus);
}

##################################################Subroutines##################################################
sub create_consensus
{
my $input_file = $_[0];
my $out_file = $_[1];
my $header = $_[2];
my %seq_hash =();
my $length =0;
my %w_hash =();
my %sum_weight =();
my %aa_weight =();

open(my $FH,"<","$dir/$input_file") or die "Couldn't read file $!\n";
{
	while(my $line = <$FH>)             #It creates a hash of array for the sequence file
		{
		chomp($line);
		#print"$line\n";
		my ($taxon,$seq) = split(/\s+/,$line);
		#print "$taxon\n";
		chomp($seq);
		my @sequences = split(//,$seq);
		for(my $j=0;$j<=$#sequences;$j++)
			{
			$seq_hash{$taxon}[$j] = $sequences[$j];
			}
			$length = scalar( @{ $seq_hash{$taxon} } );
			#print"The length of the array for $taxon is $length\n"; 
		}
}		
for(my $i=0;$i<=$length;$i++)               ##It counts different types of amino acid at a position 
{
my %aa_count =();
for my $temp(keys %seq_hash)
	{
		my $char1 = $seq_hash{$temp}[$i];
		chomp($char1);
		$aa_count{$char1}++;                  #creates a hash with aa and its count 
	}
				#print"position $i\n";
				my $diff_res = scalar(keys %aa_count);         #counts number of entries in array which is equal to different types of residue
				#print" The number of different residues at $i = $diff_res\n";
				for my $key (keys %seq_hash)                    #create a hash with position weights (1/r*s)
				{
					my $seq_char = $seq_hash{$key}[$i];
					if (exists $aa_count{$seq_char})
					{
					my $val = $aa_count{$seq_char};          #The number of times this residue occurred = no of sequences it is present
					my $weight = 1/($diff_res * $val);       # Weight of that residue as described by Henikoff 
					#print "$weight  for $seq_char position $i\n";
					$w_hash{$key}[$i] = $weight;              # create a matrix (hash of arrays) of weights for all residues in the sequence
					}
				}
}
foreach (keys %w_hash)            # get the sum of all individual aa weights. It gives the sequence weight
{
my $sum_seq = sum(@{$w_hash{$_}});
$sum_weight{$_} = $sum_seq;
#print "The sum for $_ is $sum_seq\n";
}

foreach (keys %seq_hash)            # Assign the sequence weight to individual aa in the sequence
{
 my $seq_w = $sum_weight{$_};
 #print" The weight for $_ is $seq_w\n";
 for(my $i=0;$i<=$length;$i++)
 	{
 	if( $seq_hash{$_}[$i] ne "-")
 		{
 		$aa_weight{$_}[$i] = $seq_w;        # A new hash stores sequence weight for each residue for overall sequence
 		}
 	else
 		{
 		$aa_weight{$_}[$i] = 0;
 		}
 	}
}
my @consensus = ();
for(my $i=0;$i<=$length;$i++)           # Find the maximum weight aa at the given position and identify that residue to create a consensus sequence
{
	my %aa_wc =();
	foreach (keys %seq_hash)
	{
		my $pos_char = $seq_hash{$_}[$i];
		my $pos_val = $aa_weight{$_}[$i];
		$aa_wc{$pos_char} = $aa_wc{$pos_char} + $pos_val;
	}
	foreach my $key (sort { $aa_wc{$b} <=> $aa_wc{$a} } keys %aa_wc)
	{
        #printf "%4d %s\n", $aa_wc{$key}, $key;
        push(@consensus,$key);
        last;
    }
	#print Dumper(\%aa_wc);
}
my $final_seq = join('',@consensus);

open(my $fh1,">>","$dir/$out_file") or die "Couldn't write to $!\n";
{
print $fh1 ">$header\n$final_seq";
}

#print"$final_seq\n";
#print Dumper(\%aa_weight);
#print Dumper(\%sum_weight);
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
