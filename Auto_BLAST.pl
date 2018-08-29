#/usr/bin/perl
use strict;
use warnings;
use Cwd;

# This pipeline is created to extract YouGoGlenCoco and Zelda phages ORF from a list of B type phages
# It requires: 
# <Genome_list> List of the genomes stored in BType_phage_genome folder in current working directory (Make sure the makblastdb is executed on these genomes) It should not have any ".txt" or ".fasta" in the list file 
# <Query_list> List of the query sequences stored as a single file in the query folder in current working directory 

my $working_dir = getcwd();

if ($#ARGV != 1 ) { # This part of the code is to check whether you entered all the arguments at the command line
	print "usage: perl Auto_BLAST.pl <Genome_list> <Query_list> \n"; 
	exit;
}

print "You entered the following arguments at the command line\n";
foreach my $temp(@ARGV) #prints arguments entered at the command line 
{
print"$temp\n";
}

my ($Genome_list,$Query_list) = @ARGV;
my @Query = ();
my @Btype_genome = ();

print "checking $Genome_list and $Query_list file in $working_dir\n";

if ((-s "$working_dir/$Genome_list") and (-s "$working_dir/$Query_list")) # This part of the code checks whether the genome list file and query list file is present in current working directory
	{
		print "$Genome_list and $Query_list is present in $working_dir\n";
		my @content = read_file($Genome_list);
		foreach my $temp(@content)
		{
			my $file  = $temp . ".txt";
			if (-s "$working_dir/BType_phage_genome/$file")
			{
				print "File $file is present in $working_dir\n";
				chomp($temp);
				push(@Btype_genome,$temp); 
			}
		else
			{
				die ("File $file is not found! Exiting the pipeline\n");
			}
		}
		print "Checking query files in $working_dir\n";
		my @query = read_file($Query_list);
		foreach my $temp(@query)
		{
			my $file  = $temp . ".fasta";
			if (-s "$working_dir/query/$file")
			{
				print "File $file is present in $working_dir\n";
				chomp($temp);
				push(@Query,$temp); 
			}
		else
			{
				die ("File $file is not found! Exiting the pipeline\n");
			}
		}
	}
else
	{
		die ("$Genome_list or $Query_list is not found! Exiting the pipeline\n");
	}
my $logfile = "Log.txt";
open( my $fh6 ,">>","$working_dir/$logfile");
{
foreach my $temp(@Query)# This part of the code calls tblastx  on each query and genome
{
	my $query_file = $temp . ".fasta";
	my $extract_outfile = $temp . ".extract" . ".fasta";
	foreach my $genome(@Btype_genome)
	{
		my $Blast_outfile = $temp . ".$genome" . ".out.txt";
		my $genome_file = $genome . ".txt";
		print "Calling tBLASTx on $query_file with $genome genome \n";
		system("~/Desktop/ncbi-blast-2.7.1+/bin/tblastx -query query/$query_file -db BType_phage_genome/$genome -out $Blast_outfile -outfmt 6");
		my ($db_sub_name, $db_start, $db_end,$e_val_given) = split_line($Blast_outfile); # called a function to read Blast output file to get the top hit
		if ($e_val_given <= 1e-10) # Checking significant hits
		{
		open (my $fh5, ">>", "$working_dir/$extract_outfile");
		{
		print $fh6 "$genome\t$temp\t$db_sub_name\t$db_start\t$db_end\t$e_val_given\tSignificant Hit\n";  #printing top hit on the log file
		print "$db_sub_name\t$genome_file\n";
		my $extracted_seq = extract_seq($genome_file,$db_sub_name);# Extracting the top hit
		print $fh5 ">$genome$temp\n$extracted_seq\n"; # print the extracted sequence to the output file
		}
		}
		else
		{
		print $fh6 "$genome\t$temp\t$db_sub_name\t$db_start\t$db_end\t$e_val_given\tNot significant\n"; #printing top hit on the log file
		}
		unlink("$working_dir/$Blast_outfile") or die "couldn't delete this file!";	#deleting blast output
	}
}
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

sub read_fasta_file  #This function is to read fasta files and send the sequence without any label
{
	my $file = shift;
	chomp($file);
	my @content = ();
	open(my $fh1,"<","$working_dir/$file") or die "Can't open $file\n"; #opens the input file 
	{
	while(my $line = <$fh1>)
	{
	chomp($line);
	if($line=~ /^>/)
		{
		next;
		}
	else
		{
		push(@content,$line);
		}
	}
	}
		return join("",@content);
}
sub split_line # This is to split the blast hit and return the best hit 
{
my $cur_file = shift;
my $fp;
	open($fp,"$cur_file");
	my @fp = <$fp>;
	my $line = join('',@fp);
	$line =~ tr/\n/\t/;
	my @parse = split(/\t/,$line);
    my $len= scalar(@parse);
	my @query_name = ();
	my @subject_name = ();
	my @per_idn = ();
	my @align_len = ();
	my @num_mismatch = ();
	my @num_gap_pos = ();
	my @query_start = ();
	my @query_end = ();
	my @subject_start = ();
	my @subject_end = ();
	my @e_val = ();
	my @bit_score = ();
	

	for(my $i=0;$i<=12;$i++)
	{
		for(my $j=$i;$j<=$len;($j=($j+12)))
		{
		if($i==0)
		{
		push(@query_name,$parse[$j]);
		}
		elsif($i==1)
		{
		push(@subject_name,$parse[$j]);
		}
		elsif($i==2)
		{
		push(@per_idn,$parse[$j]);
		}
		elsif($i==3)
		{
		push(@align_len,$parse[$j]);
		}
		elsif($i==4)
		{
		push(@num_mismatch,$parse[$j]);
		}
		elsif($i==5)
		{
		push(@num_gap_pos,$parse[$j]);
		}
		elsif($i==6)
		{
		push(@query_start,$parse[$j]);
		}
		elsif($i==7)
		{
		push(@query_end,$parse[$j]);
		}
		elsif($i==8)
		{
		push(@subject_start,$parse[$j]);
		}
		elsif($i==9)
		{
		push(@subject_end,$parse[$j]);
		}
		elsif($i==10)
		{
		push(@e_val,$parse[$j]);
		}
		elsif($i==11)
		{
		push(@bit_score,$parse[$j]);
		}
		}
	}
	my $db_sub_name = $subject_name[0];
	my $db_start = $subject_start[0];
	my $db_end = $subject_end[0];
	my $e_val_given = $e_val[0];
	
	return $db_sub_name, $db_start, $db_end,$e_val_given;

}

sub extract_seq # This subroutine is to extract the ORF from the top blast hit
{
my $file = $_[0];
my $scaffold = $_[1];
my $extract_seq = '';
system("sed -n \"/$scaffold/,/>/p\"  $working_dir/BType_phage_genome/$file >  $working_dir/temp.txt");
system("egrep -v \">|>\" $working_dir/temp.txt > $working_dir/temp.final.txt");
my $seq = read_fasta_file("temp.final.txt");
unlink("$working_dir/temp.txt") or die "couldn't delete this file!";
unlink("$working_dir/temp.final.txt") or die "couldn't delete this file!";
return $seq;
}