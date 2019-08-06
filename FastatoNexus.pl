#/usr/bin/perl
use strict;
use warnings;

my $dir = getcwd();


my ($input_query) = shift;
print "The arguments passed are :\n";
foreach my $temp(@ARGV) #prints arguments
{
print"$temp\n";
}

my @file_prefix =();

open(my $fh1,"<","$dir/$input_query") or die "Can't open $input_query\n"; #creates query array from file provided
{
while(my $line = <$fh1>)
{
chomp($line);
push(@file_prefix,$line);
}
}

foreach my $temp(@file_prefix)
{
	my $file_name = $temp . ".fasta";
	my @dat = ReadInFASTA($file_name);
	my @taxon = GetSeqName(@dat);
	my @seqdata = GetSeqDat(@dat);
	my $out_file = $temp . ".nex";
	print"$file_name\n";
	my $num_taxa = scalar(@taxon);
	my $numsites = length($seqdata[0]);
open(my $fh1,">>","$dir/$out_file");
	{
		
	print $fh1 "#NEXUS\n\n";
	print $fh1 "Begin data;\n";
	print $fh1 "\tDimensions ntax=$num_taxa nchar=$numsites;\n";
	print $fh1 "\tFormat datatype=protein gap=-;\n";
	print $fh1 "\tMatrix\n";
		#print 
		for(my $i=0;$i <=$#taxon; $i++)
		{
		#print"$taxon[$i]\n";
		print $fh1 "$taxon[$i]\t\t$seqdata[$i]\n";
		#last;
		}
	print $fh1 "\t;\n";
	print $fh1 "End;\n";
	}
	
	
}
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();
    my @seqName = ();
    my @seqDat = ();

    open (INFILE,"$fileloc/$infile") || die "Can't open $infile\n";

    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
	    s/\s+//g;                  # get rid of any spaces
            next if (/^$/);            # skip empty line
            s/[uU]/T/g;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . uc($_);
        }

	# checking no occurence of internal separator $sep.
	die ("ERROR: \"$sep\" is an internal separator.  Line $. of " .
	     "the input FASTA file contains this charcter. Make sure this " . 
	     "separator character is not used in your data file or modify " .
	     "variable \$sep in this script to some other character.\n")
	    if (/$sep/);

    }
    close(INFILE);

    foreach my $i (0..$#seqName) {
	$result[$i] = $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}
sub GetSeqDat {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
	@line = split (/$sep/, $i);
	push @result, $line[1];
    }

    return (@result)
}

sub GetSeqName {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
	@line = split (/$sep/, $i);
	push @result, $line[0];
    }
    return (@result)
}
