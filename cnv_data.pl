#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::MoreUtils;

my (
	$file_in, $destination, $flag_debug, $flag_help,   
	@header, %samples, @temp, $sampleID, $gene, $score, @tempID, $data		
);

my $help_message = "
This script prepares the CNV data downloaded from GDAC for oncoIMPACT.

Usage:
	cnv_data.pl [OPTIONS]

Options:
	--in = path to thresholded data file provided by GDAC *
	--destination = full path to output folder *
	--debug: prints trace to STDERR
	--help : prints this message 
	
* indicates required parameters	


Version:
	1.0

Author:
	Burton Chia - chiakhb\@gis.a-star.edu.sg
	Denis Bertrandd - bertrandd\@gis.a-star.edu.sg\n";

if ( @ARGV == 0 ) {
	print $help_message;
	exit 0;
}

GetOptions(
	"in=s" 			 => \$file_in,
	"destination=s"  => \$destination,
	"debug"          => \$flag_debug,
	"help"           => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "[CNV] input_folder: $file_in\n";
	print STDERR "[CNV] destination: $destination\n";
}

print "Processing input file. Please wait...";
# Process input file
open(IN, $file_in) or die ("File $file_in does not exists.\n");

# Read header
chomp(@header = split(/\t/, <IN>));

# Read remaining lines
while(<IN>){
	chomp(@temp = split(/\t/, $_));
	$gene = $temp[0];
	next if (index($gene, '|') != -1); # there are 159 out of 23547 genes with multiple entries (isoforms). As we need the list to be unique, these genes are dropped. 2 of these genes (CRLF2,P2RY8) are found in the cancer census database
	for (my $i=3; $i < @temp; $i++){
		$score = $temp[$i];
		next if($score == 0);
		@tempID = split(/-/, $header[$i]);
		$sampleID = "$tempID[0]-$tempID[1]-$tempID[2]-$tempID[3]";
		unless(exists $samples{$sampleID}){
			my @array = ();
			$samples{$sampleID} = \@array;
		}
		if($score > 1){
			$data = $gene . "_AMPL\t" . $score;
			push(@{$samples{$sampleID}}, $data);
		} elsif($score < -1){
			$data = $gene . "_DEL\t" . $score;
			push(@{$samples{$sampleID}}, $data);
		}
	}
}
close(IN);
print "done.\n";

print "Generating results. Please wait...";
# Generate results
foreach my $sample ( keys %samples ) {
	print STDERR "[CNV] *** SAMPLE $destination/$sample\n" if ($flag_debug);
	mkdir("$destination/$sample") if ( !-d "$destination/$sample" );

	open( NEWF, ">$destination/$sample/CNV_Data.txt" ) or die "unable to open file hander on $destination/$sample/CNV_Data.txt\n";
	print( NEWF join( "\n", @{ $samples{$sample} } ) );
	close(NEWF);
}
print "done.\n";