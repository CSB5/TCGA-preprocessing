#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ( $inFile, $outFile, $flag_debug, $flag_help );

my $help_message = "
This script trims the sample names in column headers.

Usage:
	trim_sample_names.pl [OPTIONS]

Options:
	--in = path to normalized counts matrix *
	--out = path to output file *
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
	"in=s"	 => \$inFile,
	"out=s"  => \$outFile,
	"debug"  => \$flag_debug,
	"help"   => \$flag_help
  )
  or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "[DEBUG] Input matrix: $inFile\n";
	print STDERR "[DEBUG] Output file: $outFile\n";
}

my ($header, @headers, @sample, $newHeader);

# Generating output file
open( OUT,  "> $outFile" );
open( IN, $inFile );
$newHeader = "";
chomp(@headers = split(/\t/, <IN>));
foreach $header (@headers){
	@sample = split(/-/, $header);
	$newHeader = $newHeader . "$sample[0]-$sample[1]-$sample[2]-$sample[3]" . "\t";
}
chop($newHeader);
print OUT $newHeader . "\n";
print OUT <IN>;
close(OUT);
close(IN);
