#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::MoreUtils qw/uniq/;
use File::Find;

my (
	$destination, $file_in, $flag_debug, $flag_help,
	@columns, @col_parts, $sample_name, %sample_list, $folder
);

my $help_message = "
This script prepares the SNP data downloaded from GDAC for oncoIMPACT.

Usage:
	snp_data.pl [OPTIONS]

Options:
	--in = full path to MAF file *
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
	"in=s"       	=> \$file_in,
	"destination=s" => \$destination,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "[SNV] Input parameters:\n";
	print STDERR "[SNV] maf_file: $file_in\n";
	print STDERR "[SNV] destination: $destination\n";
}

print "Processing MAF file. Please wait...";
open( FILE, "$file_in" );
print STDERR "[SNP] *** READING FILE $file_in\n" if ($flag_debug);
while (<FILE>) {
	@line = split( /\t/, $_ );
	if ( $line[0] ne "Hugo_Symbol" ) {
		@col_parts = split( /-/, $line[15] );
		$sample_name = "$col_parts[0]-$col_parts[1]-$col_parts[2]-$col_parts[3]";
		if (   $line[8] ne "Translation_Start_Site"
			#&& $line[8] ne "Splice_Site"
			&& $line[8] ne "3'Flank"
			&& $line[8] ne "5'Flank"
			&& $line[8] ne "5'UTR"
			&& $line[8] ne "RNA"
			&& $line[8] ne "Silent"
			&& $line[8] ne "Splice_Region"
			&& $line[8] ne "Intron"
			 )
		{
			if ( !exists $sample_list{$sample_name} ) {
				my @sample_columns;
				$sample_list{$sample_name} = \@sample_columns;
			}
			push( @{ $sample_list{$sample_name} }, "$line[0]\_MUT" );
		}
	}
}
print "done.\n";

print "Generating results. Please wait...";
# my $sample2; #for OV specifically
foreach my $sample ( keys %sample_list ) {
	print STDERR "[SNP] *** SAMPLE $destination/$sample\n" if ($flag_debug);
	mkdir("$destination/$sample") if ( !-d "$destination/$sample" );
#	#For OV specifically
#	if (-d "$destination/${sample}A"){
#		$sample2 = "${sample}A";
#	}
#	elsif (-d "$destination/${sample}B"){
#		$sample2 = "${sample}B";
#	} 
#	elsif (-d "$destination/${sample}C"){
#		$sample2 = "${sample}C";
#	} 
#	elsif (-d "$destination/${sample}D"){
#		$sample2 = "${sample}D";
#	} else{
#		print "Error: $sample\n";
#	}

	open( NEWF, ">$destination/temp.txt" );
	print( NEWF join( "\n", @{ $sample_list{$sample} } ) );
	close(NEWF);

	system(
		"sort $destination/temp.txt | uniq > $destination/$sample/SNP_Data.txt"
#		"sort $destination/temp.txt | uniq > $destination/$sample2/SNP_Data.txt" # For OV specifically
	);

}
system("rm $destination/temp.txt");
print "done.\n";