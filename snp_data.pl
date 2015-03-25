#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::MoreUtils qw/uniq/;
use File::Find;

my (
	$destination, @snv_folders, $flag_debug, $flag_help,
	@columns, @col_parts, $sample_name, %sample_list, $folder
);

my $help_message = "
This script prepares the SNP data downloaded from TCGA for oncoIMPACT.

Usage:
	snp_data.pl [OPTIONS]

Options:
	--input = full path to folder(s) containing SNP data *
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
	"input=s"       => \@snv_folders,
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
	print STDERR "[SNV] snv_file(s): @snv_folders\n";
	print STDERR "[SNV] destination: $destination\n";
}

foreach $folder (@snv_folders) {
	opendir( DIR, $folder );
	@files = readdir(DIR);
	foreach $file (@files){
		next if ( $file =~ m/^\./);
		open( FILE, "$folder/$file" );
		print STDERR "[SNP] *** READING FILE $folder/$file\n" if ($flag_debug);
		while (<FILE>) {
			@line = split( /\t/, $_ );
			if ( $line[0] ne "Hugo_Symbol" ) {
				@col_parts = split( /-/, $line[15] );
				$sample_name = "$col_parts[0]-$col_parts[1]-$col_parts[2]-$col_parts[3]";
				if (   $line[8] ne "Translation_Start_Site"
					&& $line[8] ne "Splice_Site"
					&& $line[8] ne "3'Flank"
					&& $line[8] ne "5'Flank"
					&& $line[8] ne "5'UTR"
					&& $line[8] ne "RNA"
					&& $line[8] ne "Silent"
					&& $line[8] ne "Splice_Region" )
				{
	
					if ( !exists $sample_list{$sample_name} ) {
						my @sample_columns;
						$sample_list{$sample_name} = \@sample_columns;
					}
					push( @{ $sample_list{$sample_name} }, "$line[0]\_MUT" );
				}
			}
		}
	}	
}

foreach my $sample ( keys %sample_list ) {
	print STDERR "[SNP] *** SAMPLE $destination/$sample\n" if ($flag_debug);
	mkdir("$destination/$sample") if ( !-d "$destination/$sample" );

	open( NEWF, ">$destination/temp.txt" );
	print( NEWF join( "\n", @{ $sample_list{$sample} } ) );
	close(NEWF);

	system(
		"sort $destination/temp.txt | uniq > $destination/$sample/SNP_Data.txt"
	);

}

system("rm $destination/temp.txt");
