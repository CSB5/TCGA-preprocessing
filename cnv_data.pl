#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::MoreUtils;

my (
	$input_folder,   $file_manifest,     $destination, $region_file,
	$file_extention, $clinical_data,     $up,          $down,
	$scripts_folder, $flag_debug,        $flag_help,   $size_threshold,
	@files,          @lines,             $line,        @columns,
	@col_parts,      $header,            $sample_name, %sample_list,
	$status,         %arrayID_to_sample, %female_sample, $patient_ID
);
$up             = 2;
$down           = -2;
$size_threshold = -1;
$status         = "";

my $help_message = "
This script prepares the CNV data downloaded from TCGA for oncoIMPACT.

Usage:
	cnv_data.pl [OPTIONS]

Options:
	--input = full path to input folder *
	--fileManifest = path to file manifest provided by TCGA *
	--destination = full path to output folder *
	--region = path to region file *
	--extention = suffix of CNV files *
	--clinical = path to clinical data *
	--minUp = minimum value to consider region up-regulated [default: 2]
	--minDown = minimum value to consider region down-regulated [default: -2] 
	--scripts = full path to scripts folder *
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
	"input=s"        => \$input_folder,
	"fileManifest=s" => \$file_manifest,
	"destination=s"  => \$destination,
	"region=s"       => \$region_file,
	"extention=s"    => \$file_extention,
	"clinical=s"     => \$clinical_data,
	"minUp:f"        => \$up,
	"minDown:f"      => \$down,
	"scripts=s"      => \$scripts_folder,
	"debug"          => \$flag_debug,
	"help"           => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "[CNV] Input parameters:\n";
	print STDERR "[CNV] input_folder: $input_folder\n";
	print STDERR "[CNV] file_manifest: $file_manifest\n";
	print STDERR "[CNV] destination: $destination\n";
	print STDERR "[CNV] region_file: $region_file\n";
	print STDERR "[CNV] file_extension: $file_extention\n";
	print STDERR "[CNV] clinical_data: $clinical_data\n";
	print STDERR "[CNV] up: $up\n";
	print STDERR "[CNV] down: $down\n";
	print STDERR "[CNV] scripts_folder: $scripts_folder\n";
}

#The sample info
%female_sample = ();
open( FILE, $clinical_data ) or die "Clinical_data \"$clinical_data\" not found";
@temp = split( /\s+/, <FILE> );
my $genderIndex = List::MoreUtils::first_index { $_ eq "gender" } @temp;
while (<FILE>) {
	@line   = split( /\t/, $_ );
	$sample = $line[0];
	$sex    = $line[$genderIndex];
	if ( $sex eq "FEMALE" ) {
		$female_sample{$sample} = 1;
	}
}
close(FILE);
print STDERR "[CNV] *** clinical_data reading completed\n" if ($flag_debug);

open( FILE, $file_manifest );
%arrayID_to_sample = ();
while (<FILE>) {
	chop $_;
	@lines = split( /\t/, $_ );
	if ( $lines[0] eq "CNV (SNP Array)" ) {
		@tmp = split( /\-/, $lines[5] );
		$sample_ID = $tmp[0] . "-" . $tmp[1] . "-" . $tmp[2] . "-" . $tmp[3];
		$arrayID_to_sample{ $lines[6] } = $sample_ID;
		print STDERR "[CNV] (arrayID_to_sample) $lines[6] :: $sample_ID\n"
		  if ($flag_debug);
	}
}
close(FILE);
print STDERR "[CNV] *** file_manifest reading completed\n" if ($flag_debug);

my $nb_removed_del = 0;
my $nb_new_ampl    = 0;

opendir( DIR, $input_folder );
@files = readdir(DIR);
foreach my $file (@files) {
	next if ( $file =~ m/^\./ );
	if ( index( $file, $file_extention ) > -1 ) {
		print STDERR "[CNV] *** READING FILE $input_folder/$file\n"
		  if ($flag_debug);
		open( FILE, "$input_folder/$file" );
		while (<FILE>) {
			$line = $_;
			@columns = split( /\t/, $line );
			if ( $columns[0] eq "Sample" ) {
				$header = $line;
				next();
			}

			$sample_name = $arrayID_to_sample{ $columns[0] . $file_extention };
			@tmp = split( /\-/, $sample_name );
			$patient_ID = $tmp[0] . "-" . $tmp[1] . "-" . $tmp[2];
			print STDERR
"[CNV] *** Sample name: $columns[0]$file_extention |$sample_name|\n"
			  if ($flag_debug);
			my $seg_mean = $columns[-1];
			my $chrom    = $columns[1];
			chomp($seg_mean);
			next() if ( index( $seg_mean, "NA" ) > -1 );
			if ( $chrom ne "Y" && $chrom ne "X" ) {

				if ( ( $seg_mean < $down ) || ( $seg_mean > $up ) ) {
					if ( !exists $sample_list{$sample_name} ) {
						my @sample_columns;
						$sample_list{$sample_name} = \@sample_columns;
					}
					push( @{ $sample_list{$sample_name} }, $line );
					print STDERR
					  "[CNV] === Pushing to array:${sample_name}::$line\n"
					  if ($flag_debug);
				}
			}
			else {
				if ( exists $female_sample{$patient_ID} ) {
					if ( $chrom eq "X" && ( $seg_mean < $down )
						|| ( $seg_mean > $up ) )
					{
						push( @{ $sample_list{$sample_name} }, $line );
					}
				}
				else {
					if (   ( $seg_mean + 0.5 < $down )
						|| ( $seg_mean + 0.5 > $up ) )
					{
						if ( $seg_mean + 0.5 > $up && $seg_mean < $up ) {
							$nb_new_ampl++;
							print STDERR
							  " ********************** NEW AMPL $line\n"
							  if ($flag_debug);
						}
						push( @{ $sample_list{$sample_name} }, $line );
					}
					if ( $seg_mean + 0.5 > $down && $seg_mean < $down ) {
						$nb_removed_del++;
						print STDERR " *********** REMOVE DEL $line\n"
						  if ($flag_debug);
					}
				}
			}
		}
		close(FILE);
	}
}

print STDERR "[CNV] Moving to next step\n" if ($flag_debug);

foreach my $sample ( keys %sample_list ) {
	print STDERR "[CNV] *** Sample:$sample\n" if ($flag_debug);
	system("mkdir -p $destination/$sample")
	  unless ( -d "$destination/$sample" );
	print STDERR "[CNV] $destination/$sample/CNV_genes.txt" if ($flag_debug);
	open( NEWF, ">", "$destination/$sample/CNV_genes.txt" );
	print( NEWF $header . join( "", @{ $sample_list{$sample} } ) );
	close(NEWF);

	system(
"perl $scripts_folder/CNV_reg2gene_conv.pl $destination/$sample/CNV_genes.txt $region_file $destination/$sample/CNV_Data.txt $size_threshold"
	);
	print STDERR
"[CNV] perl $scripts_folder/CNV_reg2gene_conv.pl $destination/$sample/CNV_genes.txt $region_file $destination/$sample/CNV_Data.txt $size_threshold"
	  if ($flag_debug);
}
