#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($input_folder, $flag_help);

my $help_message = "
This script performs cleaning of intermediate data as final preparation for oncoIMPACT.

Usage:
	merge_and_clean.pl [OPTIONS]

Options:
	--input = full path to folder containing sample folders of intermediate results *
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
	"input=s" => \$input_folder,
	"help"    => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

opendir( DIR, $input_folder );
my @the_TCGA_DIR = readdir(DIR);
close(DIR);

system("mkdir $input_folder/COMPLETE_SAMPLES")
  unless -d "$input_folder/COMPLETE_SAMPLES";
system("mkdir $input_folder/INCOMPLETE_SAMPLES")
  unless -d "$input_folder/INCOMPLETE_SAMPLES";

foreach my $dir (@the_TCGA_DIR) {
	if ( index( $dir, "TCGA" ) > -1 ) {
		$sample_dir = "$input_folder/$dir";
		$cnv_file   = "$sample_dir/CNV_Data.txt";
		$snv_file   = "$sample_dir/SNP_Data.txt";
		$expr_file  = "$sample_dir/EXPR_Data.txt";
		if (   -e $cnv_file
			&& -e $snv_file
			&& -e $expr_file )
		{

			$out_file_name = "$sample_dir/Genelist_Status.txt";
			system("cat $cnv_file $snv_file $expr_file > $out_file_name");
			system("mv $sample_dir $input_folder/COMPLETE_SAMPLES/");
		} elsif (   -e $cnv_file
			|| -e $snv_file
			|| -e $expr_file )
		{
			system("mv $sample_dir $input_folder/INCOMPLETE_SAMPLES/");
		}
	}
}
