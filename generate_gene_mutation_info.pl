#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my ($input_folder, $output_folder, $flag_help);

my $help_message = "
This script generates additional gene mutation information.

Usage:
	generate_gene_mutation_info.pl [OPTIONS]

Options:
	--inDir = full path to oncoIMPACT COMPLETE_SAMPLES folder *
	--outDir = full path to output folder *
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
	"inDir=s" => \$input_folder,
	"outDir=s" => \$output_folder,
	"help"    => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}


opendir( DIR, $input_folder );
my @complete_DIR = grep{!/^\./} readdir(DIR);
close(DIR);


# Generate list of samples with all three mutation data (SNV, CNV and EXP)
open(OUT, "> $output_folder/complete_samples_list.txt");
print OUT join("\n", @complete_DIR) . "\n";
close(OUT);


# Copy selected normals file
system("cp $input_folder/../RNA-SEQ/selected_normals.dat $output_folder");

my (%samples, %info, $geneName, $mutationType, @temp, @mutation);

foreach my $dir (@complete_DIR) {
	open( FILE, "$input_folder/$dir/Genelist_Status.txt");
	while(<FILE>){
		chomp(@temp = split(/\t/, $_));
		chomp(@mutation = split(/_/, $temp[0]));
		$geneName = $mutation[0];
		$mutationType = $mutation[1];
		
		next if($mutationType eq "UP" || $mutationType eq "DOWN");
		
		if(! exists $samples{$geneName}){
			my @sampleList = ();
			my %mutationInfo = ();
			$samples{$geneName} = \@sampleList;
			$info{$geneName} = \%mutationInfo;
		}
		
		push( @{$samples{$geneName}}, $dir);
		if(! exists $info{$geneName}->{$dir}){
			my @mutationList = ();
			$info{$geneName}->{$dir} = \@mutationList;
		}
		push( @{$info{$geneName}->{$dir}}, $mutationType );
	}
	close(FILE);
}

open (OUT, "> $output_folder/gene_mutation_frequency.txt");
print OUT "GeneName\tFrequency\tSampleInfo\n";
foreach my $gene ( sort keys %samples ) {
	print OUT "$gene\t";	# print gene name
	print OUT (scalar uniq(@{$samples{$gene}}) / scalar @complete_DIR) . "\t";	# print frequency
	foreach my $sample ( sort keys %{ $info{$gene} } ) {
		print OUT "$sample:" . join(",", @{$info{$gene}->{$sample}}) . ";";
	}
	print OUT "\n";
}
