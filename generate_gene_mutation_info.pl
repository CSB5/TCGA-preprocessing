#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use 5.010;

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

my (%all, %allInfo, %cnv, %cnvInfo, %snp, %snpInfo, $geneName, $mutationType, @temp, @mutation);

foreach my $dir (@complete_DIR) {
	open( FILE, "$input_folder/$dir/Genelist_Status.txt");
	while(<FILE>){
		chomp(@temp = split(/\t/, $_));
		chomp(@mutation = split(/_/, $temp[0]));
		$geneName = $mutation[0];
		$mutationType = $mutation[1];

		given($mutationType){
			when(["DEL", "AMPL"]){
				if(! exists $cnv{$geneName}){
					my @sampleList = ();
					my %mutationInfo = ();
					$cnv{$geneName} = \@sampleList;
					$cnvInfo{$geneName} = \%mutationInfo;
				}

				push( @{$cnv{$geneName}}, $dir);
				if(! exists $cnvInfo{$geneName}->{$dir}){
					my @mutationList = ();
					$cnvInfo{$geneName}->{$dir} = \@mutationList;
				}
				push( @{$cnvInfo{$geneName}->{$dir}}, $mutationType );
				continue;
			}
			when("MUT"){
				if(! exists $snp{$geneName}){
					my @sampleList = ();
					my %mutationInfo = ();
					$snp{$geneName} = \@sampleList;
					$snpInfo{$geneName} = \%mutationInfo;
				}

				push( @{$snp{$geneName}}, $dir);
				if(! exists $snpInfo{$geneName}->{$dir}){
					my @mutationList = ();
					$snpInfo{$geneName}->{$dir} = \@mutationList;
				}
				push( @{$snpInfo{$geneName}->{$dir}}, $mutationType );
				continue;
			}
			when(["DEL","AMPL","MUT"]){
				if(! exists $all{$geneName}){
					my @sampleList = ();
					my %mutationInfo = ();
					$all{$geneName} = \@sampleList;
					$allInfo{$geneName} = \%mutationInfo;
				}

				push( @{$all{$geneName}}, $dir);
				if(! exists $allInfo{$geneName}->{$dir}){
					my @mutationList = ();
					$allInfo{$geneName}->{$dir} = \@mutationList;
				}
				push( @{$allInfo{$geneName}->{$dir}}, $mutationType );
			}
			default{
				# do nothing
			}
		}
	}
	close(FILE);
}

# print combined mutation info
open (OUT, "> $output_folder/gene_mutation_frequency.txt");
print OUT "GeneName\tFrequency\tSampleInfo\n";
foreach my $gene ( sort keys %all ) {
	print OUT "$gene\t";	# print gene name
	print OUT (scalar uniq(@{$all{$gene}}) / scalar @complete_DIR) . "\t";	# print frequency
	foreach my $sample ( sort keys %{ $allInfo{$gene} } ) {
		print OUT "$sample:" . join(",", @{$allInfo{$gene}->{$sample}}) . ";";
	}
	print OUT "\n";
}
close(OUT);

# print CNV mutation info
open (OUT, "> $output_folder/cnv_mutation_frequency.txt");
print OUT "GeneName\tFrequency\tSampleInfo\n";
foreach my $gene ( sort keys %cnv ) {
	print OUT "$gene\t";	# print gene name
	print OUT (scalar uniq(@{$cnv{$gene}}) / scalar @complete_DIR) . "\t";	# print frequency
	foreach my $sample ( sort keys %{ $cnvInfo{$gene} } ) {
		print OUT "$sample:" . join(",", @{$cnvInfo{$gene}->{$sample}}) . ";";
	}
	print OUT "\n";
}
close(OUT);

# print SNP mutation info
open (OUT, "> $output_folder/snp_mutation_frequency.txt");
print OUT "GeneName\tFrequency\tSampleInfo\n";
foreach my $gene ( sort keys %snp ) {
	print OUT "$gene\t";	# print gene name
	print OUT (scalar uniq(@{$snp{$gene}}) / scalar @complete_DIR) . "\t";	# print frequency
	foreach my $sample ( sort keys %{ $snpInfo{$gene} } ) {
		print OUT "$sample:" . join(",", @{$snpInfo{$gene}->{$sample}}) . ";";
	}
	print OUT "\n";
}
close(OUT);
