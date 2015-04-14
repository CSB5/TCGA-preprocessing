#!/usr/bin/perl

use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my ($in_file, $cancerCensus_file, $output_folder, $flag_help);

my $help_message = "
This script generates sample-based gene mutation information.

Usage:
	generate_sample_based_gene_mutation_info.pl [OPTIONS]

Options:
	--in = gene_mutation_frequency.txt file *
	--cc = cancer census database *
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
	"in=s" => \$in_file,
	"cc=s"	=> \$cancerCensus_file,
	"outDir=s" => \$output_folder,
	"help"    => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

## Read cancer census file ##
my %cancerCensus = ();
open(FILE, $cancerCensus_file);
<FILE>;	# skip header
my @temp;
while(<FILE>){
	chomp(@temp = split(/\t/, $_));
	$cancerCensus{$temp[0]} = "";
}
close(FILE);

## Process mutation file
my (@samples, @mutation, $sample, $outFile, %ht_all, %ht_cc) = ();
open(FILE, $in_file);
<FILE>;	#skip header
while(<FILE>){
	chomp(@temp = split(/\t/, $_));
	@samples = split(/;/, $temp[2]);
	foreach $sample (@samples){
		@mutation = split(/:/, $sample);
		unless(exists $ht_all{$mutation[0]}){
			my @mt = ();
			$ht_all{$mutation[0]} = \@mt;
		}
		push( @{ $ht_all{$mutation[0]} }, $mutation[1]);
		if(exists $cancerCensus{$temp[0]}){
			unless(exists $ht_cc{$mutation[0]}){
			my @mt = ();
			$ht_cc{$mutation[0]} = \@mt;
			}
			push( @{ $ht_cc{$mutation[0]} }, $mutation[1]);
		}
	}
}
close(FILE);

## Generate results file
# Mutations per sample
$outFile = "$output_folder/mutations_per_sample-all_genes.dat";
open(OUT, ">$outFile");
print OUT "Sample\tMutations\t%_of_SNVs\n";
my ($key, @snv);
foreach $key (sort keys %ht_all){	
	@snv = grep (/MUT/, @{ $ht_all{$key} });
	print OUT $key . "\t" . scalar @{ $ht_all{$key} } . "\t" . (scalar @snv / scalar @{ $ht_all{$key} }) . "\n";
}
close(OUT);

$outFile = "$output_folder/mutations_per_sample-cc_genes.dat";
open(OUT, ">$outFile");
print OUT "Sample\tMutations\t%_of_SNVs\n";
foreach $key (sort keys %ht_cc){
	@snv = grep (/MUT/, @{ $ht_cc{$key} });
	print OUT $key . "\t" . scalar @{ $ht_cc{$key} } . "\t" . (scalar @snv / scalar @{ $ht_cc{$key} }) . "\n";
}
close(OUT);