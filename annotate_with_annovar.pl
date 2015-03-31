#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($destination, $maf_dir, $flag_debug, $flag_help, $command);
my $maf_database = "/mnt/genomeDB/annovar/humandb";

my $help_message = "
This script performs ANNOVAR annotation on TCGA MAF file.

Usage:
	annotate_with_annovar.pl [OPTIONS]

Options:
	--input = full path to folder containing SNP data *
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
	"input=s"       => \$maf_dir,
	"destination=s" => \$destination,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "[ANNOVAR] Input parameters:\n";
	print STDERR "[ANNOVAR] SNV folder: $maf_dir\n";
	print STDERR "[ANNOVAR] destination: $destination\n";
}

my @files = glob "$maf_dir/*.maf";

## convert MAF to VCF
$command = "cat $files[0] | cut -f 5,7,11,13,16 | awk '{if(\$3 != \"-\" && \$4 != \"-\") print \$1\"\\t\"\$2\"\\t.\\t\"\$3\"\\t\"\$4\"\\t88\\tPASS\\t\"\$5}' | grep -v Chrom > $destination/annovar.vcf";
print STDERR "$command\n" if $flag_debug;
system($command);

## Annotate VCF
$command = "table_annovar.pl $destination/annovar.vcf $maf_database -buildver hg19 -out $destination/annovar -remove -protocol knownGene,1000g2014oct_all,snp142,ljb26_all -operation g,f,f,f -nastring . -vcfinput";
print STDERR "$command\n" if $flag_debug;
system($command);

## Cleanup
$command = "rm -rf $destination/annovar.vcf $destination/annovar.avinput";
print STDERR "$command\n" if $flag_debug;
system($command);
