#!/usr/bin/perl

use warnings;
use Getopt::Long;

my (
	$TCGA_NORMAL_ID_FILE, $TCGA_ID_FILE, $TCGA_RNA_SEQ_DIR,
	$outDir,              @sample_order, @sample_name_order,
	%sample_to_ID, %gene_expression, $cmp_ID, $flag_help, $flag_debug
);

my $help_message = "
This script prepares the expression data downloaded from TCGA for oncoIMPACT.

Usage:
	merge_normal_cancer_RNA_seq.pl [OPTIONS]

Options:
	--normal = path to FILE_SAMPLE_MAP_NORMAL.txt *
	--id = path to FILE_SAMPLE_MAP.txt *
	--rnaseq = path to RNAseq folder *
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
	"normal=s"      => \$TCGA_NORMAL_ID_FILE,
	"id=s"          => \$TCGA_ID_FILE,
	"rnaseq=s"      => \$TCGA_RNA_SEQ_DIR,
	"destination=s" => \$outDir,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "[EXP-MNC] Input parameters:\n";
	print STDERR "[EXP-MNC] TCGA_NORMAL_ID_FILE: $TCGA_NORMAL_ID_FILE\n";
	print STDERR "[EXP-MNC] TCGA_ID_FILE: $TCGA_ID_FILE\n";
	print STDERR "[EXP-MNC] TCGA_RNA_SEQ_DIR: $TCGA_RNA_SEQ_DIR\n";
	print STDERR "[EXP-MNC] destination: $outDir\n";
}

$cmp_ID = 0;

#The normal files
print STDERR "[EXP-MNC] *** GET SAMPLE LIST\n" if ($flag_debug);
open( FILE, $TCGA_NORMAL_ID_FILE );
while (<FILE>) {
	chop $_;
	@line = split( /\t/, $_ );
	if ( @line != 0 ) {
		$path = $TCGA_RNA_SEQ_DIR . "/" . $line[0];
		if ( index( $path, "rsem.genes.results" ) != -1 ) {
			$sample = $line[1];

			push( @sample_order,      $path );
			push( @sample_name_order, $sample );
			$sample_to_ID{$path} = $cmp_ID;
			$cmp_ID++;
		}
	}
}
close(FILE);

#The TCGA data
open( FILE, $TCGA_ID_FILE );
while (<FILE>) {
	chop $_;
	@line = split( /\t/, $_ );
	if ( @line != 0 ) {
		$path = $TCGA_RNA_SEQ_DIR . "/" . $line[0];
		if ( index( $path, "rsem.genes.results" ) != -1
			&& !exists $sample_to_ID{$path} )
		{
			$sample = $line[1];
			push( @sample_order,      $path );
			push( @sample_name_order, $sample );
			$sample_to_ID{$path} = $cmp_ID;
			$cmp_ID++;
		}
	}
}
close(FILE);

#Produce the matrix
print STDERR "[EXP-MNC] *** GET GENE LIST\n" if ($flag_debug);
for ( my $i = 0 ; $i < @sample_order ; $i++ ) {
	$file      = $sample_order[$i];
	$sample_ID = $sample_to_ID{$file};
	open( FILE, $file );
	print STDERR "[EXP-MNC] (FileNotFound) File $file not found!\n" if(! -e $file);
	print STDERR "[EXP-MNC] *** READ FILE $file $sample_ID\n" if ($flag_debug);
	$cmp_gene = 0;
	while (<FILE>) {
		chop $_;
		@line  = split( /\t/, $_ );
		$gene  = $line[0];
		$count = $line[1];
		if ( $gene ne "gene_id" ) {
			$gene = convert_name($gene);
			if ( $gene ne "?" ) {
				if ( !exists $gene_expression{$gene} ) {
					my @tab = ();
					for ( $j = 0 ; $j < @sample_order ; $j++ ) {
						$tab[$j] = 0;
					}
					$gene_expression{$gene} = \@tab;
				}
				$gene_expression{$gene}->[$sample_ID] = $count;
				$cmp_gene++;
			}
		}
	}
	close(FILE);
}

print STDERR "[EXP-MNC] *** EXPORT MATRIX\n" if ($flag_debug);
open( OUT, ">$outDir/RNA_SEQ_COUNT_MATRIX.dat" );
print OUT "" . ( join( "\t", @sample_name_order ) ) . "\n";
foreach $gene ( keys %gene_expression ) {
	$count = $gene_expression{$gene};
	print OUT $gene . "\t" . ( join( "\t", @{$count} ) ) . "\n";
}
close(OUT);

sub convert_name {
	my ($name) = @_;
	@str = split( /\|/, $name );
	return $str[0];
}

