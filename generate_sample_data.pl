#!/usr/bin/perl
use warnings;
use Getopt::Long;

my (
	$dir, $outDir, $flag_help, $flag_debug
);

my $help_message = "
This script prepares the expression data downloaded from TCGA for oncoIMPACT.

Usage:
	generate_sample_data.pl [OPTIONS]

Options:
	--inDir = full path to input folder *
	--outDir = full path to output folder *
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
	"inDir=s" => \$dir,
	"outDir=s" => \$outDir,
	"debug"         => \$flag_debug,
	"help"          => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

if ($flag_debug) {
	print STDERR "Input parameters:\n";
	print STDERR "input directory: $dir\n";
	print STDERR "destination: $outDir\n";
}

opendir( DIR, $dir );
my @samples = sort(grep(!/^(\.|\.\.)$/, readdir(DIR)));
close(DIR);

my ( %cnv, %exp, %snp ) = ();

for ( my $sample = 0 ; $sample < @samples ; $sample++ ) {
	$inDir = "$dir/$samples[$sample]";

	# cnv
	open( FILE, "$inDir/CNV_Data.txt" );
	print STDERR "[cnv] opening file:$inDir/CNV_Data.txt\n" if $debug;
	while (<FILE>) {
		chomp( @temp = split( /\t/, $_ ) );
		@temp = split( /_/, $temp[0] );
		if ( !exists $cnv{ $temp[0] } ) {
			my @gene;
			$cnv{ $temp[0] } = \@gene;
			for ( $j = 0 ; $j < @samples ; $j++ ) {
				push( @{ $cnv{ $temp[0] } }, 0 );
				print STDERR "[cnv] pushing cnv[$temp[0]], 0\n" if $debug;
			}
		}

		if ( $temp[1] eq "AMPL" ) {
			$cnv{ $temp[0] }->[$sample] = 1;
			print STDERR "[cnv] pushing cnv[$temp[0]], 1\n" if $debug;
		}
		else {
			$cnv{ $temp[0] }->[$sample] = -1;
			print STDERR "[cnv] pushing cnv[$temp[0]], -1\n" if $debug;
		}
	}
	close(FILE);

	# exp
	open( FILE, "$inDir/EXPR_Data.txt" );
	print STDERR "[exp] opening file:$inDir/EXPR_Data.txt\n" if $debug;
	while (<FILE>) {
		chomp( @temp = split( /\t/, $_ ) );
		@temp2 = split( /_/, $temp[0] );
		if ( !exists $exp{ $temp2[0] } ) {
			my @gene;
			$exp{ $temp2[0] } = \@gene;
			for ( $j = 0 ; $j < @samples ; $j++ ) {
				push( @{ $exp{ $temp2[0] } }, 0 );
				print STDERR "[exp] pushing exp[$temp2[0]], 0\n" if $debug;
			}
		}
		$exp{ $temp2[0] }->[$sample] = $temp[1];
		print STDERR "[exp] pushing exp[$temp2[0]], $temp[1]\n" if $debug;
	}
	close(FILE);

	# snp
	open( FILE, "$inDir/SNP_Data.txt" );
	print STDERR "[snp] opening file:$inDir/SNP_Data.txt\n" if $debug;
	while (<FILE>) {
		chomp( @temp = split( /\t/, $_ ) );
		@temp = split( /_/, $temp[0] );
		if ( !exists $snp{ $temp[0] } ) {
			my @gene;
			$snp{ $temp[0] } = \@gene;
			for ( $j = 0 ; $j < @samples ; $j++ ) {
				push( @{ $snp{ $temp[0] } }, 0 );
				print STDERR "[snp] pushing snp[$temp[0]], 0\n" if $debug;
			}
		}
		$snp{ $temp[0] }->[$sample] = 1;
		print STDERR "[snp] pushing snp[$temp[0]], 1\n" if $debug;
	}
	close(FILE);
}
print "size of cnv:" . (keys %cnv) . "\n";
print "size of exp:" . (keys %exp) . "\n";
print "size of snp:" . (keys %snp) . "\n";

# CNV
open( OUT, "> $outDir/CNV.txt" );
print OUT "Genes";
for ( my $sample = 0 ; $sample < @samples ; $sample++ ) {
	print OUT "\t$samples[$sample]";
}
print OUT "\n";
foreach $gene ( keys %cnv ) {
	print OUT "$gene\t";
	print OUT join "\t", @{ $cnv{$gene} };
	print OUT "\n";
}
close OUT;

# SNP
open( OUT, "> $outDir/SNP.txt" );
print OUT "Genes";
for ( my $sample = 0 ; $sample < @samples ; $sample++ ) {
	print OUT "\t$samples[$sample]";
}
print OUT "\n";
foreach $gene ( keys %snp ) {
	print OUT "$gene\t";
	print OUT join "\t", @{ $snp{$gene} };
	print OUT "\n";
}
close OUT;

# EXP
open( OUT, "> $outDir/EXPR.txt" );
print OUT "Genes";
for ( my $sample = 0 ; $sample < @samples ; $sample++ ) {
	print OUT "\t$samples[$sample]";
}
print OUT "\n";
foreach $gene ( keys %exp ) {
	print OUT "$gene\t";
	print OUT join "\t", @{ $exp{$gene} };
	print OUT "\n";
}
close OUT;
