#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ( $outDir, $rootDir, @samples, $flag_debug, $flag_help, $command );

my $help_message = "
This script downloads data from GDAC via FireHose.

Usage:
	download_from_gdac.pl [OPTIONS]

Options:
	--outDir: output directory *
	--samples: acronymns of cancer samples to download (eg: blca, gbm). Multiple entries allowed, separated with \",\" *
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
	"outDir=s" => \$rootDir,
	"samples=s" => \@samples,    
	"debug"    => \$flag_debug,
	"help"     => \$flag_help
  )
  or die("Error in command line arguments.\n");

@samples = split(/,/,join(',',@samples));

if ($flag_help) {
	print $help_message;
	exit 0;
}

foreach my $sample (@samples){
	print "Currently processing sample: $sample\n";
	$outDir = "$rootDir/" . uc $sample;
	$command = "mkdir -p $outDir/GISTIC2 $outDir/MUTSIGCV";
	print STDERR "$command\n" if $flag_debug;
	system($command);
	
	print "Downloading data...";
	$command = "firehose_get -tasks mutsig gistic analyses latest " . (lc $sample) . "< $rootDir/cmd";
	print STDERR "$command\n" if $flag_debug;
	system("$command >> firehose.log");
	print "done.\n";
	
	print "Unpacking data...";
	$command = "tar -zxvf analyses__*/" . uc $sample . "/*/*MutSigNozzleReportCV.Level_4*.tar.gz -C $outDir/MUTSIGCV";
	print STDERR "$command\n" if $flag_debug;
	system($command);
	$command = "tar -zxvf analyses__*/" . uc $sample . "/*/*CopyNumber_Gistic2.Level_4*.tar.gz -C $outDir/GISTIC2";
	print STDERR "$command\n" if $flag_debug;
	system($command);
	$command = "rm -rf analyses__*";
	print STDERR "$command\n" if $flag_debug;
	system($command);	
}
