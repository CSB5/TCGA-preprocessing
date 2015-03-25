#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ( $configFile, $flag_debug, $flag_help, %config, @queue, $command );
my $qsub =
"qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -v PATH,PERL5LIB,R_LIBS_SITE,MOSEKLM_LICENSE_FILE,AUGUSTUS_CONFIG_PATH,CLASSPATH";

my $help_message = "
This script runs DEanalysis.

Usage:
	run_DEanalysis.pl [OPTIONS]

Options:
	--config = path to config file *
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
	"config=s" => \$configFile,
	"debug"    => \$flag_debug,
	"help"     => \$flag_help
  )
  or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

# Read config file
print "Reading config file. Please wait...";
Config::Simple->import_from( $configFile, \%config );
print "done.\n";

system("mkdir $config{'general.logsDir'}")
  unless ( -d "$config{'general.logsDir'}" );
open( TRACE,
	">$config{'general.logsDir'}/$config{'general.disease'}.error.log" );
@queue = ();
my $runtime = $config{'cluster.runtime'} * 60 * 60;

my (
	$expTempDir, $noNormals, $noSamples,      $i,
	$datTable,   $outPrefix, $message,        @currentJobs,
	$jobList,    $jobID,     @refreshedArray, $x,
	$y
  )
  = ();
$expTempDir = "$config{'general.analysisDir'}/RNA-SEQ";
unless ( -e $expTempDir ) {
	$command = "mkdir -p $expTempDir";
	system($command);
	print TRACE "[EXP] $command\n" if ($flag_debug);
}
unless ( -s "$expTempDir/EXCLUDED_NORMALS.txt") {
	$command = "cut -f2 $config{'general.dataDir'}/$config{'exp.sdrf'} | sort | uniq | grep -f $config{'general.dataDir'}/$config{'exp.normalList'} - | grep -vf $config{'exp.selectedNormals'} - > $expTempDir/EXCLUDED_NORMALS.txt";
	system($command);
	print TRACE "[EXP] $command\n" if ($flag_debug);
}
unless (-s "$expTempDir/EXCLUDED_SAMPLES.txt") {
	$command = "cat $expTempDir/EXCLUDED_NORMALS.txt $config{'exp.excludedTumors'} > $expTempDir/EXCLUDED_SAMPLES.txt";
	system($command);
	print TRACE "[EXP] $command\n" if ($flag_debug);
}
unless ( -s "$expTempDir/FILE_SAMPLE_MAP.txt" ) {
	$command =
"cut -f2,22 $config{'general.dataDir'}/$config{'exp.sdrf'} | awk '{print \$2\"\\t\"\$1}' | grep \"rsem.genes.results\" | sort -k2,2 | grep -vf $expTempDir/EXCLUDED_SAMPLES.txt - > $expTempDir/FILE_SAMPLE_MAP.txt";
	system($command);
	print TRACE "[EXP] $command\n" if ($flag_debug);
}
unless ( -s "$expTempDir/FILE_SAMPLE_MAP_NORMAL.txt" ) {
	$command =
"grep -f $config{'exp.selectedNormals'} $expTempDir/FILE_SAMPLE_MAP.txt | sort -k2,2 > $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt";
	system($command);
	print TRACE "[EXP] $command\n" if ($flag_debug);
}
unless ( -s "$expTempDir/RNA_SEQ_COUNT_MATRIX.dat" ) {
	$command =
"perl $config{'general.scriptsDir'}/merge_normal_cancer_RNA_seq.pl --normal $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt --id $expTempDir/FILE_SAMPLE_MAP.txt --rnaseq $config{'general.dataDir'}/$config{'exp.rnaseqDir'} --destination $expTempDir";
	system($command);
	print TRACE "[EXP] $command\n" if ($flag_debug);
}
$command = "wc -l < $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt | cut -f 1";
chomp( $noNormals = `$command` );
$command = "wc -l < $expTempDir/FILE_SAMPLE_MAP.txt | cut -f 1";
chomp( $noSamples = `$command` );
$datTable = "$expTempDir/RNA_SEQ_COUNT_MATRIX.dat";

print TRACE "[EXP] Num of normals: $noNormals\n" if ($flag_debug);
print TRACE "[EXP] Num of samples: $noSamples\n" if ($flag_debug);
print TRACE "[EXP] Dat table: $datTable\n"       if ($flag_debug);

$command =
"$config{'general.scriptsDir'}/run_DEanalysis_batch_pool.sh $expTempDir/RNA_SEQ_COUNT_MATRIX.dat $noNormals $noSamples $config{'exp.resultsDir'}";
system($command);
print TRACE "[EXP] $command\n" if ($flag_debug);

close(TRACE);