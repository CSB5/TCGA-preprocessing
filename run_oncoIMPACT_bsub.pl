#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ( $configFile, $flag_debug, $flag_help, %config, @queue, $command );

my $help_message = "
This script prepares the TCGA data for oncoIMPACT run before launching oncoIMPACT.

Usage:
	run_oncoIMPACT.pl [OPTIONS]

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
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

# Read config file
print "Reading config file. Please wait...";
Config::Simple->import_from( $configFile, \%config );
print "done.\n";

system("mkdir $config{'general.logsDir'}") unless (-d "$config{'general.logsDir'}");
open( TRACE,
	">$config{'general.logsDir'}/$config{'general.disease'}.error.log" );
@queue = ();

if ( $config{'general.flagCNV'} ) {
	print "Preparing CNV data. Please wait...";
	$command = "bsub -q $config{'cluster.queue'} -M $config{'cluster.mem'} -n 1 -W $config{'cluster.runtime'}:0 -J $config{'general.disease'}_cnv -e $config{'general.logsDir'}/$config{'general.disease'}_cnv.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_cnv.run.log perl $config{'general.scriptsDir'}/cnv_data.pl --input $config{'general.dataDir'}/$config{'cnv.cnvFolder'} --fileManifest  $config{'general.dataDir'}/$config{'general.fileManifest'} --destination $config{'general.analysisDir'} --region $config{'cnv.region'} --extention $config{'cnv.extention'} --minUp $config{'cnv.minUp'} --minDown $config{'cnv.minDown'} --scripts $config{'general.scriptsDir'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	print "done.\n";
}

if ( $config{'general.flagSNP'} ) {
	print "Preparing SNP data. Please wait...";
	$command = "bsub -q $config{'cluster.queue'} -M $config{'cluster.mem'} -n 1 -W $config{'cluster.runtime'}:0 -J $config{'general.disease'}_snp -e $config{'general.logsDir'}/$config{'general.disease'}_snp.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_snp.run.log perl $config{'general.scriptsDir'}/snp_data.pl";
	if ( ref($config{'snp.snpFolders'} ) eq 'ARRAY'){
		foreach ( @{$config{'snp.snpFolders'}} ) {
			$command = $command . " --input $config{'general.dataDir'}/$_";
		}
	} else{
		$command = $command . " --input $config{'general.dataDir'}/$config{'snp.snpFolders'}";
	}
	$command = $command . " --destination $config{'general.analysisDir'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	print "done.\n";
}

if ( $config{'general.flagEXP'} ) {
	print "Preparing Expression data. Please wait...";
	prep_exp();
	print "done.\n";
}

if ( $config{'general.flagMNC'} ) {
	print "Checking if all preparation subroutines have completed. Please wait...";
	while (@queue) {
		my $jobList = `bjobs | tail -n+2 | grep -v "^ " | cut -f1 -d " "`;
		my @refreshedQueue = ();
		foreach (@queue) {
			if ( $jobList =~ m/$_/ ) {
				$jobID = $_;
				push( @refreshedQueue, $jobID );
				print TRACE "\t$jobID still running.\n" if ($flag_debug);
			}
		}
		@queue = @refreshedQueue;
		sleep 60;
	}
	print "done.\n";

	print "Merging and cleaning data. Please wait...";
	$command = "bsub -q $config{'cluster.queue'} -M $config{'cluster.mem'} -n 1 -W $config{'cluster.runtime'}:0 -J $config{'general.disease'}_mnc -e $config{'general.logsDir'}/$config{'general.disease'}_mnc.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_mnc.run.log perl $config{'general.scriptsDir'}/merge_and_clean.pl --input $config{'general.analysisDir'}";
	submit($command);
	print "done.\n";
}

if ( $config{'general.flagONCOIMPACT'} ) {
	print "Running oncoIMPACT. Please wait...\n";
	print "Generating basic stats. Please wait...";
	my $sample_stats_dir = "$config{'general.analysisDir'}/oncoIMPACT_ready/SAMPLE_STATS";
	my $sample_stats_file = "$sample_stats_dir/basic_stats.dat";

	my $basic_stats_path = "$config{'general.scriptsDir'}/plot_basic_stats.pl";
	my $basic_stats_R_path = "$config{'general.scriptsDir'}/sample_data.R";

    system("mkdir $sample_stats_dir") unless (-d $sample_stats_dir);
    $command = "bsub -q $config{'cluster.queue'} -M $config{'cluster.mem'} -n 1 -W $config{'cluster.runtime'}:0 -J $config{'general.disease'}_cnv -e $config{'general.logsDir'}/$config{'general.disease'}_gbs.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_gbs.run.log $basic_stats_path $config{'general.dataDir'} $config{'oncoIMPACT.networkType'} > $sample_stats_file";
    submit($command);

    system("cp $basic_stats_R_path $sample_stats_dir/");
    $command = "bsub -q $config{'cluster.queue'} -M $config{'cluster.mem'} -n 1 -W $config{'cluster.runtime'}:0 -J $config{'general.disease'}_cnv -e $config{'general.logsDir'}/$config{'general.disease'}_gbs2.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_gbs2.run.log sample_data.R";
    submit($command);
    print "done.\n";
}


### Sub-routines ###
sub prep_exp {
	my (
		$expTempDir, $noNormals, $noSamples,    $i,
		$datTable,   $outPrefix, $bsub_message, @currentJobs,
		$jobList,    $jobID, @refreshedArray, $x, $y
	) = ();
	$expTempDir = "$config{'general.analysisDir'}/expTemp";
	unless ( -d $expTempDir ){
		$command = "mkdir -p $expTempDir";
		system($command);
		print TRACE "[EXP] $command\n" if ($flag_debug);		
	}
	unless ( -e "$expTempDir/FILE_SAMPLE_MAP.txt" ){
		$command = "cut -f2,22 $config{'general.dataDir'}/$config{'exp.sdrf'} | awk '{print \$2\"\\t\"\$1}' | grep \"rsem.genes.results\" | sort -k2,2 > $expTempDir/FILE_SAMPLE_MAP.txt"; 
		system($command);
		print TRACE "[EXP] $command\n" if ($flag_debug);
	}
	unless ( -e "$expTempDir/FILE_SAMPLE_MAP_NORMAL.txt"){
		$command = "grep \"rsem.genes.results\" $config{'general.dataDir'}/$config{'exp.normalList'} | cut -f4 | tr \"::\" \"\\t\" | cut -f1 | grep -f - $expTempDir/FILE_SAMPLE_MAP.txt | sort --random-sort -k2,2 | head -n $config{'exp.manNorm'} > $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt"; 
		system($command);
		print TRACE "[EXP] $command\n" if ($flag_debug);
	}
	unless ( -e "$expTempDir/RNA_SEQ_COUNT_MATRIX.dat"){
		$command = "perl $config{'general.scriptsDir'}/merge_normal_cancer_RNA_seq.pl --normal $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt --id $expTempDir/FILE_SAMPLE_MAP.txt --rnaseq $config{'general.dataDir'}/$config{'exp.rnaseqDir'} --destination $config{'general.analysisDir'}/expTemp";
		system($command);
		print TRACE "[EXP] $command\n" if ($flag_debug);
	}
	$command = "wc -l < $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt | cut -f 1";
	chomp($noNormals = `$command`);
	print TRACE "[EXP] $command\n" if ($flag_debug);
	$command = "wc -l < $expTempDir/FILE_SAMPLE_MAP.txt | cut -f 1";
	chomp($noSamples = `$command`);
	print TRACE "[EXP] $command\n" if ($flag_debug);
	$datTable = "$expTempDir/RNA_SEQ_COUNT_MATRIX.dat";

	print TRACE "[EXP] Num of normals: $noNormals\n" if ($flag_debug);
	print TRACE "[EXP] Num of samples: $noSamples\n" if ($flag_debug);
	print TRACE "[EXP] Dat table: $datTable\n"       if ($flag_debug);

	@currentJobs = ();
	for ( my $i = ( $noNormals + 1 ) ; $i <= $noSamples ; $i++ ) {
		chomp($outPrefix = `head -n 1 $datTable | cut -f $i`);
		unless( -e "$config{'general.analysisDir'}/expTemp/${outPrefix}.dat" ){
			open( OUT, ">$config{'general.analysisDir'}/expTemp/${outPrefix}.dat" );
			print TRACE "[EXP] Generating $config{'general.analysisDir'}/expTemp/${outPrefix}.dat\n";
			$command = "head -n 1 $datTable | cut -f 1-$noNormals,$i";
			print OUT `$command`;
			print TRACE "[EXP] $command\n" if ($flag_debug);
			$x = $noNormals + 1;
			$y = $i + 1;
			$command = "tail -n+2 $datTable | cut -f 1-$x,$y";		
			print OUT `$command`;
			print TRACE "[EXP] $command\n" if ($flag_debug);
			close OUT;		
		}
		unless( -e "$config{'general.analysisDir'}/expTemp/${outPrefix}.$config{'exp.extention'}" ){
			$command = "bsub -q $config{'cluster.queue'} -M $config{'cluster.mem'} -n 1 -W $config{'cluster.runtime'}:0 -J $config{'general.disease'}_$config{'exp.software'} $config{'general.scriptsDir'}/call_DE.R $config{'general.analysisDir'}/expTemp/${outPrefix}.dat $noNormals $config{'exp.extention'} $config{'exp.software'}";
			$bsub_message = `$command`;		
			print TRACE "[EXP] $command\n" if ($flag_debug);
			print TRACE "[EXP] Bsub return: $bsub_message\n" if ($flag_debug);
	
			for ( split /^/, $bsub_message ) {
				if (/Job <(\d*)> is submitted to queue <(\w*)>./) {
					push( @currentJobs, $1 );
					print TRACE "[EXP] Job $1 submitted.\n" if ($flag_debug);
				}
			}
		}
	}

	while (@currentJobs) {
		@refreshedArray = ();
		$jobList = `bjobs | tail -n+2 | grep -v "^ " | cut -f1 -d " "`;
		foreach (@currentJobs) {
			if ( $jobList =~ m/$_/ ) {
				$jobID = $_;
				push( @refreshedArray, $jobID );
				print TRACE "[EXP] $jobID still running.\n" if ($flag_debug);
			}
		}
		@currentJobs = @refreshedArray;
		sleep 60;
	}

	$command =
"bsub -q $config{'cluster.queue'} -M $config{'cluster.mem'} -n 1 -W $config{'cluster.runtime'}:0 -J $config{'general.disease'}_exp -e $config{'general.logsDir'}/$config{'general.disease'}_exp.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_exp.run.log perl $config{'general.scriptsDir'}/expression_data.pl --input $config{'general.analysisDir'}/expTemp --destination $config{'general.analysisDir'} --extention $config{'exp.extention'} --software $config{'exp.software'} --fdr $config{'exp.fdr'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	print TRACE "[EXP] $command\n" if ($flag_debug);

}    # end sub prep_exp

sub submit {
	$command = "@_";
	my $return = `$command`;
	print TRACE "[Command] $command\n"           if ($flag_debug);
	print TRACE "\tBsub return: $return\n" if ($flag_debug);
	for ( split /^/, $return ) {
		if (/Job <(\d*)> is submitted to queue <(\w*)>./) {
			push( @queue, $1 );
			print TRACE "\tJob $1 submitted.\n" if ($flag_debug);
		}
	}
}    # end sub submit
