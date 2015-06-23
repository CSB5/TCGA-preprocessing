#!/usr/bin/perl

use warnings;
use Config::Simple;
use Getopt::Long;

my ( $configFile, $flag_debug, $flag_help, %config, @queue, $command );
my $qsub = "qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -v PATH,PERL5LIB,R_LIBS_SITE,MOSEKLM_LICENSE_FILE,AUGUSTUS_CONFIG_PATH,CLASSPATH";

my $help_message = "
This script prepares the TCGA data for oncoIMPACT.

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

system("mkdir $config{'general.logsDir'}") unless (-e "$config{'general.logsDir'}");
open( TRACE,
	">$config{'general.logsDir'}/$config{'general.disease'}.error.log" );
@queue = (9999);
my $runtime = "$config{'cluster.runtime'}:0:0";

if ( $config{'general.flagCNV'} ) {
	print "Preparing CNV data. Please wait...";
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_cnv -e $config{'general.logsDir'}/$config{'general.disease'}_cnv.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_cnv.run.log $config{'general.scriptsDir'}/cnv_data.pl --in $config{'cnv.inputFile'} --destination $config{'general.analysisDir'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}

if ( $config{'general.flagSNP'} ) {
	print "Preparing SNP data. Please wait...";
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -N $config{'general.disease'}_snp -pe OpenMP 1 -e $config{'general.logsDir'}/$config{'general.disease'}_snp.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_snp.run.log $config{'general.scriptsDir'}/snp_data.pl --in $config{'snp.maf'} --destination $config{'general.analysisDir'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}

if ( $config{'general.flagEXP'} ) {
	print "Preparing Expression data. Please wait...";
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -N $config{'general.disease'}_exp -pe OpenMP 1 -e $config{'general.logsDir'}/$config{'general.disease'}_exp.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_exp.run.log $config{'general.scriptsDir'}/expression_data.pl --input $config{'exp.resultsDir'} --destination $config{'general.analysisDir'} --extention $config{'exp.extention'} --software $config{'exp.software'} --fdr $config{'exp.fdr'} --reads $config{'exp.reads'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	#prep_exp();
	
	print "Job submitted.\n";
}

if ( $config{'general.flagMNC'} ) {
	print "Merging and cleaning data. Please wait...";
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -N  $config{'general.disease'}_mnc -pe OpenMP 1 -e $config{'general.logsDir'}/$config{'general.disease'}_mnc.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_mnc.run.log -hold_jid " . join(",", @queue) . " $config{'general.scriptsDir'}/merge_and_clean.pl --input $config{'general.analysisDir'}";
	submit($command);	
	
	print "Job submitted.\n";	
}

if ( $config{'general.flagMUT'} ) {
	print "Generating gene mutation information. Please wait...";
	
	system("mkdir $config{'general.outDir'}") unless (-e "$config{'general.outDir'}");
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -N  $config{'general.disease'}_mut -pe OpenMP 1 -e $config{'general.logsDir'}/$config{'general.disease'}_mut.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_mut.run.log -hold_jid " . join(",", @queue) . " $config{'general.scriptsDir'}/generate_gene_mutation_info.pl --inDir $config{'general.analysisDir'}/COMPLETE_SAMPLES --outDir $config{'general.outDir'}";
	submit($command);	
	
	print "Job submitted.\n";	
	
	
	print "Filtering TCGA MAF file against complete_samples_list. Please wait...";
	
	my $lastID = $queue[-1];
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -N  $config{'general.disease'}_filterMAF -pe OpenMP 1 -e $config{'general.logsDir'}/$config{'general.disease'}_filterMAF.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_filterMAF.run.log -hold_jid $lastID $config{'general.scriptsDir'}/filter_maf.pl --samples $config{'general.outDir'}/complete_samples_list.txt --maf $config{'snp.maf'} --outDir $config{'general.outDir'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
	
	
	print "Generating sample-based mutation information. Please wait...";
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -N  $config{'general.disease'}_mut -pe OpenMP 1 -e $config{'general.logsDir'}/$config{'general.disease'}_mut.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_mut.run.log -hold_jid " . join(",", @queue) . " $config{'general.scriptsDir'}/generate_sample_based_gene_mutation_info.pl --in $config{'general.outDir'}/gene_mutation_frequency.txt --cc $config{'general.cancerCensus'} --outDir $config{'general.outDir'}";
	submit($command);	
	
	print "Job submitted.\n";	
}

if ( $config{'general.flagONCOIMPACT'} ) {
	print "Generating matrices for oncoIMPACT. Please wait...";
	
	system("mkdir $config{'general.outDir'}") unless (-e "$config{'general.outDir'}");
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_oncoIMPACT -e $config{'general.logsDir'}/$config{'general.disease'}_oncoIMPACT.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_oncoIMPACT.run.log -hold_jid " . join(",", @queue) . " $config{'general.scriptsDir'}/generate_sample_data.pl --inDir $config{'general.analysisDir'}/COMPLETE_SAMPLES --outDir $config{'general.outDir'}";
	$command = $command . " --debug" if ($flag_debug);
    submit($command);
    
    print "Job submitted.\n";
}

if ( $config{'general.flagDRIVERNET'} ) {
	print "Generating expression data for DriverNet. Please wait...";
	
	system("mkdir $config{'general.outDir'}") unless (-e "$config{'general.outDir'}");
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_DriverNet_exp -e $config{'general.logsDir'}/$config{'general.disease'}_DriverNet_exp.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_DriverNet_exp.run.log -hold_jid " . join(",", @queue) . " $config{'general.scriptsDir'}/generate_expression_table_for_DRIVERNET.pl --normal $config{'general.dataDir'}/$config{'exp.normalList'} --sdrf $config{'general.dataDir'}/$config{'exp.sdrf'} --rnaseq $config{'general.dataDir'}/$config{'exp.rnaseqDir'} --destination $config{'general.outDir'} --scriptsDir $config{'general.scriptsDir'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
	
	print "Generating mutation matrix for DriverNet. Please wait...";
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_DriverNet_mut -e $config{'general.logsDir'}/$config{'general.disease'}_DriverNet_mut.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_DriverNet_mut.run.log -hold_jid " . join(",", @queue) . " $config{'general.scriptsDir'}/produce_mut_matrix_driver_net.pl $config{'general.analysisDir'}/COMPLETE_SAMPLES $config{'general.outDir'}/Mutation_matrix_DRIVERNET.txt $config{'general.scriptsDir'}";
	submit($command);
	
	print "Job submitted.\n";
}

if ( $config{'general.flagANNOVAR'} ) {
	print "Running ANNOVAR annotation. Please wait...";
	
	system("mkdir $config{'general.outDir'}") unless (-e "$config{'general.outDir'}");
	
	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP 1 -N $config{'general.disease'}_ANNOVAR -e $config{'general.logsDir'}/$config{'general.disease'}_ANNOVAR.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_ANNOVAR.run.log -hold_jid " . join(",", @queue) . " $config{'general.scriptsDir'}/annotate_with_annovar.pl --in $config{'general.outDir'}/GDAC_somatic_mutations.filtered.maf --destination $config{'general.outDir'}";
	$command = $command . " --debug" if ($flag_debug);
	submit($command);
	
	print "Job submitted.\n";
}

if ( $config{'general.flagEXPR'} ) {
	print "Generating normalized expression counts. Please wait...";
	
	system("mkdir $config{'general.outDir'}") unless (-e "$config{'general.outDir'}");
	
	$command = "$config{'general.scriptsDir'}/trim_sample_names.pl --in $config{'general.analysisDir'}/RNA-SEQ/RNA_SEQ_COUNT_MATRIX.DESeq.normalized_filtered.counts.txt --out $config{'general.outDir'}/normalized_expression_matrix.txt";
	submit($command);
	
	$command = "cp $config{'general.analysisDir'}/RNA-SEQ/selected_normals.dat $config{'general.outDir'}/";
	submit($command);
	
	print "Job submitted.\n";
}



### Sub-routines ###
#sub prep_exp {
#	my (
#		$expTempDir, $noNormals, $noSamples,    $i,
#		$datTable,   $outPrefix, $message, @currentJobs,
#		$jobList,    $jobID, @refreshedArray, $x, $y
#	) = ();
#	$expTempDir = "$config{'general.analysisDir'}/expTemp";
#	unless ( -d $expTempDir ){
#		$command = "mkdir -p $expTempDir";
#		system($command);
#		print TRACE "[EXP] $command\n" if ($flag_debug);		
#	}
#	unless ( -e "$expTempDir/FILE_SAMPLE_MAP.txt" ){
#		$command = "cut -f2,22 $config{'general.dataDir'}/$config{'exp.sdrf'} | awk '{print \$2\"\\t\"\$1}' | grep \"rsem.genes.results\" | sort -k2,2 > $expTempDir/FILE_SAMPLE_MAP.txt"; 
#		system($command);
#		print TRACE "[EXP] $command\n" if ($flag_debug);
#	}
#	unless ( -e "$expTempDir/FILE_SAMPLE_MAP_NORMAL.txt"){
#		$command = "grep -f $config{'general.dataDir'}/$config{'exp.normalList'} $expTempDir/FILE_SAMPLE_MAP.txt | sort --random-sort -k2,2 | head -n $config{'exp.manNorm'} > $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt"; 
#		system($command);
#		print TRACE "[EXP] $command\n" if ($flag_debug);
#	}
#	unless ( -e "$expTempDir/RNA_SEQ_COUNT_MATRIX.dat"){
#		$command = "$config{'general.scriptsDir'}/merge_normal_cancer_RNA_seq.pl --normal $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt --id $expTempDir/FILE_SAMPLE_MAP.txt --rnaseq $config{'general.dataDir'}/$config{'exp.rnaseqDir'} --destination $config{'general.analysisDir'}/expTemp";
#		system($command);
#		print TRACE "[EXP] $command\n" if ($flag_debug);
#	}
#	$command = "wc -l < $expTempDir/FILE_SAMPLE_MAP_NORMAL.txt | cut -f 1";
#	chomp($noNormals = `$command`);
#	print TRACE "[EXP] $command\n" if ($flag_debug);
#	$command = "wc -l < $expTempDir/FILE_SAMPLE_MAP.txt | cut -f 1";
#	chomp($noSamples = `$command`);
#	print TRACE "[EXP] $command\n" if ($flag_debug);
#	$datTable = "$expTempDir/RNA_SEQ_COUNT_MATRIX.dat";
#
#	print TRACE "[EXP] Num of normals: $noNormals\n" if ($flag_debug);
#	print TRACE "[EXP] Num of samples: $noSamples\n" if ($flag_debug);
#	print TRACE "[EXP] Dat table: $datTable\n"       if ($flag_debug);
#
#	@currentJobs = ();
#	for ( my $i = ( $noNormals + 1 ) ; $i <= $noSamples ; $i++ ) {
#		chomp($outPrefix = `head -n 1 $datTable | cut -f $i`);
#		unless( -e "$config{'general.analysisDir'}/expTemp/${outPrefix}.dat" ){
#			open( OUT, ">$config{'general.analysisDir'}/expTemp/${outPrefix}.dat" );
#			print TRACE "[EXP] Generating $config{'general.analysisDir'}/expTemp/${outPrefix}.dat\n";
#			$command = "head -n 1 $datTable | cut -f 1-$noNormals,$i";
#			print OUT `$command`;
#			print TRACE "[EXP] $command\n" if ($flag_debug);
#			$x = $noNormals + 1;
#			$y = $i + 1;
#			$command = "tail -n+2 $datTable | cut -f 1-$x,$y";		
#			print OUT `$command`;
#			print TRACE "[EXP] $command\n" if ($flag_debug);
#			close OUT;		
#		}
#		unless( -e "$config{'general.analysisDir'}/expTemp/${outPrefix}.$config{'exp.extention'}" ){
#			$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_$config{'exp.software'} $config{'general.scriptsDir'}/call_DE.R $config{'general.analysisDir'}/expTemp/${outPrefix}.dat $noNormals $config{'exp.extention'} $config{'exp.software'}";
#			# $command = "$qsub -N $config{'general.disease'}_$config{'exp.software'} -e $config{'general.logsDir'}/$config{'general.disease'}_$config{'exp.software'}.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_$config{'exp.software'}.run.log $config{'general.scriptsDir'}/call_DE.R $config{'general.analysisDir'}/expTemp/${outPrefix}.dat $noNormals $config{'exp.extention'} $config{'exp.software'}";
#			$message = `$command`;		
#			print TRACE "[EXP] $command\n" if ($flag_debug);
#			push( @currentJobs, $message );
#			print TRACE "[EXP] Job $message submitted.\n" if ($flag_debug);			
#		}
#	}
#
#	while (@currentJobs) {
#		@refreshedArray = ();
#		$jobList = `qstat | tail -n+3 | cut -f3 -d " "`;
#		foreach (@currentJobs) {
#			if ( $jobList =~ m/$_/ ) {
#				$jobID = $_;
#				push( @refreshedArray, $jobID );
#				print TRACE "[EXP] $jobID still running.\n" if ($flag_debug);
#			}
#		}
#		@currentJobs = @refreshedArray;
#		sleep 60;
#	}
#	
#	sleep 60;
#
#	$command = "$qsub -l mem_free=$config{'cluster.mem'}G,h_rt=$runtime -pe OpenMP $config{'cluster.numThreads'} -N $config{'general.disease'}_exp -e $config{'general.logsDir'}/$config{'general.disease'}_exp.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_exp.run.log $config{'general.scriptsDir'}/expression_data.pl --input $config{'general.analysisDir'}/expTemp --destination $config{'general.analysisDir'} --extention $config{'exp.extention'} --software $config{'exp.software'} --fdr $config{'exp.fdr'} --reads $config{'exp.reads'}";
#	$command = "$qsub -N $config{'general.disease'}_exp -e $config{'general.logsDir'}/$config{'general.disease'}_exp.error.log -o $config{'general.logsDir'}/$config{'general.disease'}_exp.run.log $config{'general.scriptsDir'}/expression_data.pl --input $config{'general.analysisDir'}/expTemp --destination $config{'general.analysisDir'} --extention $config{'exp.extention'} --software $config{'exp.software'} --fdr $config{'exp.fdr'}";
#	$command = $command . " --debug" if ($flag_debug);
#	submit($command);
#	print TRACE "[EXP] $command\n";# if ($flag_debug);
#
#}    # end sub prep_exp


#sub checkQueue {
#	print "Checking if all preparation subroutines have completed. Please wait...";
#	while (@queue) {
#		my $jobList = `qstat | tail -n+3 | cut -f3 -d " "`;
#		my @refreshedQueue = ();
#		foreach (@queue) {
#			if ( $jobList =~ m/$_/ ) {
#				$jobID = $_;
#				push( @refreshedQueue, $jobID );
#				print TRACE "\t$jobID still running.\n" if ($flag_debug);
#			}
#		}
#		@queue = @refreshedQueue;
#		sleep 60;
#	}
#	print "done.\n";	
#}

sub submit {
	$command = "@_";
	print TRACE "[Command] $command\n"	if ($flag_debug);
	my $return = `$command`;
	chomp($return);	
	push( @queue, $return );
	print TRACE "\tJob $return submitted.\n" if ($flag_debug);
}    # end sub submit
