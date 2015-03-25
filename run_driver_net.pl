#!/usr/bin/perl

use warnings;
use Getopt::Long;

my ($data_dir, $driver_net_data, $nb_proc, $nb_random_test_by_proc, $script_dir);
$nb_proc = 20;

my $help_message = "
This script runs DriverNet.

Usage:
	run_driver_net.pl [OPTIONS]

Options:
	--dataDir = full path to oncoIMPACT_ready folder *
	--expData = full path to expression data for DriverNet *
	--numProc = number of processors to use [Default: 20]
	--scriptDir = path to oncoIMPACT scripts directory * 
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
	"dataDir=s" => \$data_dir,
	"expData=s" => \$driver_net_data,
	"numProc:i" => \$nb_proc,
	"scriptDir=s" => \$script_dir,
	"help"    => \$flag_help
) or die("Error in command line arguments.\n");

if ($flag_help) {
	print $help_message;
	exit 0;
}

$nb_random_test_by_proc = 500 / $nb_proc;

$driver_net_dir = "$data_dir/DRIVER_NET/DRIVER_NET_GENE_LIST";
system("mkdir $driver_net_dir");
	
#Construct the driver net data mutation data
$driver_net_mut_matrix = "$driver_net_dir/patMutMatrix.txt";
print STDERR "$script_dir/produce_mut_matrix_driver_net.pl $data_dir $driver_net_mut_matrix\n";
system("$script_dir/produce_mut_matrix_driver_net.pl $data_dir $driver_net_mut_matrix");
	
open (OUT_cmd, ">$driver_net_dir/cmd.txt");
	
#Construct the R file
for($i = 0; $i < $nb_proc; $i++){
    $file = "$driver_net_dir/rand_$i.R";
    open(OUT, ">$file");
    print OUT "library(DriverNet)\n";
    print OUT "data(sampleGeneNames)\n";
    print OUT "my_patMut <- read.table(\"$driver_net_mut_matrix\", header=T)\n";
    print OUT "load(\"$script_dir/influenceGraph.rda\")\n";

    #Data given in the paper
    if(index($driver_net_data, ".rda") != -1){
	print OUT "load(\"$driver_net_data\")\n";
    }
    #Use an outlier matrix computed from the data
    else{
	print OUT "patExpMatrix <- read.table(\"$driver_net_data\", header=T)\n";
	print OUT "patOutMatrix <- getPatientOutlierMatrix (as.matrix(patExpMatrix), th=2)\n";
    }

    print OUT "randomDriversResult = computeRandomizedResult(patMutMatrix=my_patMut, patOutMatrix=patOutMatrix, influenceGraph=influenceGraph, geneNameList= sampleGeneNames, outputFolder=NULL, printToConsole=FALSE,numberOfRandomTests=$nb_random_test_by_proc, weight=FALSE, purturbGraph=FALSE, purturbData=TRUE)\n";
    
    print OUT "save(randomDriversResult, file=\"$driver_net_dir/rand_$i.rda\")\n";
    close(OUT);
	    
    print OUT_cmd "/mnt/software/bin/R-3.1.0 --vanilla < $file\n";
}
close(OUT_cmd);
	
#run the analysis in parallel
$cmd = "cat $driver_net_dir/cmd.txt | xargs -I cmd --max-procs=$nb_proc bash -c cmd > /dev/null \n";
system("$cmd");

#Combine the results
$file = "$driver_net_dir/rand_comb.R";
open(OUT, ">$file");
for($i = 0; $i < $nb_proc; $i++){
    print OUT "load(\"$driver_net_dir/rand_$i.rda\")\n";
    print OUT "X_$i = randomDriversResult\n";
    print OUT "randomDriversResult_comb = c(X_0, X_1)\n" if($i == 1);
    print OUT "randomDriversResult_comb = c(randomDriversResult_comb, X_$i)\n" if($i > 1);
}
print OUT "save(randomDriversResult_comb, file=\"$driver_net_dir/rand_comb.rda\")\n";
close(OUT);
system("/mnt/software/bin/R-3.1.0 --vanilla < $file\n");

#Get the pvalue
$file = "$driver_net_dir/driver_net.R";
open(OUT, ">$file");
print OUT "library(DriverNet)\n";
print OUT "data(sampleGeneNames)\n";
print OUT "my_patMut <- read.table(\"$driver_net_mut_matrix\", header=T)\n";
print OUT "load(\"/home/bertrandd/pathways/Fine_Arts/DRIVER_NET/SOFTWARE/data/paperData/influenceGraph.rda\")\n";
#
#Data given in the paper
if(index($driver_net_data, ".rda") != -1){
    print OUT "load(\"$driver_net_data\")\n";
}
#Use an outlier matrix computed from the data
else{
    print OUT "patExpMatrix <- read.table(\"$driver_net_data\", header=T)\n";
    print OUT "patOutMatrix <- getPatientOutlierMatrix(as.matrix(patExpMatrix), th=2)\n";
    
}
#
#
print OUT "driversList = computeDrivers(my_patMut,patOutMatrix,influenceGraph,outputFolder=NULL, printToConsole=FALSE)\n";
print OUT "load(\"$driver_net_dir/rand_comb.rda\")\n";
print OUT "res = resultSummary(driversList, randomDriversResult_comb, my_patMut, influenceGraph, outputFolder=NULL, printToConsole=FALSE)\n";
print OUT "write.table(res[,1:4], file=\"$driver_net_dir/res_driver_net.dat\", sep=\"\t\", row.names=T, col.names=T, quote=FALSE)\n";
close(OUT);
system("/mnt/software/bin/R-3.1.0 --vanilla < $file\n");
    
#For the annotation
system("./PIPELINE/annotate_driver_net_results.pl $driver_net_dir/res_driver_net.dat | PIPELINE/compute_CS_concordance.pl - 5 CS >  $driver_net_dir/res_driver_concordance.dat");


#PIPELINE/run_stability_over_subsample_analysis_driver_net.pl Output_Destination/STABILITY_OVER_SUBSAMPLE_ANALYSIS /home/bertrandd/pathways/Fine_Arts/DRIVER_NET/SOFTWARE/data/paperData/GBM_data.rda \"-01\"
#PIPELINE/run_stability_over_subsample_analysis_driver_net.pl Ovarian_sample/STABILITY_OVER_SUBSAMPLE_ANALYSIS/ /home/bertrandd/pathways/Fine_Arts/DRIVER_NET/SOFTWARE/data/paperData/OV_data.rda A
#./PIPELINE/run_driver_net.pl PROSTATE_CANCER/ PROSTATE_CANCER/DRIVER_NET/DRIVER_NET_GENE_LIST/sampleExpr.dat \"\"
