#!/usr/bin/perl
use warnings;


my ($read_count_file, $control_list_file, $out_file) = @ARGV;


print STDERR " *** READ control file\n";
my %control = ();
open(FILE, $control_list_file);
while(<FILE>){
    chop $_;
    $control_list_file{$control} = 1;
}


print STDERR " *** Read read count file\n";
open(FILE, $read_count_file);
my @sample_order_full = ();
my %sample_used = ();
my $str_normal;
my $str_tumor;
my $flag_normal = 0;
my $cmp_tumor = 0;
my $cmp_normal = 0;
open(OUT, ">$out_file");

%normal_ID = ();
%tumor_ID = ();

while(<FILE>){
    chop $_;
    @line = split(/\t/, $_);
    if(@sample_order_full == 0){
	for(my $i = 0; $i < @line; $i++){
	    $sample_name = $line[$i];
	    $flag_normal = 0;
	    if($sample_name =~ /TCGA(\-|\.)(.+)(\-|\.)(.+)(\-|\.)11(A|B)(\-|\.)(.+)/){
		$flag_normal = 1;
	    }
	    
	    push(@sample_order_full, $sample_name);
	    
	    if(!$flag_normal){
		$sample_used{$sample_name} = 1;
		$str_tumor .= $sample_name."\t";
		$tumor_ID{$sample_name} = $cmp_tumor;
		$cmp_tumor++;
	    }
	    else{
		if(exit $control{$sample_name}){
		    $str_normal .= $sample_name."\t";
		    $normal_ID{$sample_name} = $cmp_normal;
		    $cmp_normal++;
		}
	    }
	}
	chop $str_tumor;
	print OUT $str_normal.$str_tumor."\n";
    }
    
    else{
	#exit(0);
	@gene_exp = ();
	for(my $i = 1; $i < @line; $i++){
	    $sample_name = $sample_order_full[$i-1];
	    if(exists $sample_used{$sample_name}){
		if(exists $tumor_ID{$sample_name}){
		    $gene_exp[$cmp_normal+$tumor_ID{$sample_name}] = $line[$i];
		}
		else{
		    $gene_exp[$normal_ID{$sample_name}] = $line[$i];
		}
	    }
	}
	print OUT $line[0]."\t".(join("\t", @gene_exp))."\n";
    }
}

close(FILE);
close(OUT);


#Prostate
#PIPELINE/filter_expression_matrix.pl PROSTATE_CANCER /mnt/pnsg10_projects/bertrandd/oncoimpact/tcga/PROSTATE_ADENOCARCINOMA/ANALYSIS/DEanalysis/RNA_SEQ_COUNT_MATRIX_8normals.dat 8 /mnt/pnsg10_projects/bertrandd/oncoimpact/tcga/PROSTATE_ADENOCARCINOMA/ANALYSIS/DEanalysis/RNA_SEQ_COUNT_MATRIX_8normals_incomplete.dat
#PIPELINE/filter_expression_matrix.pl PROSTATE_CANCER_FULL /mnt/pnsg10_projects/bertrandd/oncoimpact/tcga/PROSTATE_ADENOCARCINOMA/ANALYSIS/DEanalysis/RNA_SEQ_COUNT_MATRIX_8normals.dat 8 /mnt/pnsg10_projects/bertrandd/oncoimpact/tcga/PROSTATE_ADENOCARCINOMA/ANALYSIS/DEanalysis/RNA_SEQ_COUNT_MATRIX_8normals_full.dat

#Bladder
