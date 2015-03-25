#!/usr/bin/perl
use warnings;


my ($read_count_file, $control_list_file, $call_type, $out_dir, $onco_impact_dir) = @ARGV;

$nb_normal = `wc -l $control_list_file | awk '{print \$1}'`; chop $nb_normal;


#Construct the matrix used for the RNA-SEQ DE ANALYSIS
my $read_count_normal_file = "$out_dir/RNA_SEQ_COUNT_MATRIX_NORMAL.dat";
run_exe("./select_expression_matrix_new_set.pl $read_count_file $control_list_file $read_count_normal_file");

$all_sample_str = `head -n1`;chop $all_sample_str;
@sample_tab = split(/\t/, $all_sample_str);
$nb_sample = @sample_tab+0;

#Run the DE analaysis
#Call without pool normalization
if($call_type eq "PAIR"){
    run_exe("run_DEanalysis_batch.sh $read_count_normal_file $nb_normal $nb_sample");
}

#Call with pool normalization
if($call_type eq "POOL"){
    run_exe("call_DE_whole_matrix $read_count_normal_file $nb_normal DE DESeq");
    my $read_count_normal_normalized_file = "$out_dir/RNA_SEQ_COUNT_MATRIX_NORMAL.DESeq.normalized_counts.txt.dat";
    run_exe("/home/chiakhb/SCRIPTS/oncoIMPACT/run_DEanalysis_batch_pool.sh $read_count_normal_normalized_file $nb_normal $nb_sample");
}

#Update the oncoIMPACT data
run_exe("./PIPELINE/filter_expression_data_read_count.pl $onco_impact_dir $out_dir pool.DE_result_cqn_DEseq 1 0");
run_exe("./PIPELINE/merge_and_clean.pl");

sub run_exe{
    my ($exe) = @_;
    $run = 1;
    print STDERR $exe."\n";;
    print STDERR `$exe` if($run);
}


