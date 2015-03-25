#!/usr/bin/perl
use warnings;
#use strict;

my ($sample_data_dir, $out_file, $scripts_dir) = @ARGV;

require "$scripts_dir/Construct_network.pl";

my %gene_to_index;
my @index_to_gene;
my @connections;
construct_driver_net_network(\@index_to_gene, \%gene_to_index, \@connections, $scripts_dir);

my %sample_mutated_gene = ();
my %all_mutated_gene = ();
my @all_mutated_gene_order;

#Collect the fold change of each gene to compute the module impact
opendir(DIR, $sample_data_dir);
@the_DATA_DIR = readdir(DIR);
close(DIR);
foreach my $dir_sample (@the_DATA_DIR){
    $mutation_file_name = "$sample_data_dir/$dir_sample/Genelist_Status.txt";
    $mutation_file_name_cell = "$sample_data_dir/$dir_sample/Genelist_Status_cell.txt";
    #print STDERR " *** $mutation_file_name\n";
    if(-e $mutation_file_name || -e $mutation_file_name_cell){
	if( -e $mutation_file_name_cell){
	    $mutation_file_name = $mutation_file_name_cell;
	}
	
	my %gene_map = ();
	$sample_mutated_gene{$dir_sample} = \%gene_map;

	open(FILE, "$mutation_file_name");
	print STDERR " *** read file $mutation_file_name\n";#<STDIN>;
	#read the file to obtain the dysregulated and mutated genes
	while(<FILE>){
	    chop ($_);
	    @line = split(/\t/, $_);
	    my @parts = split(/_/,$line[0]);
	    my $gene_name = $parts[0];
	    if(exists $gene_to_index{$gene_name}){
		my $status = $parts[1];
		if (! ($status eq "UP" || $status eq "DOWN")){
		    #It is a mutation 
		    $sample_mutated_gene{$dir_sample}->{$gene_name} = 1;
		    if(!exists $all_mutated_gene{$gene_name}){
			$all_mutated_gene{$gene_name} = 1;
			push(@all_mutated_gene_order, $gene_name);
		    }
		}
	    }
	}
	close(FILE);
    }
}

#Write the matrix
open(OUT, ">$out_file");
my $res = "";
for($i = 0; $i < @all_mutated_gene_order; $i++){
    $gene_name = $all_mutated_gene_order[$i];
    $res .= $gene_name."\t";
}
chop $res; print OUT $res."\n";

foreach $sample (keys %sample_mutated_gene){
    print OUT $sample;
    $res = "";
    for($i = 0; $i < @all_mutated_gene_order; $i++){
	$gene_name = $all_mutated_gene_order[$i];
	if(exists $sample_mutated_gene{$sample}->{$gene_name}){
	    $res .= "\t1";
	}
	else{$res .= "\t0";}
    }
    print OUT $res."\n";#<STDIN>;
}

