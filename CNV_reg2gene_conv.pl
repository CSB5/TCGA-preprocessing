#!/usr/bin/perl
use warnings "all";

my $cnv_file       = $ARGV[0];
my $region_file    = $ARGV[1];
my $output_file    = $ARGV[2];
my $size_threshold = -1;
if ( defined $ARGV[3] ) {
	$size_threshold = $ARGV[3];
}

open( COORD, $cnv_file );
my %chrom;
my $nb_diff = 0;

my $size = 0;

while (<COORD>) {
	@coord = split( /\s+/, $_ );
	if ( $coord[0] ne "Sample"
		&& index( $coord[0], "#" ) == -1 )
	{

		$chrom_name = "chr$coord[1]\.fa";
		$status     = "AMPL\t$coord[-1]";
		$status     = "DEL\t$coord[-1]" if ( $coord[-1] < 0 );
		$size       = ( $coord[3] - $coord[2] );
		my @inter =
		  ( $coord[2], $coord[3], $status, $size, "$coord[2]:$coord[3]" );
		if ( $size_threshold == -1 || $size < $size_threshold ) {
			if ( !( exists $chrom{$chrom_name} ) ) {
				my @c = ();
				$chrom{$chrom_name} = \@c;
			}
			push( @{ $chrom{$chrom_name} }, \@inter );
		}
	}
}
close(COORD);

my $n    = 1;
my $f    = 1;
my $good = 1;
my $interval_pos;

open( REGION, "$region_file" );
open( NEWF, ">", "$output_file" );
my $duplicated_region;

my $flag_only_one_start_site = 1;
my %map_gene;

while (<REGION>) {
	@coord = split( /\s+/, $_ );

	print STDERR "--> $n\n" if ( $n % 10000 == 1 );
	$chrom_name = $coord[0];
	$chrom_name = $chrom_name . ".fa" if ( index( $chrom_name, ".fa" ) == -1 );
	$cdsStart = $coord[4];
	$cdsEnd = $coord[5];

	if ( exists $chrom{$chrom_name}
		&& ( $interval_pos =
			is_in_interval( $chrom{$chrom_name}, $coord[2], 1 ) ) != -1
		&& $interval_pos ==
		is_in_interval( $chrom{$chrom_name}, $coord[3], 1 ) 
		&& $cdsEnd != $cdsStart		# condition added to filter pseudo genes which have identical CDS start and end sites
		)
	{
		$duplicated_region = $chrom{$chrom_name}->[$interval_pos];
		$size              = $duplicated_region->[1] - $duplicated_region->[0];
		if ( !$flag_only_one_start_site || !( exists $map_gene{ $coord[6] } ) )
		{
			print( NEWF "$coord[6]\_$duplicated_region->[2]\n" );
		}
		$map_gene{ $coord[6] } = 1;
		$f++;
	}
	$n++;
}
$f--;
print STDERR "nb region duplicated: $n $f " . ( $f / $n ) . "\n";

#Return -1 if it not in the interval, or the interval position
sub is_in_interval {
	my ( $ptr, $pos, $flag_in ) = @_;

	my @tab_inter = @$ptr;
	my $curr      = int( ( @tab_inter + 0 ) / 2 );
	my $min       = 0;
	my $max       = @tab_inter - 1;

	while ( ( $max - $min ) > 1
		&& !(  $pos >= $tab_inter[$curr][0]
			&& $pos <= $tab_inter[ $curr + 1 ][0] ) )
	{

		if ( $pos > $tab_inter[$curr][0] ) {
			$min = $curr;
		}
		else { $max = $curr; }

		$curr = int( $min + ( $max - $min ) / 2 );
	}

	my $res = -1;

	#res is on the current interval
	$res = $curr
	  if ( $pos >= $tab_inter[$curr][0] && $pos <= $tab_inter[$curr][1] );

	#Case of we can choose between 2 intervals
	if ( $max - $min == 1 ) {
		$curr = $min
		  if ( $pos >= $tab_inter[$min][0] && $pos <= $tab_inter[$min][1] );
		$curr = $max
		  if ( $pos >= $tab_inter[$max][0] && $pos <= $tab_inter[$max][1] );
	}

	$good++ if ( $res != -1 );

	return $res;
}
