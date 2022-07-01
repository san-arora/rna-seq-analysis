#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
$| = 1;

# Grab and set all options
my %OPTIONS;
getopts('co:r:', \%OPTIONS);

die qq(
Usage:   HTSeq_combine_output.pl [OPTIONS] <group matrix>

OPTIONS:		-c		reduce the matrix by removing any feature with no counts
			-o	STR	output file name for expression matrix
			-r	STR	the name of the output report
) unless (@ARGV == 1);

my @counts;
my @features;
my %report;
my @samplenames;
my @groups;
my @files;
my $groupcount = 0;
my %grouphash;

# read in the group matrix table
open(IN, $ARGV[0]);
while(<IN>) {
	chomp;

	# skip any empty lines
	next if $_ eq "";

	# split line
	my @sample = split("\t", $_);

	# check bam file exists
	die("Error: Sample HTseq file $sample[2] do not exist, please check\n") unless -f $sample[2];

	if(! defined $grouphash{$sample[0]}) {
		$groupcount++;
		$grouphash{$sample[0]} = "G${groupcount}:$sample[0]";
	}

	push @groups, $grouphash{$sample[0]};
	push @samplenames, $sample[1];
	push @files, $sample[2];
}
close(IN);

for(my $index = 0; $index <= $#files; $index++) {
	my $input_file = $files[$index];
	my $sample = $samplenames[$index];

	# read in hiseq results
	open(IN, "$input_file") or die("$input_file: $!\n");


	$report{$sample} = "Htseq file: $input_file\n";
	my $i = 0;
	while(my $row = <IN>) {
		# store the report is an hash
		if(grep /^_*no_feature|^_*ambiguous|^_*too_low_aQual|^_*not_aligned|^_*alignment_not_unique/, $row) {
			$report{$sample} .= $row;
		} else {
			# store the counts in a matrix
			chomp $row;
			my ($feature, $value) = split "\t", $row;
			$counts[$i][$index] = $value;
			if($input_file eq $files[0]) {
				push @features, $feature;
			}
		}
		$i++;
	}
}

# print the matrix
open(MATRIX, ">$OPTIONS{o}") || die "Could Not Create Output File $OPTIONS{o}: $!\n";
print MATRIX "#\t".join("\t", @groups)."\n";
print MATRIX "#Feature\t".join("\t", @samplenames)."\n";
for(my $row = 0; $row <= $#features; $row++) {
	if(defined $OPTIONS{c}) {
		my $rowsum = 0;
		$rowsum += $_ foreach @{ $counts[$row] };
		if(!$rowsum) {
			next;
		}
	}
	print MATRIX "$features[$row]\t".join("\t", @{ $counts[$row] })."\n";
}
close(MATRIX);

# print the report
open(REPORT, ">$OPTIONS{r}") || die "Could Not Create Output File $OPTIONS{r}: $!\n";
print REPORT "$groups[$_]:$samplenames[$_]\n$report{$samplenames[$_]}\n" foreach (0..$#samplenames);
close(REPORT);



