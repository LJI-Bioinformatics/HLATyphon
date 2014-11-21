#!/share/apps/perl/perl-5.18.1/bin/perl

# Author: Alex Fu (zfu@liai.org)
# Date: April 18th., 2013
# Updated: April 28th., 2013
# Version: 1.0.0

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
use Date::Simple qw(today);


use File::Basename;
use File::Find::Rule;
use File::Find;
use File::Copy qw(move copy);
use File::Copy::Recursive qw(dirmove);
use File::Path qw(make_path remove_tree rmtree);
use File::stat;
use Storable qw(dclone);
use Getopt::Long;

my $HLA_Dir; my $OUTPUT_Dir; my $Run_ID; my $Validation_Dir;

GetOptions( "run-data-path|r=s" 	=> \$HLA_Dir,
			"output-path|o=s"		=> \$OUTPUT_Dir, 
			"run-id|id=s"       	=> \$Run_ID,
			"validation-data|v=s"   => \$Validation_Dir,
		  );

chomp $HLA_Dir; chomp $OUTPUT_Dir; chomp $Run_ID; chomp $Validation_Dir;

my $GenDX_Typing_File      = $HLA_Dir."/Sample_Locus.csv";
my $Validation_Output_File = $OUTPUT_Dir."/Validation.txt";
my $Unknown_Output_File	   = $OUTPUT_Dir."/Unknown.txt";

my $trimmed_reads_file     = $HLA_Dir."/countTrimmedReads.txt";
my %trimmed_reads_hash;
my $raw_reads_file         = $HLA_Dir."/countRawReads.txt";
my %raw_reads_hash;

print "***********************************************************************\n\n";

open (RAWREADS, $raw_reads_file) || die "Can not open $raw_reads_file.\n";
while (<RAWREADS>){
	
	chomp;
	my @each_fastq_file = split (/\s+/, $_);
	my @each_fastq_file_name = split ('_', $each_fastq_file[0]);
	$raw_reads_hash{$each_fastq_file_name[0]} = $each_fastq_file[1];
	
	print $each_fastq_file_name[0]," ===> ",$raw_reads_hash{$each_fastq_file_name[0]},"\n";
}
close RAWREADS;

print "***********************************************************************\n\n";

open (TRIMMED, $trimmed_reads_file) || die "Can not open $trimmed_reads_file.\n";
while (<TRIMMED>){
	
	chomp;
	my @each_fastq_file = split (/\s+/, $_);
	my @each_fastq_file_name = split ('_', $each_fastq_file[0]);
	$trimmed_reads_hash{$each_fastq_file_name[0]} = $each_fastq_file[1];
	
	print $each_fastq_file_name[0]," ===> ",$trimmed_reads_hash{$each_fastq_file_name[0]},"\n";	
}
close TRIMMED;

print "***********************************************************************\n\n";

open (OUTPUT,  ">", $Validation_Output_File) || die "Cannot open $Validation_Output_File.\n";
open (OUTPUT1, ">", $Unknown_Output_File)    || die "Cannot open $Unknown_Output_File.\n";
open (GENDX,   "<", $GenDX_Typing_File)      || die "Cannot open $GenDX_Typing_File.\n";

#################################################
#
# Initialize output variables.
#
#################################################

my $Sample_ID;
my $Locus;
my $trimmedReadsCounts;
my $rawReadsCounts;
my $readsFilteredPercentage;

my $Luminex_Allele;
my $Luminex_Allele_1;
my $Luminex_Allele_2;
my $Allele_Combination_1;
my $Putative_Allele_1;
my $Putative_Allele_2;
my $Allele_Combination_2;
my $Alternative_Combinations;
my $Allele_Combination_1_AF;
my $Allele_Combination_1_AF_Ratio;
my $Allele_Combination_1_Evidence;
my $Allele_Combination_1_Reads;
my $LIAI_match_Luminex;

my $Allele1_Exon2_Reads;
my $Allele1_Exon3_Reads;
my $Allele1_Exon4_Reads;
my $Allele2_Exon2_Reads;
my $Allele2_Exon3_Reads;
my $Allele2_Exon4_Reads;

my $Allele1_Exon2_Bases;
my $Allele1_Exon3_Bases;
my $Allele1_Exon4_Bases;
my $Allele2_Exon2_Bases;
my $Allele2_Exon3_Bases;
my $Allele2_Exon4_Bases;

my $Allele_Combination_Mapped_Reads;
#my $Mapped_Reads_Hash = get_mapped_reads($HLA_Dir);


#################################################
#
# Open GenDX Result.
#
#################################################

my @GENDX_Typing_Result = <GENDX>;

GENDXLOG: foreach (@GENDX_Typing_Result) {
	
	chomp; my @Each_GenDX_Outcome = split(',', $_, 2);
	$Each_GenDX_Outcome[0] =~ s/\s+//g; $Sample_ID = $Each_GenDX_Outcome[0];
	$Each_GenDX_Outcome[1] =~ s/\s+//g; $Locus     = $Each_GenDX_Outcome[1];
	my $sample_ID_Locus = $Sample_ID.$Locus;

	my $Luminex_Results_Hash = get_luminex_validation($Sample_ID, $Locus, $Validation_Dir);
	
	print "********************************************\n\n";
	print "=== Sample ID:    $Sample_ID\n";
	print "=== LOCUS:        $Locus\n";
	print "=== The sample ID and Locus is: $sample_ID_Locus.\n";
	
	if (! exists $Luminex_Results_Hash->{$sample_ID_Locus}) {
		
		print "Luminex validation did not exist.\n\n";
	}
	else {
		
		$Luminex_Allele = join(',', @{$Luminex_Results_Hash->{$sample_ID_Locus}});
		$Luminex_Allele_1 = $Luminex_Results_Hash->{$sample_ID_Locus}->[0];
		$Luminex_Allele_2 = $Luminex_Results_Hash->{$sample_ID_Locus}->[1];
	}
	print "=== LUMINEX ALLELE: $Luminex_Allele\n";
	
	#####################################################
	#
	# Get number of mapped reads.
	#
	#####################################################
				
	foreach my $Sample_Locus_Name (keys %raw_reads_hash) {
				
		if (grep /$Sample_Locus_Name/, $sample_ID_Locus) {
					
			$trimmedReadsCounts      = $trimmed_reads_hash{$Sample_Locus_Name};
			$rawReadsCounts          = $raw_reads_hash{$Sample_Locus_Name};
			$readsFilteredPercentage = $trimmedReadsCounts / $rawReadsCounts;
			$readsFilteredPercentage = sprintf ('%.4f', $readsFilteredPercentage);
			
			}
	} #end of itinerary of %raw_reads_hash
	
	#####################################################
	#
	# Get LIAI HLA Typing Results.
	#
	#####################################################
	
	my $LIAI_results = get_LIAI_results($Validation_Dir, $HLA_Dir, $Sample_ID, $Locus);
	
	
        $Allele_Combination_1			= $LIAI_results->[0];
	$Allele_Combination_2			= $LIAI_results->[1];
	$Allele_Combination_1_AF		= $LIAI_results->[2];
	$Allele_Combination_1_AF_Ratio	= $LIAI_results->[3];
	$Allele_Combination_1_Evidence	= $LIAI_results->[4];
	$Allele_Combination_1_Reads		= $LIAI_results->[5];
	
	########################################################
	#
	# Print out samples with no validation data
	#
	########################################################
	
	
	if (! exists $Luminex_Results_Hash->{$sample_ID_Locus}) {
		
		if ($Allele_Combination_1 eq "N/A") {
			
			$Putative_Allele_1 = "N/A";
			$Putative_Allele_2 = "N/A";
			$Alternative_Combinations = "N/A";
			
			$Allele1_Exon2_Reads = "N/A";
 			$Allele1_Exon3_Reads = "N/A";
 			$Allele1_Exon4_Reads = "N/A";
 			$Allele2_Exon2_Reads = "N/A";
 			$Allele2_Exon3_Reads = "N/A";
 			$Allele2_Exon4_Reads = "N/A";

 			$Allele1_Exon2_Bases = "N/A";
 			$Allele1_Exon3_Bases = "N/A";
 			$Allele1_Exon4_Bases = "N/A";
 			$Allele2_Exon2_Bases = "N/A";
 			$Allele2_Exon3_Bases = "N/A";
 			$Allele2_Exon4_Bases = "N/A";

 			$Allele_Combination_Mapped_Reads = "N/A";

		}
		elsif (grep /;/, $Allele_Combination_1)	{
			
			$Putative_Allele_1 = "N/A";
			$Putative_Allele_2 = "N/A";
			$Alternative_Combinations = $Allele_Combination_1;
			
			$Allele1_Exon2_Reads = "N/A";
 			$Allele1_Exon3_Reads = "N/A";
 			$Allele1_Exon4_Reads = "N/A";
 			$Allele2_Exon2_Reads = "N/A";
 			$Allele2_Exon3_Reads = "N/A";
 			$Allele2_Exon4_Reads = "N/A";

 			$Allele1_Exon2_Bases = "N/A";
 			$Allele1_Exon3_Bases = "N/A";
 			$Allele1_Exon4_Bases = "N/A";
 			$Allele2_Exon2_Bases = "N/A";
 			$Allele2_Exon3_Bases = "N/A";
 			$Allele2_Exon4_Bases = "N/A";

 			$Allele_Combination_Mapped_Reads = "N/A";
		}
		else {
			
			($Putative_Allele_1, $Putative_Allele_2) = split (",", $Allele_Combination_1, 2);
			$Alternative_Combinations = $Allele_Combination_2;
			
			my ($Allele1_Exon_Bases, $Allele1_FullName) = get_Exon_Bases($Putative_Allele_1, $HLA_Dir, $Sample_ID, $Locus);
			my ($Allele2_Exon_Bases, $Allele2_FullName) = get_Exon_Bases($Putative_Allele_2, $HLA_Dir, $Sample_ID, $Locus);

			my $Allele1_Exon_Reads = get_Exon_Reads($Allele1_FullName, $HLA_Dir, $Sample_ID, $Locus);
			my $Allele2_Exon_Reads = get_Exon_Reads($Allele2_FullName, $HLA_Dir, $Sample_ID, $Locus);
			
			$Allele1_Exon2_Reads = $Allele1_Exon_Reads->[0];
 			$Allele1_Exon3_Reads = $Allele1_Exon_Reads->[1];
 			$Allele1_Exon4_Reads = $Allele1_Exon_Reads->[2];
 			$Allele2_Exon2_Reads = $Allele2_Exon_Reads->[0];
 			$Allele2_Exon3_Reads = $Allele2_Exon_Reads->[1];
 			$Allele2_Exon4_Reads = $Allele2_Exon_Reads->[2];

 			$Allele1_Exon2_Bases = $Allele1_Exon_Bases->[0];
 			$Allele1_Exon3_Bases = $Allele1_Exon_Bases->[1];
 			$Allele1_Exon4_Bases = $Allele1_Exon_Bases->[2];
 			$Allele2_Exon2_Bases = $Allele2_Exon_Bases->[0];
 			$Allele2_Exon3_Bases = $Allele2_Exon_Bases->[1];
 			$Allele2_Exon4_Bases = $Allele2_Exon_Bases->[2];
			
			$Allele_Combination_Mapped_Reads = get_AlleleCombo_Reads($Putative_Allele_1, $Putative_Allele_2, $HLA_Dir, $Sample_ID, $Locus);
			
			
		}
		
		print OUTPUT1 $Run_ID,"\t";
		print OUTPUT1 $Sample_ID,"\t";
		print OUTPUT1 $Locus,"\t";
		print OUTPUT1 $trimmedReadsCounts,"\t";
		print OUTPUT1 $rawReadsCounts,"\t";
		print OUTPUT1 $readsFilteredPercentage,"\t";
		print OUTPUT1 $Allele_Combination_Mapped_Reads,"\t";
		print OUTPUT1 $Allele_Combination_1_Reads,"\t";
		
		
		print OUTPUT1 ' '.$Putative_Allele_1,"\t";
		print OUTPUT1 ' '.$Putative_Allele_2,"\t";
		
		print OUTPUT1 $Allele1_Exon2_Bases, "\t";
		print OUTPUT1 $Allele1_Exon3_Bases, "\t";
		print OUTPUT1 $Allele1_Exon4_Bases, "\t";
		print OUTPUT1 $Allele2_Exon2_Bases, "\t";
		print OUTPUT1 $Allele2_Exon3_Bases, "\t";
		print OUTPUT1 $Allele2_Exon4_Bases, "\t";

		print OUTPUT1 $Allele1_Exon2_Reads, "\t";
		print OUTPUT1 $Allele1_Exon3_Reads, "\t";
		print OUTPUT1 $Allele1_Exon4_Reads, "\t";
		print OUTPUT1 $Allele2_Exon2_Reads, "\t";
		print OUTPUT1 $Allele2_Exon3_Reads, "\t";
		print OUTPUT1 $Allele2_Exon4_Reads, "\t";
		
		
		print OUTPUT1 ' '.$Alternative_Combinations,"\t";
	
		print OUTPUT1 $Allele_Combination_1_AF,"\t";
		print OUTPUT1 $Allele_Combination_1_AF_Ratio,"\t";
		print OUTPUT1 $Allele_Combination_1_Evidence,"\n";


	}
	
	########################################################
	#
	# Print out samples with validation data
	#
	########################################################
	
	else {
		
		if ($Allele_Combination_1 eq "N/A") {
			
			$Putative_Allele_1 = "N/A";
			$Putative_Allele_2 = "N/A";
			$Alternative_Combinations = "N/A";
			$LIAI_match_Luminex = "NO_CALL";

			$Allele1_Exon2_Reads = "N/A";
 			$Allele1_Exon3_Reads = "N/A";
 			$Allele1_Exon4_Reads = "N/A";
 			$Allele2_Exon2_Reads = "N/A";
 			$Allele2_Exon3_Reads = "N/A";
 			$Allele2_Exon4_Reads = "N/A";

 			$Allele1_Exon2_Bases = "N/A";
 			$Allele1_Exon3_Bases = "N/A";
 			$Allele1_Exon4_Bases = "N/A";
 			$Allele2_Exon2_Bases = "N/A";
 			$Allele2_Exon3_Bases = "N/A";
 			$Allele2_Exon4_Bases = "N/A";

 			$Allele_Combination_Mapped_Reads = "N/A";

		}
		elsif (grep /;/, $Allele_Combination_1)	{
			
			$Putative_Allele_1 = "N/A";
			$Putative_Allele_2 = "N/A";
			$Alternative_Combinations = $Allele_Combination_1;
			$LIAI_match_Luminex = "NO_CALL";
			
			$Allele1_Exon2_Reads = "N/A";
 			$Allele1_Exon3_Reads = "N/A";
 			$Allele1_Exon4_Reads = "N/A";
 			$Allele2_Exon2_Reads = "N/A";
 			$Allele2_Exon3_Reads = "N/A";
 			$Allele2_Exon4_Reads = "N/A";

 			$Allele1_Exon2_Bases = "N/A";
 			$Allele1_Exon3_Bases = "N/A";
 			$Allele1_Exon4_Bases = "N/A";
 			$Allele2_Exon2_Bases = "N/A";
 			$Allele2_Exon3_Bases = "N/A";
 			$Allele2_Exon4_Bases = "N/A";

 			$Allele_Combination_Mapped_Reads = "N/A";
			
		}
		else {
			
			($Putative_Allele_1, $Putative_Allele_2) = split (",", $Allele_Combination_1, 2);
			$Alternative_Combinations = $Allele_Combination_2;
			$LIAI_match_Luminex = compare_with_Luminex($Allele_Combination_1, $Luminex_Allele);

			my ($Allele1_Exon_Bases, $Allele1_FullName) = get_Exon_Bases($Putative_Allele_1, $HLA_Dir, $Sample_ID, $Locus);
			my ($Allele2_Exon_Bases, $Allele2_FullName) = get_Exon_Bases($Putative_Allele_2, $HLA_Dir, $Sample_ID, $Locus);

			my $Allele1_Exon_Reads = get_Exon_Reads($Allele1_FullName, $HLA_Dir, $Sample_ID, $Locus);
			my $Allele2_Exon_Reads = get_Exon_Reads($Allele2_FullName, $HLA_Dir, $Sample_ID, $Locus);
			
			$Allele1_Exon2_Reads = $Allele1_Exon_Reads->[0];
 			$Allele1_Exon3_Reads = $Allele1_Exon_Reads->[1];
 			$Allele1_Exon4_Reads = $Allele1_Exon_Reads->[2];
 			$Allele2_Exon2_Reads = $Allele2_Exon_Reads->[0];
 			$Allele2_Exon3_Reads = $Allele2_Exon_Reads->[1];
 			$Allele2_Exon4_Reads = $Allele2_Exon_Reads->[2];

 			$Allele1_Exon2_Bases = $Allele1_Exon_Bases->[0];
 			$Allele1_Exon3_Bases = $Allele1_Exon_Bases->[1];
 			$Allele1_Exon4_Bases = $Allele1_Exon_Bases->[2];
 			$Allele2_Exon2_Bases = $Allele2_Exon_Bases->[0];
 			$Allele2_Exon3_Bases = $Allele2_Exon_Bases->[1];
 			$Allele2_Exon4_Bases = $Allele2_Exon_Bases->[2];
			
			$Allele_Combination_Mapped_Reads = get_AlleleCombo_Reads($Putative_Allele_1, $Putative_Allele_2, $HLA_Dir, $Sample_ID, $Locus);
			
		}
		

		chomp $Luminex_Allele_1;
		chomp $Luminex_Allele_2;

		print OUTPUT $Run_ID,"\t";
		print OUTPUT $Sample_ID,"\t";
		print OUTPUT $Locus,"\t";
		print OUTPUT $trimmedReadsCounts,"\t";
		print OUTPUT $rawReadsCounts,"\t";
		print OUTPUT $readsFilteredPercentage,"\t";
		print OUTPUT $Allele_Combination_Mapped_Reads,"\t";
		print OUTPUT $Allele_Combination_1_Reads,"\t";

		print OUTPUT ' '.$Luminex_Allele_1,"\t";
		print OUTPUT ' '.$Luminex_Allele_2,"\t";
		
		print OUTPUT ' '.$Putative_Allele_1,"\t";
		print OUTPUT ' '.$Putative_Allele_2,"\t";
		
		print OUTPUT $LIAI_match_Luminex,"\t";
		
		print OUTPUT $Allele1_Exon2_Bases, "\t";
		print OUTPUT $Allele1_Exon3_Bases, "\t";
		print OUTPUT $Allele1_Exon4_Bases, "\t";
		print OUTPUT $Allele2_Exon2_Bases, "\t";
		print OUTPUT $Allele2_Exon3_Bases, "\t";
		print OUTPUT $Allele2_Exon4_Bases, "\t";

		print OUTPUT $Allele1_Exon2_Reads, "\t";
		print OUTPUT $Allele1_Exon3_Reads, "\t";
		print OUTPUT $Allele1_Exon4_Reads, "\t";
		print OUTPUT $Allele2_Exon2_Reads, "\t";
		print OUTPUT $Allele2_Exon3_Reads, "\t";
		print OUTPUT $Allele2_Exon4_Reads, "\t";
		
		print OUTPUT ' '.$Alternative_Combinations,"\t";
	
		print OUTPUT $Allele_Combination_1_AF,"\t";
		print OUTPUT $Allele_Combination_1_AF_Ratio,"\t";
		print OUTPUT $Allele_Combination_1_Evidence,"\n";
		

	}
	
	print "=== Done $Sample_ID HLA_$Locus.\n";
		
}# end of each GenDX result.

close OUTPUT;
close OUTPUT1;
close GENDX;

sub get_allele_frequency {
	
	my $validation_dir = shift;
	my $alleleFreqFile = $validation_dir."/Allele_Frequencies.csv";
	
	open (FILE, $alleleFreqFile) || die "Can not open $alleleFreqFile.\n";
	my %AF_hash;
	
	while (<FILE>) {

		chomp $_; 		
		next if (!$_); 		
		$_ =~ s/^\s+$// if (grep/^\s+$/, $_);
		my @line = split(",", $_);
		$AF_hash{$line[0]} = $line[1];

	}
	
	close FILE;
	return \%AF_hash;
		
}

sub get_AlleleCombo_Reads {
	
	my $allele1 = shift;
	my $allele2 = shift;

	my $Input_Dir  	   = shift;
	my $Sample_ID      = shift;
	my $Locus_Type     = shift;
	
	my $parent_Dir = $Input_Dir."/".$Sample_ID."/".$Locus_Type."_nuc";
	my $mappedReadsFile = $parent_Dir."/Collapse_Count_Reads.txt";
	
	open (FILE, $mappedReadsFile) || die "Can not open $mappedReadsFile.\n";
	my $totalMappedReads = 0;
	
	while (<FILE>){

		chomp $_; 		
		next if (!$_); 		
		my @line = split(/\s+/, $_);
		
		if ( ($line[0] eq $allele1) and ($line[1] eq $allele2) ) {
			
			$totalMappedReads = $line[2];
			last;
		}
		
		if ( ($line[1] eq $allele1) and ($line[0] eq $allele2) ) {
			
			$totalMappedReads = $line[2];
			last;
		}
		

	}
	
	close FILE;
	return $totalMappedReads;
}

sub get_Exon_Bases {

	my $ALLELE = shift;

	my $Input_Dir  	   = shift;
	my $Sample_ID      = shift;
	my $Locus_Type     = shift;
	
	my $parent_Dir  = $Input_Dir."/".$Sample_ID."/".$Locus_Type."_nuc";
	my $exonFile    = $parent_Dir."/detail/mappedBasePerExon.txt";
	my %allele_Hash = {};
	my %name_Hash   = {};
	
	open (FILE, $exonFile) || die "Can not open $exonFile.\n";
	
	while (<FILE>) {
		
		chomp $_; next if (!$_);
		my @line = split(/\s+/, $_);
		my ($locus, $allele) = split("_", $line[0]);
		my @allele_sets = split(":", $allele);
		my $key = $allele_sets[0].":".$allele_sets[1];
		if ($key eq $ALLELE){
			
			print "Now key is $key.\n";
			print "Now ALLELE is $ALLELE.\n";

			my @value = ($line[1], $line[2], $line[3]);
			my $sum   = $line[1] + $line[2] + $line[3];
			
			print $line[1], "\n";
			print $line[2], "\n";
			print $line[3], "\n";
			print @value, "\n";
			print "Sum: $sum\n";

			$allele_Hash{$sum} = \@value;
			$name_Hash{$sum}   = $allele;
		}
	}
	
	close FILE;
	
	my @allValues = sort {$b <=> $a} keys %allele_Hash;
	
	my $maxKey = $allValues[0];
	print "Max Key: $maxKey\n";
	return ($allele_Hash{$maxKey}, $name_Hash{$maxKey});

}

sub get_Exon_Reads {

	my $ALLELE = shift;

	my $Input_Dir  	   = shift;
	my $Sample_ID      = shift;
	my $Locus_Type     = shift;
	
	my $parent_Dir  = $Input_Dir."/".$Sample_ID."/".$Locus_Type."_nuc";
	my $exonFile    = $parent_Dir."/detail/mappedReadsPerExon.txt";
	my %allele_Hash = {};
	
	open (FILE, $exonFile) || die "Can not open $exonFile.\n";
	
	while (<FILE>) {
		
		chomp $_; next if (!$_);
		my @line = split(/\s+/, $_);
		my ($locus, $allele) = split("_", $line[0]);
		#my @allele_sets = split(":", $allele);
		#my $key = $allele_sets[0].":".$allele_sets[1];
		if ($allele eq $ALLELE){
			
			print "Now key is $allele.\n";
			print "Now ALLELE is $ALLELE.\n";

			my @value = ($line[1], $line[2], $line[3]);
			my $sum   = $line[1] + $line[2] + $line[3];
			
			print $line[1], "\n";
			print $line[2], "\n";
			print $line[3], "\n";
			print @value, "\n";;
			print "Sum: $sum\n";
			
			my @value = ($line[1], $line[2], $line[3]);
			my $sum   = $line[1] + $line[2] + $line[3];
			$allele_Hash{$sum} = \@value;
			last;
		}
	}
	
	close FILE;
	
	my @allValues = sort {$b <=> $a} keys %allele_Hash;
	
	my $maxKey = $allValues[0];
	print "Max Key: $maxKey\n";
	return $allele_Hash{$maxKey};

}


sub get_mapSingleRefs_Hash {
	
	my $mapSingleRefs_File = shift;
	open (FILE, "<", $mapSingleRefs_File) || die "Can not open $mapSingleRefs_File.\n";
	my %mapSingleRefs_hash;
	
	while (<FILE>){

		chomp $_; next if (!$_);
		my @line = split(/\s+/, $_);
		my ($locus, $allele) = split("_", $line[0]);
		my @allele_sets = split(":", $allele);
		my $key = $allele_sets[0].":".$allele_sets[1];
		my $value = $line[1] + $line[2] + $line[3];
		
		next if ( (exists $mapSingleRefs_hash{$key}) and ($mapSingleRefs_hash{$key} >= $value) );
		$mapSingleRefs_hash{$key} = $value;

	}
	
	close FILE;
	return \%mapSingleRefs_hash;

}

sub get_luminex_validation {
	
	my $sample_ID  = shift;
	my $locus_type = shift;
	my $validation_dir = shift;
	my $validation_file = $validation_dir."/Overall_Validation.csv";
	
	open (FILE, $validation_file) || die "Can not open $validation_file.\n";
	my %luminex_hash;
	
	while (<FILE>){

		chomp $_; 		next if (!$_); 		$_ =~ s/^\s+$// if (grep /^\s+$/, $_);
		my @line = split(',', $_);
		my $luminex_HASH_KEY = $line[0].$line[1];
		push (@{$luminex_hash{$luminex_HASH_KEY}}, $line[2], $line[3]);

	}
	
	close FILE;
	return \%luminex_hash;
}

sub get_LIAI_results {
	
	my $validation_dir = shift;
	my $Input_Dir  	   = shift;
	my $Sample_ID      = shift;
	my $Locus_Type     = shift;
	
	my $parent_Dir = $Input_Dir."/".$Sample_ID."/".$Locus_Type."_nuc";
	my $countFile  = $parent_Dir."/Collapse_Count_Bases.txt";
	my $mapSingleRefsFile = $parent_Dir."/detail/mappedBasePerExon.txt";
	
	my @Return_LIAI_Results; my $NA = "N/A";
	
	if (-z $countFile) {
			
			push @Return_LIAI_Results, $NA,$NA,$NA,$NA,$NA,$NA;
			return(\@Return_LIAI_Results);
	}
	
	
	open (LOCIINPUT, "<", $countFile) || die "Can not open $countFile.\n";
	my @Loci_Input = <LOCIINPUT>;
	
	if ( grep /^$NA/, $Loci_Input[0] ) {
		
		push @Return_LIAI_Results, $NA,$NA,$NA,$NA,$NA,$NA;
		return(\@Return_LIAI_Results);		
	}
	
	my %mappedReads_Hash;
		
	foreach my $each_line(@Loci_Input) {
			
			chomp $each_line;
			my @each_typing = split ('\s+', $each_line);
			my $alleleCombination = $each_typing[0].",".$each_typing[1];
			my $alleleCombinationCounts = $each_typing[2];
			
			if (exists $mappedReads_Hash{$alleleCombination}) {
				
				my $old_counts = $mappedReads_Hash{$alleleCombination};
				my $new_counts = $alleleCombinationCounts;
				
				if ($old_counts > $new_counts) { $mappedReads_Hash{$alleleCombination} = $old_counts; }
				else						   { $mappedReads_Hash{$alleleCombination} = $new_counts; }
				
			}
			
			else {
				
				$mappedReads_Hash{$alleleCombination} = $alleleCombinationCounts;
			}
			
			
	}
	
	my %uniqueAlleleCombination_Hash;
	
	while (my ($uniqueAC, $uniqueAC_Counts) = each %mappedReads_Hash) {
		
		push (@{$uniqueAlleleCombination_Hash{$uniqueAC_Counts}}, $uniqueAC);
	}
	
	my @uniqueAC_MappedReads = sort {$b <=> $a} keys %uniqueAlleleCombination_Hash;
	
	#######################################################
	#
	# Initialize returned values.
	#
	#######################################################
	
	my $AC_1; 	       my $AC_2; 	      my $AC_1_AF; 
	my $AC_1_AF_Ratio; my $AC_1_Evidence; my $AC_1_Counts; 
	
	#######################################################
	#
	# Get Frequency Hash.
	#
	#######################################################

	my $alleleFrequency_Hash = get_allele_frequency($validation_dir);

	#######################################################
	#
	# Get Single Reference Mappability.
	#
	#######################################################

	my $mapSingleRefs_Hash = get_mapSingleRefs_Hash($mapSingleRefsFile);
	
	#######################################################
	#
	# No Ties.
	#
	#######################################################
	
	my @top1_allele_combination = @{$uniqueAlleleCombination_Hash{$uniqueAC_MappedReads[0]}};
	$AC_1_Counts = $uniqueAC_MappedReads[0];
	
	if ( scalar(@top1_allele_combination) == 1) {
		
		$AC_1          = $top1_allele_combination[0];
		$AC_2	       = $NA;
		$AC_1_AF   	   = "NO_TIE";
		$AC_1_AF_Ratio = "NO_TIE";
		$AC_1_Evidence = "NO_TIE";
		
		print "\n\n%%%%%%%%%%%%%%%%%%%% NO TIE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
		
	}
	
	###############################################################
	#
	# use Allele_Frequencies to determine allele combination in ties.
	#
	###############################################################
	
	
	elsif ( scalar(@top1_allele_combination) > 1 ) {
	
		my %alleleComboFrequency;
		my %alleleComboFrequency_Add;
	
		my $infinite_small = 10 ** -22;
		
		foreach my $alleleCombo (@top1_allele_combination) {
			
			my ($allele1, $allele2) = split(',', $alleleCombo, 2);
			my $allele1_FreqKey = $Locus_Type."*".$allele1;
			my $allele2_FreqKey = $Locus_Type."*".$allele2;
			
			print "\n----Now the allele1_FreqKey is $allele1_FreqKey\n";
			print "\n----Now the allele2_FreqKey is $allele2_FreqKey\n";
			
			my $allele1_Freq; my $allele2_Freq; my $alleles_Freq; my $alleles_Freq_Add;
			
			if (exists $alleleFrequency_Hash->{$allele1_FreqKey}){
				
				$allele1_Freq = $alleleFrequency_Hash->{$allele1_FreqKey};
				
				print "\n----Find allele1_FreqKey\n";
			}
			else { $allele1_Freq = 0; print "\n----Do not find allele1_FreqKey\n";}
			
			if (exists $alleleFrequency_Hash->{$allele2_FreqKey}){
				
				$allele2_Freq = $alleleFrequency_Hash->{$allele2_FreqKey};
				print "\n----Find allele2_FreqKey\n";
			}
			else { $allele2_Freq = 0; print "\n----Do not find allele2_FreqKey\n";}
			
			$alleles_Freq 	  = $allele1_Freq * $allele2_Freq;
			$alleles_Freq_Add = $allele1_Freq + $allele2_Freq;
			
			if ($alleles_Freq) { 
				push @{$alleleComboFrequency{$alleles_Freq}}, $alleleCombo;
			}
			else {
				push @{$alleleComboFrequency{$infinite_small}}, $alleleCombo;
			}
			
			if ($alleles_Freq_Add) { 
				push @{$alleleComboFrequency_Add{$alleles_Freq_Add}}, $alleleCombo;
			}
			else {
				push @{$alleleComboFrequency_Add{$infinite_small}}, $alleleCombo;
			}

		}# end of itinerary top 1 allele combination;
		
		my @frequencyValues     = sort {$b <=> $a} keys %alleleComboFrequency;
		my @frequencyValues_Add = sort {$b <=> $a} keys %alleleComboFrequency_Add;
		
		print "##############  PRINT OUT FREQUNCIES  ########################\n\n";
		
		foreach my $each_frequency (@frequencyValues) {
			
			foreach my $each_AC (@{$alleleComboFrequency{$each_frequency}}) {
				print $each_AC," ---->  ", $each_frequency,"\n";
			}
		}
		
		print "\n##############  PRINT OUT FREQUNCIES sum up ########################\n\n";
		
		foreach my $each_frequency_Add (@frequencyValues_Add) {
			
			foreach my $each_AC_Add (@{$alleleComboFrequency_Add{$each_frequency_Add}}) {
				print $each_AC_Add," ---->  ", $each_frequency_Add,"\n";
			}
		}
		
		print "\n##############  END ::: PRINT OUT FREQUNCIES ########################\n\n";
		
		if ( (scalar(@frequencyValues) == 1) && ($frequencyValues[0] == $infinite_small) ) {

			# all Allele_Frequencies = 0;
			
			if ( (scalar(@frequencyValues_Add) == 1) && ($frequencyValues_Add[0] == $infinite_small) ){
				
				print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
				my %single_mapping_HASH;
				
				foreach my $putativeAllele (@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) {
					
					my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
					$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
					print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
					print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
					print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
				}

				print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
				my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
				my @alternative_alleles = @all_alleles[1..$#all_alleles];

				$AC_1          = $all_alleles[0];
				$AC_2	       = join (';', @alternative_alleles);
				$AC_1_AF	= 0;
				$AC_1_AF_Ratio = "INFINITE";
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "SUM UP: All ties have frequncy = 0\n\n";
			
			}
		
			if ( (scalar(@frequencyValues_Add) == 1) && ($frequencyValues_Add[0] != $infinite_small) ){
			
				print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
				my %single_mapping_HASH;
				
				foreach my $putativeAllele (@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) {
					
					my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
					$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
					print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
					print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
					print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
				}

				print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
				my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
				my @alternative_alleles = @all_alleles[1..$#all_alleles];

				$AC_1          = $all_alleles[0];
				$AC_2	       = join (';', @alternative_alleles);
				$AC_1_AF   	   = $frequencyValues_Add[0];
				$AC_1_AF_Ratio = 1;
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "SUM UP: All ties have the same frequncy yet are not 0\n\n";
			
			}
		
			if ( (scalar(@frequencyValues_Add) == 2) && ($frequencyValues_Add[1] == $infinite_small) ){
			
				if ( scalar(@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) == 1 ) {

					$AC_1          = join (';', @{$alleleComboFrequency_Add{$frequencyValues_Add[0]}});
					$AC_2	       = join (';', @{$alleleComboFrequency_Add{$frequencyValues_Add[1]}});
					$AC_1_AF   	   = $frequencyValues_Add[0];
					$AC_1_AF_Ratio = "INFINITE";
					$AC_1_Evidence = "Allele_Frequencies";
			
					print "SUM UP: Two different frequncy values, one = 0.\n\n";
				}

				if ( scalar(@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) > 1 ) {	
				
					print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
					my %single_mapping_HASH;
				
					foreach my $putativeAllele (@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) {
					
						my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
						$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
						print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
						print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
						print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
					}

					print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
					my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
					my @alternative_alleles = @all_alleles[1..$#all_alleles];

					$AC_1          = $all_alleles[0];
					$AC_2	       = join (';', @alternative_alleles);
					$AC_1_AF   	   = $frequencyValues_Add[0];
					$AC_1_AF_Ratio = 1;
					$AC_1_Evidence = "Allele_Frequencies";
			
					print "SUM UP: Two different frequncy values, one = 0, the other have ties.\n\n";
				}	
			}
		
			if ( (scalar(@frequencyValues_Add) == 2) && ($frequencyValues_Add[1] != $infinite_small) ){
				
				if ( scalar(@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) == 1 ) {

					$AC_1          = join (';', @{$alleleComboFrequency_Add{$frequencyValues_Add[0]}});
					$AC_2	       = join (';', @{$alleleComboFrequency_Add{$frequencyValues_Add[1]}});
					$AC_1_AF   	   = $frequencyValues_Add[0];
					$AC_1_AF_Ratio = $frequencyValues_Add[0]/$frequencyValues_Add[1];
					$AC_1_AF_Ratio = sprintf ('%.4f', $AC_1_AF_Ratio);
					$AC_1_Evidence = "Allele_Frequencies";
			
					print "SUM UP: Two different frequncy values, none = 0.\n\n";
				}
				
				if ( scalar(@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) > 1 ) {
				
					print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
					my %single_mapping_HASH;
				
					foreach my $putativeAllele (@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) {
					
						my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
						$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
						print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
						print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
						print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
					}

					print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
					my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
					my @alternative_alleles = @all_alleles[1..$#all_alleles];

					$AC_1          = $all_alleles[0];
					$AC_2	       = join (';', @alternative_alleles);
					$AC_1_AF   	   = $frequencyValues_Add[0];
					$AC_1_AF_Ratio = 1;
					$AC_1_AF_Ratio = sprintf ('%.4f', $AC_1_AF_Ratio);
					$AC_1_Evidence = "Allele_Frequencies";
			
					print "SUM UP: Two different frequncy values, none = 0, the biggest one has ties.\n\n";				
				}
			
			}
		
			if ( scalar(@frequencyValues_Add) > 2 ) {
				
				if ( scalar(@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) == 1 ) {

					my @alternativeAlleleCombo_Add;
			
					for (my $i = 1; $i < scalar(@frequencyValues_Add); $i++) {
				
						push @alternativeAlleleCombo_Add, @{$alleleComboFrequency_Add{$frequencyValues_Add[$i]}};
					}
			
					$AC_1          = join (';', @{$alleleComboFrequency_Add{$frequencyValues_Add[0]}});
					$AC_2	       = join (';', @alternativeAlleleCombo_Add);
					$AC_1_AF   	   = $frequencyValues_Add[0];
					$AC_1_AF_Ratio = $frequencyValues_Add[0]/$frequencyValues_Add[1];
					$AC_1_AF_Ratio = sprintf ('%.4f', $AC_1_AF_Ratio);
					$AC_1_Evidence = "Allele_Frequencies";
			
					print "SUP UP: At lest three different frequncy values.\n\n";
				}
				
				if ( scalar(@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) > 1 ) {
				
					print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
					my %single_mapping_HASH;
				
					foreach my $putativeAllele (@{$alleleComboFrequency_Add{$frequencyValues_Add[0]}}) {
					
						my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
						$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
						print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
						print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
						print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
					}

					print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";

					my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
					my @alternative_alleles = @all_alleles[1..$#all_alleles];

					$AC_1          = $all_alleles[0];
					$AC_2	       = join (';', @alternative_alleles);
					$AC_1_AF   	   = $frequencyValues_Add[0];
					$AC_1_AF_Ratio = 1;
					$AC_1_AF_Ratio = sprintf ('%.4f', $AC_1_AF_Ratio);
					$AC_1_Evidence = "Allele_Frequencies";
			
					print "SUP UP: At lest three different frequncy values, the biggest value has ties.\n\n";	
				
				}
								
			}
		
			
		} #End of if ( (scalar(@frequencyValues) == 1) && ($frequencyValues[0] == $infinite_small) )
		
		if ( (scalar(@frequencyValues) == 1) && ($frequencyValues[0] != $infinite_small) ){
			
			print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
			my %single_mapping_HASH;
				
			foreach my $putativeAllele (@{$alleleComboFrequency{$frequencyValues[0]}}) {
					
				my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
				$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
				print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
				print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
				print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
				}

			print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
			my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
			my @alternative_alleles = @all_alleles[1..$#all_alleles];

			
			$AC_1          = $all_alleles[0];
			$AC_2	       = join (';', @alternative_alleles);
			$AC_1_AF   	   = $frequencyValues[0];
			$AC_1_AF_Ratio = 1;
			$AC_1_Evidence = "Allele_Frequencies";
			
			print "All ties have the same frequncy yet are not zero.\n\n";
			
		}
		
		if ( (scalar(@frequencyValues) == 2) && ($frequencyValues[1] == $infinite_small) ){
			
			if ( scalar(@{$alleleComboFrequency{$frequencyValues[0]}}) == 1 ) {
	
				$AC_1          = join (';', @{$alleleComboFrequency{$frequencyValues[0]}});
				$AC_2	       = join (';', @{$alleleComboFrequency{$frequencyValues[1]}});
				$AC_1_AF   	   = $frequencyValues[0];
				$AC_1_AF_Ratio = "INFINITE";
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "Two different frequncy values, one = 0.\n\n";
			}
			
			if ( scalar(@{$alleleComboFrequency{$frequencyValues[0]}}) > 1 ) {
				
				print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
				my %single_mapping_HASH;
				
				foreach my $putativeAllele (@{$alleleComboFrequency{$frequencyValues[0]}}) {
					
					my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
					$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
					print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
					print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
					print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
				}

				print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
			
				my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
				my @alternative_alleles = @all_alleles[1..$#all_alleles];

			
				$AC_1          = $all_alleles[0];
				$AC_2	       = join (';', @alternative_alleles);
				$AC_1_AF   	   = $frequencyValues[0];
				$AC_1_AF_Ratio = 1;
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "Two different frequncy values, one = 0, the other have ties.\n\n";	
			}
			
		}
		
		if ( (scalar(@frequencyValues) == 2) && ($frequencyValues[1] != $infinite_small) ){
			
			if ( scalar(@{$alleleComboFrequency{$frequencyValues[0]}}) == 1 ) {

				$AC_1          = join (';', @{$alleleComboFrequency{$frequencyValues[0]}});
				$AC_2	       = join (';', @{$alleleComboFrequency{$frequencyValues[1]}});
				$AC_1_AF   	   = $frequencyValues[0];
				$AC_1_AF_Ratio = $frequencyValues[0]/$frequencyValues[1];
				$AC_1_AF_Ratio = sprintf ('%.4f', $AC_1_AF_Ratio);
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "Two different frequncy values, none = 0.\n\n";
			}

			if ( scalar(@{$alleleComboFrequency{$frequencyValues[0]}}) > 1 ) {
				
				print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
				my %single_mapping_HASH;
				
				foreach my $putativeAllele (@{$alleleComboFrequency{$frequencyValues[0]}}) {
					
					my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
					$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
					print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
					print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
					print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
				}

				print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
			
				my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
				my @alternative_alleles = @all_alleles[1..$#all_alleles];

			
				$AC_1          = $all_alleles[0];
				$AC_2	       = join (';', @alternative_alleles);
				$AC_1_AF   	   = $frequencyValues[0];
				$AC_1_AF_Ratio = 1;
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "Two different frequncy values, none = 0, the biggest one have ties.\n\n";	
			}
			
		}
		
		if ( scalar(@frequencyValues) > 2 ) {
			
			if ( scalar(@{$alleleComboFrequency{$frequencyValues[0]}}) == 1 ) {
			
				my @alternativeAlleleCombo;
			
				for (my $i = 1; $i < scalar(@frequencyValues); $i++) {
				
					push @alternativeAlleleCombo, @{$alleleComboFrequency{$frequencyValues[$i]}};
				}
				
				$AC_1          = join (';', @{$alleleComboFrequency{$frequencyValues[0]}});
				$AC_2	       = join (';', @alternativeAlleleCombo);
				$AC_1_AF   	   = $frequencyValues[0];
				$AC_1_AF_Ratio = $frequencyValues[0]/$frequencyValues[1];
				$AC_1_AF_Ratio = sprintf ('%.4f', $AC_1_AF_Ratio);
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "At lest three different frequncy values.\n\n";
			}
			
			if ( scalar(@{$alleleComboFrequency{$frequencyValues[0]}}) > 1 ) {
				
				print "******** USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
				
				my %single_mapping_HASH;
				
				foreach my $putativeAllele (@{$alleleComboFrequency{$frequencyValues[0]}}) {
					
					my ($Putative_Allele_1, $Putative_Allele_2) = split (",", $putativeAllele, 2);
				
					$single_mapping_HASH{$putativeAllele} = $mapSingleRefs_Hash->{$Putative_Allele_1} + $mapSingleRefs_Hash->{$Putative_Allele_2};
					
					print "allele 1 ::: $Putative_Allele_1 ::: $mapSingleRefs_Hash->{$Putative_Allele_1} \n";
					print "allele 2 ::: $Putative_Allele_2 ::: $mapSingleRefs_Hash->{$Putative_Allele_2} \n";
					print "COMBO	::: $putativeAllele    ::: $single_mapping_HASH{$putativeAllele}     \n\n";
			
				}

				print "******** DONE ::: USING SINGLE MAPPING DATA TO DETERMINE THE GENOTYPE *********\n\n";
			
				my @all_alleles = sort { $single_mapping_HASH{$b} <=> $single_mapping_HASH{$a} } keys %single_mapping_HASH;
				my @alternative_alleles = @all_alleles[1..$#all_alleles];

			
				$AC_1          = $all_alleles[0];
				$AC_2	       = join (';', @alternative_alleles);
				$AC_1_AF   	   = $frequencyValues[0];
				$AC_1_AF_Ratio = 1;
				$AC_1_Evidence = "Allele_Frequencies";
			
				print "At lest three different frequncy values, the biggest one have ties.\n\n";	
			}
			
		}
		
	} # end of elsif ( scalar(@top1_allele_combination) > 1 ).
	
	else { die "Error: Top1 allele combo.\n"; }
	
	push @Return_LIAI_Results, $AC_1, 		   $AC_2, 		   $AC_1_AF; 
	push @Return_LIAI_Results, $AC_1_AF_Ratio, $AC_1_Evidence, $AC_1_Counts;

	close LOCIINPUT;
	
	if ( scalar(@Return_LIAI_Results) == 6 ) { return(\@Return_LIAI_Results); }
	else									 { die "Error: wrong elements in return array.\n"; }
	
}

sub compare_with_Luminex {

	my $putative_allele = shift;
	my $luminex_allele  = shift;
	$luminex_allele  =~ s/\s+//g;
	
	my $matching_score = 0;
	my ($luminex_allele_1, $luminex_allele_2) = split(",", $luminex_allele);
	my ($allele_1, $allele_2) 				  = split(",", $putative_allele);
	
	my ($allele_1_family, $allele_1_subtype) = split (':', $allele_1);
	my ($allele_2_family, $allele_2_subtype) = split (':', $allele_2);
	
	my $luminex_allele_1_agree    = 0;
	my $luminex_allele_1_disagree = 0;
	my $luminex_allele_2_agree    = 0;
	my $luminex_allele_2_disagree = 0;
	
	#################################################
	#
	# Hemizygous
	#
	#################################################
	
	if ( $luminex_allele_2 eq 'HEMI' ) {
		
		print "^^^^^ This allele combination is Hemizygous.\n";
		
		my ($luminex_allele_1_family, $luminex_allele_1_subtype) = split (':', $luminex_allele_1);
		
		if ($allele_1 ne $allele_2 ) {
			
			if ( $luminex_allele_1_subtype eq '*' ) {
				
				if ($luminex_allele_1_family eq $allele_1_family) {
					
					if ($luminex_allele_1_family eq $allele_2_family) {
							$matching_score = 1;   }
					else {	$matching_score = 0.5; }	
					
				}
				
				else {
					
					if ($luminex_allele_1_family eq $allele_2_family) {
							$matching_score = 0.5; }
					else {	$matching_score = 0;   }
					
				}
			
			} # luminex_allele_1_subtype eq '*'
			
			else {
				
				if ($luminex_allele_1 eq $allele_1) {
					
					if ($luminex_allele_1 eq $allele_2) {
							$matching_score = 1;   }
					else {	$matching_score = 0.5; }	
					
				}
				
				else {
					
					if ($luminex_allele_1 eq $allele_2) {
							$matching_score = 0.5;   }
					else {	$matching_score = 0; }
					
				}
				
			} # luminex_allele_1_subtype ne '*'	
			
		}# $allele_1 ne $allele_2
		
		elsif ( $allele_1 eq $allele_2 ) {
		
			if ( $luminex_allele_1_subtype eq '*' ) {
				
				if ($luminex_allele_1_family eq $allele_1_family) {
					
					$matching_score = 1;
					
				}
				
				else {
					
					$matching_score = 0; 
					
				}
			
			} # luminex_allele_1_subtype eq '*'
			
			else {
				
				if ($luminex_allele_1 eq $allele_1) {
					
					$matching_score = 1;
					
				}
				
				else {
					
					$matching_score = 0;
					
				}
				
			} # luminex_allele_1_subtype ne '*'
					
		} # $allele_1 eq $allele_2
		
		else { die "Allele Typing Error: Hemizygous.\n"; }
		
		return $matching_score;
		
	} # end of Hemizygous
	
	
	#################################################
	#
	# Homozygous
	#
	#################################################
	
	elsif ( ( $luminex_allele_1 eq $luminex_allele_2 ) and ( $luminex_allele_2 ne 'HEMI' ) ) {
		
		print "^^^^^ This allele combination is Homozygous.\n";
		
		my ($luminex_allele_1_family, $luminex_allele_1_subtype) = split (':', $luminex_allele_1);
		my ($luminex_allele_2_family, $luminex_allele_2_subtype) = split (':', $luminex_allele_2); 
		
		if ( ($luminex_allele_1_subtype eq '*') and ($luminex_allele_2_subtype eq '*') ) {
			
			print "^^^^^ Two *.\n";
			
			if 	  ( ($luminex_allele_1_family eq $allele_1_family) and ($luminex_allele_2_family eq $allele_2_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family eq $allele_2_family) and ($luminex_allele_2_family eq $allele_1_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family eq $allele_1_family) and ($luminex_allele_2_family ne $allele_2_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family eq $allele_2_family) and ($luminex_allele_2_family ne $allele_1_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family ne $allele_1_family) and ($luminex_allele_2_family eq $allele_2_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family ne $allele_2_family) and ($luminex_allele_2_family eq $allele_1_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family ne $allele_1_family) and ($luminex_allele_2_family ne $allele_2_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family ne $allele_2_family) and ($luminex_allele_2_family ne $allele_1_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			else{ die "Unknown homozygous type.\n"; }
		
		} # $luminex_allele_1_subtype = $luminex_allele_2_subtype = '*'
		
		elsif ( ($luminex_allele_1_subtype ne '*') and ($luminex_allele_2_subtype ne '*') ) {
			
			print "^^^^^ No *.\n";
			
			if 	  ( ($luminex_allele_1 eq $allele_1) and ($luminex_allele_2 eq $allele_2) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 eq $allele_2) and ($luminex_allele_2 eq $allele_1) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 eq $allele_1) and ($luminex_allele_2 ne $allele_2) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 eq $allele_2) and ($luminex_allele_2 ne $allele_1) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 ne $allele_1) and ($luminex_allele_2 eq $allele_2) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 ne $allele_2) and ($luminex_allele_2 eq $allele_1) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 ne $allele_1) and ($luminex_allele_2 ne $allele_2) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 ne $allele_2) and ($luminex_allele_2 ne $allele_1) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			else{ die "Unknown homozygous type.\n"; }
		
			
		} # $luminex_allele_1_subtype = $luminex_allele_2_subtype != '*'
		
		else { die "Allele Typing Error: Homozygous.\n"; }
	}
	
	#################################################
	#
	# Heterzygous
	#
	#################################################

	elsif ( ( $luminex_allele_1 ne $luminex_allele_2 ) and ( $luminex_allele_2 ne 'HEMI' ) ){
		
		print "^^^^^ This allele combination is Heterzygous.\n";
		
		my ($luminex_allele_1_family, $luminex_allele_1_subtype) = split (':', $luminex_allele_1);
		my ($luminex_allele_2_family, $luminex_allele_2_subtype) = split (':', $luminex_allele_2); 
		
		if ( ($luminex_allele_1_subtype eq '*') and ($luminex_allele_2_subtype eq '*') ) {
			
			print "^^^^^ Two *.\n";
			
			if 	  ( ($luminex_allele_1_family eq $allele_1_family) and ($luminex_allele_2_family eq $allele_2_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family eq $allele_2_family) and ($luminex_allele_2_family eq $allele_1_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family eq $allele_1_family) and ($luminex_allele_2_family ne $allele_2_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family eq $allele_2_family) and ($luminex_allele_2_family ne $allele_1_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family ne $allele_1_family) and ($luminex_allele_2_family eq $allele_2_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family ne $allele_2_family) and ($luminex_allele_2_family eq $allele_1_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family ne $allele_1_family) and ($luminex_allele_2_family ne $allele_2_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family ne $allele_2_family) and ($luminex_allele_2_family ne $allele_1_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			else{ die "Unknown heterzygous type.\n"; }
		
		} # $luminex_allele_1_subtype = $luminex_allele_2_subtype = '*'
		
		elsif ( ($luminex_allele_1_subtype eq '*') and ($luminex_allele_2_subtype ne '*') ) {
			
			print "^^^^^ * in allele 1.\n";
			
			if 	  ( ($luminex_allele_1_family eq $allele_1_family) and ($luminex_allele_2 eq $allele_2) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family eq $allele_2_family) and ($luminex_allele_2 eq $allele_1) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family eq $allele_1_family) and ($luminex_allele_2 ne $allele_2) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family eq $allele_2_family) and ($luminex_allele_2 ne $allele_1) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family ne $allele_1_family) and ($luminex_allele_2 eq $allele_2) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family ne $allele_2_family) and ($luminex_allele_2 eq $allele_1) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1_family ne $allele_1_family) and ($luminex_allele_2 ne $allele_2) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1_family ne $allele_2_family) and ($luminex_allele_2 ne $allele_1) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			else{ die "Unknown heterzygous type.\n"; }
		
		} # $luminex_allele_1_subtype = '*' != $luminex_allele_2_subtype 
		
		elsif ( ($luminex_allele_1_subtype ne '*') and ($luminex_allele_2_subtype eq '*') ) {
			
			print "^^^^^ * in allele 2.\n";
			
			if 	  ( ($luminex_allele_1 eq $allele_1) and ($luminex_allele_2_family eq $allele_2_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 eq $allele_2) and ($luminex_allele_2_family eq $allele_1_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 eq $allele_1) and ($luminex_allele_2_family ne $allele_2_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 eq $allele_2) and ($luminex_allele_2_family ne $allele_1_family) )
				{$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 ne $allele_1) and ($luminex_allele_2_family eq $allele_2_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 ne $allele_2) and ($luminex_allele_2_family eq $allele_1_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 ne $allele_1) and ($luminex_allele_2_family ne $allele_2_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 ne $allele_2) and ($luminex_allele_2_family ne $allele_1_family) )
				{$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			else{ die "Unknown heterzygous type.\n"; }
		
		} # $luminex_allele_1_subtype != $luminex_allele_2_subtype = '*'
		
		elsif ( ($luminex_allele_1_subtype ne '*') and ($luminex_allele_2_subtype ne '*') ) {
			
			print "^^^^^ No *.\n";
			
			if 	  ( ($luminex_allele_1 eq $allele_1) and ($luminex_allele_2 eq $allele_2) )
				{	$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 eq $allele_2) and ($luminex_allele_2 eq $allele_1) )
				{	$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 eq $allele_1) and ($luminex_allele_2 ne $allele_2) )
				{	$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 eq $allele_2) and ($luminex_allele_2 ne $allele_1) )
				{	$luminex_allele_1_agree = 1; $luminex_allele_1_disagree = 0;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 ne $allele_1) and ($luminex_allele_2 eq $allele_2) )
				{	$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 ne $allele_2) and ($luminex_allele_2 eq $allele_1) )
				{	$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 1; $luminex_allele_2_disagree = 0;}
			elsif ( ($luminex_allele_1 ne $allele_1) and ($luminex_allele_2 ne $allele_2) )
				{	$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			elsif ( ($luminex_allele_1 ne $allele_2) and ($luminex_allele_2 ne $allele_1) )
				{	$luminex_allele_1_agree = 0; $luminex_allele_1_disagree = 1;$luminex_allele_2_agree = 0; $luminex_allele_2_disagree = 1;}
			else{ die "Unknown homozygous type.\n"; }
		
			
		} # $luminex_allele_1_subtype != * != $luminex_allele_2_subtype 
		
		else { die "Allele Typing Error: Heterzygous.\n"; }
		
	}
	
	#################################################
	#
	# Zygosity unknown
	#
	#################################################

	else { die "Can not identify zygosity. Quit.\n"; }
		
	if ($luminex_allele_1_agree    * $luminex_allele_2_agree)    { $matching_score = 1;   }
	if ($luminex_allele_1_agree    * $luminex_allele_2_disagree) { $matching_score = 0.5; }
	if ($luminex_allele_1_disagree * $luminex_allele_2_agree)    { $matching_score = 0.5; }
	if ($luminex_allele_1_disagree * $luminex_allele_2_disagree) { $matching_score = 0;   }

	return $matching_score;
}



