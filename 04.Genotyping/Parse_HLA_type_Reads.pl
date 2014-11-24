#!/share/apps/perl/perl-5.18.1/bin/perl

# Author: Alex Fu (zfu@liai.org)
# Date: April 9th., 2014
# Version: 2.0

# This version use new names of reads count file.

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

my $HLA_Dir; my $OUTPUT_Dir;

GetOptions( "run-data-path|r=s" => \$HLA_Dir,  
            "output-path|o=s"	=> \$OUTPUT_Dir, );

chomp $HLA_Dir; chomp $OUTPUT_Dir;


my $trimmed_reads_file     = $HLA_Dir."/countTrimmedReads.txt";
my %trimmed_reads_hash;
my $raw_reads_file         = $HLA_Dir."/countRawReads.txt";
my %raw_reads_hash;
my $sample_locus_file      = $HLA_Dir."/Sample_Locus.csv";
my %HLA_type;

#######################################################################
#
# Create some hashes.
#
#######################################################################

print "***********************************************************************\n\n";

open (RAWREADS, $raw_reads_file) || die "Can not open $raw_reads_file.\n";
while (<RAWREADS>){
	
	chomp;
	my @each_fastq_file = split (/\s+/, $_);
	my @each_fastq_file_name = split ('_', $each_fastq_file[0]);
	$raw_reads_hash{$each_fastq_file_name[0]} = $each_fastq_file[1];
	
	print $each_fastq_file_name[0],"===>",$raw_reads_hash{$each_fastq_file_name[0]},"\n";
}
close RAWREADS;

print "***********************************************************************\n\n";

open (TRIMMED, $trimmed_reads_file) || die "Can not open $trimmed_reads_file.\n";
while (<TRIMMED>){
	
	chomp;
	my @each_fastq_file = split (/\s+/, $_);
	my @each_fastq_file_name = split ('_', $each_fastq_file[0]);
	$trimmed_reads_hash{$each_fastq_file_name[0]} = $each_fastq_file[1];
	
	print $each_fastq_file_name[0],"===>",$trimmed_reads_hash{$each_fastq_file_name[0]},"\n";	
}
close TRIMMED;

print "***********************************************************************\n\n";

open (SAMPLELOCUS, $sample_locus_file) || die "Can not open $sample_locus_file.\n";
while (<SAMPLELOCUS>){
	
	chomp;
	
	my ($sampleID, $Loci_Type) = split (',', $_);
	
	next if (!$sampleID);
	next if (!$Loci_Type);

	push (@{$HLA_type{$Loci_Type}}, $sampleID);
	
	print $sampleID,"===>",$Loci_Type,"\n";	
}
close SAMPLELOCUS;

print "***********************************************************************\n\n";

#######################################################################
#
# Print header of overall detail file.
#
#######################################################################

my $overall_detail_file  = $OUTPUT_Dir."/Overall_Detail_Reads.txt";

open (OVERALL, "> $overall_detail_file") || die "Can not open $overall_detail_file.\n";

print OVERALL "Sample_ID\t";
print OVERALL "Locus\t";
print OVERALL "Amount_of_Trimmed_Reads\t";
print OVERALL "Amount_of_Raw_Reads\t";
print OVERALL "Trimmed_Reads_Percentage\t";
print OVERALL "Allele_1\t";
print OVERALL "Allele_2\t";
print OVERALL "Allele_Combinations_Mapped_Reads\n";


#######################################################################
#
# Check each sample per locus.
#
#######################################################################


while (my ($Loci, $Loci_Samples) = each %HLA_type) {
	
	print "******************************************\n";
	print "**\n";
	print "Now dealing with Loci $Loci:\n";
	print "**\n";
	print "******************************************\n\n";
		
	LOCI: foreach my $sample_ID (sort {$a cmp $b} @{$Loci_Samples}) {
		
		print "##### Now parse sample: $sample_ID\n";
		print $HLA_Dir, "\n";
		print $sample_ID, "\n";
		print $Loci, "\n";
		
		chomp $Loci;
		chomp $sample_ID;
		
		my $parent_Dir = $HLA_Dir."/".$sample_ID."/".$Loci."_nuc";
		print $parent_Dir, "\n";
		
		if (! -d $parent_Dir) {
			#system("mkdir -p $parent_Dir");
		        print $parent_Dir, " does not exist!\n";
		}	 

		
		my $countFile          = $parent_Dir."/AlleleCombo_MappedReads.txt";
		my $truncatedCountFile = $parent_Dir."/AlleleCombo_MappedReads.txt";
		
		if (! -e $countFile) { $countFile = $truncatedCountFile; }
		
		my $sortedCountFile    = $parent_Dir."/Sorted_Count_Reads.txt";
		my $collapseCountFile  = $parent_Dir."/Collapse_Count_Reads.txt";
		
		my $trimmedReadsCounts;
		my $rawReadsCounts;
		my $readsFilteredPercentage;
		
		open (SORTED,   ">", $sortedCountFile)   || die "Can not open $sortedCountFile.\n";
	        open (COLLAPSE, ">", $collapseCountFile) || die "Can not open $collapseCountFile.\n";
	    
	        my $sample_ID_Locus = $sample_ID.$Loci;
		print "##### The sample ID and Locus is: $sample_ID_Locus.\n";
		
		foreach my $Sample_Locus_Name (keys %raw_reads_hash) {
				
			if (grep /$Sample_Locus_Name/, $sample_ID_Locus) {
					
					$trimmedReadsCounts      = $trimmed_reads_hash{$Sample_Locus_Name};
					$rawReadsCounts          = $raw_reads_hash{$Sample_Locus_Name};
					$readsFilteredPercentage = $trimmedReadsCounts / $rawReadsCounts;
					$readsFilteredPercentage = sprintf ('%.4f', $readsFilteredPercentage);
			
			}
		} #end of reads number hash;
		
		#######################################################################
		#
		# If count_2.txt file did not exist.
		#
		#######################################################################
		
		if ( (! -e $countFile) || (-z $countFile) ){
			
			print OVERALL  "$sample_ID\t";
			print OVERALL  "$Loci\t";
			print OVERALL  $trimmedReadsCounts,"\t";
			print OVERALL  $rawReadsCounts,"\t";	
			print OVERALL  $readsFilteredPercentage, "\t";
			print OVERALL  "N/A\tN/A\tN/A\n";
			print SORTED   "N/A\tN/A\tN/A\n";
			print COLLAPSE "N/A\tN/A\tN/A\tN/A\tN/A\tN/A\n";
			
			close SORTED; close COLLAPSE;
			next LOCI;
		}
		
		#######################################################################
		#
		# If count_2.txt file do exist.
		#
		#######################################################################
		
		open (LOCIINPUT, "<", $countFile) || die "Can not open $countFile.\n";
		my @Loci_Input = <LOCIINPUT>;
		my %mappedReads_Hash;
		
		foreach my $each_line(@Loci_Input) {
			
			chomp $each_line;
			my ($each_allele, $allele_count) = split ('\s+', $each_line, 2);
			push (@{$mappedReads_Hash{$allele_count}}, $each_allele);
			
		}
		
		LOG: foreach my $COUNTS (sort {$b <=> $a} keys %mappedReads_Hash) {
			
				my %sortedCountHash;
				my %collapseCountHash;
				
		ALLELE: foreach my $alleles (@{$mappedReads_Hash{$COUNTS}}) {
					
				my @allele_combo = split ('_', $alleles);	
					
				my @allele_1_Label = split (':', $allele_combo[1]);
	        		my @allele_2_Label = split (':', $allele_combo[3]);
	        		
	        		# Initilize the keys and values of %sortedCountHash and %collapseCountHash
	        		
	        		my $sortedCountHash_KEY; 
	        		my $sortedCountHash_allele_1; 
	        		my $sortedCountHash_allele_2;
	        		my $collapseCountHash_KEY; 
	        		my $collapseCountHash_allele_1; 
	        		my $collapseCountHash_allele_2;
	        		
	        		# Initilize the keys and values of %sortedCountHash and %collapseCountHash
	        		
	        		my $allele_1_Set1; my $allele_1_Set2; my $allele_1_Set3; my $allele_1_Set4;
	        		my $allele_2_Set1; my $allele_2_Set2; my $allele_2_Set3; my $allele_2_Set4;
	        		
	        		$allele_1_Set1 = '1'.$allele_1_Label[0];
	        		$allele_2_Set1 = '1'.$allele_2_Label[0];
	        		
	        		if ( length($allele_1_Label[1]) == 3 ) 	    { $allele_1_Set2 = '0'.$allele_1_Label[1];  }
	        		elsif ( length($allele_1_Label[1]) == 2 ) 	{ $allele_1_Set2 = '00'.$allele_1_Label[1]; }
	        		else								  	    { $allele_1_Set2 = $allele_1_Label[1];      }
	        		
	        		if ( length($allele_2_Label[1]) == 3 ) 	 { $allele_2_Set2 = '0'.$allele_2_Label[1]; }
	        		elsif ( length($allele_2_Label[1]) == 2 ){ $allele_2_Set2 = '00'.$allele_2_Label[1]; }
	        		else								  	 { $allele_2_Set2 = $allele_2_Label[1];     }
	        		
	        		if ( scalar(@allele_1_Label) <= 2 )   	 { $allele_1_Set3 = '0000';     			   }
	        		elsif ( length($allele_1_Label[2]) == 3 ){ $allele_1_Set3 = '0'.$allele_1_Label[2]; }
	        		elsif ( length($allele_1_Label[2]) == 2 ){ $allele_1_Set3 = '00'.$allele_1_Label[2]; }
	        		else									 { $allele_1_Set3 = $allele_1_Label[2];     }
	        		
	        		if ( scalar(@allele_2_Label) <= 2 )   	 { $allele_2_Set3 = '0000';     			   }
	        		elsif ( length($allele_2_Label[2]) == 3 ){ $allele_2_Set3 = '0'.$allele_2_Label[2]; }
	        		elsif ( length($allele_2_Label[2]) == 2 ){ $allele_2_Set3 = '00'.$allele_2_Label[2]; }
	        		else									 { $allele_2_Set3 = $allele_2_Label[2];     }
	        		
	        		if ( scalar(@allele_1_Label) <= 3 )   	 { $allele_1_Set4 = '0000';     			   }
	        		elsif ( length($allele_1_Label[3]) == 3 ){ $allele_1_Set4 = '0'.$allele_1_Label[3]; }
	        		elsif ( length($allele_1_Label[3]) == 2 ){ $allele_1_Set4 = '00'.$allele_1_Label[3]; }
	        		else									 { $allele_1_Set4 = $allele_1_Label[3];     }
	        		
	        		if ( scalar(@allele_2_Label) <= 3 )   	 { $allele_2_Set4 = '0000';     			   }
	        		elsif ( length($allele_2_Label[3]) == 3 ){ $allele_2_Set4 = '0'.$allele_2_Label[3]; }
	        		elsif ( length($allele_2_Label[3]) == 2 ){ $allele_2_Set4 = '00'.$allele_2_Label[3]; }
	        		else									 { $allele_2_Set4 = $allele_2_Label[3];     }
	        		
	        		# Assign the keys and values of %sortedCountHash and %collapseCountHash
	        		
	        		if ( $allele_1_Set1 < $allele_2_Set1 ) {
	        			
	        			$sortedCountHash_KEY = $allele_1_Set1.$allele_1_Set2.$allele_1_Set3.$allele_1_Set4;
	        			$sortedCountHash_KEY = $sortedCountHash_KEY.'.';
	        			$sortedCountHash_KEY = $sortedCountHash_KEY.$allele_2_Set1.$allele_2_Set2.$allele_2_Set3.$allele_2_Set4;
	        			
	        			$collapseCountHash_KEY = $allele_1_Set1.$allele_1_Set2.'.'.$allele_2_Set1.$allele_2_Set2;
	        			
	        			$sortedCountHash_allele_1 = $allele_combo[1]; 
	        		        $sortedCountHash_allele_2 = $allele_combo[3];
	        		
	        		        $collapseCountHash_allele_1 = $allele_1_Label[0].":".$allele_1_Label[1]; 
	        		        $collapseCountHash_allele_2 = $allele_2_Label[0].":".$allele_2_Label[1];
	        			
	        		}
	        		
	        		elsif ( $allele_1_Set1 > $allele_2_Set1 ) {
	        		
	        			$sortedCountHash_KEY = $allele_2_Set1.$allele_2_Set2.$allele_2_Set3.$allele_2_Set4;
	        			$sortedCountHash_KEY = $sortedCountHash_KEY.'.';
	        			$sortedCountHash_KEY = $sortedCountHash_KEY.$allele_1_Set1.$allele_1_Set2.$allele_1_Set3.$allele_1_Set4;
	        			
	        			$collapseCountHash_KEY = $allele_2_Set1.$allele_2_Set2.'.'.$allele_1_Set1.$allele_1_Set2;
	        			
	        			$sortedCountHash_allele_1 = $allele_combo[3]; 
	        		        $sortedCountHash_allele_2 = $allele_combo[1];
	        		
	        		        $collapseCountHash_allele_1 = $allele_2_Label[0].":".$allele_2_Label[1]; 
	        		        $collapseCountHash_allele_2 = $allele_1_Label[0].":".$allele_1_Label[1];
	        			
	        		}
	        		
	        		elsif ( $allele_1_Set1 == $allele_2_Set1 ) {
	        			
	        			if ( $allele_1_Set2 <= $allele_2_Set2 ) {
	        				
	        				$sortedCountHash_KEY = $allele_1_Set1.$allele_1_Set2.$allele_1_Set3.$allele_1_Set4;
	        				$sortedCountHash_KEY = $sortedCountHash_KEY.'.';
	        				$sortedCountHash_KEY = $sortedCountHash_KEY.$allele_2_Set1.$allele_2_Set2.$allele_2_Set3.$allele_2_Set4;
	        			
	        				$collapseCountHash_KEY = $allele_1_Set1.$allele_1_Set2.'.'.$allele_2_Set1.$allele_2_Set2;
	        			
	        				$sortedCountHash_allele_1 = $allele_combo[1]; 
	        		    	$sortedCountHash_allele_2 = $allele_combo[3];
	        		
	        		    	$collapseCountHash_allele_1 = $allele_1_Label[0].":".$allele_1_Label[1]; 
	        		    	$collapseCountHash_allele_2 = $allele_2_Label[0].":".$allele_2_Label[1];
	        		        				
	        			}
	        			
	        			elsif ( $allele_1_Set2 > $allele_2_Set2 ) {
	        			
	        				$sortedCountHash_KEY = $allele_2_Set1.$allele_2_Set2.$allele_2_Set3.$allele_2_Set4;
	        				$sortedCountHash_KEY = $sortedCountHash_KEY.'.';
	        				$sortedCountHash_KEY = $sortedCountHash_KEY.$allele_1_Set1.$allele_1_Set2.$allele_1_Set3.$allele_1_Set4;
	        			
	        				$collapseCountHash_KEY = $allele_2_Set1.$allele_2_Set2.'.'.$allele_1_Set1.$allele_1_Set2;
	        			
	        				$sortedCountHash_allele_1 = $allele_combo[3]; 
	        		    	$sortedCountHash_allele_2 = $allele_combo[1];
	        		
	        		    	$collapseCountHash_allele_1 = $allele_2_Label[0].":".$allele_2_Label[1]; 
	        		    	$collapseCountHash_allele_2 = $allele_1_Label[0].":".$allele_1_Label[1];
	        			
	        			}
	        			
	        			else { die "allele_1_Set2 and allele_2_Set2 have error.\n"; }
	        		
	        		}# end of elsif ( $allele_1_Set1 == $allele_2_Set1 )
	        		
	        		else { die "allele_1_Set1 and allele_2_Set1 have error.\n"; }
	        		
	        		##############################
	        		#
	        		# Put key-value pair to hash.
	        		#
	        		##############################
	        		
	        		if (! exists $sortedCountHash{$sortedCountHash_KEY}) {
	        			
	        			push @{$sortedCountHash{$sortedCountHash_KEY}}, $sortedCountHash_allele_1;
	        			push @{$sortedCountHash{$sortedCountHash_KEY}}, $sortedCountHash_allele_2;
	        			
	        		}
	        		
	        		if (! exists $collapseCountHash{$collapseCountHash_KEY}) {
	        			
	        			push @{$collapseCountHash{$collapseCountHash_KEY}}, $collapseCountHash_allele_1;
	        			push @{$collapseCountHash{$collapseCountHash_KEY}}, $collapseCountHash_allele_2;
	        			
	        		}
	        		
				} # end of ALLELE;
			
			foreach my $printoutSortedCountKey (sort {$a <=> $b} keys %sortedCountHash) {
				
				print "Sorted Key: $printoutSortedCountKey\n";
				
				print OVERALL "$sample_ID\t";
				print OVERALL "$Loci\t";
				print OVERALL $trimmedReadsCounts,"\t";
				print OVERALL $rawReadsCounts,"\t";	
				print OVERALL $readsFilteredPercentage, "\t";
				print OVERALL " ".$sortedCountHash{$printoutSortedCountKey}->[0],"\t";
				print OVERALL " ".$sortedCountHash{$printoutSortedCountKey}->[1],"\t";
				print OVERALL $COUNTS,"\n";
			    
			        print SORTED $sortedCountHash{$printoutSortedCountKey}->[0],"\t";
				print SORTED $sortedCountHash{$printoutSortedCountKey}->[1],"\t";
				print SORTED $COUNTS,"\n";
			   
			}
			
			foreach my $printoutCollapseCountKey (sort {$a <=> $b} keys %collapseCountHash) {
				
				print "Collapse Key: $printoutCollapseCountKey\n";
				    
			        print COLLAPSE $collapseCountHash{$printoutCollapseCountKey}->[0],"\t";
				print COLLAPSE $collapseCountHash{$printoutCollapseCountKey}->[1],"\t";
				print COLLAPSE $COUNTS,"\t";
				print COLLAPSE $trimmedReadsCounts,"\t";
				print COLLAPSE $rawReadsCounts,"\t";
				print COLLAPSE $readsFilteredPercentage,"\n";
			
			}
		
		} # end of LOG: each $COUNTS.
		
	close SORTED; close COLLAPSE; close LOCIINPUT;	
		
	}# LOCI END: end of each sample.
	
	print "******************************************\n";
	print "**\n";
	print "HLA_$Loci was done.\n";
	print "**\n";
	print "******************************************\n\n";
	
}# end of each Loci type: A, B, C and etc.

close OVERALL;
