#!/share/apps/perl/perl-5.18.1/bin/perl

# Author: Alex Fu (zfu@liai.org)
# Date: April 23th., 2013
# Updated: April 23th., 2013
# Version: 1.0.0

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
use File::Basename;
use File::Find::Rule;
use File::Find;
use File::Copy qw(move copy);
use File::Copy::Recursive qw(dirmove);
use File::Path qw(make_path remove_tree rmtree);
use File::stat;
use Storable qw(dclone);
use Getopt::Long;

my $HLA_Dir; my %Sample_ID_Hash; my $equipment = "M00688";

GetOptions( "run-data-path|r=s" => \$HLA_Dir,
			"equipment-ID|e=s"  => \$equipment,
		  );

chomp $HLA_Dir;
my $removed_Dir = $HLA_Dir."/emptyTrimmedReadsRemoved";
my $merged_Dir  = $HLA_Dir."/emptyTrimmedReadsRemoved/merged";

#system("mkdir -p $removed_Dir");
#system("mkdir -p $merged_Dir");

system("id");
system("echo $removed_Dir");
system("echo $merged_Dir");
mkdir("$removed_Dir") or die "Cannot be created because: $!";
mkdir("$merged_Dir") or die "Cannot be created because: $!";

opendir(my $HLA_DIR_HUNDLE, $HLA_Dir) || die "Can not open folder $HLA_Dir.\n";


FILE: while (readdir $HLA_DIR_HUNDLE) {

	next if ($_ eq '.'); next if ($_ eq '..');
	my $fileName = $_; chomp $fileName;
	next if (!grep/.fastq$/, $fileName);
	next if (grep/single/, $fileName);
	next if (grep/Undetermined/, $fileName);
    
    print "\n\n\n******************************************************\n";
    print "Now select: $HLA_Dir/$fileName\n";
    
    my @sample_name = split("_", $fileName);
    my $sample_ID = $sample_name[0];
    print "The sample ID is: $sample_ID\n";
        
    next FILE if (exists $Sample_ID_Hash{$sample_ID});
    
    my $forward_reads; my $reverse_reads;
    my $merged_reads   = $merged_Dir."/".$sample_ID."_merged.fastq";
    my $new_forward_reads = $removed_Dir."/".$sample_ID."_R1_001.fastq";
    my $new_reverse_reads = $removed_Dir."/".$sample_ID."_R2_001.fastq";
    
    if (grep/_R1_/, $fileName) {
    	
    	$forward_reads = $HLA_Dir."/".$fileName;
    	$reverse_reads = $forward_reads;
    	$reverse_reads =~ s/_R1_/_R2_/g; 
    }
    
    if (grep/_R2_/, $fileName) {
    		
    	$reverse_reads = $HLA_Dir."/".$fileName;
    	$forward_reads = $reverse_reads;
    	$forward_reads =~ s/_R2_/_R1_/g; 
    }
    
    print "*** Forward: $forward_reads\n";
    print "*** Reverse: $reverse_reads\n";
    
    my %forward_reads_ID_Hash; my %reverse_reads_ID_Hash; 
    
    open (FORWARD, "< $forward_reads") || die "Can not open $forward_reads.\n";
    my @all_forward_reads = <FORWARD>;
    for (my $i = 0; $i < scalar(@all_forward_reads); $i++) {
    	
    	if (grep /$equipment/, $all_forward_reads[$i]) {
    		
    		my ($main_info, $sub_info) = split (/\s+/, $all_forward_reads[$i], 2);
    		my @reads_ID = split (':', $main_info);
    		
    		#print "Reads ID: $all_forward_reads[$i]\n";
    		
    		my $x; my $y;
    		
    		if (length($reads_ID[5]) == 5) {$x = $reads_ID[5];}
    		if (length($reads_ID[5]) == 4) {$x = "0".$reads_ID[5];}
    		if (length($reads_ID[5]) == 3) {$x = "00".$reads_ID[5];}
    		if (length($reads_ID[5]) == 2) {$x = "000".$reads_ID[5];}
    		if (length($reads_ID[5]) == 1) {$x = "0000".$reads_ID[5];}
    		
    		if (length($reads_ID[6]) == 5) {$y = $reads_ID[6];}
    		if (length($reads_ID[6]) == 4) {$y = "0".$reads_ID[6];}
    		if (length($reads_ID[6]) == 3) {$y = "00".$reads_ID[6];}
    		if (length($reads_ID[6]) == 2) {$y = "000".$reads_ID[6];}
    		if (length($reads_ID[6]) == 1) {$y = "0000".$reads_ID[6];}
    		
    		my $coordinates = $reads_ID[4].".".$x.$y;
    		
    		#print "Forward coordinates: $coordinates\n";
    		
    		if (exists $forward_reads_ID_Hash{$coordinates}) {
    			
    			die "forward coordinates error!\n";
    		}
    		else {
    			
    			push @{$forward_reads_ID_Hash{$coordinates}}, $all_forward_reads[$i];
    			push @{$forward_reads_ID_Hash{$coordinates}}, $all_forward_reads[$i + 1];
    			push @{$forward_reads_ID_Hash{$coordinates}}, $all_forward_reads[$i + 2];
    			push @{$forward_reads_ID_Hash{$coordinates}}, $all_forward_reads[$i + 3];
    			
    			
    		}
    		
    	}# if sequence ID.
    	
    }# end of each forward read.
    close FORWARD;
    
    open (REVERSE, "< $reverse_reads") || die "Can not open $reverse_reads.\n";
    my @all_reverse_reads = <REVERSE>;
    for (my $i = 0; $i < scalar(@all_reverse_reads); $i++) {
    	
    	if (grep /$equipment/, $all_reverse_reads[$i]) {
    		
    		my ($main_info, $sub_info) = split (/\s+/, $all_reverse_reads[$i], 2);
    		my @reads_ID = split (':', $main_info);
    		
    		#print "Reads ID: $all_reverse_reads[$i]\n";
    		
    		my $x; my $y;
    		
    		if (length($reads_ID[5]) == 5) {$x = $reads_ID[5];}
    		if (length($reads_ID[5]) == 4) {$x = "0".$reads_ID[5];}
    		if (length($reads_ID[5]) == 3) {$x = "00".$reads_ID[5];}
    		if (length($reads_ID[5]) == 2) {$x = "000".$reads_ID[5];}
    		if (length($reads_ID[5]) == 1) {$x = "0000".$reads_ID[5];}
    		
    		if (length($reads_ID[6]) == 5) {$y = $reads_ID[6];}
    		if (length($reads_ID[6]) == 4) {$y = "0".$reads_ID[6];}
    		if (length($reads_ID[6]) == 3) {$y = "00".$reads_ID[6];}
    		if (length($reads_ID[6]) == 2) {$y = "000".$reads_ID[6];}
    		if (length($reads_ID[6]) == 1) {$y = "0000".$reads_ID[6];}
    		
    		my $coordinates = $reads_ID[4].".".$x.$y;
    		
    		#print "Reverse coordinates: $coordinates\n";
    		
    		if (exists $reverse_reads_ID_Hash{$coordinates}) {
    			
    			die "Reverse coordinates error!\n";
    		}
    		else {
    			
    			push @{$reverse_reads_ID_Hash{$coordinates}}, $all_reverse_reads[$i];
    			push @{$reverse_reads_ID_Hash{$coordinates}}, $all_reverse_reads[$i + 1];
    			push @{$reverse_reads_ID_Hash{$coordinates}}, $all_reverse_reads[$i + 2];
    			push @{$reverse_reads_ID_Hash{$coordinates}}, $all_reverse_reads[$i + 3];
    			
    			
    		}
    		
    	}# if sequence ID.
    	
    }# end of each forward read.
    close REVERSE;
    
    my @sorted_forward_reads = sort {$b <=> $a} keys %forward_reads_ID_Hash;
    my @sorted_reverse_reads = sort {$b <=> $a} keys %reverse_reads_ID_Hash;
    
    my $forward_counter = scalar(@sorted_forward_reads) - 1;
    my $reverse_counter = scalar(@sorted_reverse_reads) - 1;
    my $total_reads = scalar(@sorted_forward_reads) + scalar(@sorted_reverse_reads) + 2;
    
    print "*** Number of forward reads: ", $forward_counter + 1, "\n";
    print "*** Number of reverse reads: ", $reverse_counter + 1, "\n";
    print "*** Number of total reads  : $total_reads\n";
    
    
    open (MERGED, "> $merged_reads")      || die "Can not open $merged_reads.\n";
    open (NEWFR,  "> $new_forward_reads") || die "Can not open $new_forward_reads.\n";
    open (NEWRR,  "> $new_reverse_reads") || die "Can not open $new_reverse_reads.\n";
    
    foreach my $readsDataKeys (@sorted_forward_reads) {
    	
    	my $forwardReadsLength = length($forward_reads_ID_Hash{$readsDataKeys}->[1]);
    	my $reverseReadsLength = length($reverse_reads_ID_Hash{$readsDataKeys}->[1]);
    	
    	if ( ($forwardReadsLength > 49) && ($reverseReadsLength > 49) ) {
    		
    		print NEWFR $forward_reads_ID_Hash{$readsDataKeys}->[0];
    		print NEWFR $forward_reads_ID_Hash{$readsDataKeys}->[1];
    		print NEWFR $forward_reads_ID_Hash{$readsDataKeys}->[2];
    		print NEWFR $forward_reads_ID_Hash{$readsDataKeys}->[3];
    		
    		print NEWRR $reverse_reads_ID_Hash{$readsDataKeys}->[0];
    		print NEWRR $reverse_reads_ID_Hash{$readsDataKeys}->[1];
    		print NEWRR $reverse_reads_ID_Hash{$readsDataKeys}->[2];
    		print NEWRR $reverse_reads_ID_Hash{$readsDataKeys}->[3];
    		
    		print MERGED $forward_reads_ID_Hash{$readsDataKeys}->[0];
    		print MERGED $forward_reads_ID_Hash{$readsDataKeys}->[1];
    		print MERGED $forward_reads_ID_Hash{$readsDataKeys}->[2];
    		print MERGED $forward_reads_ID_Hash{$readsDataKeys}->[3];
    		
    		print MERGED $reverse_reads_ID_Hash{$readsDataKeys}->[0];
    		print MERGED $reverse_reads_ID_Hash{$readsDataKeys}->[1];
    		print MERGED $reverse_reads_ID_Hash{$readsDataKeys}->[2];
    		print MERGED $reverse_reads_ID_Hash{$readsDataKeys}->[3];
    		
    	}
    	
    	
    }
    
    close MERGED;
    close NEWFR;
    close NEWRR;
    
    
    if (! -e $merged_reads) { die "File $merged_reads did not exist.\n"; }
    else { 
    	
    	system ("echo 'Totol lines in OLD forward reads file'");
    	system ("wc -l $forward_reads");
    	system ("echo 'Totol lines in OLD reverse reads file'");
    	system ("wc -l $reverse_reads");
    	system ("echo 'Totol lines in NEW forward reads file'");
    	system ("wc -l $new_forward_reads");
    	system ("echo 'Totol lines in NEW reverse reads file'");
    	system ("wc -l $new_reverse_reads");
    	system ("echo 'Totol lines in merged reads file'");
    	system ("wc -l $merged_reads");
    	print "*** DONE: $merged_reads.\n"; }
    
}# end of open HLA sequence folder.
