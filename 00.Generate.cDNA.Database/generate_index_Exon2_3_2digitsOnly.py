import sys
import operator
import re
import os
import glob

input_index_dir = "/Bioinformatics/Users/zfu/HLA_Typing/HLA_cDNA_Database/IMGT.Release.3.12.0/Index/"
exonBoundary_dir = "/Bioinformatics/Users/zfu/HLA_Typing/HLA_cDNA_Database/IMGT.Release.3.12.0/ExonBoundary/"
output_index_dir = "/Bioinformatics/Users/zfu/HLA_Typing/HLA_cDNA_Database/IMGT.Release.3.12.0/Index.Test.Single/Exon_2_3/"

Loci_List = ['A', 'B', 'C', 'DRB1', 'DPA1', 'DPB1', 'DQA1', 'DQB1']
#Loci_List = ['DPA1', 'DPB1', 'DQA1', 'DQB1']

for locus in Loci_List:
    
    if not os.path.exists(output_index_dir + locus + "_nuc"):
        os.makedirs(output_index_dir + locus + "_nuc") 
    
    input_id_FH = open(input_index_dir + locus + "_nuc/" + locus + "_nuc_id.txt", "r")
    input_seq_FH = open(input_index_dir + locus + "_nuc/" + locus + "_nuc_seq.txt", "r")
    output_id_FH = open(output_index_dir + locus + "_nuc/" + locus + "_nuc_id.txt", "w")
    output_seq_FH = open(output_index_dir + locus + "_nuc/" + locus + "_nuc_seq.txt", "w")
    
    exonBoundary_dict = {}
    with open(exonBoundary_dir + locus + "_nuc_EX.txt", "r") as input_exonBoundary:
        for line in input_exonBoundary:
            (ID, Boundary) = re.split("\s+", line.strip("\n"), 1)
            Boundary_List = re.split(",", Boundary.strip(","))
            exonBoundary_dict[ID] = Boundary_List
    
    allele_ID_list = input_id_FH.readlines()
    allele_sequence_list = input_seq_FH.readlines()
    
    for i in range(0, len(allele_ID_list)):
        
        allele_name_list = re.split(":", allele_ID_list[i].strip("\n"))
        
        if (len(allele_name_list[1]) < 3):
            output_id_FH.write(allele_ID_list[i])
            exon_start = int(exonBoundary_dict[allele_ID_list[i].strip("\n")][1]) - 1
            exon_end   = int(exonBoundary_dict[allele_ID_list[i].strip("\n")][3]) - 1
            output_seq_FH.write(allele_sequence_list[i][exon_start:exon_end] + "\n")
        else:
            print "Locus: " + locus + " " + allele_ID_list[i] + "\n"
            print len(allele_name_list[1])
            print "\n"
            print "Locus: " + locus + " " + allele_name_list[1] + "\n"

    input_id_FH.close()
    input_seq_FH.close()
    output_id_FH.close()
    output_seq_FH.close()
    
    
    
    