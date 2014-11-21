import sys
import operator
import re
import os
import glob

data_dir = sys.argv[1]
NUMBER_OF_EXONS = sys.argv[2]
index_dir = "/Bioinformatics/Users/zfu/HLA_Typing/HLA_cDNA_Database/IMGT.Release.3.12.0/Index"
exonBoundary_dir = "/Bioinformatics/Users/zfu/HLA_Typing/HLA_cDNA_Database/IMGT.Release.3.12.0/ExonBoundary"
os.chdir(data_dir)
file_list = glob.glob("*.txt")

for file in file_list:
    
    exonBoundary_dict = {}
    cDNASequence_dict = {}
    
    fileNoSuffix = re.sub(r'.txt', r'', file)
    locus = re.sub(r'Converted_', r'', fileNoSuffix)
    locus_index_dir = index_dir + "/" + locus + "_nuc"
    
    if not os.path.exists(locus_index_dir):
        os.makedirs(locus_index_dir)
    
    exonBoundaryFileName = locus + "_nuc_EX.txt"
    alleleIDFileName = locus + "_nuc_id.txt"
    allelecDNAFileName = locus + "_nuc_seq.txt"
    
    exonBoundaryPathName = os.path.join(exonBoundary_dir, exonBoundaryFileName)
    alleleIDPathName     = os.path.join(locus_index_dir, alleleIDFileName)
    allelecDNAPathName   = os.path.join(locus_index_dir, allelecDNAFileName)
    
    exonBoundaryPathName_FH = open(exonBoundaryPathName, "w")
    alleleIDPathName_FH     = open(alleleIDPathName, "w")
    allelecDNAPathName_FH   = open(allelecDNAPathName, "w")
    
    with open(data_dir + "/" + file, "r") as input:
        for line in input:
            
            (allele_ID, sequence) = re.split(r'\s+', line.strip("\n"), 1)
            
            print("The original ID: " + allele_ID + "\n")
            
            allele_ID = re.sub(r'HLA-', r'', allele_ID)
            
            print("After remove HLA-: " + allele_ID + "\n")
            
            allele_ID = re.sub(r'\*', r'_', allele_ID)
            
            print("After remove *: " + allele_ID + "\n")
            
            exons_list = []
            exons_list = re.split(r'\|', sequence)
            
            print("Number of exons: " + str(len(exons_list)) + "\n")
            
            if ( len(exons_list) != int(NUMBER_OF_EXONS) ):
                raise ValueError("The number of exons is not equal to: " + str(NUMBER_OF_EXONS))
            
            if (not allele_ID in exonBoundary_dict):
                exonBoundary_dict[allele_ID] = []
            
            if (not allele_ID in cDNASequence_dict):
                cDNASequence_dict[allele_ID] = []
            
            for exon in exons_list:
                
                exon = re.sub(r'\*', r'', exon)
                exon = re.sub(r'\.', r'', exon)
                
                exon_length = len(exon)
                
                exonBoundary_dict[allele_ID].append(exon_length)
                cDNASequence_dict[allele_ID].append(exon)
        
    
    print("The locus file is: " + data_dir + "/" + file + "\n")
    print("Number of Alleles: " + str(len(exonBoundary_dict)) + "\n")
    
    for alleles in exonBoundary_dict:
        
        alleleIDPathName_FH.write(alleles + "\n")
        
        cDNA_Sequence = "".join(cDNASequence_dict[alleles])
        allelecDNAPathName_FH.write(cDNA_Sequence + "\n")
        
        exonBoundaryStartPosition = 1
        exonBoundaryPathName_FH.write(alleles + "\t" + str(exonBoundaryStartPosition) + ",")
        
        for length in exonBoundary_dict[alleles]:
            exonBoundaryStartPosition = exonBoundaryStartPosition + length
            exonBoundaryPathName_FH.write(str(exonBoundaryStartPosition) + ",")
        
        exonBoundaryPathName_FH.write("\n")
    
    alleleIDPathName_FH.close()
    exonBoundaryPathName_FH.close()
    allelecDNAPathName_FH.close()
    
    print("DONE: The locus file is: " + data_dir + "/" + file + "\n\n")

print("DONE: " + data_dir + "\n\n")
                
            
            
            
            
    
    
    