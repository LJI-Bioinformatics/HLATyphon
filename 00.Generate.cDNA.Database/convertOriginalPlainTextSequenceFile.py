import sys
import operator
import re
import os
import glob

original_sequence_file_dir = sys.argv[1]

os.chdir(original_sequence_file_dir)
file_list = glob.glob("*.txt")

for original_sequence_file in file_list:
    
    sequenceDirName = original_sequence_file_dir
    sequenceFileName = original_sequence_file
    #(sequenceDirName, sequenceFileName) = os.path.split(original_sequence_file)
    convertedSequenceFileName = "Converted_" + sequenceFileName;
    convertedSequencePathName = os.path.join(sequenceDirName, convertedSequenceFileName)
    
    Allele_Dictionary = {}
    
    convertedSequenceFile_FH = open(convertedSequencePathName, "w")
    
    with open(original_sequence_file, "r") as input:
    # read file containing locus allele IDs.
        for line in input:
            modified_line = re.sub(r'^\s+', r'', line.strip("\n"))
            match = re.match(r'HLA', modified_line)
            if (match is not None):
                (allele, sequence) = re.split("\s+", modified_line, 1)
                modified_sequence = re.sub(r'\s+', r'', sequence)
                if (not allele in Allele_Dictionary):
                    Allele_Dictionary[allele] = []
                Allele_Dictionary[allele].append(modified_sequence)
    
    for each_allele in Allele_Dictionary:
        
        convertedSequenceFile_FH.write(each_allele + "\t" + "".join(Allele_Dictionary[each_allele]) + "\n")
        
    convertedSequenceFile_FH.close()
 
            
            