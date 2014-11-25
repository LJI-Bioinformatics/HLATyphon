import sys
import operator
import re
import os
import glob

exonBoundary_dir = "/Bioinformatics/Users/zfu/HLA_Typing/HLA_cDNA_Database/IMGT.Release.3.12.0/ExonBoundary"
output_dir       = "/Bioinformatics/Users/zfu/HLA_Typing/HLA_cDNA_Database/IMGT.Release.3.12.0/ExonBoundary_Exon2"
os.chdir(exonBoundary_dir)
file_list = glob.glob("*.txt")

for file in file_list:
    
    output_file_FH = open(output_dir + "/" + file, "w")
    
    with open (exonBoundary_dir + "/" + file, "r") as INPUT_FILE:
        for line in INPUT_FILE:
            (ID, Boundary) = re.split("\s+", line.strip("\n"), 1)
            Boundary_List = re.split(",", Boundary.strip(","))
            Exon2_Boundary_List = [int(elem) - int(Boundary_List[1]) + 1 for elem in Boundary_List]
            output_file_FH.write(ID + "\t" + ",".join([str(boundary) for boundary in Exon2_Boundary_List]) + "\n")
    
    output_file_FH.close()