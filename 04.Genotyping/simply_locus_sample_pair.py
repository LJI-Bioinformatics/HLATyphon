import sys,os,glob,re

# Note that the dictionary of loci should be modified by the experimental run Fastq file name.

dict_loci = {"A":    "A", 
             "B":    "B", 
             "C":    "C",
             "DPA":  "DPA1",
             "DPB":  "DPB1",
             "DQA":  "DQA1",
             "DQB":  "DQB1",
             "DRB1": "DRB1"
            }

data_dir    = sys.argv[1]
output_dir  = sys.argv[2]

os.chdir(data_dir)
file_list = glob.glob("*.fastq")

name_list    = []

for file in file_list:
    
    sample_name = ""
    
    file_name_list = file.split("_", 1)
    
    if ("A" in file_name_list[0] and "DP" not in file_name_list[0] and "DQ" not in file_name_list[0]):
        sample_name = file_name_list[0].replace("A",",A")
    if ("B" in file_name_list[0] and "DP" not in file_name_list[0] and "DQ" not in file_name_list[0]):
        sample_name = file_name_list[0].replace("B",",B")
    if ("C" in file_name_list[0]):
        sample_name = file_name_list[0].replace("C",",C")
    if ("DPA" in file_name_list[0]):
        sample_name = file_name_list[0].replace("DPA",",DPA1")
    if ("DPB" in file_name_list[0]):
        sample_name = file_name_list[0].replace("DPB",",DPB1")
    if ("DQA" in file_name_list[0]):
        sample_name = file_name_list[0].replace("DQA",",DQA1")
    if ("DQB" in file_name_list[0]):
        sample_name = file_name_list[0].replace("DQB",",DQB1")
    if ("DRB1" in file_name_list[0]):
        sample_name = file_name_list[0].replace("DRB1",",DRB1")
    
    name_list.append(sample_name)
        
name_set = set(name_list)

with open(output_dir + "/Sample_Locus.csv", "w") as FILE:
    FILE.write("\n".join(name_set))
    
    
        
    
    
    
    
        
    
    
    
    
    
