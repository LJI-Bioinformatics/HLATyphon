import sys,os,glob,re,shutil
import math

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
output_file = sys.argv[2]
os.chdir(data_dir)
file_list = glob.glob("*.fastq")

output_FILE_HANDLE = open(output_file, "w")

for file in file_list:
    for loci in dict_loci:
        m = re.search('\S+(?='+loci+'_)' , file)
        if(loci == "A" and m): 
            if(("DP" not in file) and ("DQ" not in file)):
                output_FILE_HANDLE.write(m.group(0) + "," + dict_loci[loci] + "\n")
        if(loci == "B" and m):
            if(("DP" not in file) and ("DQ" not in file) and ("DR" not in file)):
		output_FILE_HANDLE.write(m.group(0) + "," + dict_loci[loci] + "\n")
        if(loci != "A" and loci != "B" and m):
		output_FILE_HANDLE.write(m.group(0) + "," + dict_loci[loci] + "\n")
           
output_FILE_HANDLE.close()





