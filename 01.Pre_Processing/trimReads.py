import sys,os,glob,subprocess

data_dir = "/Bioinformatics/NGS_raw_data/LIAI/MiSeq/HLA_Typing_Runs/Run52/original"

os.chdir(data_dir)
file_list = glob.glob("*.fastq")

for file in file_list:
    
    ## run external perl script
    if("_R1_" in file):
        file1 = file
        file2 = file.replace("_R1_", "_R2_")
        scmd = "/Bioinformatics/Users/zfu/HLA_Typing/src/trim.pl --type=1 --offset=0 --qual-threshold=28  --qual-type=0 --length-threshold=50" + \
               " --pair1=" + file1 + " --pair2=" + file2 + \
               " --outpair1=/Bioinformatics/NGS_raw_data/LIAI/MiSeq/HLA_Typing_Runs/Run52/trimmed/" + file1 + " --outpair2=/Bioinformatics/NGS_raw_data/LIAI/MiSeq/HLA_Typing_Runs/Run52/trimmed/" + file2 + \
               " --single=/Bioinformatics/NGS_raw_data/LIAI/MiSeq/HLA_Typing_Runs/Run52/trimmed/single.fastq"
        subprocess.call(scmd, shell=True)
