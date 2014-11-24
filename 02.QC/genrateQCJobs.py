import sys,os,glob,subprocess

data_dir = sys.argv[1]

original_data_dir   = data_dir + "/original"
trimmed_data_dir    = data_dir + "/trimmed/emptyTrimmedReadsRemoved"
proj_dir            = data_dir + "/QC"
original_output_dir = proj_dir + "/original"
trimmed_output_dir  = proj_dir + "/trimmed"

if os.path.exists(proj_dir):
	print "Exists: " + proj_dir + "\n"
else:
	os.system("mkdir " + proj_dir)

if os.path.exists(original_output_dir):
	print "Exists: " + original_output_dir + "\n"
else:
	os.system("mkdir " + original_output_dir)

if os.path.exists(trimmed_output_dir):
	print "Exists: " + trimmed_output_dir + "\n"
else:
	os.system("mkdir " + trimmed_output_dir)

os.chdir(original_data_dir)
original_file_list = glob.glob("*.fastq")
os.chdir(trimmed_data_dir)
trimmed_file_list = glob.glob("*.fastq")

ORIGINAL_JOB_FILE = proj_dir + "/OriginalSeqJobs.txt"
TRIMMED_JOB_FILE  = proj_dir + "/TrimmedSeqJobs.txt"

ORIGINAL_FH = open(ORIGINAL_JOB_FILE, "w")
TRIMMED_FH  = open(TRIMMED_JOB_FILE, "w")

for file in original_file_list:
	inputfile = original_data_dir + "/" + file
	ORIGINAL_FH.write("/share/apps/FastQC/fastqc " + inputfile + " -o " + original_output_dir + "\n")

for file in trimmed_file_list:
	inputfile = trimmed_data_dir + "/" + file
	TRIMMED_FH.write("/share/apps/FastQC/fastqc " + inputfile + " -o " + trimmed_output_dir + "\n")

ORIGINAL_FH.close()
TRIMMED_FH.close()

ORIGINAL_QSUB_FH = open(proj_dir + "/OriginalSeqJobs.QSUB.sh", "w")
TRIMMED_QSUB_FH  = open(proj_dir + "/TrimmedSeqJobs.QSUB.sh", "w")

ORIGINAL_QSUB_FH.write("#!/bin/bash\n\n")
ORIGINAL_QSUB_FH.write("#PBS -N QC_Original" + "\n")
ORIGINAL_QSUB_FH.write("#PBS -o " + proj_dir + "/QC_Original.out" + "\n")
ORIGINAL_QSUB_FH.write("#PBS -e " + proj_dir + "/QC_Original.err" + "\n")
ORIGINAL_QSUB_FH.write("#PBS -l walltime=1:00:00\n")
ORIGINAL_QSUB_FH.write("#PBS -q default\n")
ORIGINAL_QSUB_FH.write("#PBS -l mem=9600MB\n")
ORIGINAL_QSUB_FH.write("#PBS -m abe\n")
ORIGINAL_QSUB_FH.write("#PBS -M liai.hpc.jobs@gmail.com\n")
ORIGINAL_QSUB_FH.write("#PBS -l nodes=1:ppn=8\n\n")
ORIGINAL_QSUB_FH.write("/share/apps/parallel/bin/parallel -j 0 -a " + ORIGINAL_JOB_FILE + "\n\n")
ORIGINAL_QSUB_FH.write("exit\n")

TRIMMED_QSUB_FH.write("#!/bin/bash\n\n")
TRIMMED_QSUB_FH.write("#PBS -N QC_TRIMMED" + "\n")
TRIMMED_QSUB_FH.write("#PBS -o " + proj_dir + "/QC_TRIMMED.out" + "\n")
TRIMMED_QSUB_FH.write("#PBS -e " + proj_dir + "/QC_TRIMMED.err" + "\n")
TRIMMED_QSUB_FH.write("#PBS -l walltime=1:00:00\n")
TRIMMED_QSUB_FH.write("#PBS -q default\n")
TRIMMED_QSUB_FH.write("#PBS -l mem=9600MB\n")
TRIMMED_QSUB_FH.write("#PBS -m abe\n")
TRIMMED_QSUB_FH.write("#PBS -M liai.hpc.jobs@gmail.com\n")
TRIMMED_QSUB_FH.write("#PBS -l nodes=1:ppn=8\n\n")
TRIMMED_QSUB_FH.write("/share/apps/parallel/bin/parallel -j 0 -a " + TRIMMED_JOB_FILE + "\n\n")
TRIMMED_QSUB_FH.write("exit\n")

ORIGINAL_QSUB_FH.close()
TRIMMED_QSUB_FH.close()











