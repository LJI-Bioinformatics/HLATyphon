#!/bin/bash

#PBS -N MERGE.FASTQ
#PBS -o /Bioinformatics/Users/zfu/HLA_Typing/src/HLA_Typing_Parsing_Codes/MERGE.FASTQ.out
#PBS -e /Bioinformatics/Users/zfu/HLA_Typing/src/HLA_Typing_Parsing_Codes/MERGE.FASTQ.err
#PBS -l walltime=1:00:00
#PBS -q default
#PBS -l mem=4800MB
#PBS -m abe
#PBS -M liai.hpc.jobs@gmail.com
#PBS -l nodes=1:ppn=4

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo ' '
echo ' '


sg fs_bioinformatics -c "/share/apps/parallel/bin/parallel -j 0 -a /Bioinformatics/Users/zfu/HLA_Typing/src/HLA_Typing_Parsing_Codes/parallelMerge_CMD.txt"

exit
