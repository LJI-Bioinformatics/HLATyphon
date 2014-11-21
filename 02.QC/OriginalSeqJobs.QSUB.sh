#!/bin/bash

#PBS -N QC_Original
#PBS -o /Bioinformatics/Users/zfu/HLA_Typing/QC/Reanalyze.32Samples/QC_Original.out
#PBS -e /Bioinformatics/Users/zfu/HLA_Typing/QC/Reanalyze.32Samples/QC_Original.err
#PBS -l walltime=4:00:00
#PBS -q default
#PBS -l mem=9600MB
#PBS -m abe
#PBS -M liai.hpc.jobs@gmail.com
#PBS -l nodes=1:ppn=8

/share/apps/parallel/bin/parallel -j 0 -a /Bioinformatics/Users/zfu/HLA_Typing/QC/Reanalyze.32Samples/OriginalSeqJobs.txt

exit
