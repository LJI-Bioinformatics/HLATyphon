HLATyphon
=========

High-Fidelity HLA typing pipeline from next-generation sequencing data

Authors: Dr. Zheng Fu and Dr. Jason Greenbaum  
Date: Nov. 2014  
Version: 1.0  

# Overview

HLATyphon is a novel NGS data analysis pipeline allowing user to make HLA alleles genotyping quickly and efficiently. It supports all cDNA template of HLA locus A, B, C, DPA1, DPB1, DQA1, DQB1 and DRB1 in IMGT/HLA database. The input of HLATyphon could be either genomic DNA or RNA paired-end short reads. Currently the internal HLA allele reference sequences of HLATyphon were from IMGT/HLA databese release 3.17.0 including 12,010 alleles. 

# Requirments

Since most codes of HLATyphon were writen by python and perl and it was parallelized for grid-computing. Thus the following packages should be installed on your computer before you performing HLATyphon.

- Python 2.7
- Perl 5.16.3
- Bowtie2 Version 2.2.2 (not Bowtie1 since it doesn't support 'soft clip')
- Portable Batch System
- GNU parallel

# Contact

For more information about HLATyphon please contact zfu at liai dot org.









