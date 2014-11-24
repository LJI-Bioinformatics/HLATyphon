#!/bin/bash
mkdir -p ${PROJ}/detail
mkdir -p ${PROJ}/SAM_FILES
mkdir -p ${SCPROJ}/SAM_FILES
${python} ${PROJ}/countMappedReads.py ${id_file} ${seq_file} ${exon_bd} ${PROJ} ${SCPROJ} ${primer_start} ${primer_end} ${EQUIPMENT}
sort -k2,2n ${PROJ}/count_2.txt | tail -10 > ${PROJ}/HLA_type.txt


