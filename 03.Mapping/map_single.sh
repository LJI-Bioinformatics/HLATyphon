line=`awk ' END { print NR } ' ${id_file}`
mkdir -p ${PROJ}/map_single
mkdir -p ${SCPROJ}/map_single
rm ${PROJ}/map_single/*
rm ${SCPROJ}/map_single/*
cp ${fq_first} ${SCPROJ}/fq_first
cp ${fq_second} ${SCPROJ}/fq_second

if [ $debug = 0 ]
then
    export BOWTIE2_INDEXES=${bw2_index_dir}/${bw2_index_pre}
    for((j=1; j<=$line; j++));
    do
        id=`sed -n "${j}p" < ${id_file}`
        seq=`sed -n "${j}p" < ${seq_file}`
        # algin reads
        ${bw2_align} --local -N 0 -L 25 -i S,1,0.23 --mp 1000,1000 --np 1000 --rdg 1000,1000 --rfg 1000,1000 --score-min L,100,0 --fr --no-discordant -x ${bw2}/${id} --ignore-quals -q -1 ${SCPROJ}/fq_first -2 ${SCPROJ}/fq_second -S ${SCPROJ}/map_single/${id}.bowtie2.sam
        echo ${id} " aligned"
        if [ $paired = 1 ]
        then 
            ${samtools} view -h -F 12 -S ${SCPROJ}/map_single/${id}.bowtie2.sam > ${SCPROJ}/map_single/${id}.report.bowtie2.sam
        fi
        if [ $paired = 0 ]
        then
            ${samtools} view -h -F 4 -S ${SCPROJ}/map_single/${id}.bowtie2.sam > ${SCPROJ}/map_single/${id}.report.bowtie2.sam
        fi
	rm ${SCPROJ}/map_single/${id}.bowtie2.sam
    done
fi
