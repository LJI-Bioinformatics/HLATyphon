# map to individual reference sequence
if [ $map_single = 1 ]
then
    source ${PROJ}/map_single.sh
fi

# filter out unreasonable individual reference sequence
if [ $select_single = 1 ]
then
    source ${PROJ}/select_single.sh
fi

# remove the scratch directory
rm -rf ${SCPROJ}
rm -rf ${PROJ}/map_comb
rm -rf ${PROJ}/map_single
rm -rf ${PROJ}/bw2_index_remain


