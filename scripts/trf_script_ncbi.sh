#!/bin/bash

DIR='/mnt/data/satelome/users/mpopova/ncbi'
OUTPUT_DIR='/mnt/data/satelome/users/mpopova/ncbi_trf_results'
TRF='/home/mpopova/soft/trf_scripts/trf_search.py'

for i in $(cat ${DIR}/files_name.txt); do
        gzip -dkn ${DIR}/${i}.fna.gz;
        python ${TRF} -i ${DIR}/${i}.fna -o ${OUTPUT_DIR} -p tandem -r ${i}_results.yaml -t 40;
        rm ${DIR}/${i}.fna;
done
