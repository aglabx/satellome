#!/bin/bash

DIR='/mnt/data/satelome/users/mpopova/dnazoo'
QUAST='/home/mpopova/miniconda3/envs/quast_env/bin/quast'

for i in $(cat ${DIR}/files_name.txt); do
        python ${QUAST} -o ${DIR}/quast_results/${i} -t 80 -l ${i} -e --large ${DIR}/${i}.fasta.gz ;
done
