#!/bin/bash

DIR='/mnt/data/satelome/users/mpopova/dnazoo'
BUSCO='/mnt/data/satelome/users/mpopova/busco_database'
for i in $(cat ${DIR}/files_name.txt); do
	gzip -dkn ${DIR}/${i}.fasta.gz;
        busco -o ${i} -m genome -c 120 -l ${BUSCO}/carnivora_odb10 -i ${DIR}/${i}.fasta --offline;
	rm ${DIR}/${i}.fasta;
done

