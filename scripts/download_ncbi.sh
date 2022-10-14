#!/bin/bash

DIR='/mnt/data/satelome/users/mpopova'
OUTPUT_DIR='/mnt/data/satelome/users/mpopova/ncbi'

while IFS= read -r line; do
	wget -O ${OUTPUT_DIR}/${line};
done < ${DIR}/to_download.txt
