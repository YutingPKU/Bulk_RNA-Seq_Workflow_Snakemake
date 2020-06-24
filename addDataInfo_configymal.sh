#!/bin/bash

while read -r line
do
	echo "  $line:" >> config.yaml
	echo "    - data/$line/${line}_R1.fq.gz" >> config.yaml
	echo "    - data/$line/${line}_R2.fq.gz" >> config.yaml
done < sample.ls
