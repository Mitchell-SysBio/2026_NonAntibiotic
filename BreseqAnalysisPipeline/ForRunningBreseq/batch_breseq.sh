#!/bin/bash

for name in */;
do
cd ${name}
pwd

breseq -r ../MG1655_NC_00913.gbk *R1_001.fastq.gz *R2_001.fastq.gz
cd ../

done
