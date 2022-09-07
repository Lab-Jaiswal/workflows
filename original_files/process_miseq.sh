#!/bin/bash

seq_path=$1

cd $seq_path
module load bcl2fastq2
bcl2fastq -o ./fastq -p 8
