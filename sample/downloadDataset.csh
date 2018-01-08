#!/bin/sh
# Sample dataset from the Redundans project ( https://github.com/Gabaldonlab/redundans )
# This provides a useful test of the workflow without long runtimes.
wget https://github.com/Gabaldonlab/redundans/raw/master/test/ref.fa -O 600-ref.fa
wget https://github.com/Gabaldonlab/redundans/raw/master/test/600_1.fq.gz -O - | gunzip -c > 600_1.fastq
wget https://github.com/Gabaldonlab/redundans/raw/master/test/600_2.fq.gz -O - | gunzip -c > 600_2.fastq


