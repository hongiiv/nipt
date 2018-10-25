#!/bin/bash


#bedtools makewindows -g hg38.chrom.sizes -w 100000 > bins
#bedtools coverage -abam HG00650.mapped.ILLUMINA.bwa.CHS.low_coverage.20101123.bam -b bins > cov.bins
python nipt.py

