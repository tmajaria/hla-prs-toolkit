#!/usr/bin/bash

echo "Running script to generate allele calls from PLINK data..."
python 1_plink2call.py --vcf demo/demo_cohort.vcf.gz --mapping demo/mapping.txt --plink /Users/tmajaria/Documents/src/plink/plink_mac_20190617/plink
echo "Running script to generate scores from allele calls..."
python 2_cat2scores.py --cat demo/demo_cohort_cat.txt --score demo/scorefile.txt 
