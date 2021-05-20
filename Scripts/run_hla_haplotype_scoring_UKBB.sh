#!/bin/bash

#$ -N ukbb_sharp_t1d_grs67_hla_may032021
#$ -l h_rt=24:00:00 ####### amount of time (hr:min:sec)
#$ -l h_vmem=16G #### amount of memory. default 2G (i think)
#$ -cwd ####### always have this: use current working directory
#$ -j y ####### always have this: join out and error outputs
source /broad/software/scripts/useuse
reuse Python-3.9

SAVEDIR=/humgen/florezlab/users/tmajaria/ukbb_t1d_grs/results

cd $SAVEDIR

RSIDFILE=/humgen/florezlab/users/tmajaria/ukbb_t1d_grs/PGS000024_unique_rsids.txt
QCTOOL=/humgen/florezlab/users/tmajaria/opt/qctool_v2.0.6-CentOS_Linux7.3.1611-x86_64/qctool
PLINK2=/humgen/florezlab/users/tmajaria/opt/plink2
PLINK1=/humgen/florezlab/users/tmajaria/opt/plink

# extract variants to new bgen file
${QCTOOL} \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr6_v3.bgen \
	-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
	-incl-rsids ${RSIDFILE} \
	-og ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.bgen \
	-os ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.sample

# convert to plink2 format
${PLINK2} \
	--bgen ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.bgen ref-first \
	--sample ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.sample \
	--make-pgen \
	--out ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021
​
# convert to plink1 format
${PLINK2} \
	--pfile ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021 \
	--make-bed \
	--rm-dup force-first \
	--out ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021

# make sure we have a stable copy of the fam file
cp ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.fam ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.STABLE.fam

# generate 10k sample subsets (running all UKBB samples at once takes tons of memory)
split -l 10000 ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.STABLE.fam ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.sample_subset_

# create list of sample subset files to loop through
printf '%s\n' ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.sample_subset_* > hlahaplo_sample_files.txt

# run haplotype scoring on each sample subset, remove intermediate files
while read samplefile; do
	awk '{ print $1, $2 }' ${samplefile} > ${samplefile}_samplelist.txt
	$PLINK1 --bfile ${SAVEDIR}/ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021 --keep ${samplefile}_samplelist.txt --make-bed --out ${SAVEDIR}/${samplefile}
	python3 /humgen/florezlab/users/tmajaria/opt/hla-prs-toolkit/Scripts/1_plink2call.py \
		--bfile ${SAVEDIR}/${samplefile} \
		--mapping /humgen/florezlab/users/tmajaria/opt/hla-prs-toolkit/Snplists/T1D_GRS67/mapping_TOPMED.txt \
		--plink $PLINK1 \
		--out_dir ${SAVEDIR}

	python3 /humgen/florezlab/users/tmajaria/opt/hla-prs-toolkit/Scripts/2_cat2scores.py \
		--cat ${SAVEDIR}/${samplefile}_cat.txt \
		--score /humgen/florezlab/users/tmajaria/opt/hla-prs-toolkit/Snplists/T1D_GRS67/scorefile.txt \
		--out_dir ${SAVEDIR}

	rm ${samplefile}_samplelist.txt ${samplefile}.bed ${samplefile}.fam ${samplefile}.bim ${samplefile}_dosage_anot.txt ${samplefile}_count.txt ${samplefile}_table.txt ${samplefile}_cat.txt ${samplefile}
	
done <hlahaplo_sample_files.txt
