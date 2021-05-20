#!/bin/bash
#$ -N ukbb_sharp_t1d_grs67_may032021
#$ -pe smp 4
#$ -l h_vmem=16G
#$ -cwd
#$ -S /bin/bash
#$ -o $JOB_NAME_$TASK_ID.out
#$ -e $JOB_NAME_$TASK_ID.err

INPUTFILES=($(ls -d /broad/ukbb/imputed_v3/ukb_imp_chr*_v3.bgen))
INPUTFILENAME="${INPUTFILES[$SGE_TASK_ID - 1]}"
OUTPREF1=${INPUTFILENAME##*/}
OUTPREF=${OUTPREF1%.*}

SAVEDIR=/humgen/florezlab/users/tmajaria/ukbb_gels/results
cd $SAVEDIR

source /broad/software/scripts/useuse
reuse GCC-5.2
QCTOOL=/humgen/florezlab/users/tmajaria/opt/qctool_v2.0.6-CentOS_Linux7.3.1611-x86_64/qctool
PLINK2=/humgen/florezlab/users/tmajaria/opt/plink2

# extract variants to new bgen file
${QCTOOL} \
	-g ${INPUTFILENAME} \
	-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
	-incl-rsids /humgen/florezlab/users/tmajaria/ukbb_t1d_grs/PGS000024_unique_rsids.txt \
	-og ${SAVEDIR}/${OUTPREF}.sharp_t1d_grs67_nohaplo.may032021.bgen \
	-os ${SAVEDIR}/${OUTPREF}.sharp_t1d_grs67_nohaplo.may032021.sample

# convert to plink2 format
${PLINK2} \
	--bgen ${SAVEDIR}/${OUTPREF}.sharp_t1d_grs67_nohaplo.may032021.bgen ref-first \
	--sample ${SAVEDIR}/${OUTPREF}.sharp_t1d_grs67_nohaplo.may032021.sample \
	--make-pgen \
	--out ${SAVEDIR}/${OUTPREF}.sharp_t1d_grs67_nohaplo.may032021
​
# score
${PLINK2} \
	--pfile ${SAVEDIR}/${OUTPREF}.sharp_t1d_grs67_nohaplo.may032021 \
	--score /humgen/florezlab/users/tmajaria/ukbb_t1d_grs/PGS000024_ukb_imputed_v3.score_file.txt 1 2 3 header-read ignore-dup-ids cols=+scoresums \
	--out ${SAVEDIR}/${OUTPREF}.sharp_t1d_grs67_nohaplo.may032021

