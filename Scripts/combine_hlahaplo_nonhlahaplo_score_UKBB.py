import pandas as pd
import os
import numpy as np
from functools import reduce

# paths & outfiles
res_dir = '/humgen/florezlab/users/tmajaria/ukbb_t1d_grs/results'
score_out = res_dir+'/ukb_imputed_v3.sharp_t1d_grs67.full_score.may132021.csv'

# gather hla & non-hla haplo score files
hla_files = [f for f in os.listdir(res_dir) if f.startswith('ukb_imputed_v3.chr6.sharp_t1d_grs67_hlahaplo_variants.may032021.sample_subset_') and f.endswith('_cat_Scored.txt')]
nonhla_files = [f for f in os.listdir(res_dir) if f.startswith('ukb_imp_chr') and f.endswith('_v3.sharp_t1d_grs67_nohaplo.may032021.sscore')]

# load hla haplo-only scores
hla_df = pd.concat([pd.read_csv(res_dir+'/'+f, sep="\t") for f in hla_files])
hla_df.columns = ['FID', 'IID', 'sharp_t1d_grs67_hlahaplo_score']

# load non-hla haplo scores
nonhla_dfs = list()
for f in nonhla_files:
	chr = f.split('_')[2]
	df = pd.read_csv(res_dir+'/'+f, sep='\t')
	df.columns = ['FID', 'IID', 'sharp_t1d_grs67_nonhla_ALLELE_CT_'+chr, 'sharp_t1d_grs67_nonhla_NAMED_ALLELE_DOSAGE_SUM_'+chr, 'sharp_t1d_grs67_nonhla_score_AVG_'+chr, 'sharp_t1d_grs67_nonhla_score_SUM_'+chr]
	nonhla_dfs.append(df)

# merge chromosome-specific non-hla haplo scores, add all chromosome-specific scores to get one non-hla haplo score
nonhla_df = reduce(lambda left, right: pd.merge(left, right, on=['FID', 'IID'], how='outer'), nonhla_dfs)
nonhla_df['sharp_t1d_grs67_nonhla_score'] = nonhla_df.loc[:,[col for col in nonhla_df.columns if col.startswith('sharp_t1d_grs67_nonhla_score_SUM_')]].sum(axis=1)

# merge in hla haplo score
score_df = pd.merge(hla_df, nonhla_df, on=['FID', 'IID'], how='outer')

# combine for full prs
score_df['sharp_t1d_grs67_score'] = score_df.sharp_t1d_grs67_hlahaplo_score + score_df.sharp_t1d_grs67_nonhla_score
score_df.loc[:,['FID', 'IID', 'sharp_t1d_grs67_hlahaplo_score', 'sharp_t1d_grs67_nonhla_score', 'sharp_t1d_grs67_score']].to_csv(score_out, index = False)

# write it out
score_df.to_csv(score_out, index = False)