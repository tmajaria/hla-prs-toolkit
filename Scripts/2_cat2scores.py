#####################################################################
##  cat2scores - Table of DIPLOTYPES to SCORE values
#####################################################################

""" cat2scores.py

	Usage:
		cat2scores.py --cat <catfile> --score <scorefile> [--out_dir <out_dir>]

	Options:
		--cat	: Show this help message
		--score	: Training feature set
		--out_dir	: optional location to write files

"""	
from docopt import docopt
import re
import pandas as pd
import numpy as np
import subprocess
import sys
import os

#####################################################################
##  Main
#####################################################################
def main(docopt_args):
	catfile=docopt_args["<catfile>"]
	scorefile=docopt_args["<scorefile>"]
	if docopt_args['--out_dir']:
		out_dir=docopt_args['<out_dir>']
	else:
		out_dir = os.path.dirname(catfile)
	prefix1=str(re.sub('(.txt)','',catfile)).split('/')[-1]
	prefix2=str(re.sub('(.txt)','',scorefile)).split('/')[-1]

	#Load the files
	scores=pd.read_csv(scorefile, header=0, sep="\s+|\t+|\s+\t+|\t+\s+", engine='python')
	table_cat=pd.read_csv(catfile, header=0, sep="\s+|\t+|\s+\t+|\t+\s+", engine='python')
	score_out=table_cat[['FID','IID']].copy()
	score_out['SCORE']=0

	print("Matching combinations of alleles to scores...")
	#Match using only the diplotypes
	tmp=scores.dropna()
	for i in range(len(table_cat)):
		for j in range(len(tmp)):
			if ((table_cat.loc[table_cat.index[i],'GENO1'] == tmp.loc[tmp.index[j],'ALLELE1']) and
				(table_cat.loc[table_cat.index[i],'GENO2'] == tmp.loc[tmp.index[j],'ALLELE2']) or
				(table_cat.loc[table_cat.index[i],'GENO2'] == tmp.loc[tmp.index[j],'ALLELE1'] and 
				table_cat.loc[table_cat.index[i],'GENO1'] == tmp.loc[tmp.index[j],'ALLELE2'])):
					score_out.loc[score_out.index[i],'SCORE'] = tmp.loc[tmp.index[j],'BETA']

	print("Check for and sum NULL labelled additive alleles...")
	#Match remaining using additive haplotypes if exists
	tmp=scores[scores.isna().any(axis=1)]
	if tmp.empty==False:
		for i in range(len(table_cat)):
			if score_out.loc[score_out.index[i],'SCORE']==0:
				for col in ['GENO1','GENO2']:
					for j in range(len(tmp)):
						if (table_cat.loc[table_cat.index[i],col] ==  tmp.loc[tmp.index[j],'ALLELE2']):
							score_out.loc[score_out.index[i],'SCORE'] += tmp.loc[tmp.index[j],'BETA']
		print("Success! Matched the NULL labelled additive alleles")
	else:
		print("None detected - skipping...")

	score_out.to_csv(os.path.join(out_dir, prefix1+"_Scored.txt"), header=True, index=False, sep="\t")
	print("Success! Created "+prefix1+"_Scored.txt")
	print("Finished.")

if __name__ == "__main__":
	args = docopt(__doc__)
	main(args)
