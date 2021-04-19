#####################################################################
##  plink2call - PLINK hard genotypes to HLA allele calls
#####################################################################
## Modified by TDM Apr. 19 2021

""" plink2call.py

	Usage:
		plink2call.py (--vcf | --bfile) <file> --mapping <mapping> --plink <plink>

	Options:
		--bfile	: plink file set input
		--vcf	: vcf input
		--mapping	: Training feature set
		--plink	: location of plink executable (Tim local: /Users/tmajaria/Documents/src/plink/plink_mac_20190617/plink)

"""	
from docopt import docopt
import pandas as pd
import numpy as np
import subprocess
import sys
import os

#####################################################################
##  Main
#####################################################################
def main(docopt_args):
	mapping=docopt_args["<mapping>"]
	plink=docopt_args["<plink>"]
	file=docopt_args['<file>']

	#PLINK Location
	#plink="/gpfs/mrc0/projects/Research_Project-MRC158833/programs/plink/plink --silent"

	#Generate frequencies
	print("Generating frequencies with PLINK...")
	if docopt_args['--bfile']:
		command=plink+" --bfile "+file+" --freq --out "+file
		file_pref=file
	else:
		file_pref=file.replace(".vcf", "").replace(".gz", "")
		command=plink+" --vcf "+file+" --freq --out "+file_pref
		
	subprocess.run(command,shell=True)
	freq=pd.read_csv(file_pref+".frq",header=0, delim_whitespace=True)
	ma=freq[["SNP","A1"]]

	#Append MAF and tagged allele to SNP matrix
	print("Matching SNP allele to HLA allele...")
	tags=pd.read_csv(mapping, header=0, sep="\s+|\t+|\s+\t+|\t+\s+", engine='python')
	vmap=pd.merge(ma,tags)

	#Use plink score function to generate a table of allele dosages	
	print("Generating table of allele dosages...")
	tmp=file_pref+"_temp.tmp"
	for i in range(len(vmap)):
		scorefile=vmap[['SNP','A1']]
		scorefile['VAL']=1
		scorefile.iloc[[i]].to_csv(tmp,header=False, index=False, sep="\t")
		if docopt_args['--bfile']:
			command=plink+" --bfile "+file+" --score "+tmp+" no-mean-imputation"
		else:
			command=plink+" --vcf "+file+" --score "+tmp+" no-mean-imputation"
		subprocess.run(command,shell=True)
		temp=pd.read_csv("plink.profile",header=0, sep="\s+|\t+|\s+\t+|\t+\s+", 
									engine='python')
		temp['SCORE']=temp['SCORE']*2

		#Label the correct allele
		allele=str(vmap.loc[i,'ALLELE'])
		if i==0:
			table_hla=temp[['FID','IID','SCORE']]
			table_hla=table_hla.rename(columns = {'SCORE': allele})
		if i>0:
			table_hla[allele]=temp['SCORE']
	#Clean up and out
	command="rm "+file+".log "+file+".frq "+file+".nosex "+tmp+" 2> /dev/null"
	subprocess.run(command,shell=True)
	command="rm plink.profile plink.log plink.nosex 2> /dev/null"
	subprocess.run(command,shell=True)
	#table_hla.to_csv(bfile+"_dosage_anot.txt", header=True, index=False, sep="\t")
	#print("Success! Created "+bfile+"_dosage_anot.txt")

	#Sum rows and check
	print("Checking row sums for excess alleles...")
	table_hla['COUNT']=table_hla.iloc[:,3:].sum(axis=1).to_frame()
	table_count=table_hla[['FID','IID','COUNT']]
	table_count.to_csv(file_pref+"_count.txt", header=True, index=False, sep="\t")
	print("Created count file for failure checking "+file_pref+"_count.txt")
	table_clean=table_hla[(table_hla['COUNT'] < 3)]
	del table_clean['COUNT']
	table_clean.to_csv(file_pref+"_table.txt", header=True, index=False, sep="\t")
	print("Success! Created cleaned table "+file_pref+"_clean.txt")

	#Create categorised list
	print("Creating categorical list of genotypes per sample...")
	table_cat=table_clean[['FID','IID']].copy()
	table_cat['GENO1']="X"
	table_cat['GENO2']="X"
	for i in range(len(table_clean)):
		for j in range(2,len(table_clean.columns)):
			col=str(table_clean.columns[j])
			if table_clean.iloc[i,j]==2:
				col=table_clean.columns[j]
				table_cat.loc[table_cat.index[i],'GENO1']=col
				table_cat.loc[table_cat.index[i],'GENO2']=col
			elif table_clean.iloc[i,j]==1:
				if table_cat.loc[table_cat.index[i],'GENO1']=="X":
					table_cat.loc[table_cat.index[i],'GENO1']=col
				else:
					table_cat.loc[table_cat.index[i],'GENO2']=col
			
	table_cat.to_csv(file_pref+"_cat.txt", header=True, index=False, sep="\t")
	print("Success! Created "+file_pref+"_cat.txt")

	print("Finished.")

if __name__ == "__main__":
	args = docopt(__doc__)
	main(args)
