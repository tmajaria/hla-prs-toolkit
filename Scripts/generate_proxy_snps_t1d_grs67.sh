#!/bin/bash
# EUR LD proxies for T1D GRS SNPs

cd /Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d

# populations to get proxies from
POP='CEU%2BTSI%2BFIN%2BGBR%2BIBS'

# token to use in ld link api
TOKEN=get_your_own_token_online
RSIDFILE=/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/PGS000024_rsids.txt

# extract all RSIDs
python - << EOF
import pandas as pd
prs_file = '/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/PGS000024.txt.gz'
out_file = '/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/PGS000024_rsids.txt'
prs = pd.read_csv(prs_file, skiprows=9, sep = "\t")[['rsID']].apply(lambda x: x.str.split('_x_').explode())
prs.to_csv(out_file, sep = '\n', header = False, index = False)
EOF

# call ldlink api to get proxies
while read rsid; do
  curl -k -X GET 'https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var='${rsid}'&pop='${POP}'&r2_d=r2&window=500000&token='${TOKEN} > /Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/eu_proxy/${rsid}_ldproxy_EUR.txt
done <$RSIDFILE

# process variants, output snp->proxy map and full list of snps & proxies to be used for querying bravo
python - << EOF
import pandas as pd
import os

base_dir = '/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/eu_proxy'
proxy_dfs = dict()
for f in [base_dir + '/' + f for f in os.listdir(base_dir) if f.endswith('ldproxy_EUR.txt')]:
	rs = f.replace(base_dir+'/','').replace('_ldproxy_EUR.txt', '')
	# proxy_dfs[rs] = pd.read_csv(f, sep = "\t").iloc[0:10,]
	proxy_dfs[rs] = pd.read_csv(f, sep = "\t")
	proxy_dfs[rs][['index_var']] = rs

proxy_var = pd.concat(proxy_dfs)[['RS_Number','index_var','Dprime','R2','Correlated_Alleles']]
proxy_var = proxy_var.loc[proxy_var.RS_Number != ".",:]
proxy_var = proxy_var.loc[proxy_var.R2 >= 0.6,:]
proxy_var.columns = [ 'EUR_proxy', 'Sharp_t1d_snp', 'Dprime', 'R2','Correlated_Alleles']
proxy_var.iloc[:,[1,0,2,3,4]].to_csv('/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/sharp_t1d_EUR_proxy.all_variants.map.txt', sep = '\t', index=None)

all_var = pd.concat([proxy_var.Sharp_t1d_snp, proxy_var.EUR_proxy]).unique()
pd.DataFrame(all_var).to_csv('/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/sharp_t1d_EUR_proxy.all_variants.rsids_only.txt', sep = '\t', index=None, header=None)
EOF

# query bravo (bravo has to be setup first, see: https://bravo.sph.umich.edu/freeze5/hg38/help)
BRAVOFILE=/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/sharp_t1d_EUR_proxy.all_variants.rsids_only.txt
while read rsid; do
	curl -H "Authorization: Bearer `/Users/tmajaria/Documents/src/bravo print-access-token`" "https://bravo.sph.umich.edu/freeze5/hg38/api/v1/variant?variant_id="$rsid > eu_proxy/${rsid}_topmed_bravo_response.json
done <$BRAVOFILE

# add the topmed ids to the snps&proxies
python - << EOF
import pandas as pd
import os
import json

def correlated_allele_weight(row):
	if row.rsID == row.EUR_proxy:
		return(row.effect_allele)
	else:
		index_effect_allele = row.effect_allele
		corr_alleles = row.Correlated_Alleles.split(',')
		proxy_effect_allele = ""
		if "-" in row.Correlated_Alleles:
			base_to_add = index_effect_allele.replace([c for c in corr_alleles if "-" not in c.split('=')[0]][0].split('=')[0], '')
			corr_alleles = [base_to_add+c.replace('-','') for c in corr_alleles]
		for alleles in corr_alleles:
			als = alleles.split('=')
			if (als[0] == index_effect_allele) or (als[0] == '-') : 
				proxy_effect_allele = als[1]
		if proxy_effect_allele == "":
			return(None)
		else:
			return(proxy_effect_allele)

proxy_map = pd.read_csv('/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/sharp_t1d_EUR_proxy.all_variants.map.txt', sep = '\t')
base_dir = '/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/eu_proxy'
prs_file = '/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/PGS000024.txt.gz'
bravo_dfs = dict()
for fn in [base_dir + '/' + f for f in os.listdir(base_dir) if f.endswith('topmed_bravo_response.json')]:
	rs = fn.replace(base_dir+'/','').replace('_topmed_bravo_response.json', '')
	with open(fn, 'r') as f:
		data = json.load(f)

	if 'data' in data.keys(): 
		bravo_dfs[rs] = pd.DataFrame(data['data'])
		bravo_dfs[rs][['og_rs']] = rs

bravo_var = pd.concat(bravo_dfs)
bravo_var['SNP'] = bravo_var.rsids.apply(lambda x: x[0])
bravo_var['MARKER'] = bravo_var.variant_id.apply(lambda x: 'chr'+x.replace('-',':'))
bravo_var = bravo_var.loc[bravo_var.SNP.isin(pd.concat([proxy_map.EUR_proxy, proxy_map.Sharp_t1d_snp])),['SNP','MARKER']]
pd.DataFrame(bravo_var.MARKER.unique()).to_csv('/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/sharp_t1d_EUR_proxy.all_variants.var_ids_only.txt', sep = '\t', index=None, header=None)

bravo_proxy = pd.merge(proxy_map, bravo_var, left_on='EUR_proxy', right_on='SNP', how='left').loc[:,['Sharp_t1d_snp','EUR_proxy','MARKER','Dprime','R2','Correlated_Alleles']]
bravo_proxy = bravo_proxy.loc[pd.notnull(bravo_proxy.MARKER),:]
bravo_proxy.to_csv('/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/sharp_t1d_EUR_proxy.all_variants.rsid_var_id.map.txt', sep = '\t', index=None)

prs = pd.read_csv(prs_file, skiprows=9, sep = "\t")
# prs = prs.loc[~prs.is_haplotype,['rsID','effect_allele','effect_weight','OR','is_haplotype']]
prs = prs.loc[:,['rsID','effect_allele','effect_weight','OR','is_haplotype']]
score_df = pd.merge(prs, bravo_proxy, left_on = 'rsID', right_on = 'Sharp_t1d_snp', how = 'right')

score_df[['proxy_effect_allele']] = score_df.apply(correlated_allele_weight, axis=1)
score_df = score_df.loc[score_df.proxy_effect_allele != None,['Sharp_t1d_snp','EUR_proxy','MARKER','proxy_effect_allele','effect_weight','Dprime','R2','OR','is_haplotype']]
score_df.columns = ['Sharp_t1d_snp','EUR_proxy_snp','TOPMED_marker','effect_allele','effect_weight','Dprime','R2','OR','is_haplotype']
score_df.to_csv('/Users/tmajaria/Documents/projects/biobanks/mgb/data/prs/sharp_t1d/sharp_t1d_EUR_proxy.all_variants.scores_raw.txt', sep = '\t', index=None)

EOF
