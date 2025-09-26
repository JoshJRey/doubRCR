#!/usr/bin/env python

### Import packages
import numpy as np
import allel
import pandas as pd
from argparse import ArgumentParser

### Set up our command line argument parser
parser = ArgumentParser()
parser.add_argument("sim")
parser.add_argument("type")

### Set variables from the command line
args = parser.parse_args()
sim = args.sim
type = args.type

### Load in vcf
vcf_path = "/rds/general/user/jjr18/results/" + str(type) + "/" + str(sim) + ".vcf"
callset= allel.read_vcf(vcf_path)

### Access genotype data
gt = allel.GenotypeChunkedArray(callset['calldata/GT'])

alt = callset['variants/ALT'][:].astype('U')
ref = callset['variants/REF'][:].astype('U')

### Generate allele count array
ac_f = gt.count_alleles()

### Calculate the number of genomes in the VCF
ngenomes = 2*np.count_nonzero(callset['samples'])

### Identify every double-singleton site (one count for two different minor alleles), and what bases are present at each site
double_singles = ((ac_f[:,0]==(ngenomes-2)) & (ac_f[:,1]==1) & (ac_f[:,2]==1) & (ac_f[:,3]==0))

alt_singles = alt[double_singles]
ref_singles = ref[double_singles]


### Count double-singletons for each possible base combination
ds_counts = []
for a in np.unique(ref_singles):
    base = (alt_singles[:,:2])[ref_singles==a]
    for b in ['A', 'C', 'T', 'G']:
        for c in ['A', 'C', 'T', 'G']:
            both = (base==b) + (base==c)
            summ = np.sum(both, axis = 1)
            if np.count_nonzero((summ==2))!=0:
                ds_counts.append(np.count_nonzero((summ==2)))
                

ds_counts = np.array(ds_counts)
# Remove duplicates
ds_counts = ds_counts[[0,1,3,6,7,9,12,13,15,18,19,21]]


### Convert counts into proportion of reference base sites, for each reference base
ref_counts = np.unique(ref, return_counts=True)[1]
ds_props = np.zeros_like(ds_counts, dtype=float)
for i in range(len(ref_counts)):
    ds_props[i*3:(i+1)*3] = ds_counts[i*3:(i+1)*3] / ref_counts[i]
    

ref_ds = ['A','A','A','C','C','C','G','G','G','T','T','T']
alt_ds_1 = ['C','C','T','A','A','T','A','A','C','A','A','C']
alt_ds_2 = ['T','G','G','T','G','G','C','T','T','C','G','G']

ds_frame = pd.DataFrame({'ref': ref_ds, 'alt1': alt_ds_1, 'alt2': alt_ds_2, 'props': ds_props})

### Calculate the expected number of doubletons for each possible nucleotide change
doub_alt = ['C','G','T','A','G','T','A','C','T','A','C','G']
doub_preds = pd.DataFrame({'ref': ref_ds, 'alt':doub_alt, 'prob': np.zeros(12), 'expected': np.zeros(12)})

for i in range(len(doub_preds)):
    refcheck = (doub_preds.iloc[i,0] == ds_frame['ref'])
    ref_frame = ds_frame[refcheck]
    altcheck = (doub_preds.iloc[i,1] == ref_frame['alt1'])|(doub_preds.iloc[i,1] == ref_frame['alt2'])
    alters = ref_frame[altcheck]
    non_alters = (ref_frame[-altcheck])
    doub_preds.iloc[i,2] = ((alters.iloc[0,3])*(alters.iloc[1,3])/(non_alters.iloc[0,3]))

for i in range(len(ref_counts)):
    doub_preds.iloc[i*3:(i+1)*3,3] = doub_preds.iloc[i*3:(i+1)*3,2] * ref_counts[i]
    
### Compare overall number of doubletons to expected number, to give prediction of proportion that are non-recurrent,
### and then the probability that a given one is non-recurrent   
expect_doubs = sum(doub_preds.iloc[:,3])
final_prop = sum(doub_preds.iloc[:,3])/np.count_nonzero((ac_f[:,0]==(ngenomes-2)) & (ac_f[:,1]==2) & (ac_f[:,2]==0) & (ac_f[:,3]==0))
final_prob = final_prop/2
print(final_prob)

### Write results to file
filename = "/rds/general/user/jjr18/results/" + str(type) + "/tri_pro_"+ str(sim) +".txt"

doubletons = np.count_nonzero((ac_f[:,0]==(ngenomes-2)) & (ac_f[:,1]==2) & (ac_f[:,2]==0) & (ac_f[:,3]==0))

with open(filename, "w") as file:
    # Write the values
    file.write(f"{expect_doubs}\t{doubletons}\n")
