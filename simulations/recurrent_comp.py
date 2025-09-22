#!/usr/bin/env python3

### Import packages
import msprime
import tskit
import numpy as np
import pandas as pd
import allel

### Below code can be run to compare the true proportion of recurrent doubletons from a tree sequence
### to the expected proportions from the number of tri-allelic singletons (sometimes referred to here as double-singletons).
### In it's present state, this is set up to be run in a Jupyter notebook, but could easily be adapted to run as a script.

hky_model = msprime.HKY(
    equilibrium_frequencies=[0.3, 0.2, 0.2, 0.3]  # Frequencies for A, C, G, T
)

ts_mutated_hky = msprime.sim_mutations(ts, rate=8.125e-8, random_seed=1, model=hky_model)

print("Pi = ", ts_mutated_hky.diversity())

variants = ts_mutated_hky.variants()
genotypes = np.array([variant.genotypes for variant in variants])

genotypes_3d = genotypes.reshape(genotypes.shape[0], -1, 2)

# Convert the reshaped genotype matrix to an Allel GenotypeArray
genotype_array = allel.GenotypeArray(genotypes_3d)

ngenomes = 2*genotype_array.n_samples

allele_counts = genotype_array.count_alleles()


if (np.count_nonzero(allele_counts.allelism()==4))>0:
    doubletons = ((allele_counts[:,0]==(ngenomes-2)) & (allele_counts[:,1]==2) & (allele_counts[:,2]==0) & (allele_counts[:,3]==0))
else:
    doubletons = ((allele_counts[:,0]==(ngenomes-2)) & (allele_counts[:,1]==2) & (allele_counts[:,2]==0))



print("No. of doubletons = ", np.count_nonzero(doubletons))
num_doubs = np.count_nonzero(doubletons)
      
same_base_sites = []
different_base_sites = []

tree_iter = ts_mutated_hky.trees()
tree = next(tree_iter)
tree.seek(0)

for site in ts_mutated_hky.sites():
    # Advance tree if needed
    while site.position >= tree.interval.right:
        tree = next(tree_iter)

    if len(site.mutations) == 2:
        derived_states = []

        for mut in site.mutations:
            count = sum(1 for _ in tree.samples(mut.node))
            if count == 1:
                derived_states.append(mut.derived_state)

        if len(derived_states) == 2:
            pos = int(site.position)
            if derived_states[0] == derived_states[1]:
                same_base_sites.append(pos)
            else:
                different_base_sites.append(pos)

                
print("No. recurrent doubletons = ", len(same_base_sites))
print((len(same_base_sites)/np.count_nonzero(doubletons))*100)



ts_mutated_hky.dump("/Users/josh/phd/sims/example_instant_change.ts")

!tskit vcf example_instant_change.ts > example_instant_change.vcf

callset = allel.read_vcf("/Users/josh/phd/sims/example_instant_change.vcf")

ac_f = allele_counts

alt = callset['variants/ALT'][:].astype('U')
ref = callset['variants/REF'][:].astype('U')

ngenomes = 2*genotype_array.n_samples

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
ds_counts = ds_counts[[0,1,3,6,7,9,12,13,15,18,19,21]]
print(ds_counts)

### Convert counts into proportion of reference base sites, for each reference base
ref_counts = np.unique(ref, return_counts=True)[1]
ds_props = np.zeros_like(ds_counts, dtype=float)
for i in range(len(ref_counts)):
    ds_props[i*3:(i+1)*3] = ds_counts[i*3:(i+1)*3] / ref_counts[i]
    
ref_ds = ['A','A','A','C','C','C','G','G','G','T','T','T']
alt_ds_1 = ['C','C','T','A','A','T','A','A','C','A','A','C']
alt_ds_2 = ['T','G','G','T','G','G','C','T','T','C','G','G']

print(ds_props)

ds_frame = pd.DataFrame({'ref': ref_ds, 'alt1': alt_ds_1, 'alt2': alt_ds_2, 'props': ds_props})

### Calculate the expected number of doubletons for each possible base change
doub_alt = ['C','G','T','A','G','T','A','C','T','A','C','G']
doub_preds = pd.DataFrame({'ref': ref_ds, 'alt':doub_alt, 'prob': np.zeros(12), 'expected': np.zeros(12)})

for i in range(len(doub_preds)):
    refcheck = (doub_preds.iloc[i,0] == ds_frame['ref'])
    ref_frame = ds_frame[refcheck]
    altcheck = (doub_preds.iloc[i,1] == ref_frame['alt1'])|(doub_preds.iloc[i,1] == ref_frame['alt2'])
    alters = ref_frame[altcheck]
    non_alters = (ref_frame[-altcheck])
    doub_preds.iloc[i,2] = ((alters.iloc[0,3])*(alters.iloc[1,3])/(2*non_alters.iloc[0,3]))

for i in range(len(ref_counts)):
    doub_preds.iloc[i*3:(i+1)*3,3] = doub_preds.iloc[i*3:(i+1)*3,2] * ref_counts[i]
    
### Compare overall number of doubletons to expected number, to give prediction of proportion that are non-recurrent,
### and then the probability that a given one is non-recurrent   

doubletons = ((ac_f[:,0]==(ngenomes-2)) & (ac_f[:,1]==2) & (ac_f[:,2]==0) & (allele_counts[:,3]==0))

final_prop = sum(doub_preds.iloc[:,3])/np.count_nonzero((ac_f[:,0]==(ngenomes-2)) & (ac_f[:,1]==2) & (ac_f[:,2]==0) & (allele_counts[:,3]==0))
final_prob = final_prop #/2
print(sum(doub_preds.iloc[:,3]))
print(final_prob*100)