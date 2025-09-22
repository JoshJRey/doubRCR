Readme for scripts used for simulations of doubleton haplotypes.

All scripts are set with some parameters (population sizes, seed, genome length) that can be varied from the command line. Other parameters (recombination rate, number of samples, ancestry model, time of expansion (Expansion only), number of demes (Island only), migration rate (Island only)) are not initially set to be varied on the command line, but could be easily modified if desired.

Below is an overview of the steps required to run simulations to test the effectiveness of our approach to identifying RCR doubletons.

Generate simulated genomes:
1. Run simulation (here using msprime, although one could also use SLiM).

2. Overlay mutations.

3. Convert tree sequence to VCF.

Extract information on doubletons and their haplotypes:
4. Identify recurrent doubleton sites within tree sequence.

5. Identify doubleton haplotypes within VCF.

6. Identify other low frequency variants within VCF.

7. Identify and count occurrences of linked rare variants on doubleton haplotypes.

Extract information on random pairs:
8. Identify haplotypes around randomly sampled pairs of individuals at doubleton sites.

9. Identify and count occurrences of linked rare variants on random pair haplotypes.

Identify putative RCR doubletons
10. Estimate recurrent rate from triallelic singletons.

11. Combine data from each simulation.

11. Estimate probability of being RCR.

