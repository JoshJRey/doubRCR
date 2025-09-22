# doubRCR
Scripts for the identification of doubletons who are reciprocal closest relatives (RCR) from population genomic data.


##Simulations

This directory contains scripts used for simulations to test the effectiveness of our approach to identifying RCR doubletons. Much of this process is the same 
as for empirical data. More information can be found in the readme in /simulations.

Below is an overview of the steps required to run simulations to test the effectiveness of our approach to identifying RCR doubletons. More detailed documentation found in /simulations.

1. Run simulation (here using msprime, although one could also use SLiM).

2. Overlay mutations.

3. Convert tree sequence to VCF.

4. Identify recurrent sites within tree sequence.

5. Identify doubleton haplotypes within VCF.

6. Identify other low frequency variants within VCF.

7. Identify and count occurrences of linked rare variants on doubleton haplotypes.

8. Identify 
