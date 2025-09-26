Readme file for scripts for identification of doubletons, doubleton haplotypes, linked rare variants, and random position haplotypes.

Many of the scripts here use those published at https://github.com/mathii/f2 for their paper "Demography and the age of rare variants". I would recommend first cloning their repo, and then adding/replacing the following that we have provided. Please ensure that if using these to also cite their work.

Following files are present:

doubleton_haplofinder.sh - A modified version of a Mathieson and McVean script for doubleton haplotype identification.

random_haplofinder.sh - A modified version of a Mathieson and McVean script for random pair haplotype identification.

fn_finder.sh - A slightly modified version of a Mathieson and McVean script, designed to identify variants between a lower and upper count in the population and determine the individuals sharing them.

counting.R - A script to count the number of variants of a given allele count that are found on each doubleton haplotype.

random_counting.R - A mirror of the above, but designed to do this for random pair haplotypes.

random_pair_generator.R - A small script that, given a file with the positions of doubletons and the individuals who share them, will randomly sample other pairs of individuals at those positions that do not carry the doubleton.


Basic pipeline:

1. Run doubleton_haplofinder to find doubleton haplotypes in chosen VCF.

2. Run fn_finder to find other rare variants, up to desired limit.

3. Run fn_finder to find linked rare variants on doubleton haplotypes.

4. Run random_pair_generator and random_haplofinder to find matched haplotypes around random pairs of individuals at doubleton sites.

5. Run random_counting to find linked rare variants on random pair haplotypes.

6. Run tri_rate to find the expected rate of recurrent mutation





