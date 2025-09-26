#!/bin/bash -vx

# Modified version of Mathieson and McVean (2014) script.
# Finds haplotypes between any two individuals at any given site in genome.

##############################################################################################################
##Command line inputs

CHR=$1  ##The chromosome (or chromosome arm) we are working with
path_to_seq_data=$2  ##The location of our sequence data, in vcf.gz format. This must must must be biallelic and have no missing data, or the scripts will not work.
nbp=$3  ##Length of the chromosome (in base pairs)
ncores=$4 ##Number of cores to use, for parallelisation of haplotype finding script
type=$5 ##Used if performing replicates (i.e. for simulations) - gives a categorisation to the results folder
##############################################################################################################

# Edit parameters here. 
RES_DIR=/rds/general/user/jjr18/results/${type}/${CHR}

# Where are the recombination maps, in impute format
HM2_MAP=/rds/general/user/jjr18/home/f2/setup/Ag_map_${CHR}_real.txt   ##Recombination rate map. Example will be provided to demonstrate format.

# Where is the code - this point to the directory you downloaded from github
CODE_DIR=/rds/general/user/jjr18/home/f2  ##mathii/f2 repository

##############################################################################################################

##Below, we are setting up the directory structure used by the scripts. In the TH (haplotypes) folder, we will store the adjusted genotype files, as well as files with the location of any f1 or f2 snps.
##The results folder will contain the results (obvs).
TH=${RES_DIR}/haplotypes
RD=${RES_DIR}/results
CD=${CODE_DIR}


##Here, we are piping any verbose outcomes of these scripts into a log file. Specifically, we're assigning a variable to be a file name, and then the exec line tells outputs to pipe into that file.
LOG=${RD}/log_randoms.txt
exec > ${LOG} 2>&1

##############################################################################################################

# set max number of file descriptors - sometimes if working with big VCFs default is too low
ulimit -n 60000


# 1) Find f2 haplotypes
## Here, we are running a script that allows us to identify our haplotypes. When given a set of focal mutations, this will identify our haplotypes surrounding them, and then generate a table with 
## the summary. More details in the script file. The R flags --vanilla --quiet --slave are used to pipe our scripts into a suitable R session.
## IMPORTANT NOTE - this script is very slow, as we are doing many pairwise comparisons. Use ncores to try and parallelise the process if running slow. Note that this will produce duplicated haplotypes if they extend into multiple chunks, for more info see 
## f2/scripts/haplotypes_from_f2.R. If want to parallelise further, I have developed a different script which allows for more splitting, please reach out if desired. 

R --vanilla --quiet --slave --args ${CD} ${TH} ${RD}/random_haplotypes_${CHR}.txt \
    ${HM2_MAP} 1 1 1 ${ncores} < ${CD}/scripts/haplotypes_from_random.R
gzip -f ${RD}/random_haplotypes_${CHR}.txt


