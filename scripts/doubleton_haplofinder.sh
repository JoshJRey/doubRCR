#!/bin/bash -v

### This is a slightly modified version of the script from Mathieson and McVean (2014), available at https://github.com/mathii/f2

##############################################################################################################
##Command line inputs

CHR=$1  ##The chromosome (or chromosome arm) we are working with
path_to_seq_data=$2  ##The location of our sequence data, in vcf.gz format. This must must must be biallelic with no missing data, or the scripts will not work.
nbp=$3  ##Length of the chromosome (in base pairs)
ne=$4 ##Effective population size
mu=$5 ##Mutation rate
type=$6 ##Used if performing replicates (i.e. for simulations) - gives a categorisation to the results folder
ncores=$7 ##Number of cores to use, for parallelisation of haplotype finding script

##Results directory

##############################################################################################################

# Edit parameters here.

# Where do you want the results to go?
RES_DIR=/rds/general/user/jjr18/results/${type}/${CHR} ##This will be the directory that all of your results go into, so name accordingly

# Where are the recombination maps, in impute format
HM2_MAP=/rds/general/user/jjr18/home/mig/haplos/setup/Ag_map_{CHR}.txt   ##Recombination rate map. Example will be provided to demonstrate format.

# Where is the code - this point to the directory you downloaded from github (##mathii/f2 repository)
CODE_DIR=/rds/general/user/jjr18/home/f2

# two colum tab delim file - col1=sample name, col2=population  - example will be provided
PANEL=/rds/general/user/jjr18/home/f2/setup/${RES_DIR}_panel.txt


##############################################################################################################

##Below, we are setting up the directory structure used by the scripts. In the TH (haplotypes) folder, we will store the adjusted genotype files, as well as files with the location of any f1 or f2 snps.
##The results file will contain the results (obvs).
RD=${RES_DIR}/results
CD=${RES_DIR}

for dir in ${WD} ${TH}/by_sample ${RD}
do
    mkdir -p ${dir}
done

##Here, we are piping any verbose outcomes of these scripts into a log file. Specifically, we're assigning a variable to be a file name, and then the exec line tells outputs to pipe into that file.
LOG=${RD}/log_haplo_setup.txt
exec > ${LOG} 2>&1

##############################################################################################################

# 1) Parse data into correct format
# 1.1) Extract sample names
##Self explanatory, but these are the samples you want to work with. So, if you want to work with a subset, first exclude those from your panel file.
cut -f1 ${PANEL} > ${TH}/samples.txt

# 1.2) Sequence data
# set max number of file descriptors - sometimes if working with big VCFs default is too low
ulimit -n 3000


# Keep only the samples we want - can skip this step if you know that you
# are keeping all the samples and sites in the chip data. 
## Scripts do not work if we have missing data, so we exlcude sites with missing data using --max-missing 1. This has also been modified from the original to now work with vcftools 0.1.14. If you want 
## only a subset of samples, be sure to make sure your samples.txt has that subset. Effectively, this step is just excluding missing sites, and tidying the vcffile, before piping it into a new one 
## using --recode and --stdout

vcftools --gzvcf ${path_to_seq_data} --keep ${TH}/samples.txt \
    --max-missing 1 \
    --min-alleles 2 --max-alleles 2 \
    --out ${TH}/chr${CHR}.tmp --recode --stdout | gzip -c \
    > ${TH}/chr${CHR}.tmp.vcf.gz

# Convert to flat format
## By flat format, they mean with genotype calls combined (i.e. 1/1 would be 2) - we do this as our data is unphased. Note that the columns called by awk ($2, i=10 etc) have been altered, as 
## Ag vcf files have a different structure to human ones. This script basically says "take the location, and then the number either side of each call, and nothing else!)

zgrep -v "^#" ${TH}/chr${CHR}.tmp.vcf.gz \
   | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
   | gzip -c > ${TH}/all.chr${CHR}.haps.gz

# Split up by sample. 
##Largely self explanatory, but now we are taking these "flattened" haplotypes, and separating them so that each individual has a file detailing their own haplotype, and only theirs. This is uses
##a script contained in f2/scripts. The outputs of this are stored in results/your_folder_name/haplotypes/by_sample, and should have files for each individual, with a series of 0, 1, and 2s.

python ${CD}/scripts/haps_to_gt_bysample.py -h ${TH}/all.chr${CHR}.haps.gz \
    -o ${TH}/by_sample/ -s ${TH}/samples.txt


# Cleanup
##For neatness, we can clean out temporary files. If debugging, you may want to comment these lines out to keep the intermediate files.
rm ${TH}/chr${CHR}.tmp.log
rm ${TH}/chr${CHR}.tmp.vcf.gz
rm ${TH}/all.chr${CHR}.haps.gz

# 1.3) Sequence data - extract singletons and doubletons
## This calls on a script in f2/scripts
## If not using singleton or doubleton data, no need to run this code!

for n in 1 2
do
   vcftools --gzvcf ${TH}/chr${CHR}.tmp.vcf.gz  --out ${TH}/chr${CHR}.f${n}.log.tmp \
        --recode --remove-indels --max-missing 1 --stdout --mac $n --max-mac $n \
       | grep -v "^#" \
       | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
       | gzip -c > ${TH}/chr${CHR}.f${n}.haps.gz

   python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/chr${CHR}.f${n}.haps.gz \
       -o ${TH}/pos.idx.f${n}.gz -n $n > ${TH}/pos.idx.f${n}.tmp.log 
done

# 2) Find f2 haplotypes
## Here, we are running a script that allows us to identify our haplotypes. When given a set of focal mutations, this will identify our haplotypes surrounding them, and then generate a table with 
## the summary. More details in the script file. The R flags --vanilla --quiet --slave are used to pipe our scripts into a suitable R session.
## VERY IMPORTANT NOTE - this script is very slow, as we are doing many pairwise comparisons. The 32 here is the number of cores used, to try and parallelise the process. This is the best trade off
## available on the HPC regarding available cores and script length. Note that this will produce duplicated haplotypes if they extend into multiple chunks, for more info see 
## f2/scripts/haplotypes_from_f2.R. If want to parallelise further, I have developed a different script which allows for more splitting. 

R --vanilla --quiet --slave --args ${CD} ${TH} ${RD}/f2_haplotypes_${CHR}.txt \
    ${HM2_MAP} 1 1 1 ${ncores} < ${CD}/scripts/haplotypes_from_f2.R
gzip -f ${RD}/f2_haplotypes_${CHR}.txt 




