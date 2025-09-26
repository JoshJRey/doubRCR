#!/bin/bash -vx

# Modified version of Mathieson and McVean (2014) script. Finds positions of alleles of given count, and individuals carrying them. 				

##############################################################################################################

CHR=$1 ##Chromosome (or chromosome arm) we are working with
path_to_seq_data=$2 ##The location of our sequence data, in vcf.gz format. This must must must be biallelic and have no missing data, or the scripts will not work.
lower=$3 ##Lower bound of allele count
upper=$4 ##Upper bound of allele count
type=$5 ##Used if performing replicates (i.e. for simulations) - gives a categorisation to the results folder
##############################################################################################################

# Edit parameters here. 

# Where do you want the simulations to go?
RES_DIR=/rds/general/user/jjr18/results/${type}/${CHR}
# Where is the code - this point to the directory you downloaded from github
CODE_DIR=/rds/general/user/jjr18/home/f2

##############################################################################################################

TH=${RES_DIR}/haplotypes
RD=${RES_DIR}/results
CD=${RES_DIR}


LOG=${RD}/log_fn.txt
exec > ${LOG} 2>&1
# assume that doubleton_haplofinder has already been run.

for n in $(seq ${lower} ${upper})   
do
   vcftools --gzvcf ${path_to_seq_data}  --out ${TH}/chr${CHR}.f${n}.log.tmp \
       --recode --stdout --mac $n --max-mac $n --max-missing 1 \
       --min-alleles 2 --max-alleles 2 | grep -v "^#" \
       | awk '{{printf "%d",$2};for(i=10; i<=NF; i++){printf "\t%d\t%d",substr($i,1,1),substr($i,3,1)};printf "\n"}'  \
       | gzip -c > ${TH}/chr${CHR}.f${n}.haps.gz

   python ${CD}/scripts/calculate_fn_sites.py -h ${TH}/chr${CHR}.f${n}.haps.gz \
       -o ${TH}/pos.idx.f${n}.gz -n $n > ${TH}/pos.idx.f${n}.tmp.log 

done

