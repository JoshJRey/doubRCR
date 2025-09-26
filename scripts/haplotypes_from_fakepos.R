# find the haplotypes using the f2 variants.
# args should be 1) the directory with the haplotypes in (f1, f2 and by sample) 2) The output file and 3) The map file/rate
# 4) Singleton power, 5) mutation rate
## Singleton power is not used here.
## I have added a new arguement to this file called batch, to allow for further parallelisation. Each batch will split the chromosome into that many fragments, which in turn are split into as many
## chunks as specified by n.cores. Each split is run over a separate core, to allow for parallel processing. Each batch needs to be submitted as a separate task on the HPC, so you'll have several batches running 
## at once to speed things up. If not wanted, just set batch to 1. Make sure splits and n.cores are the same! Splits refers to spliting of the batches, was set up when a batch didn't complete in the time allowance.

args <- commandArgs(TRUE)

n.cores <- 1
f2.filename <- "/pos.idx.random.gz"
if((length(args) %in% c(7,8,9))){
  code.dir <- args[1]
  hap.root <- args[2]
  out.file <- args[3]
  map.file <- args[4]
  verbose <- as.numeric(args[5])
  batch <- as.numeric(args[6])
  splits <- as.numeric(args[7])
  if(length(args)==8){
      n.cores <- as.numeric(args[8])
  }
    if(length(args)==9){
      f2.filename <- args[9]
  }

} else{
  stop("Need to specify 5 or 6 or 7 arguments")
}
## Here, we are sourcing the code from the files indicated by include.R
source(paste(code.dir, "/libs/include.R", sep=""))

## Here, we load in the files with the locations of f1 and f2 variants
f1.file <- paste(hap.root, "/pos.idx.f1.gz", sep="")
f2.file <- paste(hap.root, f2.filename, sep="")


## Here, we're loading in the pos file, which details the locations of all of the variants in our sample.
pos.file <- paste(hap.root, "/by_sample/pos.gz", sep="")
by.sample.gt.root <- paste(hap.root, "/by_sample", sep="")
samples <- scan(paste(hap.root, "/samples.txt", sep=""), quiet=FALSE, what="")
sim.pop.map <- rep("ALL", length(samples))
names(sim.pop.map) <- samples


print(pos.file)
print(by.sample.gt.root)
print(f2.file)
print(f1.file)

## Below, we can either run on a single core, or paralellise. If we parallelise, we'll split the chromosome into separate chunks to scan.

if(n.cores==1){
    haps <- find.haplotypes.from.f2(f1.file, f2.file, pos.file, by.sample.gt.root, sim.pop.map, map.file, verbose=verbose)
}else{
    require(parallel)
    pos <- scan(pos.file, quiet=FALSE)
    pos <- range(pos)+c(-1e6, 1e6)
    breaks <- seq(pos[1], pos[2], length.out=(n.cores+1))
    poss <- list()
    for(i in (1+((batch-1)*(n.cores/splits))):((n.cores/splits)*batch)){
        poss[[i]] <- breaks[c(i,i+1)]
    }
    
    haps <- mclapply(poss, find.haplotypes.from.f2.vector, mc.cores=n.cores, f1.file=f1.file, f2.file=f2.file, pos.file=pos.file, by.sample.gt.root=by.sample.gt.root, pop.map=sim.pop.map,  map.file=map.file, verbose=TRUE)
    haps <- do.call(rbind, haps)
}



# Just remove the ones with zero length.

cat(paste0("Removed ", sum(haps$map.len==0), " haplotypes with length 0\n"))
haps <- haps[haps$map.len>0,]

write.table(haps, out.file, row.names=FALSE, col.names=TRUE, quote=FALSE)
