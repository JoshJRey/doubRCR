### Script set up for remote running of linked rare variant finding for random pair haplotypes. 

#First, the packages.
require(tidyr)
library(plyr)
library(dplyr)

### Set variables from command line
args <- commandArgs(TRUE)
lower <- as.numeric(args[1]) ## Lower bound of allele count
upper <- as.numeric(args[2]) ## Upper bound of allele count
chr <- args[3] ## Chromosome (or simulation number)
n.cores <- as.numeric(args[4]) ## Number of cores to use
type <- args[5] ## Used for further subsetting of resuts (predominantly for simulations)


load.fn <- function(fn.file){
  if (file.info(fn.file)$size == 132) {
    warning(paste("Empty file:", fn.file))
    return(data.frame(pos = numeric(0), ID1 = integer(0), ID2 = integer(0)))
  }

  fn <- read.table(fn.file, as.is = TRUE, header = FALSE)

  if (nrow(fn) == 0) {
    warning(paste("File has no lines:", fn.file))
    return(data.frame(pos = numeric(0), ID1 = integer(0), ID2 = integer(0)))
  }

  idx.range <- 2:length(fn)
  fn[, idx.range] <- fn[, idx.range] + 1
  colnames(fn) <- c("pos", paste("ID", idx.range - 1, sep = ""))
  fn[, idx.range] <- ceiling(fn[, idx.range] / 2)
  return(fn)
}


hap.file <- paste0("/rds/general/user/jjr18/results/", type, "/", chr, "/results/random_haplotypes_", chr, ".txt.gz")
haplotypes <- read.table(hap.file, header = T)

pos.file <- paste0("/rds/general/user/jjr18/results/", type, "/", chr, "/haplotypes/by_sample/pos.gz")
### Set out counting functions


fn.count.vector <- function(pos.lim){
  return(fn.count(pos.lim=pos.lim))
}
  
fn.count <- function(pos.lim){
  haplotypes <- haplotypes[haplotypes$pos>=pos.lim[1] & haplotypes$pos<pos.lim[2],]
  haplotypes <- haplotypes[order(haplotypes$ID1,haplotypes$ID2),]
  variant_count <- data.frame(haplotypes$pos)
  names(variant_count)[1] <- "pos"
  for(i in lower:upper){
    variants <- load.fn(paste0("/rds/general/user/jjr18/results/", type, "/", chr, "/haplotypes/pos.idx.f", i, ".gz"))
    
    variant_count[[paste0("f", i)]] <- rep(0, nrow(variant_count))
    ## load in the required variant frequency and set up a new column with it's name
    for(j in 1:NROW(haplotypes)){
      cat(paste("\r(", j, "/", NROW(haplotypes),"), ", i, sep=""))
      iden1 <- variants %>%
        filter_at(vars(starts_with("ID")), any_vars(. == haplotypes$ID1[j]))%>%
        filter_at(vars(starts_with("ID")), any_vars(. == haplotypes$ID2[j])) %>%
        filter_at(vars(starts_with("pos")), any_vars(haplotypes$hap.left[j] < .)) %>%
        filter_at(vars(starts_with("pos")), any_vars(haplotypes$hap.right[j] > .))
      count <- NROW(iden1)
      variant_count[[paste0("f",i)]][j] <- count
    }
  
  }
  return(variant_count)
}

###Split positions, and then run

if(n.cores==1){
  cat("Error")
}else{
    require(parallel)
    pos <- scan(pos.file, quiet=TRUE)
    pos <- range(pos)+c(-1e3, 1e3)
    breaks <- seq(pos[1], pos[2], length.out=n.cores+1)
    poss <- list()
    for(i in 1:n.cores){
        poss[[i]] <- breaks[c(i,i+1)]
    }
    haps <- mclapply(poss, fn.count.vector, mc.cores=n.cores)
    haps <- do.call(rbind, haps)
}


  
write.table(haps, paste0("/rds/general/user/jjr18/results/", type, "/", chr, "/results/random_", chr, "_f", lower, "_f", upper, ".txt"), quote = F, row.names = F)

