args <- commandArgs(TRUE)
file.loc <- args[1]
nsamp <- as.integer(args[2])

setwd(file.loc)
posfile <- read.table('pos.idx.f2.gz')
chroms <- nsamp*2


for(i in 1:nrow(posfile)){
  repeat {
    id1 <- sample(0:(chroms-2),1)
    # exit if the condition is met
    if (id1 != posfile$V2[i]) break
  }
  posfile$V2[i] <- id1
  repeat {
    id2 <- sample(1:(chroms-1),1)
    # exit if the condition is met
    same <- (id2 != posfile$V3[i])
    bigger <- id2 > posfile$V2[i]
    if (same && bigger == T) break
  }
  posfile$V3[i] <- id2
}

write.table(posfile, file = "pos.idx.random.gz", row.names = F, sep = "\t", quote = F, col.names = F)
