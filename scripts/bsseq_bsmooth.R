#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(bsseq)
library(tidyverse)
library(dplyr)
library(broom)
library(tidyr)
library(parallel)
library(BiocParallel)

library(future)
library(future.apply)
plan("multiprocess", workers = 30)
options(future.globals.maxSize = 47000 * 1024^2)

bis_files = read.delim(args[2])

filenames=bis_files[,4]
samplenames=bis_files[,1]

chr=args[1]

fun1=function(x,y,z){
  tt = read.bismark(files=x , colData = DataFrame(row.names = y),
                    strandCollapse=T, rmZeroCov=T, verbose=T,
                    BPPARAM = MulticoreParam(workers = 20),
                    loci=GRanges(z, IRanges(1:249000000, width=1)))
  return(tt)
}
print("Finished uploading the data")
tt1=mapply(fun1, filenames,samplenames,chr)
print("Transformed the data with mapply")
names(tt1)=samplenames

BS_samples = combineList(tt1)
print("Finished combineList the data")
pData(BS_samples)$Type=bis_files$group

saveRDS(BS_samples, file=args[3])
print(paste("BS_samples is saved",chr))

BS_samples_fit <- BSmooth(
  BSseq = BS_samples, BPPARAM = MulticoreParam(workers = 30),
  verbose = TRUE)
print("Finished BSmoothing the data")
saveRDS(BS_samples_fit, file=args[4])
print(paste("BS_samples_fit is saved",chr))

avg.cov=round(colMeans(getCoverage(BS_samples)), 1)

write.table(avg.cov,file=args[5],quote=F,sep=",",col.names = T, row.names = F)
