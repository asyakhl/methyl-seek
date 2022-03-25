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

BS_samples=readRDS(args[1])
BS_samples_fit=readRDS(args[2])
bis_files = read.delim(args[3])

pData(BS_samples_fit)$group=bis_files$group
pData(BS_samples)$group=bis_files$group

BS.cov <- getCoverage(BS_samples_fit)

comp=args[4]
a=unlist(strsplit(comp,"vs"))[1]
b=unlist(strsplit(comp,"vs"))[2]
keepLoci.ex <- which(rowSums(BS.cov[, which(BS_samples$group ==a)] >= 2) >= 2 &
                       rowSums(BS.cov[, which(BS_samples$group == b)] >= 2) >= 2)
length(keepLoci.ex)

BS_samples_fit <- BS_samples_fit[keepLoci.ex,]
i=args[5]
if (dim(BS_samples_fit)[1]==0){
  print(cat(a "vs" b "chr",i, "empty BS_samples_fit"))
}

BS_samples_fit_tstat <- BSmooth.tstat(BS_samples_fit,
                                      group1 = which(BS_samples_fit$group==a),
                                      group2 = which(BS_samples_fit$group==b),
                                      estimate.var = "group2",
                                      local.correct = TRUE,
                                      verbose = TRUE)

dmrs0 <- dmrFinder(BS_samples_fit_tstat, qcutoff = c(0.025, 0.975))

dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

if (dim(dmrs)[1]==0){
  print(cat(a "vs" b "chr",i, " no dmrs"))
}
write.table(dmrs,file=args[6],
            quote=F,sep="\t",col.names = T, row.names = F)
