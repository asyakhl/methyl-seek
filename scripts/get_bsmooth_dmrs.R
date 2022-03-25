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
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 20000 * 1024^2)

BS_samples=args[1]
BS_samples_fit=args[2]
bis_files = read.delim(args[3])

comp=args[4]
a=unlist(strsplit(comp,"vs"))[1]
b=unlist(strsplit(comp,"vs"))[2]

count=0
l1=list()
for (i in seq_along((1:22, "X","M","Y"))){
  samples=readRDS(BS_samples[i])
  fit=readRDS(BS_samples_fit[i])
  pData(fit)$group=bis_files$group
  pData(samples)$group=bis_files$group

  BS.cov <- getCoverage(fit)

  keepLoci.ex <- which(rowSums(BS.cov[, which(samples$group ==a)] >= 2) >= 2 &
                       rowSums(BS.cov[, which(samples$group == b)] >= 2) >= 2)
  length(keepLoci.ex)

  fit <- fit[keepLoci.ex,]
  if (dim(fit)[1]==0){
    print(cat(a "vs" b "chr",i, "empty fit"))
    next
  }

  fit_tstat <- BSmooth.tstat(fit,
                                      group1 = which(fit$group==a),
				      group2 = which(fit$group==b),
				      estimate.var = "group2",
                                      local.correct = TRUE,
                                      verbose = TRUE)

  dmrs0 <- dmrFinder(fit_tstat, qcutoff = c(0.025, 0.975))

  dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

  if (dim(dmrs)[1]==0){
    print(cat(a "vs" b "chr",i, " no dmrs"))
  }else{
    count=count+1
    l1[[count]]=dmrs
  }
}

df=bind_rows(l1)
write.table(df,file=args[5],
            quote=F,sep="\t",col.names = T, row.names = F)
