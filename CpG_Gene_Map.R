#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
chr.id <- as.numeric(slurm_arrayid)

library(data.table)
library(optparse)
###Input option list###
option_list <- list(
  make_option("--trait", type = "character", default = FALSE, action = "store", help = "Tissue"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

trait <- opt$trait

setwd(paste0("/gpfs/research/chongwu/austinking/ACAT_project/summary-statistics/",trait))
save_datdir = "/gpfs/research/chongwu/austinking/ACAT_project/MWAS/phenotype"
save_datdir = paste(save_datdir,"/",trait,"/",sep="")

##read in CpG sites; MWAS_result.rda has already been run through ReIVpred step 1 and 2
load(paste(trait,"-MWAS_result.rda",sep=""))

final.dat = readRDS("/gpfs/research/chongwu/austinking/ACAT_project/MWAS/gencode.v26.hg19.genes_chr1_22.rds")
final.dat = as.data.frame(final.dat)

final.dat$CpG = NA

final.dat = final.dat[final.dat[,"chr"]==chr.id,]


dat2 = dat[dat[,5]==paste("chr",chr.id,sep=""),]


##Map CpG sites to ensembl gene ids
for(indx in 1:dim(final.dat)[1]) {
  tmp.dat = dat2[dat2$P0 >= (final.dat[indx,"startbp"] - 500) & dat2$P0 <= (final.dat[indx,"endbp"] + 500),]
  if(dim(tmp.dat)[1] > 0) {
    final.dat[indx,"CpG"] = paste(tmp.dat[,"gene"],collapse=";") 
  }
}
saveRDS(final.dat,file = paste(save_datdir,"processed_gene_CpG_inf_CHR",chr.id,".rds",sep=""))



