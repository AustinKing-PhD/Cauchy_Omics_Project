#!/usr/bin/env Rscript
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
#chr.id <- as.numeric(slurm_arrayid)

library(data.table)
library(optparse)
###Input option list###
option_list <- list(
  make_option("--trait", type = "character", default = FALSE, action = "store", help = "Tissue"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

trait <- opt$trait

setwd(paste0("/gpfs/research/chongwu/austinking/ACAT_project/MWAS/phenotype/",trait))
save_datdir = "/gpfs/research/chongwu/austinking/ACAT_project/MWAS/phenotype"
save_datdir = paste(save_datdir,"/",trait,"/",sep="")

final.dat = NULL
for(chr.id in 1:22) {
  tmp = readRDS(paste("processed_gene_CpG_inf_CHR",chr.id,".rds",sep=""))
  final.dat = rbind(final.dat,tmp)
}

saveRDS(final.dat,file = paste(save_datdir,"processed_gene_CpG_inf.rds",sep = ""))

#Add gene ids to dat and wght
library(tidyr)
library(dplyr)

setwd(paste0("/gpfs/research/chongwu/austinking/ACAT_project/summary-statistics/",trait))
load(paste(trait,"-MWAS_result.rda",sep=""))

final.dat2 = final.dat[,c("gene_id","CpG")]

final.dat2_sep = final.dat2 %>% separate_rows(CpG, sep = ";")
tmp.final.dat = inner_join(dat,final.dat2_sep, by = c("gene"="CpG"))
tmp.final.dat = tmp.final.dat[,c("gene","#SNPs","Zscore","pValue","CHR","P0","gene_id")]
dat = tmp.final.dat

ID = tmp.final.dat$gene
wght = wght[wght[,"gene"] %in% ID,]

#Add Gene IDs to SNPs
wght = as.data.frame(wght)
wght = inner_join(wght,final.dat2_sep, by = c("gene"="CpG"))
wght = as.matrix(wght)

save(dat,wght, file = paste(save_datdir,"/MWAS_result_mapped.rda",sep=""))
