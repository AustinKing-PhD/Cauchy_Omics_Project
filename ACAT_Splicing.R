library(data.table)
library(optparse)

##ACAT function from Yaowu Liu github##
ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
  Pvals<-as.matrix(Pvals)
  if (is.check){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
      stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
      stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(colSums(Pvals==0)>=1)
    is.one<-(colSums(Pvals==1)>=1)
    if (sum((is.zero+is.one)==2)>0){
      stop("Cannot have both 0 and 1 p-values in the same column!")
    }
    
    if (sum(is.zero)>0){
      warning("There are p-values that are exactly 0!")
    }
    if (sum(is.one)>0){
      warning("There are p-values that are exactly 1!")
    }
    
  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights)){
    is.weights.null<-TRUE
  }else{
    is.weights.null<-FALSE
    weights<-as.matrix(weights)
    if (sum(dim(weights)!=dim(Pvals))>0){
      stop("The dimensions of weights and Pvals must be the same!")
    }else if (is.check & (sum(weights<0)>0)){
      stop("All the weights must be nonnegative!")
    }else{
      w.sum<-colSums(weights)
      if (sum(w.sum<=0)>0){
        stop("At least one weight should be positive in each column!")
      }else{
        for (j in 1:ncol(weights)){
          weights[,j]<-weights[,j]/w.sum[j]
        }
      }
    }
    
  }
  
  #### check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small<-(Pvals<1e-15)
  if (is.weights.null){
    Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-1/Pvals[is.small]/pi
    cct.stat<-colMeans(Pvals)
  }else{
    Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
    cct.stat<-colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval<-pcauchy(cct.stat,lower.tail = F)
  return(pval)
}



###Input option list###
option_list <- list(
  make_option("--trait", type = "character", default = FALSE, action = "store", help = "Tissue"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

trait <- opt$trait

setwd(paste0("/gpfs/research/chongwu/austinking/ACAT_project/Splicing/phenotype/",trait))
save_datdir = "/gpfs/research/chongwu/austinking/ACAT_project/Splicing/phenotype"
save_datdir = paste(save_datdir,"/",trait,"/",sep="")

##read in Splicing data
load("Splicing_result_mapped.rda")

#check input number of unique genes
gene.list = unique(dat$gene_id)
length(gene.list)

##ACAT each gene individually
dat2 = NULL
for (w in 1:length(gene.list)) {
  gene.id = gene.list[w] #get tissue gene ID
  tmp2 = dat[dat$gene_id==gene.id,] 
  tmp2$ACAT_pValue = ACAT(tmp2$pValue,weights = NULL)
  dat2 = rbind(dat2,tmp2)
}

#Check number of unique genes after ACAT
gene.list2 = unique(dat2$gene_id)
length(gene.list2)

#Filter wght file 
dat = dat2
wght = wght[wght[,1] %in% dat[,1],]

#check number of SNPs post filter
nrow(wght)

save(dat,wght, file = paste(save_datdir,"/Splicing_ACAT_result.rda",sep=""))
