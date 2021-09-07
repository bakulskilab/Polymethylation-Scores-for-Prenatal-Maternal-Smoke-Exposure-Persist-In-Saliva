#libraries
library(sjlabelled)
library(tidyr)
library(dplyr)
library(purrr)
library(furrr)
library(stringr)
library(sva)


#directory based programming
data.dir="/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/"
output.dir=paste0(data.dir, 'CreatedData/')
#variables 
var.filter=5000

#load in betaqc matrix
betaqc<-readRDS(file=paste0(data.dir, 'OGData/', "noob_filtered.rds"))
file.path(paste0(data.dir, 'CreatedData/'), list.files(paste0(data.dir, 'CreatedData/')))[(str_detect(file.path(paste0(data.dir, 'CreatedData/'), list.files(paste0(data.dir, 'CreatedData/'))), 'NoSmk'))]%>%map(load, envir=.GlobalEnv)

pheno<-rbind(child_nosmoke, teen_nosmoke) # dataset with the most exclusion i.e. squaring dataset off for all models (excluding kids who smoke in figure)

pheno.ancestry=pheno%>% group_by(ancestry)%>%nest()

#function to apply to each ancestry group
my.sva=function(x, vf=1000,...){
  bqc <-betaqc[, colnames(betaqc) %in% x$MethID]
  try(if(nrow(x)!=ncol(bqc)) stop("# obs in betaqc & pheno do not match"))
  x=x[match(colnames(bqc), x$MethID), ]
  try(if(identical(colnames(bqc), remove_label(x$MethID))==FALSE) stop('STOP: Sample names do not match'))
  #make model matrix including variable of interest as factor variable & adjustment variables 
  mod=model.matrix(~as.factor(smkPreg_binary), data=x)
  #make null model containing only adjustment variables
  mod0=model.matrix(~1, data=x)
  svobj=sva(bqc, mod, mod0, vfilter=vf)
  sv=data.frame(svobj$sv)
  colnames(sv)=paste('sv', 1:ncol(sv), sep='')
  sv=cbind(x%>%select(id, idnum, MethID), sv)
  return(sv)
}

#apply to each ancestry group
sv.ancestry=pheno.ancestry%>%mutate(sv=map(data, my.sva, vf=5000))

#
save(sv.ancestry, file=paste0(output.dir, 'svAncestry5k.Rds'))
print('Finished running SVA at Variance Filter of 5000 at different ancestries ')