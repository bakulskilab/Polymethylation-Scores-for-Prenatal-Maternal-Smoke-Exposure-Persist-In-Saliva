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

#make a matrix of model data 
#filter beta matrix to only those in analytic subset 
betaqc <-betaqc[, colnames(betaqc) %in% pheno$MethID]
#check that dim of pheno and betaqc match 
try(if(nrow(modeldata)!=ncol(betaqc)) stop("# obs in betaqc & pheno do not match"))

#want to make sure that the order of pheno and betaqc match
pheno=pheno[match(colnames(betaqc), pheno$MethID), ]
try(if(identical(colnames(betaqc), pheno$MethID)==FALSE) stop('STOP: Sample names do not match'))

#make model matrix including variable of interest as factor variable & adjustment variables 
mod=model.matrix(~as.factor(smkPreg_binary), data=pheno)

#make null model containing only adjustment variables
mod0=model.matrix(~1, data=pheno)

#run with reasonable variance filter (5k)
svobj=sva(betaqc, mod, mod0, vfilter=var.filter)
sv.v5k=data.frame(svobj$sv); colnames(sv.v5k)=paste('sv', 1:ncol(sv.v5k), sep='')
sv.v5k=cbind(pheno%>%select(id, idnum, MethID), sv.v5k)
save(sv.v5k, file=paste0(output.dir, 'sv5k.Rds'))
print('Finished running SVA at Variance Filter of 5000')

#iterate over possible variance filters
possible.n=c(1000, 5000, 10000, 50000, 100000)

my_num.sv=function(x, y='be'){
  n=num.sv(betaqc, mod, vfilter=x, method=y)
  my.n=c(paste0('variance filter=', x, 'w/ ', y, 'method'), n)
  return(my.n)
}

plan(multicore, workers=5)
sv.n.be=possible.n %>% future_map(my_num.sv)
save(sv.n.be, file=paste0(output.dir, 'nSVbe.Rdata'))

plan(multicore, workers=5)
sv.n.leek=possible.n %>% future_map(my_num.sv, y='leek')
save(sv.n.leek, file=paste0(output.dir, 'nSVleek.Rdata'))

#run with no variance filter 
svobj=sva(betaqc, mod, mod0)
sv.vNone=data.frame(svobj$sv); colnames(sv.vNone)=paste('sv', 1:ncol(sv.vNone))
sv.vNone=cbind(pheno%>%select(id, idnum, MethID), sv.vNone)
save(sv.vNone, file=paste0(output.dir, 'svNone.Rds'))
print('Finished running SVA with no variance filter')


