#ewas
library(dplyr)
library(furrr) 
library(purrr)
library(lme4) 
library(lmerTest)
library(matrixStats)
library('variancePartition')
library('edgeR')
library('BiocParallel')

datadir="/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data"
outputdir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Output"

#variance trim of beta matrix
betaqc<-readRDS(file=file.path(datadir, 'OGData', "betaqc.rds"))
varCut<-0.0002
betaqc<-betaqc[rowVars(betaqc)>varCut,]
dim(betaqc)

#filter to only samples in final analytic subsample
load(file=file.path(datadir, 'CreatedData', 'completeCaseSVA.Rdata'))
betaqc<-betaqc[, colnames(betaqc) %in% modeldata$MethID]

#match rows to columns
modeldata=modeldata[match(colnames(betaqc), modeldata$MethID),]
#rownames of modeldata 
rownames(modeldata)=modeldata$MethID

#functions for  EWAS 
model_lmm<-function(y, local=F, model='1'){
  #create new variable for PCS
  if(local==F){modeldata=modeldata %>% mutate(PC1=global_PC1, PC2=global_PC2)}
  if(local==T){modeldata=modeldata %>% mutate(PC1=local_PC1, PC2=local_PC2)}
  if(model=='1'){my_model=suppressMessages(suppressWarnings(lmer(y~smkPreg_binary+PC1+PC2+
                                                    cm1bsex+cm1inpov+ChildAgeComposite+Leukocytes_saliva+(1|id), data=modeldata, REML=F)))}
  if(model=='2'){my_model=suppressMessages(suppressWarnings(lmer(y~smkPreg_binary+PC1+PC2+
                                                                   cm1bsex+cm1inpov+ChildAgeComposite+Leukocytes_saliva+
                                                                   m1g2_YesNoPreg+m1g3_YesNoPreg+(1|id), data=modeldata, REML=F)))}
  if(model=='3'){my_model=suppressMessages(suppressWarnings(lmer(y~smkPreg_binary+PC1+PC2+
                                                          cm1bsex+cm1inpov+ChildAgeComposite+Leukocytes_saliva+
                                                          m1g2_YesNoPreg+m1g3_YesNoPreg+PostnatalMaternalSmokingAny+
                                                          SmkAtVisitPastmonth1|id), data=modeldata, REML=F))}
  my_tidy=broom.mixed::tidy(my_model)
  return(my_tidy)
}




# do longitduinal EWAS in all children then compare TOP SNPS across ancestries and times, 
#only need to run models three times, once for each model (covariates)

plan(multicore, workers=2)
mixed_models_1=as.data.frame(t(betaqc))%>%future_map_dfr(model_lmm, .id='CPG', model='1')%>%mutate(modeltype='Base model')
saveRDS(mixed_models_1, file = file.path(datadir, 'CreatedData', 'lme_ewas_base.RDS'))

mixed_models_2=as.data.frame(t(betaqc))%>%future_map_dfr(model_lmm, .id='CPG', model='2')%>%mutate(modeltype='Prenatal exposure model')
mixed_models_3=as.data.frame(t(betaqc))%>%future_map_dfr(model_lmm, .id='CPG', model='3')%>%mutate(modeltype='Postnatal smoke exposure model')

mixed_models_lme=rbind(mixed_models_1, mixed_models_2, mixed_models_3)

saveRDS(mixed_models_lme, file=file.path(datadir, 'CreatedData', 'lme_ewas.RDS'))

rm(list=c('mixed_models_1', 'mixed_models_2', 'mixed_models_3', 'mixed_models_lme'))

#dream method #
#see  https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html
param = SnowParam(2, "SOCK")
register(param)

# The variable to be tested must be a fixed effect
form <-  ~smkPreg_binary+global_PC1+global_PC2+cm1bsex+cm1inpov+ChildAgeComposite+Leukocytes_saliva+(1|id)

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(betaqc, form, modeldata )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, modeldata )
saveRDS(fitmm, file=file.path(datadir, 'CreatedData', 'dream_ewas.RDS'))




