#ewas
library(dplyr)
library(furrr) 
library(purrr)
library(lme4) 
library(lmerTest)
library(matrixStats)
library(limma)
library('variancePartition')
library('edgeR')
library('BiocParallel')

datadir="/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data"
outputdir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/CreatedData/EWAS"

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

#function to add annotation
add_annotations<-function(data, type='lm'){
  #add annotations
  if(type=='lme'){rownames(data)=data$CPG}
  Locations <- Locations[rownames(data), ]
  Islands.UCSC <- Islands.UCSC[rownames(data),]
  Other <- Other[rownames(data),]
  # add columns of interest from annotation to annotated version
  data.annotated <- data
  data.annotated$chr <- Locations$chr
  data.annotated$pos <- Locations$pos
  data.annotated$island <- Islands.UCSC$Relation_to_Island
  data.annotated$gene <- Other$UCSC_RefGene_Name
  data.annotated$regulatory_feature <- Other$Regulatory_Feature_Group
  return(data.annotated)
}
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)
data(Islands.UCSC)
data(Other)

##################cross-sectional EWAS###############
#########child & teen#################
#filter betaqc to child samples only 
child_nosmoke=modeldata%>%filter(childteen=='Age 9')
teen_nosmoke=modeldata%>%filter(childteen=='Age 15')
#models
mod1C=model.matrix(~factor(child_nosmoke$smkPreg_binary)+child_nosmoke$global_PC1+child_nosmoke$global_PC2+
                    factor(child_nosmoke$cm1bsex)+child_nosmoke$cm1inpov+child_nosmoke$ChildAgeComposite+
                    child_nosmoke$Leukocytes_saliva)

mod1T=model.matrix(~factor(teen_nosmoke$smkPreg_binary)+teen_nosmoke$global_PC1+teen_nosmoke$global_PC2+
                  factor(teen_nosmoke$cm1bsex)+teen_nosmoke$cm1inpov+teen_nosmoke$ChildAgeComposite+
                    teen_nosmoke$Leukocytes_saliva)

getCrossEwas=function(df, model, prefix){
  beta_new=betaqc[, colnames(betaqc) %in% df$MethID]
  df = df[match(colnames(beta_new), df$MethID), ]
  if(dim(beta_new)[2]!=dim(df)[1]){stop('sample #s dont match')}
  if(identical(colnames(beta_new), df$MethID)==FALSE){stop('sample names dont match')}
  out=lmFit(beta_new, model)
  out=eBayes(out)
  saveRDS(out, file=file.path(outputdir, paste0(prefix, '_Model.RDS')))
  ss.hits<-topTable(out, coef=2, number=nrow(beta_new), confint=T)
  ss.hits<-add_annotations(ss.hits)
  saveRDS(ss.hits, file=file.path(outputdir, paste0(prefix, '_TopHits.RDS')))
}

pmap(list('data'=list(child_nosmoke, teen_nosmoke), 'model'=list(mod1C, mod1T), 'prefix'=list('child1', 'teen1')), ~getCrossEwas(..1, ..2, ..3))


####################Longitudinal EWAS################

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
                                                          SmkAtVisitPastmonth+(1|id), data=modeldata, REML=F)))}
  my_tidy=broom.mixed::tidy(my_model, conf.int=T)
  return(my_tidy)
}


# do longitudinal EWAS in base model, all children then compare TOP SNPS across ancestries and times, w/ vs w/out controls
#for additional covariates

plan(multicore, workers=3)
mixed_models_1=as.data.frame(t(betaqc))%>%future_map_dfr(model_lmm, .id='CPG', model='1')%>%mutate(modeltype='Base model')
saveRDS(mixed_models_1, file = file.path(outputdir, 'lme_ewas_base.RDS'))
lme_smk = mixed_models_1 %>% filter(term=='smkPreg_binaryYes')%>%mutate(p.bh=p.adjust(p.value, 'hochberg'))
lme_smk=add_annotations(lme_smk, type='lme')
saveRDS(lme_smk, file=file.path(outputdir, 'lme_smoking.RDS'))
#mixed_models_2=as.data.frame(t(betaqc))%>%future_map_dfr(model_lmm, .id='CPG', model='2')%>%mutate(modeltype='Prenatal exposure model')
#mixed_models_3=as.data.frame(t(betaqc))%>%future_map_dfr(model_lmm, .id='CPG', model='3')%>%mutate(modeltype='Postnatal smoke exposure model')
#mixed_models_lme=rbind(mixed_models_1, mixed_models_2, mixed_models_3)
#saveRDS(mixed_models_lme, file=file.path(datadir, 'CreatedData', 'lme_ewas.RDS'))
#rm(list=c('mixed_models_1', 'mixed_models_2', 'mixed_models_3', 'mixed_models_lme'))

#gap hunter
allgaps=gaphunter(betaqc)
saveRDS(allgaps, file=file.path(outputdir, 'gaphunter_gaps.RDS'))

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
saveRDS(fitmm, file=file.path(outputdir, 'dream_ewas.RDS'))




