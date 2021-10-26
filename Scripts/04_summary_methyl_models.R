###################################Libraries
library(tidyr)
library(stringr)
library(dplyr)
library(compareGroups)
library(kableExtra)
library(haven)
library(sjlabelled)
library(matchmaker)
library(EpiDISH)
library("PerformanceAnalytics")
library(DescTools)
library(purrr)
library(broom)
library(broom.mixed)
library(ggplot2)
library(multcomp)
library(ggpubr)
library(cowplot)
library(geepack)
library(lme4)
library(pROC)
library(lmerTest)
library(stringr)
library(sva)
#############
#investigate why clocks have 1 fewer observation in child models
##################################Directories & data
datadir="/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/CreatedData"
outputdir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Output"
#read data 
load(paste0(datadir, '/completeCasemethyl.Rdata'))
modeldata=completecase%>%droplevels()%>%copy_labels(completecase)

#set as factors 
modeldata=modeldata %>% mutate(childteen=factor(childteen, levels=c('Age 9', 'Age 15')))
modeldata$PostnatalMaternalSmokingAny<-factor(modeldata$PostnatalMaternalSmokingAny, 
                                              levels=c("No maternal smoking at age 1 and 5", 
                                                       "Maternal smoking at age 1 or age 5"))
modeldata$SmkAtVisitPastmonth<-factor(modeldata$SmkAtVisitPastmonth, 
                                      levels=c("No smoking", 
                                               "Less than pack a day", 
                                               "Pack or more a day", "Missing"))
#modeldata$slide<-factor(modeldata$slide)
modeldata$Sample_Plate<-factor(modeldata$Sample_Plate)

#sv data - no variance filter, whole subset
sv.vNone=readRDS(paste0(datadir, '/svNone.Rds'))
is_sv=function(x) str_replace_all(x, " ", "")
sv.vNone=sv.vNone%>%rename_with(is_sv, contains('sv'))
modeldata=left_join(modeldata, sv.vNone)%>%copy_labels(modeldata)
set_label(modeldata$Sample_Plate)='Batch'
set_label(modeldata$PostnatalMaternalSmokingAny)='Any postnatal maternal smoking (Age 1 & 5)'
set_label(modeldata$SmkAtVisitPastmonth)='Primary care giver smoke past month dose'
save(modeldata, file=file.path(datadir, file='completeCaseSVA.Rdata'))

#################################Set variables for modeling
#outcomes & outcome_labels

y_vector<-c("globalmethylation", "pediatric", "anynewborn_center", "SSnewbornCT_center", "SSolder_center", "cg05575921")
outcome_labels=c('Global methylation', 'Pediatric clock',  
                 'Any smoking polymethylation score (newborns)', 
                 'Sustained smoking polymethylation score (newborns, cell-type controlled)', 
                 'Sustained smoking polymethylation score (older children)', 
                 'AHHR: cg05575921')
names(outcome_labels)=y_vector


#model variables (excepting PCs, add in individually for global & ancestry specific models)
base_model_vars<-"~smkPreg_binary+cm1bsex+cm1inpov+Leukocytes_saliva+Sample_Plate+ChildAgeComposite"
prenatal_model_vars<-paste0(base_model_vars, "+m1g2_YesNoPreg+m1g3_YesNoPreg")
secondhand_model_vars<-paste0(prenatal_model_vars, '+PostnatalMaternalSmokingAny+SmkAtVisitPastmonth')
interaction_model_vars<-paste0(secondhand_model_vars, '+ChildAgeComposite:smkPreg_binary')
#pc variables
global_pcs='+global_PC1+global_PC2'
local_pcs='+local_PC1+local_PC2'

#surrogate variable models
surrogate_model_vars<-paste0('~smkPreg_binary+', 
                             paste(colnames(modeldata)[str_detect(colnames(modeldata), 'sv')], collapse='+'))

predictors=c('base model'=base_model_vars, 
             'prenatal exposures model'=prenatal_model_vars, 
             'secondhand smoke exposure model'=secondhand_model_vars)

predictors_interaction=c(predictors, 'interaction model'=interaction_model_vars)

model_labels=c('base model', 'prenatal exposures model', 'secondhand smoke exposure model', 'surrogate variable models')

model_labels_interaction=c('base model', 'prenatal exposures model', 'secondhand smoke exposure model', 'interaction model (age*smoke)', 'surrogate variable models')


#stratification variables 
##############################
age_labels=c('Age 9', 'Age 15')
age_vector=modeldata$childteen%>%unique()
names(age_labels)=age_vector
ancestry_vector=modeldata$ancestry%>%unique()

#################################Source functions#######################################
model<-function(df, y, x){
  lm(formula(paste(y, x)), data=df)
}

model_reversed<-function(df, y, x){
  mylogit <- glm(paste(y, x), data = df, family = "binomial")
}

model_reversed_safe=safely(model_reversed)

model_gee<-function(df, y, x, corstr){
  geeglm(formula(paste(y, x)), data=df, family=gaussian(link='identity'), id=id, corstr = corstr)
}

model_lmm<-function(df, y, x){
  xnew=paste0(x, '+', '(1|id)')
  lmer(formula(paste(y, xnew)), data=df, REML=F)
}

my_clean_argnames=function(df, my_preds, childteen=T, model_labs=model_labels){
  names(model_labs)=my_preds
  df=df %>% mutate(predictors=recode_factor(predictors, !!!model_labs),
                outcome=recode_factor(outcome, !!!outcome_labels))
  if(childteen==T){df=df%>%mutate(age=recode_factor(childteen, !!!age_labels))} 
  return(df)
}
  
#################################Cross sectional models###########################
#########################Global
#set predictor variables & arguments
global_predictors=c(paste0(predictors, global_pcs), surrogate_model_vars)
args=list('childteen'=age_vector, 'outcome'=y_vector, 'predictors'=global_predictors)%>%cross_df()

#run models
global_models=modeldata%>%group_by(childteen)%>%nest()%>%left_join(args)%>%
  mutate(model=pmap(list(data, outcome, predictors), model), 
         tidied=map(model, tidy, conf.int=T), 
         glanced=map(model, glance))

#clean argument names for plotting
global_models=my_clean_argnames(global_models, global_predictors)

#########################Local
#arguments
local_predictors=c(paste0(predictors, local_pcs), surrogate_model_vars)
args=list('childteen'=c('Age 9', 'Age 15'), 'outcome'=y_vector, 'predictors'=local_predictors, 'ancestry'=ancestry_vector)%>%cross_df()

#run models
local_models=modeldata%>%group_by(childteen, ancestry)%>%nest()%>%left_join(args)%>%
  mutate(model=pmap(list(data, outcome, predictors), model), 
         tidied=map(model, tidy, conf.int=T), 
         glanced=map(model, glance))

#clean argument names for plotting
local_models=my_clean_argnames(local_models, local_predictors)

#################################Longitudinal models##################
#########################Global
#set predictor variables & arguments
#include a model with an interaction term b/t age and prenatal smoke exposure
#to determine if there is a age-slope difference b/t exposed & not exposed
global_predictors=c(paste0(predictors_interaction, global_pcs), surrogate_model_vars)
args=list('outcome'=y_vector, 'predictors'=global_predictors)%>%cross_df()

#run models
longitudinal_global_models=args %>% 
    mutate(model=pmap(list(y=outcome, x=predictors), model_lmm, df=modeldata), 
           tidied=map(model, tidy, conf.int=T), 
           glanced=map(model, glance))

#clean argument names
longitudinal_global_models=my_clean_argnames(longitudinal_global_models, global_predictors, childteen = F, model_labs = model_labels_interaction)
#########################Local
#set predictor variables and arguments
local_predictors=c(paste0(predictors_interaction, local_pcs), surrogate_model_vars)
args=list('outcome'=y_vector, 'predictors'=local_predictors, 'ancestry'=ancestry_vector)%>%cross_df()

#run models
longitudinal_local_models=modeldata%>%group_by(ancestry)%>%nest()%>%left_join(args)%>%
  mutate(model=pmap(list(data, outcome, predictors), model_lmm), 
         tidied=map(model, tidy, conf.int=T), 
         glanced=map(model, glance))

#clean argument names
longitudinal_local_models=my_clean_argnames(longitudinal_local_models, local_predictors, childteen = F, model_labs = model_labels_interaction)


################################Reverse models for roc curve########################
#add a numeric for yes/no smoking for logistic models
modeldata=modeldata%>%mutate(smkPreg_binaryN=case_when(smkPreg_binary=='Yes'~1, smkPreg_binary=='No'~0))

#global predictors, cross sectional
roc_predictors=c('', gsub('~smkPreg_binary', '', paste0(base_model_vars, global_pcs)))
roc_outcomes=c('', y_vector, 'cg05549655', 'cg22132788', 'SSnewborn_center')
args=list('childteen'=age_vector, 'methylation'=roc_outcomes, 'covariates'=roc_predictors)%>%cross_df()%>%mutate(predictors=paste0(methylation, covariates))%>%
  filter(predictors!='')


methyl_names=data.frame(name_methyl=c('No methylation', 'Global methylation', 'Pediatric clock', 
                                      'Any smoking (newborn)', 'Sustained smoking w/ cell type correction (newborn)', 
                                      'Sustained smoking (older children)',
                                      'AHRR: ch05575921', "CYP1A1: cg05549655", "MYO1G: cg22132788", 
                                      'Sustained smoking (newborn)'), 
                        methylation=roc_outcomes)

global_roc=modeldata%>%
  group_by(childteen)%>%nest()%>%left_join(args)%>%
  mutate(model=pmap(list(data, 'smkPreg_binaryN~', predictors), model_reversed), 
         tidied=map(model, tidy, conf.int=T), 
         augmented=map(model, augment, type.predict = "response"), 
         roc=map(augmented, ~roc(.x$smkPreg_binaryN, .x$.fitted)), 
         test=map(roc, ~data.frame('sens'=.x$sensitivities, spec=.x$specificities)))%>%
  left_join(methyl_names)%>%
  mutate(modeltype=case_when(covariates==''~name_methyl, 
                             covariates!=''~paste0('Base model+', name_methyl)))%>%
  mutate(auc=map_dbl(roc, ~as.numeric(gsub('Area under the curve', '', .x$auc))))

#local predictors, crosssectional 
roc_predictors=c('', gsub('~smkPreg_binary', '', paste0(base_model_vars, local_pcs)))
args=list('childteen'=age_vector, 'methylation'=c('', y_vector), 'covariates'=roc_predictors)%>%cross_df()%>%mutate(predictors=paste0(methylation, covariates))%>%filter(predictors!='')


local_roc=modeldata%>%group_by(childteen, ancestry)%>%nest()%>%left_join(args)%>%
  mutate(model=pmap(list(data, 'smkPreg_binaryN~', predictors), model_reversed), 
         tidied=map(model, tidy, conf.int=T), 
         augmented=map(model, augment, type.predict = "response"), 
         roc=map(augmented, ~roc(.x$smkPreg_binaryN, .x$.fitted)), 
         test=map(roc, ~data.frame('sens'=.x$sensitivities, spec=.x$specificities)))%>%
  left_join(methyl_names)%>%
  mutate(modeltype=case_when(covariates==''~name_methyl, 
                             covariates!=''~paste0('Base model+', name_methyl)))%>%
  mutate(auc=map_dbl(roc, ~as.numeric(gsub('Area under the curve', '', .x$auc))))



#global predictors, longitudinal 
modeldata_wide = modeldata%>%
  mutate(childteen=gsub(' ', '_', childteen))%>%
  dplyr::select(any_of(c('idnum', 'ancestry', 'childteen', 'global_PC1', 
                         'global_PC2', 'local_PC1', 'local_PC2', 
                         'smkPreg_binaryN',  y_vector, 
                         str_split(secondhand_model_vars,  '\\+')[[1]][-1])))%>%
  pivot_wider(names_from = childteen, values_from=any_of(c(y_vector, 'Leukocytes_saliva', 'Epithelial.cells_saliva', 'Sample_Plate', 'SmkAtVisitPastmonth')))

#arguments
long_predictors=c('', '+cm1bsex+cm1inpov+global_PC1+global_PC2+Leukocytes_saliva_Age_9+Sample_Plate_Age_9+Leukocytes_saliva_Age_15+Sample_Plate_Age_15')
args=list('methylation'=paste0(y_vector, '_Age_9+', y_vector, '_Age_15') , 'covariates'=long_predictors)%>%cross_df()%>%mutate(predictors=paste0(methylation, covariates))

long_global_roc=args %>% 
  mutate(model=map(predictors, model_reversed, df=modeldata_wide, y='smkPreg_binaryN~'), 
         tidied=map(model, tidy, conf.int=T), 
         augmented=map(model, augment, type.predict = "response"), 
         roc=map(augmented, ~roc(.x$smkPreg_binaryN, .x$.fitted)), 
         test=map(roc, ~data.frame('sens'=.x$sensitivities, spec=.x$specificities)))

save(list=c('long_global_roc', 'local_roc', 'global_roc', 'longitudinal_local_models', 'longitudinal_global_models', 'local_models', 'global_models'), file=paste0(datadir, '/methylation_summarymodels.Rdata'))
