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
library(lmerTest)
library(stringr)
library(sva)
##################################Directories & data
data.dir="/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/CreatedData"
outputdir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Output"
#read data 
file.path(data.dir, list.files(data.dir))[(str_detect(file.path(data.dir, list.files(data.dir)), 'NoSmk'))]%>%
  map(load, envir=.GlobalEnv)
#bind data
modeldata<-rbind(child_nosmoke, teen_nosmoke)
#label data 
set_label(modeldata)<-get_label(child_nosmoke)
#set as factors 
modeldata$PostnatalMaternalSmokingAny<-factor(modeldata$PostnatalMaternalSmokingAny, 
                                              levels=c("No maternal smoking at age 1 and age 5", 
                                                       "Maternal smoking at age 1 or age 5"))
modeldata$SmkAtVisitPastmonth<-factor(modeldata$SmkAtVisitPastmonth, 
                                      levels=c("No smoking", 
                                               "Less than pack a day", 
                                               "Pack or more a day", "Missing"))
modeldata$slide<-factor(modeldata$slide)
modeldata$Sample_Plate<-factor(modeldata$Sample_Plate)
#sv data - may change this depending on conversation w/ kelly
load(paste0(datadir, 'svNone.Rds'))
is_sv=function(x) str_replace_all(x, " ", "")
sv.vNone=sv.vNone%>%rename_with(is_sv, contains('sv'))
modeldata=left_join(modeldata, sv.vNone)

#################################Set variables for modeling
#outcomes & outcome_labels
y_vector<-c("globalmethylation", "pediatric", "grim", "anynewborn_center", "SSnewbornCT_center", "SSolder_center")
outcome_labels=c('Global methylation', 'Pediatric clock', 'GRIM clock',  
                 'Any smoking polymethylation score (newborns)', 
                 'Sustained smoking polymethylation score (newborns, cell-type controlled)', 
                 'Sustained smoking polymethylation score (older children)')
names(outcome_labels)=y_vector

#model variables (excepting PCs, add in individually for global & ancestry specific models)
base_model_vars<-"~smkPreg_binary+cm1bsex+cm1inpov+ChildAgeComposite+Epi+IC+Sample_Plate"
prenatal_model_vars<-paste0(base_model_vars, "+m1g2_YesNoPreg+m1g3_YesNoPreg")
secondhand_model_vars<-paste0(prenatal_model_vars, '+PostnatalMaternalSmokingAny+SmkAtVisitPastmonth')
#pc variables
global_pcs='+global_PC1+global_PC2'
local_pcs='+local_PC1+local_PC2'

#surrogate variable models
surrogate_model_vars<-paste0('~smkPreg_binary+', 
                             paste(colnames(modeldata)[str_detect(colnames(modeldata), 'sv')], collapse='+'))

predictors=c('base model'=base_model_vars, 
             'prenatal exposures model'=prenatal_model_vars, 
             'secondhand smoke exposure model'=secondhand_model_vars)

model_labels=c('base model', 'prenatal exposures model', 'secondhand smoke exposure model', 'surrogate variable models')



#stratification variables 
##############################
age_labels=c('Age 9', 'Age 15')
age_vector=modeldata$childteen%>%unique()
names(age_labels)=age_vector
ancestry_vector=modeldata$ancestry%>%unique()

#################################Source functions
model<-function(df, y, x){
  lm(formula(paste(y, x)), data=df)
}

model_gee<-function(df, y, x, corstr){
  geeglm(formula(paste(y, x)), data=df, family=gaussian(link='identity'), id=id, corstr = corstr)
}

model_lmm<-function(df, y, x){
  xnew=paste0(x, '+', '(1|id)')
  lmer(formula(paste(y, xnew)), data=df, REML=F)
}

my_clean_argnames=function(df, my_preds, childteen=T){
  names(model_labels)=my_preds
  df=df %>% mutate(predictors=recode_factor(predictors, !!!model_labels),
                outcome=recode_factor(outcome, !!!outcome_labels))
  if(childteen==T){df=df%>%mutate(age=recode_factor(childteen, !!!age_labels))} 
  return(df)
}

my_forest_plot<-function(df, vars){
  df %>% ggplot(aes(x=estimate, xmin=conf.low, xmax=conf.high, color=childteen))
}
  
#################################Cross sectional models
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
local_predictors=paste0(predictors, local_pcs)
args=list('childteen'=c('C', 'T'), 'outcome'=y_vector, 'predictors'=local_predictors, 'ancestry'=ancestry_vector)%>%cross_df()

#run models
local_models=modeldata%>%group_by(childteen, ancestry)%>%nest()%>%left_join(args)%>%
  mutate(model=pmap(list(data, outcome, local_predictors), model), 
         tidied=map(model, tidy, conf.int=T), 
         glanced=map(model, glance))

#clean argument names for plotting
local_models=my_clean_argnames(local_models, local_predictors)

#################################Longitudinal models
#########################Global
#set predictor variables & arguments
global_predictors=paste0(predictors, global_pcs, '+ChildAgeComposite')
args=list('outcome'=y_vector, 'predictors'=global_predictors)%>%cross_df()

#run models
longitudinal_global_models=args %>% 
    mutate(model=pmap(list(y=outcome, x=predictors), model_lmm, df=modeldata), 
           tidied=map(model, tidy, conf.int=T), 
           glanced=map(model, glance))

#clean argument names
longitudinal_global_models=my_clean_argnames(longitudinal_global_models, global_predictors, childteen = F)
#########################Local
#set predictor variables and arguments
local_predictors=paste0(predictors, local_pcs, '+ChildAgeComposite')
args=list('outcome'=y_vector, 'predictors'=local_predictors, 'ancestry'=ancestry_vector)%>%cross_df()

#run models
longitudinal_local_models=modeldata%>%group_by(ancestry)%>%nest()%>%left_join(args)%>%
  mutate(model=pmap(list(data, outcome, predictors), model_lmm), 
         tidied=map(model, tidy, conf.int=T), 
         glanced=map(model, glance))

#clean argument names
longitudinal_local_models=my_clean_argnames(longitudinal_local_models, local_predictors, childteen = F)



