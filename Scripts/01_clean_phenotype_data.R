
#####################libraries, directories, data###########################

#libraries 
library(tidyr)
library(stringr)
library(dplyr)
library(purrr)
library(tidyverse)
library(dplyr)
library(sjlabelled)
library(here)

#directory based coding
datadir<-gsub('Code', 'Data/', here::here())
outputdir<-gsub('Code', 'Output/', here::here())
codedir<-paste0(here::here(), '/')

#read in data 
pheno<-readRDS(paste0(datadir, "OGData/fullpheno.rds"))
FFmeta<-read.csv(paste0(datadir, "FFMetadata_v07.csv"))

############################relabel data###################################

#remove offending characters from variable labels in phenotype data
attributes(pheno)$var.labels<-gsub("\xad|\x92", "'",  attributes(pheno)$var.labels)
#adding in a variable/value label for id
attributes(pheno)$var.labels[which(colnames(pheno)=='id')]='ID'
attributes(pheno)$label.table[which(colnames(pheno)=='id')]=NA
#label data using sj labelled - variable names
pheno=set_label(pheno, attributes(pheno)$var.labels)
#label data using sj labelled - value labels 
pheno=set_labels(pheno, labels=attributes(pheno)$label.table)
#identify categorical variables
FF_factors=FFmeta%>%
  filter(type%in% c('Binary', 'Ordered Categorical', 'Unordered Categorical'))%>%
  pull(old_name)
FF_factors=c(FF_factors, 'm1city')
#recode binary, categorical variables so they have label names as values (nicer for table and plotting)
pheno=pheno %>% dplyr::mutate(across(any_of(FF_factors), ~sjlabelled::as_label(., drop.levels=TRUE)))

#dplyr::select variables
idvars<-c("id", "ff_id")
demovars<-c("m1city", "cm1age", "cm1bsex", "cm1ethrace", "cm1edu", "cm1hhinc", "cm1inpov", "m1h3", "m1h3a", "m1h3b", "m1i4", "m1i4a", "m1i4b", "cf1age", "cf1ethrace", "f1h3", "f1h3a", "f1h3b", "ch5agem","ch6yagem", "c6yagey", "cm5b_age", "cm5b_ageyrs", "cp6yagem", "cp6yagey", "hv5_agem", "hv6_yagem", "hv6_yagey")
smokingvars<-c("m1g4", "f1g4", "f2j5", "f2j5a", "f2j7", "f2j7a", "f3j31", "f3j32", "f4j18", "f4j19", "f5g17", "f5g18", "k5f1k", "k5f1l", "k6d40", "k6d41", "k6d41", "k6d42", "k6d43", "k6d45", "k6d46", "k6d47", "m2j5", "m2j5a", "m2j7", "m2j7a", "m4j18", "m4j19", "m5g17", "m5g18", "n5f17", "n5f18", "p3a22", "p3a23", "p3a23a", "p3a23b", "p3a24", "p4a22", "p4a23", "p5h15", "p5h15b", "p5q3cr", "p6h74", "p6h75", "p6h76", "p6h77", "p6h78")
prenataldrugusevar<-c("m1g2", "m1g3", "m1g5", "m1g6")
FF_labeled=pheno %>%dplyr::select(any_of(c(idvars, demovars, smokingvars, prenataldrugusevar)))
#apply NA codes
my_na_codes<-c(-10:-1, "-1 Refuse", "-2 Don't know", "-3 Missing", "-9 Not in wave", "-7 N/A")
na_codes <- function(x, ...) {
  x[x %in% c(...)] <- NA
  if(is.factor(x)==TRUE){x<-droplevels(x)}
  x
}
FF_labeled=FF_labeled%>%mutate(across(everything(), ~na_codes(., my_na_codes)))%>%copy_labels(FF_labeled)

############################pdqc data######################################

#As of February 11 switch to February 2022 Data Freeze, which includes some fixes to sex mismatches due to incorrect
#coding at birth interview; only clean samples which are not technical replicates are included in this file according to Jonah 
pdqc_all<-read.csv(paste0(datadir, 'OGData/Feb2022DataFreeze/', "pd_analytic_freeze1.csv"))
#Flag samples w/ incorrect sex at birth variable 
pdqc_all = pdqc_all %>% 
          left_join(FF_labeled %>% 
                      mutate(sex_old=stringr::word(cm1bsex), idnum=as.integer(id))%>% 
                      dplyr::select(idnum, sex_old))%>%
          mutate(sex_flag=ifelse(cm1bsex!=sex_old, 'recode_sex', 'correct_sex'), 
                 new_sex=cm1bsex)
sex_recodes=pdqc_all %>% filter(sex_flag=='recode_sex')%>%pull(idnum)

#add batch, slide and plate data
batchData<-read.csv(file.path(datadir, 'Methylation_450K_array_batch_information.csv'))
pdqc_clean <- pdqc_all %>% left_join(batchData)

############################create new variables###################################

#create following variables in FF_labeled:
#Methyldata - indicator for in analysis subset (have methylation 450k data) or no
#smkPregbinary Yes/no - any smoking during pregnancy
#m1g2_YesNoPreg and m1g3_YesNoPreg: yes/no alcohol and drugs (respectively)
#PostnatalMaternal smoking any: if mom smoke at either 1 or 5 - yes, if missing at either - missing, if no at both - no
#Postnatal maternal smoking does: if mom>1p/d at either 1 or 5 , if mom <1p/d at either 1 or 5, if no smoking at both, if missing + no smoking --> missing
FF_labeled=FF_labeled %>% 
  mutate(Methyldata=ifelse(id %in% pdqc_clean$idnum, 'Analysis subset', 'Not in analysis subset'), 
         smkPreg_binary=case_when(m1g4 %in% c("1 2+pk/d","2 1<pk<2","3 <1pk/d" ) ~ 'Yes', 
                                  m1g4 == '4 None'~'No'))%>%
  mutate(across(c('m1g2', 'm1g3'), 
                ~case_when(.x=='5 Never'~'No', is.na(.x)~ NA_character_, TRUE ~ 'Yes'), 
                .names="{.col}_YesNoPreg"))%>%
  mutate(PostnatalMaternalSmokingAny=case_when(m2j5 == '2 No' & m4j18=='2 No' ~ 'No maternal smoking at age 1 and 5',
                                               m2j5=='1 Yes' | m4j18=='1 Yes'~ 'Maternal smoking at age 1 or age 5',
                                               (is.na(m2j5) | m2j5=='2 No') & (is.na(m4j18) | m4j18=='2 No') ~ 'Missing'), 
         PostnatalMaternalSmokingDose=case_when(m2j5a %in% c("2 1pk/d", "3 1.5pk/d", "4 2pk/d", "5 >2pk/d") |
                                                m4j19 %in% c("2 About pack/day", "3 Pack and half/day", "4 2 packs/day", "5 More 2 packs/day")~
                                                  'Pack/day or greater when child aged 1 or 5', 
                                                m2j5a %in% c("1 <1/2 pk/d") | 
                                                m4j19 %in% c("1 Less half pack/day")~
                                                  "Less than pack/d when child aged 1 or 5", 
                                                m2j5a=='-6 Skip' & m4j19== '-6 Skip'~'No smoking when child age 1 or 5' ), 
         cm1bsex=factor(case_when(id %in% sex_recodes & cm1bsex=='1 Boy'~'2 Girl', 
                                  id %in% sex_recodes & cm1bsex=='2 Girl'~'1 Boy', 
                                  TRUE ~ as.character(cm1bsex))))%>%
  copy_labels(., FF_labeled)

#ancestry
#read in local and global pc data
id_list=list.files(path=paste0(datadir, "OGData/pcs"), pattern='*idnums*', full.names=TRUE)
ancestry_ids=id_list%>%map_dfr(read.table, .id='ancestry')


pc_files=list.files(path=paste0(datadir, "OGData/pcs"), pattern='FFCWDNA_.*.csv', full.names=TRUE)
names(pc_files)=gsub('.csv', '', list.files(path=paste0(datadir, "OGData/pcs"), pattern='FFCWDNA_.*.csv'))
names(pc_files)=paste0(str_to_title(str_match(names(pc_files), '.*_(.*)_n.*')[,2]), ' ancestry')
localpc=pc_files%>%map_dfr(read.csv, .id='ancestry')%>%rename_with(~str_c('local_', .), contains('PC'))
localpc$ancestry=case_when(localpc$ancestry=='Hispanic ancestry'~'Admixed ancestry - Latin heritage', 
                           TRUE ~ localpc$ancestry)
allPC<-read.table(paste0(datadir, "OGData/pcs/pcs_n50.eigenvec"), col.names=c("idnum", 'idnum2', paste("global_PC", 1:50, sep='')))

old_pc=list.files(path=paste0(datadir, "OGData/pcs"), pattern='.*[n|c].csv', full.names=TRUE)
names(old_pc)=gsub('.csv', '', list.files(path=paste0(datadir, "OGData/pcs"), pattern='FFCWDNA_.*.csv'))
names(old_pc)=paste0(str_to_title(str_match(names(old_pc), '.*_(.*)_n.*')[,2]), ' ancestry')
local_old=old_pc%>%map_dfr(read.csv, .id='ancestry')%>%rename_with(~str_c('local_', .), contains('PC'))
checkPC=left_join(localpc, local_old%>%mutate(old_ancestry=ancestry)%>%dplyr::select(idnum, old_ancestry))

#individuals in each ancestry specific pc group can be labeled as that ancestry for categorical ancestry variable 
FF_labeled=FF_labeled %>% left_join(localpc%>%mutate(id=as.character(idnum))%>%dplyr::select(id, ancestry, contains('local_PC')))%>%
  left_join(allPC%>%mutate(id=as.character(idnum))%>%dplyr::select(id, contains('global_PC')))%>%
  mutate(ancestry=case_when(is.na(ancestry) & is.na(global_PC1) & Methyldata=='Analysis subset' ~'Missing PC data',
                            is.na(ancestry) & Methyldata!='Analysis subset'~'Not in analysis subset',
                            !is.na(ancestry)~ancestry))%>%copy_labels(., FF_labeled)

MissingPCdata<-FF_labeled %>% filter(ancestry=="Missing PC data")%>%dplyr::select(id)
save(MissingPCdata, file = paste0(datadir, "CreatedData/MissingPC.Rdata"))
#save FFlabeled
save(FF_labeled, file = paste0(datadir, "CreatedData/FFlabeled.Rdata"))

#######################filter to analysis subset and join with pdqc###################################
myFF=FF_labeled%>%filter(Methyldata=='Analysis subset')%>%
  left_join(pdqc_clean%>% 
              dplyr::select(MethID, childteen, idnum, Sample_Plate, Slide, Array, new_sex, sex_flag)%>%
              mutate(id=as.character(idnum)))%>%
  copy_labels(., FF_labeled)

#check sex recodes
#myFF%>% filter(sex_flag=='recode_sex')%>% xtabs(~cm1bsex+new_sex, data=.)
########################make visit specific variables###################################
#child age
#mom/pgc smoking amount past month
myFF<-myFF %>% 
  mutate(ChildAgeComposite=case_when(childteen=='C' & !is.na(hv5_agem)~hv5_agem, 
                                     childteen=='C' & is.na(hv5_agem)~cm5b_age, 
                                     childteen=='T' & !is.na(hv6_yagem)~hv6_yagem, 
                                     childteen=='T' & is.na(hv6_yagem)~cp6yagem))%>%
  mutate(SmkAtVisitPastmonth=case_when(
    (childteen=='C' & ( as.numeric(substr(m5g18, 1, 1))>1 | as.numeric(substr(n5f18, 1, 1))>1 )) |
      (childteen=='T' & as.numeric(substr(p6h77, 1, 1))>2) ~ "Pack or more a day", 
    (childteen=='C' & ( as.numeric(substr(m5g18, 1, 1))==1 | as.numeric(substr(n5f18, 1, 1))==1) ) | 
      (childteen=='T' & (as.numeric(substr(p6h77, 1, 1))%in% c(1, 2)))  ~ "Less than pack a day", 
    (childteen=='C' & ( is.na(m5g18) & (is.na(n5f18)) )) |
      (childteen=='T' & is.na(p6h77))~ 'Missing', 
    childteen=='C' & ((m5g17=='2 no' & n5f17 %in% c(NA_character_, '2 no')) | (n5f17=='2 no' & m5g17%in%c('2 no', NA_character_)))|
    childteen=='T' & p6h77=='-6 Skip'~ 'No smoking'))
#throws NA warning because of as.numeric in these statements - sow rite a check to make sure there aren't NAs here 
if(myFF%>%filter(is.na(SmkAtVisitPastmonth))%>%nrow()!=0){stop('NA in SmkAtVisitPastmonth')}
myFF=myFF %>% mutate(city_binary=case_when(m1city %in% c('4 Detroit', '15 Chicago', '17 Toledo')~'Detroit, Chicago or Toledo', 
                                           !(m1city %in% c('4 Detroit', '15 Chicago', '17 Toledo'))~'Not Detroit, Chicago or Toledo'))

################################################################################add labels
set_label(myFF)=c(get_label(myFF)[1:70], 'Has 450K Illumina chip data', 'Maternal prenatal smoking', 
                  'Maternal prenatal alcohol use', 'Maternal prenatal any drug use', 
                  'Any postnatal maternal smoking when child age 1 or 5', 'Postnatal maternal smoking dose when child age 1 or 5', 
                  'Ancestry categorization from child principal components of genetic data', 
                  paste0('Within ancestry strata principal component', 1:20), 
                  paste0('Within all samples principal component', 1:50), 
                  'Methylation data ID', 'Visit', 'ID', 'Plate', 'Slide', 'Array',
                  'sex_fromjonah', 'sex recode flag',
                  'Child age at visit', 'Maternal/primary care giver smoking in month prior to visit (pks/day)', 
                  'Oversampled cities')

#######################complete case dataset###################################


basemodelvar<-c("smkPreg_binary", "ancestry", "cm1bsex", "cm1inpov", "ChildAgeComposite")
basemodeldata<-myFF %>% filter_at(all_of(basemodelvar), all_vars(!(. %in%c("Missing", "Missing PC data"))& !is.na(.)))

prenatalexposurevar<-c(basemodelvar, "m1g2_YesNoPreg", "m1g3_YesNoPreg")
prenatalexdata<-basemodeldata %>%filter_at(all_of(prenatalexposurevar), all_vars(!(. %in% c("Missing")) & !is.na(.)))

secondhandsmokevars<-c(prenatalexposurevar, "PostnatalMaternalSmokingAny", "SmkAtVisitPastmonth")
secondhandsmkdata<-prenatalexdata%>%filter_at(all_of(secondhandsmokevars), all_vars(!(. %in% c("Missing")) & !is.na(.)))

child_smoke<-secondhandsmkdata %>% filter(k5f1l!= "2 no" & childteen=='C')
#only children who report never cigarette use at age 9 OR were missing never cigarette at age 9 but reported no smoking at age 15
child_nosmoke<-secondhandsmkdata %>% filter((k5f1l=="2 no" | (is.na(k5f1l) & k6d40=='2 No')) & childteen=='C')
teen_smoke<-secondhandsmkdata %>% filter(k6d40!="2 No" & childteen=='T')
#only children who report never smoking at age 15 AND were either missing or reported not smoking at age 9
teen_nosmoke<-secondhandsmkdata %>% filter(k6d40=="2 No" & (k5f1l=='2 no' |  is.na(k5f1l)) & childteen=='T')
completecase<-rbind(child_nosmoke, teen_nosmoke)%>%copy_labels(basemodeldata)%>%copy_labels(myFF)

fathersmkdata<-completecase %>% filter(f1g4!="Missing" & !is.na(f1g4))




#####################################################################################Save all data 
save.image(file=paste0(datadir, '/CreatedData/allPhenoData.Rdata'))
#####################################################################################Save my specific data for modeling
save(completecase, file=paste0(datadir, '/CreatedData/completeCasepheno.Rdata'))

