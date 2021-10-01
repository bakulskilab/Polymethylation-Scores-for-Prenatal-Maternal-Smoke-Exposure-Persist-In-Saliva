
#####################libraries, directories, data###########################

#libraries 
library(tidyr)
library(stringr)
library(dplyr)
library(purrr)
library(tidyselect)
library(sjlabelled)

#directory based coding
datadir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/"
john_datadir<-"/nfs/turbo/bakulski1/People/johndou/Fragile_Families/ProcessedFiles/jd"
outputdir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Output/"
codedir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Code/"

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
#recode binary, categorical variables so they have label names as values (nicer for table and plotting)
pheno=pheno %>% dplyr::mutate(across(any_of(FF_factors), ~sjlabelled::as_label(., drop.levels=TRUE)))

#select variables
idvars<-c("id", "ff_id")
demovars<-c("cm1age", "cm1bsex", "cm1ethrace", "cm1edu", "cm1hhinc", "cm1inpov", "m1h3", "m1h3a", "m1h3b", "m1i4", "m1i4a", "m1i4b", "cf1age", "cf1ethrace", "f1h3", "f1h3a", "f1h3b", "ch5agem","ch6yagem", "c6yagey", "cm5b_age", "cm5b_ageyrs", "cp6yagem", "cp6yagey", "hv5_agem", "hv6_yagem", "hv6_yagey")
smokingvars<-c("m1g4", "f1g4", "f2j5", "f2j5a", "f2j7", "f2j7a", "f3j31", "f3j32", "f4j18", "f4j19", "f5g17", "f5g18", "k5f1k", "k5f1l", "k6d40", "k6d41", "k6d41", "k6d42", "k6d43", "k6d45", "k6d46", "k6d47", "m2j5", "m2j5a", "m2j7", "m2j7a", "m4j18", "m4j19", "m5g17", "m5g18", "n5f17", "n5f18", "p3a22", "p3a23", "p3a23a", "p3a23b", "p3a24", "p4a22", "p4a23", "p5h15", "p5h15b", "p5q3cr", "p6h74", "p6h75", "p6h76", "p6h77", "p6h78")
prenataldrugusevar<-c("m1g2", "m1g3", "m1g5", "m1g6")
FF_labeled=pheno %>%dplyr::select(any_of(c(idvars, demovars, smokingvars, prenataldrugusevar)))
#apply NA codes
my_na_codes<-c(-10:-1, "-1 Refuse", "-2 Don't know", "-3 Missing", "-9 Not in wave", "-7 N/A")
my_na_codes_numeric=
na_codes <- function(x, ...) {
  x[x %in% c(...)] <- NA
  if(is.factor(x)==TRUE){x<-droplevels(x)}
  x
}
FF_labeled=FF_labeled%>%mutate(across(everything(), ~na_codes(., my_na_codes)))%>%copy_labels(FF_labeled)

############################pdqc data######################################

#Note as of september 24, 2021, check with Jonah about final analytic subset sex outliers etc
pdqc_all<-readRDS(paste0(datadir, 'OGData/', "pd_qc.rds"))
#filter samples w/ >10 probe fail or sex !=predicted sex
pdqc<-pdqc_all %>% 
  filter(probe_fail_10==0 & sex==predicted_sex)%>%
  #add a new id for identifying technical replicates
  mutate(newID=paste(idnum, childteen, sep='_'))

#identify among technical replicates the sample with higher pct probe failure
reps_remove=pdqc%>%group_by(newID)%>%filter(n()>1)%>% #filter to technical replicates
  filter(probe_fail_pct==max(probe_fail_pct))%>% #filter to high pct probe fail
  pull(MethID)

#remove replicates with a higher pct probe failure & any samples identified as mother ids
pdqc_clean<-pdqc %>% filter(!(MethID %in% reps_remove) & childteen!='M')

#add batch, slide and plate data
batchData<-read.csv(file.path(datadir, 'Methylation_450K_array_batch_information.csv'))
pdqc_clean <- pdqc_clean %>% left_join(batchData)

############################create new variables###################################

#create following variables in FF_labeled:
#Methyldata - indicator for in analysis subset (have methylation 450k data) or no
#smkPregbinary Yes/no - any smoking during pregnancy
#m1g2_YesNoPreg and m1g3_YesNoPreg: yes/no alcohol and drugs (respectively)
#PostnatalMaternal smoking any: if mom smoke at either 1 or 5 - yes, if missing at either - missing, if no at both - no
#Postnatal maternal smoking does: if mom>1p/d at either 1 or 5 , if mom <1p/d at either 1 or 5, if no smoking at both, if missing + no smoking --> missing
FF_labeled=FF_labeled %>% 
  mutate(Methyldata=ifelse(id %in% pdqc_clean$id, 'Analysis subset', 'Not in analysis subset'), 
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
                                                m2j5a=='-6 Skip' & m4j19== '-6 Skip'~'No smoking when child age 1 or 5' ))%>%copy_labels(., FF_labeled)

#ancestry
#read in local and global pc data
pc_files=list.files(path=paste0(datadir, "OGData/pcs"), pattern='*.csv', full.names=TRUE)
names(pc_files)=gsub('.csv', '', list.files(path=paste0(datadir, "OGData/pcs"), pattern='*.csv'))
localpc=pc_files%>%map_dfr(read.csv, .id='ancestry')%>%rename_with(~str_c('local_', .), contains('PC'))
allPC<-read.table(paste0(datadir, "OGData/pcs/global_pcs.txt"), col.names=c("X", "idnum", paste("global_PC", 1:20, sep='')))

#individuals in each ancestry specific pc group can be labeled as that ancestry for categorical ancestry variable 
FF_labeled=FF_labeled %>% left_join(localpc%>%mutate(id=as.character(idnum))%>%select(id, ancestry, contains('local_PC')))%>%
  left_join(allPC%>%mutate(id=as.character(idnum))%>%select(id, contains('global_PC')))%>%
  mutate(ancestry=case_when(is.na(ancestry) & is.na(global_PC1) & Methyldata=='Analysis subset' ~'Missing PC data',
                            is.na(ancestry) & Methyldata!='Analysis subset'~'Not in analysis subset',
                            !is.na(ancestry)~ancestry))%>%copy_labels(., FF_labeled)

MissingPCdata<-FF_labeled %>% filter(ancestry=="Missing PC data")%>%select(id)
save(MissingPCdata, file = paste0(datadir, "CreatedData/MissingPC.Rdata"))

#######################filter to analysis subset and join with pdqc###################################
myFF=FF_labeled%>%filter(Methyldata=='Analysis subset')%>%
  left_join(pdqc_clean%>% 
              select(MethID, childteen, idnum, Sample_Plate, Slide, Array)%>%
              mutate(id=as.character(idnum)))%>%
  copy_labels(., FF_labeled)


#######################make visit specific variables###################################
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

#######################complete case dataset###################################

methylCohort<-pdqc_all %>% filter(childteen!='M')

basemodelvar<-c("smkPreg_binary", "ancestry", "cm1bsex", "cm1inpov", "ChildAgeComposite")
basemodeldata<-myFF %>% filter_at(all_of(basemodelvar), all_vars(!(. %in%c("Missing", "Missing PC data"))& !is.na(.)))

prenatalexposurevar<-c(basemodelvar, "m1g2_YesNoPreg", "m1g3_YesNoPreg")
prenatalexdata<-basemodeldata %>%filter_at(all_of(prenatalexposurevar), all_vars(!(. %in% c("Missing")) & !is.na(.)))

secondhandsmokevars<-c(prenatalexposurevar, "PostnatalMaternalSmokingAny", "SmkAtVisitPastmonth")
secondhandsmkdata<-prenatalexdata%>%filter_at(all_of(secondhandsmokevars), all_vars(!(. %in% c("Missing")) & !is.na(.)))

child_smoke<-secondhandsmkdata %>% filter(k5f1l!= "2 no" & childteen=='C')
child_nosmoke<-secondhandsmkdata %>% filter(k5f1l=="2 no" & childteen=='C')
teen_smoke<-secondhandsmkdata %>% filter(k6d40!="2 No" & childteen=='T')
teen_nosmoke<-secondhandsmkdata %>% filter(k6d40=="2 No" & childteen=='T')
completecase<-rbind(child_nosmoke, teen_nosmoke)

fathersmkdata<-completecase %>% filter(f1g4!="Missing" & !is.na(f1g4))

#####################################################################################Save all data 
save(file=paste0(datadir, '/CreatedData/allPhenoData.Rdata'))
#####################################################################################Save my specific data for modeling
save(completecase, file=paste0(datadir, '/CreatedData/completeCasepheno.Rdata'))