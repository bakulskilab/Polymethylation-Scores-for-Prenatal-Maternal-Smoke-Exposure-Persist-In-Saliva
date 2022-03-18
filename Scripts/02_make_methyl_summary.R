library(here)
#directory based coding
datadir<-gsub('Code', 'Data/', here::here())
codedir<-paste0(here::here(), '/')
clock_coefs=file.path(paste0(datadir, 'OGData'), 'clock_coefs')


####################################################################library
library(EpiDISH)
library(ewastools)
library(tidyr)
library(purrr)
library(tibble)
library(sjlabelled)
library(wateRmelon)
library(data.table)
library(DunedinPoAm)
library(knitr)
library(dplyr)
library(DNAsmokeR)
####################################################################source functions
source(file.path(codedir, 'UsefulCode', 'makePolyEpiScores.R'))

#######################################chose apriori cpgs, if any
aprioriCG<-c("cg05575921", "cg04180046", "cg05549655", "cg14179389", "cg22132788")

################################read in methylation data
betaqc<-readRDS(file=paste0(datadir, 'OGData/Feb2022DataFreeze/', "beta_analytic_freeze1.rds"))
if(nrow(betaqc)!=437588){stop("Probe number not what thought, stop and check betaqc")}
###########################################cell type proportions 
celltypes<-epidish(beta.m = betaqc, ref.m = centEpiFibIC.m, method = "RPC")
#RPC method w/ centEpiFibIC
estF_FF<-as.data.frame(celltypes$estF)%>%
  tibble::rownames_to_column(var='MethID')%>%
  rename_with(~paste0(., '_RPC'), .cols=c('Epi', 'Fib', 'IC'))
#Middleton method
estL_FF<-estimateLC(betaqc, 'saliva', constrain=TRUE)%>%
  mutate(MethID=colnames(betaqc))%>%
  rename_with(~paste0(., '_saliva'), .cols=c('Leukocytes', 'Epithelial.cells'))

cell_types=left_join(estF_FF, estL_FF)

################################################global methylation
globalmethy<-colMeans(betaqc, na.rm=T)
globalmethdf<-data.frame("MethID"=names(globalmethy), "globalmethylation"=globalmethy*100)
#add annotations for stratified methylation mean calculations
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)
data(Islands.UCSC)
data(Other)
Locations <- Locations[rownames(betaqc), ]
Islands.UCSC <- Islands.UCSC[rownames(betaqc),]
Other <- Other[rownames(betaqc),]
#list of regions
isles=Islands.UCSC%>%as.data.frame()%>%mutate(Relation_to_Island=gsub('._', '', Relation_to_Island))%>%rownames_to_column('cpgs')%>%group_by(Relation_to_Island)%>%summarise(cpgs=list(cpgs))
#mean methylation per region
methyMean_regions=lapply(isles$cpgs, function(x) as.data.frame(colMeans(betaqc[x, ], na.rm=T)))
methyMean_regions=lapply(methyMean_regions, function(x) x*100)
methyMean_regions=methyMean_regions%>% bind_cols()%>%rownames_to_column('MethID')
colnames(methyMean_regions)=c('MethID', paste0(isles$Relation_to_Island, '_mean'))
globalmethdf=left_join(globalmethdf, methyMean_regions)

#################################################apriori cpgs
aprioriCGdf<-as.data.frame(t(betaqc[aprioriCG, ]))
aprioriCGdf=aprioriCGdf*100
aprioriCGdf$MethID = colnames(betaqc)

#################################################polymethylation scores from joubert & reese
cpgs<-read.csv(paste0(datadir, 'PMC4833289_SuppFile2CpGs.csv'), stringsAsFactors = F, header=F, skip=2)
#make headers
header1<-scan(paste0(datadir, 'PMC4833289_SuppFile2CpGs.csv'), nlines=1, what=character(), sep=',')
header1[7:26]<-c(rep("SSnewborn", 5), rep("SSnewbornCT", 5), rep("SSolder", 5), rep("anynewborn", 5))
header2<-paste(header1, scan(paste0(datadir, 'PMC4833289_SuppFile2CpGs.csv'), nlines=1, skip=1, what=character(), sep=','), sep='_')
names(cpgs)<-header2

#four meta-EWASes: SS_newborn (sustained smoking in newborns); SS_newbornCT (sustained smoking in newborns controlling for cell type) SS_older (sustained smoking in older children); any_newborn (any smoking in newborns)
#pivot_longer to get single column of coefs and a column for the analysis type
#and group into a list by analysis type
cpgs<-cpgs %>% pivot_longer(cols = c(matches("_Coef|_SE|_P|_FDR|_Direction")), names_to = c("set", ".value"), names_pattern="(.+)_(.+)")%>%rename_with(~gsub('_', '', .x))%>%group_split(set)
cpgs<-cpgs%>%setNames(lapply(cpgs, function(x){x$set %>% unique()}))

#recapitulate richmond 568 and 19 CpG scores in our data
bonferroni_threshold=1.07613e-07
cpgs$richmond19=cpgs$SSolder%>%filter(P<bonferroni_threshold)
cpgs$richmond568=cpgs$SSnewborn%>%filter(P<bonferroni_threshold)
#check
try(if(nrow(cpgs$richmond19)!=19 | nrow(cpgs$richmond568)!=568) stop('check richmond scores'))

# add reese coefficients 
reese_cpgs=read.csv(paste0(datadir, 'reese_cpgs.csv'), header=F)
names(reese_cpgs)=c('CpG', 'Coef', 'Gene')
reese_cpgs$set='reeseCordblood'
reese_cpgs=tbl_df(reese_cpgs)
cpgs=as.list(cpgs)
cpgs$reese=reese_cpgs

#make scores - no transform, zscore and centermean see makePoliEpiScore for function
scores=makePolyEpiScore(m=betaqc, b=cpgs)
scores_z=makePolyEpiScore(m=betaqc, b=cpgs, transformopt = 'zscore')
scores_center=makePolyEpiScore(m=betaqc, b=cpgs, transformopt = 'meancenter')

#add Meth ID
polyscores=list('notransform'=scores, 'zscore'=scores_z, 'center'=scores_center)
polyscores<-lapply(polyscores, function(x) x%>% mutate(MethID=colnames(betaqc)))
save(polyscores, file = file.path(datadir, 'CreatedData', 'polymethylationscores.RDS'))
polyscores_wide<-lapply(seq_along(polyscores), function(i) polyscores[[i]] %>% setNames(c(paste0(colnames(polyscores[[i]])[1:(ncol(polyscores$notransform)-1)], '_', names(polyscores)[i]), 'MethID')))%>%purrr::reduce(left_join, by='MethID')

  
#################################################smokeR scores from elastic net for Rauschert et al 
library(DNAsmokeR)
beta_short=betaqc[rownames(betaqc)[rownames(betaqc) %in% DNAsmokeR:::I450K$CpG],]
beta_short_t=as.data.frame(t(beta_short))
I450K=DNAsmokeR:::I450K
#rewrite smokeScore to specify dplyr select 
smokeRFunction=function(data=data, ARRAY=c('450k', 'EPIC'), class=c('class', 'prob')){
  SCORE <- NULL
  score <- NULL
  cpg   <- NULL
  availCpGs <- names(data)[names(data) %in% DNAsmokeR:::I450K$CpG]
  for (i in availCpGs) {
    CPG <- as.numeric(I450K %>%
                        filter(CpG %in% i) %>%
                        dplyr::select(Coefficient)
    )
    score <- cbind(score , as.numeric(unlist((CPG * data[,i]))))
    cpg <- c(cpg, i)
  }
  score         <- as.data.frame(score)
  names(score)  <- cpg
  smokScore   <- rowSums(score) + as.numeric(I450K[1,2])
  normalized  <- as.numeric((smokScore - min(smokScore,na.rm=T))/(max(smokScore, na.rm=T) - min(smokScore, na.rm=T)))
  predictSMOK <- factor(ifelse(normalized >0.5, "smoke_exp", "not_exp"))
  enscore=data.frame('MethID'=rownames(data), 'enscore_probs'=normalized, 'enscore_0.5'=predictSMOK)
  return(enscore)
}
rauschert_elastic_net=smokeRFunction(beta_short_t, '450k', 'prob')
#double check that the coefficients are different from those provided in Joubert
cpg_compare=map(cpgs, ~.x %>% left_join(I450K %>%dplyr::rename(en_coef=Coefficient))) # not the same 

#################################################epigenetic clocks feb  
if(nrow(betaqc)==437588){
  clocks<-read.csv(file=paste0(datadir, "OGData/allclocks_bothstudies_withcells.csv"), header = T)
  }
#the below code will run if the # CpGs if not using precalculated clocks if you switch to some other beta matrix, this code will recreate 
# the clocks for you, except levine and GRIM clock which requires different code. 
if(nrow(betaqc)!=437588){
pdqc_all<-readRDS(paste0(datadir, 'OGData/', "pd_qc.rds"))
pd <- data.table(pdqc_all%>%filter(MethID %in% colnames(betaqc)))[, .(MethID, idnum, childteen)]
if(ncol(betaqc)!=nrow(pd)){stop('Observations not equivalent')}
#horvath
pd[, horvath := agep(betaqc)]
#skin blood clock
skinblood <- fread(paste0(clock_coefs, "skin_blood_clock.csv"))
skinbloodc <- skinblood[, Coef] %>% c
names(skinbloodc) <- skinblood[, ID]
pd[, skinblood := agep(betaqc, coeff = skinbloodc)]
#pediatric
ped <- fread(paste0(clock_coefs, "kid_clock_coefs.csv"))
pedc <- ped[, Coef]
names(pedc) <- ped[, ID]
pd[, pediatric := agep(betaqc, coeff = pedc)]
#pace of aging 
pd[, c("poam38", "poam45") := PoAmProjector(betaqc)]
clocks=pd%>%dplyr::select(-idnum)
}
################################################join methyl data together 
methyldata=left_join(cell_types, globalmethdf)%>%
  left_join(aprioriCGdf)%>%
  left_join(polyscores_wide)%>%
  left_join(clocks)%>%
  left_join(rauschert_elastic_net)

################################################labels 

methyl_labels=c('Methylation data ID', 'Epithelial cell proportion', 'Fibroblast cel proportion', 'Immune cell proportion',
                'Immune cell proportion (saliva)','Epithelial cell proportion (saliva)', 
                'Global methylation', 'Island methylation', 'Open Sea methylation', 'Shelf methylation', 'Shore methylation',
                'AHRR: cg05575921', 'MYO1G: cg04180046', 'CYP1A1: cg05549655', 'GFI1: cg14179389', "MYO1G: cg22132788")

pms_labels=c('6074 site polymethylation score: coefficients for any smoking from Joubert newborn cordblood meta-analysis - PMC4833289 (no transform)', 
             '6074 site polymethylation score: coefficients for sustained smoking from Joubert newborn cordblood meta-analysis - PMC4833289 (no transform)', 
             '6074 site polymethylation score: coefficients for sustained smoking from Joubert newborn cordblood meta-analysis, w/ cell-type control - PMC4833289 (no transform)', 
            '6074 site polymethylation score: coefficients for sustained smoking from Joubert older children peripheral blood meta-analysis - PMC4833289 (no transform)', 
             '19 site polymethylation score (Richmond): coefficients for sustained smoking from Joubert older children peripheral blood meta-analysis - PMC4833289 (no transform)',
            '568 site polymethylation score (Richmond): coefficients for sustained smoking from Joubert newborn cordblood meta-analysis - PMC4833289 (no transform)',
             '28 site polymethylation score (Reese): coefficients for sustained smoking from newborn cordblood LASSO regression - PMC5381987 (no transform)') 
             

elastic_net_labels=c('Methylation data ID', '204 site polymethylation score classification probability (Rauschert):  coefficients for sustained smoking from older children peripheral blood elastic net regression', 'Elastic Net (Rauschert et al.) classification at p>0.5')

methyl_labels_myclocks=c(methyl_labels, 
                         pms_labels, gsub('no transform', 'z-score standardized', pms_labels), gsub('no transform', 'mean-centered', pms_labels), 
                         #paste0(pms_labels, '/n (coefficients) & z-score standardized (score)'), gsub('no transform', '/n z-score standardized (coefficients) & z-score standardized (score)', pms_labels), gsub('no transform', '/n mean-centered (coefficients) & z-score standardized (score)', pms_labels),
                          'Visit','Horvath clock', 'SkinBlood clock', 'Pediatric clock', 'PoAm clock 38', 'PoAm clock 45', 
                         elastic_net_labels[2:3])

methyl_labels_jonahclocks=c(methyl_labels, 
                            pms_labels, gsub('no transform', 'z-score standardized', pms_labels), gsub('no transform', 'mean-centered', pms_labels), 
                            #paste0(pms_labels, '/n (coefficients) & z-score standardized (score)'), gsub('no transform', '/n z-score standardized (coefficients) & z-score standardized (score)', pms_labels), gsub('no transform', '/n mean-centered (coefficients) & z-score standardized (score)', pms_labels),
                            'ID', 'Visit', 'ChildAge', 'Horvath clock', 'SkinBlood clock', 'Pediatric clock', ' PhenoAge (Levine clock)', 'PoAm clock 38', 'PoAm clock 45', 
                            'DNAmGDF_15', 'DNAmB2M', 'DNAmCystatin_C', 'DNAmTIMP_1', 'DNAmadm', "DNAmpai_1", 'DNAmleptin', 'GRIM clock (pack-years)', 
                            "Plasma Blast", "CD8pCD28nCD45RAn", "CD8_naive", 'Hannum clock', 'GRIM clock', 'BIOAge4HAStatic',
                            'PCHorvath1', 'PCHorvath2', 'PCHannum', 'PCPhenoAge', 'PCGrimAge', 'Leukocytes', 'Epithelial cells', 'Epi', 'Fib', 'IC', 
                            'GR', 'NK', 'B', 'CD4', 'CD8', 'MO', 'dupe', 'study',
                            elastic_net_labels[2:3])

if(nrow(betaqc)==437588){
  set_label(methyldata)=methyl_labels_jonahclocks
}
if(nrow(betaqc)!=437588){
  set_label(methyldata)=methyl_labels_myclocks
}



#################################################checks
summary(methyldata)

#load in complete case data set
load(paste0(datadir, '/CreatedData/completeCasepheno.Rdata'))
completecase=left_join(completecase, methyldata)%>%
  mutate(childteen=case_when(childteen=='C'~'Age 9',childteen=='T'~'Age 15'))
ncolumns=ncol(completecase)

#now zscore the polymethylation scores once we have the complete case dataset
completecase=completecase%>%
  mutate(across( anynewborn_notransform:reese_center, list(scale=scale)))

set_label(completecase)=c(get_label(completecase)[1:ncolumns], paste0(pms_labels, '/n (coefficients) & z-score standardized (score)'), gsub('no transform', '/n z-score standardized (coefficients) & z-score standardized (score)', pms_labels), gsub('no transform', '/n mean-centered (coefficients) & z-score standardized (score)', pms_labels))

set_label(completecase$childteen)='Child age at saliva collection'

#################################################save 
save(methyldata, file=paste0(datadir, '/CreatedData/allmethyl.Rdata'))
save(completecase, file=paste0(datadir, '/CreatedData/completeCasemethyl.Rdata'))
