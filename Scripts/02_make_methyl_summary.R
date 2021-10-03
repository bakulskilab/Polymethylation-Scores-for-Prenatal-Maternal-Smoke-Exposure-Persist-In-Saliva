#directory based coding
datadir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/"
codedir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Code/"
clock_coefs="/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/OGData/clock_coefs/"

####################################################################source functionslibrary
library(dplyr)
library(EpiDISH)
library(ewastools)
library(tidyr)
library(purrr)
library(sjlabelled)
library(wateRmelon)
library(data.table)
library(DunedinPoAm)
library(knitr)
####################################################################source functions
source(file.path(codedir, 'UsefulCode', 'makePolyEpiScores.R'))

#######################################chose apriori cpgs, if any
aprioriCG<-c("cg05575921", "cg04180046", "cg05549655", "cg14179389")

################################read in methylation data
betaqc<-readRDS(file=paste0(datadir, 'OGData/', "noob_filtered.rds"))

###########################################cell type proportions 
celltypes<-epidish(beta.m = betaqc, ref.m = centEpiFibIC.m, method = "RPC")
#RPC method w/ centEpiFibIC
estF_FF<-as.data.frame(celltypes$estF)%>%
  tibble::rownames_to_column(var='MethID')%>%
  rename_with(~paste0(., '_RPC'), .cols=c('Epi', 'Fib', 'IC'))
#Middleton method
estL_FF<-estimateLC(betaqc, 'saliva')%>%
  mutate(MethID=colnames(betaqc))%>%
  rename_with(~paste0(., '_saliva'), .cols=c('Leukocytes', 'Epithelial.cells'))

cell_types=left_join(estF_FF, estL_FF)

################################################global methylation
globalmethy<-colMeans(betaqc, na.rm=T)
globalmethdf<-as.data.frame(cbind("MethID"=names(globalmethy), "globalmethylation"=globalmethy))

#################################################apriori cpgs
aprioriCGdf<-as.data.frame(t(betaqc[aprioriCG, ]))
aprioriCGdf$MethID = colnames(betaqc)

#################################################polymethylation scores
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

#make scores - no transform, zscore and centermean see makePoliEpiScore for function
scores=makePolyEpiScore(m=betaqc, b=cpgs)
scores_z=makePolyEpiScore(m=betaqc, b=cpgs, transformopt = 'zscore')
scores_center=makePolyEpiScore(m=betaqc, b=cpgs, transformopt = 'meancenter')

#add Meth ID
polyscores=list('notransform'=scores, 'zscore'=scores_z, 'center'=scores_center)
polyscores<-lapply(polyscores, function(x) x%>% mutate(MethID=colnames(betaqc)))
save(polyscores, file = file.path(datadir, 'CreatedData', 'polymethylationscores.RDS'))
polyscores_wide<-lapply(seq_along(polyscores), function(i) polyscores[[i]] %>% setNames(c(paste0(colnames(polyscores[[i]])[1:4], '_', names(polyscores)[i]), 'MethID')))%>%reduce(left_join, by='MethID')

#################################################epigenetic clocks
###################THIS NEEDS TO BE CHANGED -- clocks weren't calculated on johns betaqc but on jonahs, 
#have different probe filter/sample filter sets
#I've been working with Johns betaqc but using jonahs clocks. 
#clocks<-read.csv(file=paste0(datadir, "OGData/ffcw_n1776_8clocks.csv"), header = T)
#colnames(clocks)[1]<-"MethID"
#Note that there are 3 IDS missing here from Jonah's clocks because of the sex filtration step 
#i.e. John dropped sex discordant pairs, jonah dropped sex outliers. so there wer 3 sex outliers here in 
#johns betaqc that weren't in jonahs. But I think all these clocks need to be rerun 
#since they were created using jonahs betaqc but I have calculated the pms, sva & the global methylation 
#using johns betaqc 
#options: either recalculate the clocks using johns betaqc or switch back to jonahs betaqc, fitler out discordant pairs and use 
#that only 
#pdqc 
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

################################################join methyl data together 
methyldata=left_join(cell_types, globalmethdf)%>%
  left_join(aprioriCGdf)%>%
  left_join(polyscores_wide)%>%
  left_join(pd)
  #left_join(clocks)

################################################labels 
pms_labels=c('Polymethylation score: coefficients for any smoking from newborn cordblood meta-analysis (no transform)', 
             'Polymethylation score: coefficients for sustained smoking from newborn cordblood meta-analysis (no transform)', 
             'Polymethylation score: coefficients for sustained smoking from newborn cordblood meta-analysis, w/ cell-type control (no transform)', 
              'Polymethylation score: coeffcients for sustained smoking from older children peripheral blood meta-analysis (no transform)')

set_label(methyldata)=c('Methylation data ID', 'Epithelial cell proportion', 'Fibroblast cel proportion', 'Immune cell proportion',
                        'Leukocytes proportion (saliva)','Epithelial cell proportion (saliva)', 
                        'Global methylation', 'AHHR: cg05575921', 'MYO1G: cg04180046', 'CYP1A1: cg05549655', 'GFI1: cg14179389', 
                        pms_labels, gsub('no transform', 'z-score standardized', pms_labels), gsub('no transform', 'mean-centered', pms_labels), 
                        'ID', 'Visit',
                        'Horvath clock', 'SkinBlood clock', 'Pediatric clock', 'PoAm clock 38', 'PoAm clock 45')


#################################################checks
summary(methyldata)
#note the three missing clock ids when left joining to clocks. Also see above. needs to make a decision here about what to do 

#load in complete case analysis
load(paste0(datadir, '/CreatedData/completeCasepheno.Rdata'))
completecase=left_join(completecase, methyldata)


#################################################save 
save(completecase, file=paste0(datadir, '/CreatedData/completeCasemethyl.Rdata'))
