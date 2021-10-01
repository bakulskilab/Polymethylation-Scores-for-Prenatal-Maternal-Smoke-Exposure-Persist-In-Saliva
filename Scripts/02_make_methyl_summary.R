#directory based coding
datadir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/"
codedir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Code/"

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
clocks<-read.csv(file=paste0(datadir, "OGData/ffcw_n1776_8clocks.csv"), header = T)
colnames(clocks)[1]<-"MethID"


################################################join methyl data together 


#################################################checks


################################################labels 


#################################################save 
