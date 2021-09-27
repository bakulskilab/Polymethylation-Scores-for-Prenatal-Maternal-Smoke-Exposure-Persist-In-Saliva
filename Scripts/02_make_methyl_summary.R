#libraries 
library(tidyr)
library(stringr)
library(dplyr)
library(purrr)
library(tidyselect)
library(broom)
library(compareGroups)
library(kableExtra)
library(haven)
library(sjlabelled)
library(matchmaker)
library(EpiDISH)
library("PerformanceAnalytics")
library(DescTools)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggcorrplot)
library(ggalluvial)
library(cowplot)
library(DiagrammeR)
library(sva)

#directory based coding
datadir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Data/"
john_datadir<-"/nfs/turbo/bakulski1/People/johndou/Fragile_Families/ProcessedFiles/jd"
outputdir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Output/"
codedir<-"/nfs/turbo/bakulski1/People/blostein/FF_methylation/Code/"

#
#polymethylation scores, global methylation data, clocks, 
#cell-type proportions, surrogate variable analysis