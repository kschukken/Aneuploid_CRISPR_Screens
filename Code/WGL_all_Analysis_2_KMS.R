### CRISPR Screen data analysis of Whole Genome Library data for Paired cells ####
# All WGL virus screens combined
# Author: Klaske M. Schukken
# October 2023



### INTRO AND LIBRARIES ####
## 231006

### DEFINED FUNCTIONS
### visualizing gene dropout per cell line and OVERALL/COMBINED DATA
### remove positive and negative controls
## look at test sgRNAs: non-coding location dropouts
##        manually set some positive control dropouts


## Folders 
## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening and proteomics datasets: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"


### libraries
#mle analysis data
library(ggpubr)
library(reshape2)
library(tidyverse)
library(plyr)
library(ggplot2) 
library("viridis")  # color packet
library("cowplot") 
library('gprofiler2')
library(xlsx)
library(readxl)
library(readr)

###GSEA DESeq2
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("clusterProfiler", version = "3.19")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
#BiocManager::install("limma")
library(DESeq2)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)
library("limma")


## COLOR CHOICES: 

# Most dots: "black" 
# Top hits: "red" 
# Top hit (less significant): "red4"

# MitoElectronChainnAssembly             color="#F8766D"
# MitoTranslationTranscription           color="gold2"
# RNA_Processing_rRNA                    color="deepskyblue4"
# Ribosomal_Genes_noMT$Approved.symbol   color="#C77CFF"
# Proteosome_Genes$Approved.symbol       color="chartreuse1"
# RNA_Processing_Spliceosome             color="deepskyblue1"




### GET DATA #####
# MAGECK MLE data
# The ‘p-value’ is calculated by randomly permuting sgRNA labels.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6862721/


## Run function per cell line: 
#  HCT116
setwd(SchukkenData)

HCT116.Cont1.mle <- read.delim2("mle_HCT116_Cont1.gene_summary.txt", 
                                dec=".", header = TRUE)
#HCT116.Ts13.mle<- read.delim2("mle_HCT116_Ts13.gene_summary.txt", 
#                              dec=".", header = TRUE)
HCT116.Cont1.mle$Gene<- gsub("Non-Targeting_Control-N","neg ", HCT116.Cont1.mle$Gene) 
#HCT116.Ts13.mle$Gene<- gsub("Non-Targeting_Control-N","neg ", HCT116.Ts13.mle$Gene) 

HCT116.Control.guide <- read_tsv("Ryan-Count_HCT116_Cont1_i.f.head.txt")
median(HCT116.Control.guide$Cont1.initial)  #1046
median(HCT116.Control.guide$Cont1.final)    #1015



#  DLD1
setwd(SchukkenData)

DLD1.Cont1_NoMerge.mle <- read.delim2("mle_Cont.gene_summary_noMergeControl.txt", 
                                      dec=".", header = TRUE)
DLD1.Ts13_NoMerge.mle<- read.delim2("mle_Ts13.gene_summary_noMergeControl.txt", 
                                    dec=".", header = TRUE)
DLD1.Cont1_NoMerge.mle$Gene<- gsub("Non-Targeting_Control-N","neg ", DLD1.Cont1_NoMerge.mle$Gene) 
DLD1.Ts13_NoMerge.mle$Gene<- gsub("Non-Targeting_Control-N","neg ", DLD1.Ts13_NoMerge.mle$Gene) 


DLD1.Cont1.guide <- read.csv("count-DLD1_Cont-i.f.csv")
DLD1.Ts13.guide <- read.csv("count-DLD1_Ts13-i.f.csv")
median(DLD1.Cont1.guide$Cont.initial)  #1035
median(DLD1.Cont1.guide$Cont.final)    #1109
median(DLD1.Ts13.guide$Ts13.initial)   #1113
median(DLD1.Ts13.guide$Ts13.final)     #1108




#  Vaco432
setwd(SchukkenData)

Vaco432.Cont1.noMerge.mle <- read.delim2("mle_Vaco-Cont_unique-controls.gene_summary.txt", 
                                         dec=".", header = TRUE)
Vaco432.Ts13.noMerge.mle<- read.delim2("mle_Vaco-Ts13_unique-controls.gene_summary.txt", 
                                       dec=".", header = TRUE)
Vaco432.Cont1.noMerge.mle$Gene<- gsub("Non-Targeting Control","neg", Vaco432.Cont1.noMerge.mle$Gene) 
Vaco432.Ts13.noMerge.mle$Gene<- gsub("Non-Targeting Control","neg", Vaco432.Ts13.noMerge.mle$Gene) 


Vaco432.Cont1.guide <- read.csv("count-Vaco-Cont-i.f.csv")
Vaco432.Ts13.guide <- read.csv("count-Vaco-Ts13-i.f.csv")
median(Vaco432.Cont1.guide$Vaco_Control_D0)  #860
median(Vaco432.Cont1.guide$Vaco_Control_D10) #929
median(Vaco432.Ts13.guide$Vaco_Ts13_D0)      #966
median(Vaco432.Ts13.guide$Vaco_Ts13_D10)     #862




### A2780  cells
setwd(SchukkenData)

A2780.Cont1.mle <- read.delim2("mle_count-A2780-Cont1-d0-d10.gene_summary.txt", 
                               dec=".", header = TRUE)
A2780.Di1q.mle<- read.delim2("mle_count-A2780-1qDi.gene_summary.txt", 
                             dec=".", header = TRUE)
# Heatmap A2780 
A2780.Cont1.mle$Gene<- gsub("neg_", "neg ", A2780.Cont1.mle$Gene) 
A2780.Di1q.mle$Gene<- gsub("neg_", "neg ", A2780.Di1q.mle$Gene) 


# A2780 cells: guide RNA info
A2780.Cont1.guide <- read.csv("count-A2780-Cont1-d0.d10.csv", 
                              dec=".", header = TRUE)
A2780.Di1q.guide <- read.csv("count-A2780-1qDi-d0.d10.csv", 
                             dec=".", header = TRUE)
median(A2780.Cont1.guide$A2780_Cont1_D0) #432
median(A2780.Cont1.guide$A2780_Cont1_D10) #485
median(A2780.Di1q.guide$A2780_1qDi_D0) #657
median(A2780.Di1q.guide$A2780_1qDi_D10) #674




### A2058 cells +/- 1q Disomy:
setwd(SchukkenData)

A2058.Cont1.mle <- read.delim2("mle_count-A2058WT.gene_summary.txt", 
                               dec=".", header = TRUE)
A2058.Di1q.mle<- read.delim2("mle_count-A20581qDi.gene_summary.txt", 
                             dec=".", header = TRUE)

A2058.Cont1.mle$Gene<- gsub("neg_", "neg ", A2058.Cont1.mle$Gene) 
A2058.Di1q.mle$Gene<- gsub("neg_", "neg ", A2058.Di1q.mle$Gene) 


# A2058 cells: guide RNA info
A2058.Cont1.guide <- read.csv("count-A2058_WT.csv", 
                              dec=".", header = TRUE) #mean counts 860 for D=0 and 949 for D=10
A2058.Di1q.guide <- read.csv("count-A2058_1qDi.csv", 
                             dec=".", header = TRUE) #mean counts 835 for D=0 and 848 for D=10
median(A2058.Cont1.guide$A2058_WT_D.0)   #727
median(A2058.Cont1.guide$A2058_WT_D.10)  #785
median(A2058.Di1q.guide$A2058_1qDi_D.0)  #706
median(A2058.Di1q.guide$A2058_1qDi_D.10) #686




### A2058 7pDi: 
setwd(SchukkenData)

A2058.7pDi.mle <- read.delim2("mle_count-A2058_7qDi.gene_summary.txt", 
                               dec=".", header = TRUE)
A2058.7pDi.mle$Gene<- gsub("neg_", "neg ", A2058.7pDi.mle$Gene) 


# guide RNA info
A2058.7pDi.guide <- read.csv("count-A2058_7qDi.csv", dec=".", header = TRUE) 
median(A2058.7pDi.guide$A2058_7qDi_D.0) #1150
median(A2058.7pDi.guide$A2058_7qDi_D.10) #1148




### AGS  cells
setwd(SchukkenData)

AGS.Cont1.mle <- read.delim2("mle_count-AGS_Control.gene_summary.txt", 
                               dec=".", header = TRUE)
AGS.Di1q.mle<- read.delim2("mle_count-AGS_1qDi.gene_summary.txt", 
                             dec=".", header = TRUE)
# Heatmap AGS 
AGS.Cont1.mle$Gene<- gsub("neg_", "neg ", AGS.Cont1.mle$Gene) 
AGS.Di1q.mle$Gene<- gsub("neg_", "neg ", AGS.Di1q.mle$Gene) 


# AGS cells: guide RNA info
AGS.Cont1.guide <- read.csv("count-AGS_Control.csv", 
                              dec=".", header = TRUE)
AGS.Di1q.guide <- read.csv("count-AGS_1qDi.csv",       # AGS Di1q has high D=0 coverage. ~1500 reads/guide. others have ~500 reads per guide
                             dec=".", header = TRUE)
median(AGS.Cont1.guide$AGS_Control_D0) #426.  
median(AGS.Cont1.guide$AGS_Control_D10) #519
median(AGS.Di1q.guide$AGS_1qDi_D0) #1954. 
median(AGS.Di1q.guide$AGS_1qDi_D10) # 487  




### MCF10A  cells
setwd(SchukkenData)

MCF10A.Cont1.mle <- read.delim2("mle_count-MCF10A_WT.gene_summary.txt", 
                             dec=".", header = TRUE)
MCF10A.Di1q.mle<- read.delim2("mle_count-MCF10A_1qDi.gene_summary.txt", 
                           dec=".", header = TRUE)
# Heatmap MCF10A 
MCF10A.Cont1.mle$Gene<- gsub("neg_", "neg ", MCF10A.Cont1.mle$Gene) 
MCF10A.Di1q.mle$Gene<- gsub("neg_", "neg ", MCF10A.Di1q.mle$Gene) 

# MCF10A cells: guide RNA info
MCF10A.Cont1.guide <- read.csv("count-MCF10A_WT.csv", 
                            dec=".", header = TRUE)
MCF10A.Di1q.guide <- read.csv("count-MCF10A_1qDi.csv",       
                           dec=".", header = TRUE)
median(MCF10A.Cont1.guide$MCF10A_WT_D0) # 854
median(MCF10A.Cont1.guide$MACF10A_WT_D10) # 932
median(MCF10A.Di1q.guide$MCF10A_1qDi_c32_D0) # 1156
median(MCF10A.Di1q.guide$MCF10A_1qDi_c32_D10) # 901





#### Rename all Excel "Date" gene names

# Make function to rename "Sept-1" genes and "SEPT1" to their updated names "SEPTIN1" etc. for all 
#   genes who's names are misinterpretted as dates in Excel. 
# Excel also sometimes turns "dates" into a random string of numbers. so they need to converted back to gene names. 
RenameExcelGenes<- function(GeneMLEFile){
  GeneMLEFile$Gene<- as.character(GeneMLEFile$Gene)
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '1-Dec'] <- "DELEC1"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'DEC1'] <- "DELEC1"
  
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '15-Sep'] <- "SELENOF"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '1-Sep'] <- "SEPTIN1"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '2-Sep'] <- "SEPTIN2"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '3-Sep'] <- "SEPTIN3"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '4-Sep'] <- "SEPTIN4"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '5-Sep'] <- "SEPTIN5"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '6-Sep'] <- "SEPTIN6"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '7-Sep'] <- "SEPTIN7"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '8-Sep'] <- "SEPTIN8"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '9-Sep'] <- "SEPTIN9"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '10-Sep'] <- "SEPTIN10"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '11-Sep'] <- "SEPTIN11"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '12-Sep'] <- "SEPTIN12"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '14-Sep'] <- "SEPTIN14"
  
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEP15'] <- "SELENOF"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT1'] <- "SEPTIN1"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT2'] <- "SEPTIN2"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT3'] <- "SEPTIN3"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT4'] <- "SEPTIN4"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT5'] <- "SEPTIN5"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT6'] <- "SEPTIN6"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT7'] <- "SEPTIN7"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT8'] <- "SEPTIN8"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT9'] <- "SEPTIN9"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT10'] <- "SEPTIN10"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT11'] <- "SEPTIN11"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT12'] <- "SEPTIN12"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'SEPT14'] <- "SEPTIN14"
  
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '1-Mar'] <- "MARCHF1"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '2-Mar'] <- "MARCHF2"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '3-Mar'] <- "MARCHF3"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '4-Mar'] <- "MARCHF4"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '5-Mar'] <- "MARCHF5"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '6-Mar'] <- "MARCHF6"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '7-Mar'] <- "MARCHF7"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '8-Mar'] <- "MARCHF8"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '9-Mar'] <- "MARCHF9"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '10-Mar'] <- "MARCHF10"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '11-Mar'] <- "MARCHF11"
  
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH1'] <- "MARCHF1"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH2'] <- "MARCHF2"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH3'] <- "MARCHF3"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH4'] <- "MARCHF4"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH5'] <- "MARCHF5"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH6'] <- "MARCHF6"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH7'] <- "MARCHF7"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH8'] <- "MARCHF8"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH9'] <- "MARCHF9"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH10'] <- "MARCHF10"
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == 'MARCH11'] <- "MARCHF11"
  
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43525'] <- 'MARCHF1'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43534'] <- 'MARCHF10'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43535'] <- 'MARCHF11'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43526'] <- 'MARCHF2'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43527'] <- 'MARCHF3'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43528'] <- 'MARCHF4'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43529'] <- 'MARCHF5'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43530'] <- 'MARCHF6'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43531'] <- 'MARCHF7'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43532'] <- 'MARCHF8'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '43533'] <- 'MARCHF9'
  
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44986'] <- 'MARCHF1'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44987'] <- 'MARCHF2'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44988'] <- 'MARCHF3'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44989'] <- 'MARCHF4'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44990'] <- 'MARCHF5'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44991'] <- 'MARCHF6'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44992'] <- 'MARCHF7'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44993'] <- 'MARCHF8'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44994'] <- 'MARCHF9'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44995'] <- 'MARCHF10'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '44996'] <- 'MARCHF11' 
  
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45170'] <- 'SEPTIN1' #
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45171'] <- 'SEPTIN2' #
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45172'] <- 'SEPTIN3'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45173'] <- 'SEPTIN4'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45174'] <- 'SEPTIN5'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45175'] <- 'SEPTIN6'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45176'] <- 'SEPTIN7'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45177'] <- 'SEPTIN8'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45178'] <- 'SEPTIN9'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45179'] <- 'SEPTIN10'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45180'] <- 'SEPTIN11'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45181'] <- 'SEPTIN12'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45183'] <- 'SEPTIN14' #
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45184'] <- 'SELENOF'
  GeneMLEFile["Gene"][GeneMLEFile["Gene"]  == '45261'] <- 'DELEC1'
  
  return(GeneMLEFile)
} 

Vaco432.Cont1.mle<- RenameExcelGenes(Vaco432.Cont1.noMerge.mle)
Vaco432.Ts13.mle<- RenameExcelGenes(Vaco432.Ts13.noMerge.mle)
DLD1.Cont1.mle <- RenameExcelGenes(DLD1.Cont1_NoMerge.mle)
DLD1.Ts13.mle<- RenameExcelGenes(DLD1.Ts13_NoMerge.mle)
HCT116.Cont1.mle <- RenameExcelGenes(HCT116.Cont1.mle)
#HCT116.Ts13.mle<- RenameExcelGenes(HCT116.Ts13.mle). ## Removed due to karyotype sowing Ts13 loss at endpoint

A2780.Cont1.mle<- RenameExcelGenes(A2780.Cont1.mle)
A2780.Di1q.mle<- RenameExcelGenes(A2780.Di1q.mle)

A2058.Cont1.mle<- RenameExcelGenes(A2058.Cont1.mle)
A2058.Di1q.mle<- RenameExcelGenes(A2058.Di1q.mle)
A2058.7pDi.mle<- RenameExcelGenes(A2058.7pDi.mle)

AGS.Cont1.mle<- RenameExcelGenes(AGS.Cont1.mle)
AGS.Di1q.mle<- RenameExcelGenes(AGS.Di1q.mle)

MCF10A.Cont1.mle<- RenameExcelGenes(MCF10A.Cont1.mle)
MCF10A.Di1q.mle<- RenameExcelGenes(MCF10A.Di1q.mle)




## Get rid of .guide files after measuring counts
# They take up a lot of space and I don't use them after this point in script
HCT116.Control.guide<-NA
Vaco432.Cont1.guide<-NA
Vaco432.Ts13.guide<-NA
DLD1.Cont1.guide<-NA
DLD1.Ts13.guide<-NA
MCF10A.Cont1.guide<-NA
MCF10A.Di1q.guide<-NA
A2780.Cont1.guide<-NA
A2780.Di1q.guide<-NA
A2058.7pDi.guide<-NA
A2058.Cont1.guide<-NA
A2058.Di1q.guide<-NA
AGS.Cont1.guide<-NA
AGS.Di1q.guide<-NA





## DepMap: gene info & Chrm number
setwd(Dependency)
Gene_Info<- read_tsv("mart_export-2.txt")
Protein_Info3<-read_csv("Protein_location_info.csv") # data from BioMart

## downloaded human gene list with chromosome number and other IDs attached. 
## downloaded from ensemble. http://www.ensembl.org/biomart/martview/3857ff2676e33f2f2b9cff3d7fa0c8c1 
GeneList1q<- unique(subset(Protein_Info3, ChrmArm=="1q")$`Approved Symbol`)
GeneListMT<- unique(subset(Protein_Info3, Chromosome=="mitochondria")$`Approved Symbol`)
GeneList7p<- unique(subset(Protein_Info3, ChrmArm=="7p")$`Approved Symbol`)
GeneList13<- unique(subset(Protein_Info3, ChrmArm=="13q")$`Approved Symbol`)
GeneList20q<- unique(subset(Protein_Info3, ChrmArm=="20q")$`Approved Symbol`)



# Model.csv
# Sample info about DepMap cell lines
# Metadata for all of DepMap’s cancer models/cell lines. A full description of each column is available in the DepMap Release README file.
#Current DepMap Release data, including CRISPR Screens, PRISM Drug Screens, Copy Number, Mutation, Expression, and Fusions
#DepMap, Broad (2024). DepMap 24Q2 Public. Figshare+. Dataset. https://doi.org/10.25452/figshare.plus.25880521.v1
# Downloaded November 8, 2024
# 24Q2
Cell_Info<- read.csv("Model.csv", header=TRUE)



### CRISPRGeneEffect.csv
# Gene Effect scores derived from CRISPR knockout screens published by Broad’s Achilles and Sanger’s SCORE projects.
# Negative scores imply cell growth inhibition and/or death following gene knockout. Scores are normalized such that nonessential genes have a median score of 0 and independently identified common essentials have a median score of -1.
# Gene Effect scores were inferenced by Chronos ( https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02540-7 )
# Integration of the Broad and Sanger datasets was performed as described in https://doi.org/10.1038/s41467-021-21898-7, except that quantile normalization was not performed.
# Current DepMap Release data, including CRISPR Screens, PRISM Drug Screens, Copy Number, Mutation, Expression, and Fusions
# DepMap, Broad (2024). DepMap 24Q2 Public. Figshare+. Dataset. https://doi.org/10.25452/figshare.plus.25880521.v1

# Downloaded November 8, 2024
# 24Q2
CRISPR_Gene_all<- read.csv("CRISPRGeneEffect.csv", header=TRUE)
colnames(CRISPR_Gene_all)[which(names(CRISPR_Gene_all) == "X")] <- "ModelID"

## To look at reliable data only, I set a minimum cell line number to >200
# Removed 512 samples that had only 23-24 cell lines screened
keep_cols <- colSums(!is.na(CRISPR_Gene_all)) >= 200
CRISPR_Gene <- CRISPR_Gene_all[, keep_cols]



#### Aneuploidy score data from DepMap 
# aneuploidy_scores.csv
# Aneuploidy data per cell line.  Including total aneuploidy score, and cell line mean ploidy and +/- Genome doubling
# Downloaded November 9, 2023, from DepMap CCLE
# May 2023 data
Aneuploidy_Scores<- read.csv("aneuploidy_scores.csv", header=TRUE)
Aneuploidy_Scores$Score.divPloidy<- Aneuploidy_Scores$Aneuploidy.score/ round(Aneuploidy_Scores$Ploidy, digits=0)



# CCLE_segment_cn.csv
# Segment level copy number data 
# - DepMap_ID 
# - Chromosome 
# - Start (bp start of the segment) 
# - End (bp end of the segment) 
# - Num_Probes (the number of targeting probes that make up this segment) 
# - Segment_Mean (relative copy ratio for that segment) 
# - amplification status (+,-,0)
ChrmCN<-read.csv("CCLE_segment_cn.csv", header=TRUE)


## arm based aneuploidy and Copy Number per cell line
## download October 14, 2022
Arm_CN<- read.csv("Arm-level_CNAs.csv", header=TRUE)






# Get pathways: 
Proteosome_Genes<- read.csv("group-690.csv", header=TRUE)

Olfactory_Genes<- read.csv("group-141.csv", skip = 1 , header=TRUE)
Ribosomal_Genes<- read.csv("group-1054.csv", skip = 1 , header=TRUE) #HGNC Ribosomal Proteins
Ribosomal_Genes_noMT<- subset(Ribosomal_Genes, Group %in% c("L ribosomal proteins", "S ribosomal proteins"))

RNA_Processing<- read_tsv("GO_term_summary_20250720_134600.txt")
RNA_Processing1<-unique(toupper(RNA_Processing$Symbol))
RNA_Processing_Spliceosome<- unique(toupper(subset(RNA_Processing, `Annotated Term` %in% 
                                                     c("mRNA splicing, via spliceosome", 
                                                       "regulation of mRNA splicing, via spliceosome") )$Symbol ) )
RNA_Processing_rRNA<- unique(toupper(subset(RNA_Processing, `Annotated Term` %in% 
                                              c("rRNA processing", "regulation of rRNA processing"))$Symbol))


MitoTranslationTranscription <- c("TRMT10C","HSD17B10","PRORP", 
                                  "POLRMT", "TFAM", "TFB2M",
                                  "MRPL58", "MTRFR", "MTRF1", "MTRF1L", 
                                  "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10", "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16", "MRPL17", 
                                  "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22", "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30", "MRPL32", "MRPL33", "MRPL34", "MRPL35", "MRPL36", 
                                  "MRPL37", "MRPL38", "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44", "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50", 
                                  "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL58", "MRPL57", "GADD45GIP1", "MRPS30", "MRPS18A", "MRPS2", "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10", 
                                  "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17", "MRPS18C", "MRPS21", "MRPS22", "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", 
                                  "MRPS28", "DAP3", "MRPS31", "MRPS33", "MRPS34", "MRPS35", "CHCHD1", "AURKAIP1", "PTCD3", "MRPS18B", 
                                  "FASTK", "FASTKD1", "FASTKD2", "FASTKD3", "TBRG4", "FASTKD5")
MitoElectronChainnAssembly <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2", "NDUFAB1","NDUFA1",
                                "NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12",
                                "NDUFA13","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
                                "NDUFB11","NDUFC1","NDUFC2","NDUFS4","NDUFS5","NDUFS6","NDUFV3", 
                                "SDHA", "SDHB", "SDHC", "SDHD", 
                                "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB","CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11", 
                                "COX4I1","COX4I2","COX5A","COX5B","COX6A1", "COX6A2","COX6B1","COX6B2","COX6C", "COX7A1","COX7A2","COX7B","COX7B2","COX7C",
                                "COX8A","COX8C","MT-CO1","MT-CO2","MT-CO3", 
                                "ATP5F1A", "ATP5F1B", "ATP5F1C",  "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5MJ", "ATP5MK", 
                                "MT-ATP6", "MT-ATP8", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP5IF1", 
                                "ATPAF1","ATPAF2","BCS1L","CHCHD7","CMC1","CMC2", "COA1",  "COA3", "COA4","COA5", "COA6",  "COA7", "COA8",
                                "COX7A2L","COX10", "COX11", "COX14","COX15", "COX16",  "COX17", "COX18", "COX19", "COX20",
                                "DMAC1", "DMAC2", "FDXR", "FDX2", "FMC1", "FOXRED1","HCCS","HIGD1A","HIGD2A","LRPPRC", "LYRM7",
                                "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6","NDUFAF7",  "NDUFAF8","NUBPL",
                                "OCIAD2", "OXA1L", "PET100","PET117","PNKD", "SCO1", "SCO2",  "SDHAF1","SDHAF2","SDHAF3", "SDHAF4", "SFXN4",
                                "SMIM20", "SURF1","TACO1", "TIMMDC1","TMEM70","TMEM126A", "TMEM177","TMEM186","TMEM223","TMEM242","TTC19", 
                                "UQCC1","UQCC2","UQCC3", "UQCC4","UQCC5","UQCC6", 
                                "ACAD9","COA1","ECSIT","NDUFAF1","TMEM126B","TMEM186", 
                                "BOLA3", "FDXR", "FXN", "GLRX5", "IBA57", "ISCA1", "ISCA2", "ISCU", "LYRM4", "NFS1", "NFU1", "NUBPL") # HGNC groups




### PLOT RAW CONTROLS #####
# 220613
# Make heatmap for control RNA dropout. heatmap beta scores
# DEFINE FUNCTION
# get count file (from Mageck) as data frames. 

## Define function and positive/negative control RNAs.
PosCont<- c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3") 
NegCont<- c("neg ") #if you use grepl, you can find all genes containing "neg"

##Define function: for plotting control guide dropout
# Use independent datasets of mle data. plot them together. 

PlotHeatmap.CtrlRaw<-function(DropoutData = list("A2780 Control 1" =A2780.Cont1.mle, 
                                              "A2780 1q Di"=A2780.Di1q.mle), 
                           Title= "Control guide Beta-scores",
                           Pos.Cont=PosCont, Neg.Cont=c("neg")){
  ### Step 1: Get beta scores per cell line
  Pcompare<- DropoutData[[1]][, c(1,3)]
  for (i in 2:length(DropoutData)){
    Pcompare<- merge(Pcompare, DropoutData[[i]][, c(1,3)], by.x="Gene", by.y= "Gene", all=TRUE)
  }
  
  ### Step 2: Get control guides from all cell lines
  ## Subset Positive and Negative controls. 
  Poscntvalues<- data.frame()
  Negcntvalues<- data.frame()
  for (i in Pos.Cont) {# find positive control genes
    intermediate <- Pcompare[which(Pcompare$Gene==i),] #Find exact match of poscont label in $Gene
    Poscntvalues <- rbind(Poscntvalues, intermediate)
  }
  
  for (i in Neg.Cont) {# find negative control genes
    intermediate <-Pcompare[which(grepl(i, Pcompare$Gene)),] #Find all rows containing  neg values
    Negcntvalues<- rbind(Negcntvalues, intermediate)
  }
  Allcntgds<- rbind(Poscntvalues, Negcntvalues[c(1:6,8:11),])  #Use this if you have non-merged control Negative values
  #Allcntgds<- rbind(Poscntvalues, Negcntvalues) #use this if you have the merged control negative values. 
  
  ##Step 3: Plot
  ## plot positve and negative control beta scores
  colnames(Allcntgds)<- c("Gene", names(DropoutData))
  dfm <- melt(Allcntgds, id.vars = 1)
  dfm$Gene <- factor(dfm$Gene, level=unique(dfm$Gene)) #make Gene names a factor, not alphabetical
  
  plot.controls<- ggplot(dfm, aes(x = Gene, y = variable, fill = value)) + 
    geom_tile()+
    ggtitle(Title)+    
    xlab("Gene")+
    ylab("Cell line") +
    scale_fill_gradient2(low="Red", mid="white", high="Blue", midpoint=0, 
                         name="Beta-Score") +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust=1))
  
  
  plot.controls
}


## Run function per cell line: Not quantile normalized: 
PlotHeatmap.CtrlRaw(DropoutData= list("HCT116 Control"= HCT116.Cont1.mle, 
                                   #"HCT116 Ts13"= HCT116.Ts13.mle, 
                                   "DLD1 Control" = DLD1.Cont1_NoMerge.mle, 
                                   "DLD1 Ts13"= DLD1.Ts13_NoMerge.mle, 
                                   "Vaco432 Control" = Vaco432.Cont1.noMerge.mle, 
                                   "Vaco432 Ts13"= Vaco432.Ts13.noMerge.mle,
                                   "A2780 Control" =A2780.Cont1.mle, 
                                   "A2780 1qDisomy"= A2780.Di1q.mle, 
                                   "A2058 Control" =A2058.Cont1.mle, 
                                   "A2058 1qDisomy"= A2058.Di1q.mle, 
                                   "A2058 7pDisomy"= A2058.7pDi.mle,
                                   "AGS Control"= AGS.Cont1.mle, 
                                   "AGS 1qDisomy"= AGS.Di1q.mle,
                                   "MCF10A Control"= MCF10A.Cont1.mle, 
                                   "MCF10A 1qDisomy"= MCF10A.Di1q.mle
), Title="Control guide Beta-scores")
# plot.heatmap.WGL.controls.pdf
# 8x5 






### Quantile Normalization #####
## Quantile normalize the dropouts across all WGL samples
## essentially rank order them, then take the average dropout per rank, give that number per gene. 
## then run data with the new rank-order-averaged dropout scores. 

### Below function will quantile normalize your data
### !!! BUT!!! you need merge control and aneuploid data, and ONLY give function the BETA SCORE COLLUMNS!!!! 
## ALSO! Dataframe row names need to be gene names! 
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final) 
}



## to quantile normalize, merge Cont3 and Di1q data together, set gene name as row name, and look only at BETA SCORES

WGL.ALL<- merge(DLD1.Cont1.mle[,c(1,3)], DLD1.Ts13.mle[,c(1,3)], by.x="Gene", by.y="Gene")
WGL.ALL<- merge(WGL.ALL, HCT116.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
WGL.ALL<- merge(WGL.ALL, Vaco432.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene") 
WGL.ALL<- merge(WGL.ALL, Vaco432.Ts13.mle[,c(1,3)], by.x="Gene", by.y="Gene")  
WGL.ALL<- merge(WGL.ALL, A2780.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
WGL.ALL<- merge(WGL.ALL, A2780.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")
WGL.ALL<- merge(WGL.ALL, A2058.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
WGL.ALL<- merge(WGL.ALL, A2058.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")
WGL.ALL<- merge(WGL.ALL, A2058.7pDi.mle[,c(1,3)], by.x="Gene", by.y="Gene")
WGL.ALL<- merge(WGL.ALL, AGS.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene") 
WGL.ALL<- merge(WGL.ALL, AGS.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")  
WGL.ALL<- merge(WGL.ALL, MCF10A.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene") 
WGL.ALL<- merge(WGL.ALL, MCF10A.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")  
colnames(WGL.ALL)<- c("Gene", 
                                 "DLD1.Cont.Beta", "DLD1.Ts13.Beta", 
                                 "HCT116.Cont.Beta", #"HCT116.Ts13.Beta", 
                                 "Vaco432.Cont.Beta", "Vaco432.Ts13.Beta", 
                                 "A2780.Cont.Beta", "A2780.Di1q.Beta", 
                                 "A2058.Cont.Beta", "A2058.Di1q.Beta", "A2058.7pDi.Beta", 
                                 "AGS.Cont.Beta", "AGS.Di1q.Beta", 
                                 "MCF10A.Cont.Beta", "MCF10A.Di1q.Beta")


WGL.ALL.QN<- quantile_normalisation(WGL.ALL[,c(2:length(WGL.ALL))])
WGL.ALL.QN<- as.data.frame(WGL.ALL.QN)
WGL.ALL.QN$Gene<- WGL.ALL$Gene
head(WGL.ALL.QN)
WGL.ALL.QN$Gene<- as.character(WGL.ALL.QN$Gene)

quantile(WGL.ALL.QN[,1])
quantile(WGL.ALL.QN[,2])


# Save datasets
setwd(ResultsFile)
#write.csv(WGL.ALL.QN, "CRISPRScreens_BetaScores_WGL_QuantileNorm_Jun2024.csv")
#WGL.ALL.QN<- read.csv("CRISPRScreens_BetaScores_WGL_QuantileNorm_Jun2024.csv")
#WGL.ALL.QN<- WGL.ALL.QN[,-1]
#write.csv(WGL.ALL, "CRISPRScreens_BetaScores_WGL_Raw.csv")
# raw data and normalized WGL data


# Remove negative controls then look at top hits
WGL.ALL.QN_noNeg<- WGL.ALL.QN[-grep("neg ", WGL.ALL.QN$Gene), ] # get rid of negative controls


### Plot quantile normalized (QN) positive controls
# this time 
PlotHeatmap.Cont3<-function(DropoutData = list("HCT116 Control 1" =HCT116.Cont1.MegaA, 
                                               "HCT116 Control 2"=HCT116.Cont2.MegaA, 
                                               "HCT116 Ts8"=HCT116.Ts8.MegaA), 
                            Title= "Control guide Beta-scores",
                            Pos.Cont=PosCont, Neg.Cont=c("neg")){
  ### Step 1: Get beta scores per cell line
  Pcompare<- DropoutData[[1]]
  for (i in 2:length(DropoutData)){
    Pcompare<- merge(Pcompare, DropoutData[[i]], by.x="Gene", by.y= "Gene", all=TRUE)
  }
  
  ### Step 2: Get control guides from all cell lines
  ## Subset Positive and Negative controls. 
  Poscntvalues<- data.frame()
  Negcntvalues<- data.frame()
  for (i in Pos.Cont) {# find positive control genes
    intermediate <- Pcompare[which(Pcompare$Gene==i),] #Find exact match of poscont label in $Gene
    Poscntvalues <- rbind(Poscntvalues, intermediate)
  }
  
  for (i in Neg.Cont) {# find negative control genes
    intermediate <-Pcompare[which(grepl(i, Pcompare$Gene)),] #Find all rows containing  neg values
    Negcntvalues<- rbind(Negcntvalues, intermediate)
  }
  Allcntgds<- rbind(Poscntvalues, Negcntvalues[c(1:3),])
  
  ##Step 3: Plot
  ## plot positve and negative control beta scores
  colnames(Allcntgds)<- c("Gene", names(DropoutData))
  dfm <- melt(Allcntgds, id.vars = 1)
  dfm$Gene <- factor(dfm$Gene, level=unique(dfm$Gene)) #make Gene names a factor, not alphabetical
  
  plot.controls<- ggplot(dfm, aes(x = Gene, y = variable, fill = value)) + 
    geom_tile()+
    ggtitle(Title)+    
    xlab("Gene")+
    ylab("Cell line") +
    scale_fill_gradient2(low="Red", mid="white", high="Blue", midpoint=0, 
                         name="Beta-Score") +
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  
  plot.controls
}

PlotHeatmap.Cont3(DropoutData= list("DLD1 Control"= WGL.ALL.QN[,c(1,15)], 
                                    "DLD1 Ts13"= WGL.ALL.QN[,c(2,15)], 
                                    "HCT116 Control"= WGL.ALL.QN[,c(3,15)], 
                                    "Vaco432 Control"= WGL.ALL.QN[,c(4,15)], 
                                    "Vaco432 Ts13"= WGL.ALL.QN[,c(5,15)], 
                                    
                                    "A2780 Control"= WGL.ALL.QN[,c(6,15)], 
                                    "A2780 Disomy 1q"= WGL.ALL.QN[,c(7,15)], 
                                    "A2058 Control"= WGL.ALL.QN[,c(8,15)], 
                                    "A2058 Disomy 1q"= WGL.ALL.QN[,c(9,15)],
                                    "A2058 Disomy 7p"= WGL.ALL.QN[,c(10,15)], 
                                    "AGS Control"= WGL.ALL.QN[,c(11,15)], 
                                    "AGS Disomy 1q"= WGL.ALL.QN[,c(12,15)], 
                                    
                                    "MCF10A Control"= WGL.ALL.QN[,c(13,15)],
                                    "MCF10A Disomy 1q"= WGL.ALL.QN[,c(14,15)]), 
                  
                  Title="Control guide Beta-scores")
# plot.WGL.heatmap.controls.QN.pdf







### PRINCIPAL COMPONENT ANALYSIS (PCA):  WGL only #####

## Step 3:  make correlation matrix

Corr_matrix_WGL<-cor(WGL.ALL.QN_noNeg[,1:14])


WGL.pca<-princomp(Corr_matrix_WGL)
summary(WGL.pca)


#Importance of components:
#Comp.1    Comp.2    Comp.3     Comp.4     Comp.5     Comp.6     Comp.7     Comp.8     Comp.9    Comp.10
#Standard deviation     0.2264439 0.1393549 0.1077642 0.08585825 0.07287769 0.06195177 0.05312774 0.04830766 0.04348637 0.03809844
#Proportion of Variance 0.4677452 0.1771465 0.1059344 0.06724381 0.04844817 0.03501029 0.02574725 0.02128729 0.01725021 0.01324044
#Cumulative Proportion  0.4677452 0.6448917 0.7508261 0.81806987 0.86651804 0.90152833 0.92727557 0.94856286 0.96581307 0.97905352


## Step 4: plot PCA

WGL.pca2<-as.data.frame(WGL.pca$scores)
WGL.pca2$Aneuploid<-c(FALSE, TRUE, 
                      FALSE,
                      FALSE, TRUE, 
                      TRUE, FALSE, 
                      TRUE, FALSE, FALSE, 
                      TRUE, FALSE, 
                      TRUE, FALSE)
WGL.pca2$CellLine<-c("DLD1","DLD1",
                     "HCT116",
                     "Vaco432", "Vaco432", 
                     "A2780", "A2780",
                     "A2058", "A2058","A2058", 
                     "AGS", "AGS", 
                     "MCF10A", "MCF10A")


ggplot(WGL.pca2, aes(Comp.1, Comp.2, color = CellLine, shape = Aneuploid))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values=c("black", "Red", "Blue", "Green", "Grey50", "Purple", "Orange"))
# Plot.PCA.WGL.Comp1.Comp2_CellLine.pdf
# 1,3,4, etc. lots of cell line specific stuff

#Supplementary Figure 1B


### PLOT WGL QN: Scatterplot, Ranked #####
## Plot all three Beta Deltas: 
# WGL.ALL.QN_noNeg

WGL.ALL.QN.pvalue<- WGL.ALL.QN_noNeg


# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell (An - EU)
# Where A2780 and A2058 1qDi are "more euploid" than the WT
# negative means more toxic to aneuploid/trisomy cells: 
WGL.ALL.QN.pvalue$meanDiff<- (  #(WGL.ALL.QN.pvalue$HCT116.Ts13.Beta - WGL.ALL.QN.pvalue$HCT116.Cont.Beta) 
                                 (WGL.ALL.QN.pvalue$DLD1.Ts13.Beta - WGL.ALL.QN.pvalue$DLD1.Cont.Beta) 
                                + (WGL.ALL.QN.pvalue$Vaco432.Ts13.Beta - WGL.ALL.QN.pvalue$Vaco432.Cont.Beta) 
                                + (WGL.ALL.QN.pvalue$A2780.Cont.Beta - WGL.ALL.QN.pvalue$A2780.Di1q.Beta) 
                                + (WGL.ALL.QN.pvalue$A2058.Cont.Beta - WGL.ALL.QN.pvalue$A2058.Di1q.Beta) 
                                + (WGL.ALL.QN.pvalue$A2058.Cont.Beta - WGL.ALL.QN.pvalue$A2058.7pDi.Beta) 
                                + (WGL.ALL.QN.pvalue$AGS.Cont.Beta   - WGL.ALL.QN.pvalue$AGS.Di1q.Beta) 
                                + (WGL.ALL.QN.pvalue$MCF10A.Cont.Beta   - WGL.ALL.QN.pvalue$MCF10A.Di1q.Beta) 
)/7
WGL.ALL.QN.pvalue<- WGL.ALL.QN.pvalue[order(WGL.ALL.QN.pvalue$meanDiff),]

# get average beta score to see how toxic gene knockout is. 
WGL.ALL.QN.pvalue$MeanBeta<-( # WGL.ALL.QN.pvalue$HCT116.Ts13.Beta + WGL.ALL.QN.pvalue$HCT116.Cont.Beta 
                                WGL.ALL.QN.pvalue$DLD1.Ts13.Beta + WGL.ALL.QN.pvalue$DLD1.Cont.Beta
                               + WGL.ALL.QN.pvalue$Vaco432.Ts13.Beta + WGL.ALL.QN.pvalue$Vaco432.Cont.Beta
                               + WGL.ALL.QN.pvalue$A2780.Cont.Beta + WGL.ALL.QN.pvalue$A2780.Di1q.Beta 
                               + WGL.ALL.QN.pvalue$A2058.Cont.Beta + WGL.ALL.QN.pvalue$A2058.Di1q.Beta
                               + WGL.ALL.QN.pvalue$A2058.Cont.Beta 
                               + WGL.ALL.QN.pvalue$AGS.Cont.Beta + WGL.ALL.QN.pvalue$AGS.Di1q.Beta
                               + WGL.ALL.QN.pvalue$MCF10A.Cont.Beta   + WGL.ALL.QN.pvalue$MCF10A.Di1q.Beta
)/13


# Now add p-values to difference between all aneuploid and all euploid
WGL.ALL.QN.pvalue$An.minus.Eu.P<- NA
for (i in 1:length(WGL.ALL.QN.pvalue$Gene)){
  test<- t.test(as.numeric(WGL.ALL.QN.pvalue[i,c(2,5,6,8,8,11,13)]), 
                as.numeric(WGL.ALL.QN.pvalue[i,c(1,4,7,9,10,12,14)]), paired=TRUE) # t.test (aneuploid, euploid)
  WGL.ALL.QN.pvalue$An.minus.Eu.P[i]<-test$p.value
} #Paired t-test
WGL.ALL.QN.pvalue<-  WGL.ALL.QN.pvalue[order(WGL.ALL.QN.pvalue$An.minus.Eu.P),]


## Plot difference v pvalue (scatter plot)
WGL_TopHits<- subset(WGL.ALL.QN.pvalue, An.minus.Eu.P<0.05 & meanDiff < -0.2)$Gene
ggplot(WGL.ALL.QN.pvalue, aes(y=-log(An.minus.Eu.P,2), x=meanDiff, 
                                color=An.minus.Eu.P<0.05 & meanDiff< -0.2))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values=c("Black", "red"), name= "Significant")+ 
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.BetaDelta.pValue.WGL_QN.pdf
# Figure XX


# more toxic to trisomy: 
subset(WGL.ALL.QN.pvalue, meanDiff < -0.2 & An.minus.Eu.P<0.01)$Gene
# "WAPAL"      "CUL3"       "NUPR1L"     "EIF4A1"     "NDUFV1"     "GPR75-ASB3" "PPRC1"      "DNA2"       "TLK2"       "LIG1"      
# "MRPL43"     "LARP4B"     "TFB1M"      "IQGAP3"     "ANKRD49"    "ZBED6CL"    "CCDC144NL"  "UFL1"       "GPC5"       "PDGFA"     
# "LRRC6"      "EXO1"       "LRPPRC"     "AKNA"       "C1orf109"   "ZBTB5"      "REXO2"      "TBL1XR1"    "DDX10"      "RPA2"      
# "NDUFA1"     "CADM3"      "ADAM32"     "HOMER1"     "CHST6"      "MRPS18A"    "NPIPB6"     "ARHGAP21"   "ABHD16B"    "RGAG1"     
# "BCAN"       "SPATA5"     "TBC1D3K"   


#Less essential to trisomy:
subset(WGL.ALL.QN.pvalue, meanDiff > 0.2 & An.minus.Eu.P<0.01)$Gene
# "CTU2"     "FAM104A"  "ANKRD13B" "EGR4"     "SATB1"    "PRAMEF12" "CHRNA6"   "RASSF6"   "COG5"     "ACOXL"    "OR5H1"    "ZNF443"   "C4orf3"  
#  "TIPRL"    "GPR26"    "GFPT1"    "ZNF408"   "ZFP69"    "IL23R"    "AVEN"     "RNPEPL1"  "SNX15"    "FANCD2"  



## Plot Mitochondrial metabolism 
ggplot(WGL.ALL.QN.pvalue, aes(y=-log(An.minus.Eu.P,2), x=meanDiff))+
  geom_point()+
  #geom_point(data = subset(WGL.ALL.QN.pvalue, Gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="Purple") +
  #geom_point(data = subset(WGL.ALL.QN.pvalue, Gene %in% MitoElectronChainnAssembly), color = "#F8766D")+
  #geom_point(data = subset(WGL.ALL.QN.pvalue, Gene %in% MitoTranslationTranscription), color = "gold2")+
  geom_density_2d(color = "white", bins = 10)+
  theme_classic()+
  xlab("Difference in Beta score (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.BetaDelta.pValue.WGL_QN_Mito.ETC.pdf




# Plot RANK
WGL.ALL.QN.pvalue$AneuploidRank<- rank(WGL.ALL.QN.pvalue$meanDiff)
length(WGL.ALL.QN.pvalue$AneuploidRank)*0.1

ggplot(WGL.ALL.QN.pvalue, aes(x=AneuploidRank, y=meanDiff, color= AneuploidRank < 1912))+
  geom_point()+
  geom_point(data = subset(WGL.ALL.QN.pvalue, AneuploidRank < 1912))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  ylab("Difference in Beta")+
  xlab("WGL Sheltzer: Rank")
# 5x3
# plot.WGLSheltzer.Rank_Top10.pdf
# Figure XX





### Look up gene of interest: 
subset(WGL.ALL.QN.pvalue, Gene == "UBE2H")

HighlightTheseGenes<- c("UBE2H")

ggplot(WGL.ALL.QN.pvalue, aes(x=AneuploidRank, y=meanDiff, color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(WGL.ALL.QN.pvalue, Gene %in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  ylab("Beta delta (aneuploid - euploid)")+
  xlab("WGL Sheltzer: Rank")
# 5x3
# plot.WGLSheltzer.Rank_DHODH.pdf


## Plot gene of interest in difference v pvalue (scatter plot)
HighlightTheseGenes<- c("FOSL1", "UBE2H") 

ggplot(WGL.ALL.QN.pvalue, aes(y=-log(An.minus.Eu.P,2), x=meanDiff, color = Gene %in% HighlightTheseGenes))+
  geom_point()+ 
  geom_point(data = subset(WGL.ALL.QN.pvalue, Gene%in% HighlightTheseGenes))+
  #geom_hline(yintercept=0)+
  #geom_vline(xintercept=0)+
  geom_density_2d(color = "white", bins = 9)+
  scale_color_manual(values=c("Black", "Red"), name="Highlight")+ 
  theme_classic()+
  xlab("Difference in Beta score\n Aneuploid - Euploid")+
  ylab("-log2(p-value)")
# 5x4
# plot.WGL.BetaDelta.pValue_QN_TopHits_UBE2H_parnters.pdf


subset(WGL.ALL.QN.pvalue, Gene%in% HighlightTheseGenes)


# plot paired samples WGL data: 
GeneOfInterest<- "FOSL1"
df<- data.frame(Control=as.numeric(subset(WGL.ALL.QN.pvalue, Gene==GeneOfInterest)[,c(1,4,7,9,10,12,14)]), 
                Aneuploid=as.numeric(subset(WGL.ALL.QN.pvalue, Gene==GeneOfInterest)[,c(2,5,6,8,8,11,13)]) )
ggpaired(df, cond1 = "Control", cond2 = "Aneuploid", 
         title=GeneOfInterest, fill="condition")+ 
  ylab("Quantile Norm. Beta")+
  scale_fill_manual(values=c("darkolivegreen1", "skyblue"))
# 4x4 
# plot.paired.WGL_FOSL1.pdf 

t.test(x=as.numeric(subset(WGL.ALL.QN.pvalue, Gene==GeneOfInterest)[,c(1,4,7,9,10,12,14)]), 
       y= as.numeric(subset(WGL.ALL.QN.pvalue, Gene==GeneOfInterest)[,c(2,5,6,8,8,11,13)]) ,
       paired=TRUE )

subset(WGL.ALL.QN.pvalue, Gene==GeneOfInterest)





## Write datasets: 
setwd(ResultsFile)
#write.csv(WGL.ALL.QN.pvalue, "Sheltzer_WGL.ALL.QuantileNorm.pvalue.difference.csv")
# WGL.ALL.QN.pvalue<- read.csv("Sheltzer_WGL.ALL.QuantileNorm.pvalue.difference.csv")
#    WGL.ALL.QN.pvalue<- WGL.ALL.QN.pvalue[,-1]
# Dataset XX


## Write top 10% ranked genes:
# write.csv(subset(WGL.ALL.QN.pvalue, AneuploidRank < 1912)$Gene, "Sheltzer.WGL.RANK.top10.csv")
# Dataset XX



    ##### Gene Ontology: GSEA gProfiler WGL #####
# GSEA 


# I want the difference in log2 fold change 
GSEAListWGL <- WGL.ALL.QN.pvalue$meanDiff
names(GSEAListWGL) <- WGL.ALL.QN.pvalue$Gene# name the vector
GSEAListWGL<-na.omit(GSEAListWGL)# omit any NA values 
GSEAListWGL = sort(GSEAListWGL, decreasing = TRUE) # sort the list in increasing order (required for clusterProfiler)


GSEA.WGL <- gseGO(geneList=GSEAListWGL, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db", 
                 pAdjustMethod = "none")


require(DOSE)
dotplot(GSEA.WGL, showCategory=10, split=".sign") + facet_grid(.~.sign)
# plot.GSEA.WGL.AneuploidDependency.pdf 
# Sup Figure 2A
# 10x10 

GSEA.WGL$Description[1:40]
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
# "Suppressed" and genes at bottom of list, negative beta delta, more essential to aneuploid cells: 
# Pathways more essential to aneuploid cells (more toxic when knocked out)

gseaplot(GSEA.WGL, title = GSEA.WGL$Description[9], geneSetID = 9) # nucleobase-containing small molecule metabolic process 
# plot.GSEA.WGL.SmallMetabolicProcess.pdf
gseaplot(GSEA.WGL, title = GSEA.WGL$Description[3], geneSetID = 3) # Mitochondrial envelope
# plot.GSEA.WGL.MitoEnvelope.pdf
gseaplot(GSEA.WGL, title = GSEA.WGL$Description[4], geneSetID = 4) # Mitochondrial membrane
# plot.GSEA.WGL.MitoMembrane.pdf
gseaplot(GSEA.WGL, title = GSEA.WGL$Description[1], geneSetID = 1) # peptide metabolic process
# plot.GSEA.WGL.peptideMetabolicProcess.pdf
gseaplot(GSEA.WGL, title = GSEA.WGL$Description[5], geneSetID = 5) # ribonucleoprotein complex
# plot.GSEA.WGL.Ribosome.pdf
gseaplot(GSEA.WGL, title = GSEA.WGL$Description[8], geneSetID = 8) # translation
# plot.GSEA.WGL.Translation.pdf
gseaplot(GSEA.WGL, title = GSEA.WGL$Description[21], geneSetID = 21) # ribose phosphate metabolic process
#plot.GSEA.WGL.ribosephosphatemetabolicprocess.pdf
#5x4
# Sup Figure 2B

# Write GSEA results
setwd(ResultsFile)
write.csv(GSEA.WGL, "WGL.GSEA.csv")
# Dataset e





# G-Profiler: 
g.WGL.Aneu.Toxic<- gost(
  subset(WGL.ALL.QN.pvalue, meanDiff < -0.2)$Gene,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = FALSE,
  user_threshold = 0.05,
  correction_method = c("g_SCS"),
  domain_scope = "custom",
  custom_bg = as.character(WGL.ALL.QN.pvalue$Gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)

g.WGL.Aneu.Toxic_result<-g.WGL.Aneu.Toxic$result
g.WGL.Aneu.Toxic_result<-g.WGL.Aneu.Toxic_result[order(g.WGL.Aneu.Toxic_result$p_value),]

setwd(ResultsFile)
# write.csv(g.WGL.Aneu.Toxic_result[,1:13], "gProfiler.WGL.AneuploidyToxic.Sheltzer_2.csv")


##Toxic to aneuploid cells Sheltzer
P.AnTox<-g.WGL.Aneu.Toxic_result[order(g.WGL.Aneu.Toxic_result$p_value),]
P.AnTox<- subset(P.AnTox, !source %in% c("HPA", "TF"))
P.AnTox$term_name[1:200]
#P.AnTox<-P.AnTox[c(1,2, 3,5,9,10,11,21, 4,7,12,16, 15,42,102),]
P.AnTox<-P.AnTox[c(1,5, 2,3,7,9,14,20, 6,11,16,21, 15,47,129),]
P.AnTox$term_name<- factor(P.AnTox$term_name, levels= P.AnTox$term_name)
P.AnTox$Termtype<-c( "Protein complex","Protein complex",
                     
                     "Mitochondrial Metabolism","Mitochondrial Metabolism","Mitochondrial Metabolism",
                     "Mitochondrial Metabolism","Mitochondrial Metabolism","Mitochondrial Metabolism", 
                     
                     "Ribosome","Ribosome","Ribosome","Ribosome",
                     "RNA processing", "RNA processing", "RNA processing"
                   )

ggplot(P.AnTox, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key Biological terms")+
  ggtitle ("Biological terms enriched: Aneuploid toxic (Beta-Score)")+
  theme_classic()+
  coord_flip()
## 
# plot.gProfiler.sheltzer.WGL_4.pdf
# Figure 2C
# 10x5



# Write g profiler results
setwd(ResultsFile)
g.WGL.Aneu.Toxic_result_fixed <- as.data.frame(lapply(g.WGL.Aneu.Toxic_result, function(col) {
  if (is.list(col)) {
    sapply(col, function(x) paste(x, collapse = ", "))
  } else {
    col
  }
}))
# write.csv(g.WGL.Aneu.Toxic_result_fixed, "gProfiler.Sheltzer.WGL_AneuploidPathways.csv", row.names = FALSE)







### DepMap aneuploidy Score ####

##  Find genes that are more toxic to aneuploid cells in DepMap data:

# Merge aneuploidy score and DepMap CRISPR data: 
CRISPR_Gene_AneuScore<- merge(Aneuploidy_Scores, CRISPR_Gene_all, by.x="DepMap_ID", by.y="ModelID")


#Plot a significant correlation gene, BIRC6: 
ggplot(CRISPR_Gene_AneuScore, aes(y= BIRC6..57448., 
                                  x=Score.divPloidy))+
  geom_point()+
  theme_classic()+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score: aneuploid arms")+
  ylab("CRISPR effect score")

cor.test(CRISPR_Gene_AneuScore$BIRC6..57448., CRISPR_Gene_AneuScore$Score.divPloidy, method = "pearson")
# cor =  0.2328689
# p-value= 1.408e-09





### Now calculate and plot gene beta scores with aneuploidy score

CRISPR_DepMap_AneuScore<-data.frame(GeneID= character(), 
                                    Gene = character(), 
                                    Mean.Beta = numeric(),
                                    SampleNum = numeric(),
                                    
                                    Aneu.CRISPR.Corr.P= numeric(),
                                    Aneu.CRISPR.Corr.Coef= numeric(), 
                                    Aneu.divPloidy.CRISPR.Corr.P= numeric(),
                                    Aneu.divPloidy.CRISPR.Corr.Coef= numeric()
)

for (i in 6:length(CRISPR_Gene_AneuScore)){
  x<- cor.test(CRISPR_Gene_AneuScore[,i], CRISPR_Gene_AneuScore$Score.divPloidy, method = "pearson")
  y<- cor.test(CRISPR_Gene_AneuScore[,i], CRISPR_Gene_AneuScore$Aneuploidy.score, method = "pearson")
  
  CRISPR_DepMap_AneuScore<- rbind(CRISPR_DepMap_AneuScore, data.frame(GeneID= colnames(CRISPR_Gene_AneuScore)[i],
                                                                      Gene = sub("[..].*", "", colnames(CRISPR_Gene_AneuScore)[i] ),
                                                                      Mean.Beta = mean(na.omit(CRISPR_Gene_AneuScore[,i])),
                                                                      SampleNum= length(CRISPR_Gene_AneuScore[,i][!is.na(CRISPR_Gene_AneuScore[,i])]), 
                                                                      
                                                                      Aneu.CRISPR.Corr.P= y$p.value, 
                                                                      Aneu.CRISPR.Corr.Coef= y$estimate,
                                                                      
                                                                      Aneu.divPloidy.CRISPR.Corr.P= x$p.value, 
                                                                      Aneu.divPloidy.CRISPR.Corr.Coef = x$estimate ) ) 
} 

# Length of Genes analyzed: 18443 Genes
#CRISPR_DepMap_AneuScore$Gene<- sub("[..].*", "", CRISPR_DepMap_AneuScore$GeneID)
length(CRISPR_DepMap_AneuScore$Gene)
CRISPR_DepMap_AneuScore<- CRISPR_DepMap_AneuScore[order(CRISPR_DepMap_AneuScore$Aneu.divPloidy.CRISPR.Corr.Coef),]
CRISPR_DepMap_AneuScore[1:20,]


CRISPR_DepMap_AneuScore$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_AneuScore$GeneID)){
  if (CRISPR_DepMap_AneuScore$Aneu.divPloidy.CRISPR.Corr.P[i]< (0.05/length(CRISPR_DepMap_AneuScore$Aneu.divPloidy.CRISPR.Corr.P)) && CRISPR_DepMap_AneuScore$Aneu.divPloidy.CRISPR.Corr.Coef[i]< -0.1){
    CRISPR_DepMap_AneuScore$SigDiff[i]<- "TRUE"
  }
}

subset(CRISPR_DepMap_AneuScore, SigDiff=="TRUE")
CRISPR_DepMap_AneuScore$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_AneuScore$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_AneuScore$Gene<- sub("[.]", "-", CRISPR_DepMap_AneuScore$Gene) # replace single period with -. 
head(CRISPR_DepMap_AneuScore)

length(CRISPR_DepMap_AneuScore$Gene)# 18443


ggplot(CRISPR_DepMap_AneuScore, aes(y=-log(Aneu.divPloidy.CRISPR.Corr.P,2), 
                                    x=Aneu.divPloidy.CRISPR.Corr.Coef, 
                                    color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_AneuScore, SigDiff=="TRUE"))+
  theme_classic()+
  xlab("Difference in Gene Beta Score (DepMap): \n Aneuploidy Score")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.1")
# 5x4
subset(CRISPR_DepMap_AneuScore, SigDiff=="TRUE")$Gene



# Plot DepMap RANK 
CRISPR_DepMap_AneuRank<- CRISPR_DepMap_AneuScore
CRISPR_DepMap_AneuRank$AneuploidRank<- rank(CRISPR_DepMap_AneuRank$Aneu.divPloidy.CRISPR.Corr.Coef)
length(CRISPR_DepMap_AneuRank$AneuploidRank)*0.1

ggplot(CRISPR_DepMap_AneuRank, aes(x=AneuploidRank, y=Aneu.divPloidy.CRISPR.Corr.Coef, 
                                     color= AneuploidRank < 1845))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_AneuRank, AneuploidRank < 1845))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  ylab("Correlation with aneuploidy score/ploidy")+
  xlab("DepMap: Rank")
# 5x3
# plot.DepMap.Rank_Top10.pdf
# Figure XX

# write.csv(subset(CRISPR_DepMap_AneuRank, AneuploidRank < 1845)$Gene, "DepMap.Corr.RANK.top10.csv")
# write.csv(CRISPR_DepMap_AneuRank, "DepMap.Corr.RANK.csv")
# Supplementary data ranked lists








### DepMap CRISPR data: Compare DepMap to WGL Paired #####

## Filter for genes with a minimum of 200 cell line data. to get reliable correlation only: 
table(CRISPR_DepMap_AneuScore$SampleNum)
#Count: 13    14   224   232   604   605   658   659   660 
# Number:7   505    15   129   156   426     1    97 17106 
## To look at significant, reliable correlations only, I set a minimum cell line number to >200
CRISPR_DepMap_AneuScore<- subset(CRISPR_DepMap_AneuScore, SampleNum > 200)
length(CRISPR_DepMap_AneuScore$GeneID) 
#17,931 genes after filtering



##### Correlation between WGL our data and depmap Aneuploidy data

Depmap.Sheltzer.WGL.Diff1<- merge(x= CRISPR_DepMap_AneuScore, y= WGL.ALL.QN.pvalue, by.x="Gene", by.y = "Gene") #17,237 genes

Depmap.Sheltzer.WGL.Diff1$Gene2<- Depmap.Sheltzer.WGL.Diff1$Gene

length(anti_join(CRISPR_DepMap_AneuScore, WGL.ALL.QN.pvalue)$Gene) # 896 genes in Depmap not joined
length(anti_join(WGL.ALL.QN.pvalue, CRISPR_DepMap_AneuScore)$Gene) # 2077 genes in BROAD dataset not joined, no neg genes

AntiJoinDepMapCRISPR<- anti_join(WGL.ALL.QN.pvalue, CRISPR_DepMap_AneuScore) # Find genes in my dataset that didn't merge

#Change old gene names (BROAD) to newer gene names (used in DepMap) using limma library: 

AntiJoinDepMapCRISPR$Gene2<- AntiJoinDepMapCRISPR$Gene
for (i in 1:length(AntiJoinDepMapCRISPR$Gene)){
  if (length(alias2Symbol(AntiJoinDepMapCRISPR$Gene[i], species = "Hs", expand.symbols = FALSE) )>0){ # If gene can be converted. 
    AntiJoinDepMapCRISPR$Gene2[i]<-alias2Symbol(AntiJoinDepMapCRISPR$Gene[i], species = "Hs", expand.symbols = FALSE) #convert to new gene name
  }
} 
#Is this the best way to run this code? no. but it works. It takes a while to run though... 
# 2512 genes have new names, of the 2553 genes in anti-join list. 
# Gene2 is updated gene symbol according to limmis::alias2symbol

subset(AntiJoinDepMapCRISPR, Gene2 == "SGF29")  # Test! CCDC101 should be converted to SGF29
length(AntiJoinDepMapCRISPR$Gene2) # 2077


#Join DepMap with my data, subset genes with updated gene names: 
MergeX2<- merge(x=CRISPR_DepMap_AneuScore[,1:9], y=AntiJoinDepMapCRISPR, by.x="Gene", by.y="Gene2")   # 28


# Change Gene.y to Gene2 to match previous dataset
# Gene2 should be updated (DepMap) gene names. not old (BROAD institute CRISPR sgRNA) names
colnames(MergeX2)
colnames(MergeX2)[24]<- c("Gene2")

# Now merge merged DepMap-MyScreen datasets: 
Depmap.Sheltzer.WGL.Diff<- rbind(Depmap.Sheltzer.WGL.Diff1, MergeX2)

length(Depmap.Sheltzer.WGL.Diff$Gene) #17, 880 

missingGenes<- anti_join(CRISPR_DepMap_AneuScore[,1:9], Depmap.Sheltzer.WGL.Diff, join_by("Gene"=="Gene"))$Gene
alias2Symbol(missingGenes, species = "Hs", expand.symbols = FALSE)


# Depmap.Sheltzer.WGL.Diff
# Write this dataset! 
# write.csv(Depmap.Sheltzer.WGL.Diff, "Depmap.Corr.Sheltzer.WGL.Diff.csv")
# Depmap.Sheltzer.WGL.Diff<- read_csv("Depmap.Corr.Sheltzer.WGL.Diff.csv")


### Correlation between DepMap aneuploidy correlation and Paired screened: 
cor.test(Depmap.Sheltzer.WGL.Diff$Aneu.divPloidy.CRISPR.Corr.Coef, 
    Depmap.Sheltzer.WGL.Diff$meanDiff, method = "pearson")
# correlation r= -0.02052225
# correlation p-value = 0.006065


###

HighlightTheseGenes<- c("TFAM", "UQCRC1", "UQCRC2")

# Plot depmap data only: 
ggplot(Depmap.Sheltzer.WGL.Diff, aes(x=Aneu.divPloidy.CRISPR.Corr.Coef, 
                                              y=MeanBeta))+
  geom_point() +
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% HighlightTheseGenes), color="red")+
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoTranslationTranscription ), color="gold2") + 
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  theme_classic() +
  geom_hline(yintercept =0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(color = "white", bins = 9) +
  xlab("Correlation Coefficient (Depmap):\n Aneuploidy/ploidy & CRISPR Effect Score") +
  ylab("Mean Beta Score") +
  scale_color_manual(values= c("Black", "Red"), name="Genes")
# 5x4
# Plot.CRISPR.Depmap_Coeff.diff_Mito.ETC.pdf
subset(Depmap.Sheltzer.WGL.Diff, Gene %in% HighlightTheseGenes )





# plot top hits: 
ggplot(Depmap.Sheltzer.WGL.Diff, aes(y=Aneu.divPloidy.CRISPR.Corr.Coef, x=meanDiff))+
  geom_point() +
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, # MeanBeta> -1.5 & MeanBeta< -0.1 &
                           meanDiff < -0.2  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & An.minus.Eu.P <0.1), color = "red4") + 
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, # MeanBeta> -1.5 & MeanBeta< -0.1 &
                           meanDiff < -0.2  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & An.minus.Eu.P <0.05), color = "red") + 
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene == "HMGCS1"), color = "blue") + 
  theme_classic() + 
  geom_hline(yintercept =0) + 
  geom_vline(xintercept = 0) + 
  geom_density_2d(color = "white", bins = 9) + 
  ylab("Correlation (DepMap):\n Aneuploidy & CRISPR Effect Score") +
  xlab("Difference in Beta Score (Sheltzer): \n Aneuploid - Euploid") 
# 4x4
# Plot.CRISPR.DepmapCor.ploidy_vs_SheltzerBetaDelta.WGL_TopHits.pdf
# Figure XX


#plot pathways DepMap & WGL 
ggplot(Depmap.Sheltzer.WGL.Diff, aes(y=Aneu.divPloidy.CRISPR.Corr.Coef, x=meanDiff))+
  geom_point() +
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  theme_classic() + 
  geom_hline(yintercept =0) + 
  geom_vline(xintercept = 0) + 
  geom_density_2d(color = "white", bins = 9) + 
  ylab("Correlation (DepMap):\n Aneuploidy & CRISPR Effect Score") +
  xlab("Difference in Beta Score (Sheltzer): \n Aneuploid - Euploid") 
# 4x4
# Plot.CRISPR.DepmapCor.ploidy_vs_SheltzerBetaDelta.WGL_RiboMito.pdf
# Figure 2G


#plot pathways DepMap & WGL 
ggplot(Depmap.Sheltzer.WGL.Diff, aes(y=Aneu.divPloidy.CRISPR.Corr.Coef, x=meanDiff))+
  geom_point() +
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  #geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  geom_point(data = subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  theme_classic() + 
  geom_hline(yintercept =0) + 
  geom_vline(xintercept = 0) + 
  geom_density_2d(color = "white", bins = 9) + 
  ylab("Correlation (DepMap):\n Aneuploidy & CRISPR Effect Score") +
  xlab("Difference in Beta Score (Sheltzer): \n Aneuploid - Euploid") 
# 4x4
# Plot.CRISPR.DepmapCor.ploidy_vs_SheltzerBetaDelta.WGL_ProSpli.pdf
# Figure 2G


# Depmap correlation score per pathway
# Do t-test between genes in pathway and genes not in pathway
t.test(subset(Depmap.Sheltzer.WGL.Diff, ! Gene %in% Proteosome_Genes$Approved.symbol)$Aneu.divPloidy.CRISPR.Corr.Coef, 
       subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Proteosome_Genes$Approved.symbol)$Aneu.divPloidy.CRISPR.Corr.Coef ) #5.061e-06
t.test(subset(Depmap.Sheltzer.WGL.Diff, ! Gene %in% RNA_Processing_Spliceosome)$Aneu.divPloidy.CRISPR.Corr.Coef, 
       subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_Spliceosome)$Aneu.divPloidy.CRISPR.Corr.Coef ) # 6.151e-08
x<-t.test(subset(Depmap.Sheltzer.WGL.Diff, ! Gene %in% MitoElectronChainnAssembly)$Aneu.divPloidy.CRISPR.Corr.Coef, 
          subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoElectronChainnAssembly)$Aneu.divPloidy.CRISPR.Corr.Coef ) # 2.031645e-25
x<- t.test(subset(Depmap.Sheltzer.WGL.Diff, ! Gene %in% MitoTranslationTranscription)$Aneu.divPloidy.CRISPR.Corr.Coef, 
           subset(Depmap.Sheltzer.WGL.Diff, Gene %in% MitoTranslationTranscription)$Aneu.divPloidy.CRISPR.Corr.Coef ) #1.923696e-31
x<- t.test(subset(Depmap.Sheltzer.WGL.Diff, ! Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$Aneu.divPloidy.CRISPR.Corr.Coef, 
           subset(Depmap.Sheltzer.WGL.Diff, Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$Aneu.divPloidy.CRISPR.Corr.Coef ) #1.904744e-20
x<- t.test(subset(Depmap.Sheltzer.WGL.Diff, ! Gene %in% RNA_Processing_rRNA)$Aneu.divPloidy.CRISPR.Corr.Coef, 
           subset(Depmap.Sheltzer.WGL.Diff, Gene %in% RNA_Processing_rRNA)$Aneu.divPloidy.CRISPR.Corr.Coef ) # 1.182352e-20



# WGL Sheltzer lab. paired difference 
# Do t-test between genes in pathway and genes not in pathway
t.test(subset(WGL.ALL.QN.pvalue, ! Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.pvalue, Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff ) # 0.0006921
t.test(subset(WGL.ALL.QN.pvalue, ! Gene %in% RNA_Processing_Spliceosome)$meanDiff, 
       subset(WGL.ALL.QN.pvalue, Gene %in% RNA_Processing_Spliceosome)$meanDiff ) # 0.001303
t.test(subset(WGL.ALL.QN.pvalue, ! Gene %in% MitoElectronChainnAssembly)$meanDiff, 
       subset(WGL.ALL.QN.pvalue, Gene %in% MitoElectronChainnAssembly)$meanDiff ) # 4.258e-14
t.test(subset(WGL.ALL.QN.pvalue, ! Gene %in% MitoTranslationTranscription)$meanDiff, 
       subset(WGL.ALL.QN.pvalue, Gene %in% MitoTranslationTranscription)$meanDiff ) # 1.578e-15
t.test(subset(WGL.ALL.QN.pvalue, ! Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.pvalue, Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff ) #0.002644
t.test(subset(WGL.ALL.QN.pvalue, ! Gene %in% RNA_Processing_rRNA)$meanDiff, 
       subset(WGL.ALL.QN.pvalue, Gene %in% RNA_Processing_rRNA)$meanDiff ) # 6.054e-06



subset(Depmap.Sheltzer.WGL.Diff, Gene %in% HighlightTheseGenes )$Gene



# Gene of interest
subset(Depmap.Sheltzer.WGL.Diff, Gene == "UBE2H")






    ##### 13 Trisomy specific ####
WGL.ALL.QN_noNeg

WGL.ALL.QN.Ts13<- WGL.ALL.QN_noNeg[, c(1,2,4,5,15)]
# CDC16 

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell (An - EU)
# Where A2780 and A2058 1qDi are "more euploid" than the WT
# negative means more toxic to aneuploid/trisomy cells: 
WGL.ALL.QN.Ts13$meanDiff<- ( (WGL.ALL.QN.Ts13$DLD1.Ts13.Beta - WGL.ALL.QN.Ts13$DLD1.Cont.Beta) 
                             + (WGL.ALL.QN.Ts13$Vaco432.Ts13.Beta - WGL.ALL.QN.Ts13$Vaco432.Cont.Beta) 
)/2
WGL.ALL.QN.Ts13<- WGL.ALL.QN.Ts13[order(WGL.ALL.QN.Ts13$meanDiff),]

WGL.ALL.QN.Ts13$MeanBeta<-(  WGL.ALL.QN.Ts13$DLD1.Ts13.Beta + WGL.ALL.QN.Ts13$DLD1.Cont.Beta
                             + WGL.ALL.QN.Ts13$Vaco432.Ts13.Beta + WGL.ALL.QN.Ts13$Vaco432.Cont.Beta
)/4



## Plot difference v average Beta (scatter plot)
WGL_TopHits<- subset(WGL.ALL.QN.Ts13, meanDiff < -0.4 & MeanBeta < -0.2)$Gene
ggplot(WGL.ALL.QN.Ts13, aes(y=MeanBeta, x=meanDiff, 
                            color= meanDiff< -0.4 & MeanBeta < -0.2))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values=c("Black", "red"), name= "Significant")+ 
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)")+
  ylab("Mean Beta")
# 5x4
# plot.BetaDelta.pValue.WGL_QN.pdf



# more toxic to trisomy 
subset(WGL.ALL.QN.Ts13, meanDiff< -0.6 & MeanBeta < -0.5)$Gene
#  [1] "POLR3C"         "XPO1"           "POLR2L"         "VCP"            "HMGCS1"         "FAM96B"         "EIF3D"          "INTS1"          "PSMA6"         
#[10] "C7orf26"        "NDC80"          "NXT1"           "MED12"          "POLR3B"         "PPIL2"          "POLR2J"         "RBBP5"          "DDX10"         
#[19] "FBL"            "ECD"            "SART3"          "KIF11"          "POLD2"          "SLC25A28"       "RUVBL1"         "FLII"           "VPRBP"         
#[28] "NUTF2"          "GPN2"           "UTP15"          "ORAOV1"         "RPL17-C18orf32" "AGK"            "HSPD1"          "LRPPRC"         "MCM4"          
#[37] "POLR2D"  



## Plot difference v average beta score (scatter plot)
HighlightTheseGenes<- c("UBE2H")
HighlightTheseGenes<- subset(WGL.ALL.QN.Ts13, Gene %in% GeneList13)$Gene
HighlightTheseGenes<- subset(WGL.ALL.QN.Ts13, meanDiff < -0.4)$Gene

ggplot(WGL.ALL.QN.Ts13, aes(y=MeanBeta, x=meanDiff, 
                            color=Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_density_2d(color = "white", bins = 10)+
  geom_point(data = subset(WGL.ALL.QN.Ts13, Gene%in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"), name="Highlight")+ 
  theme_classic()+
  xlab("Difference in Beta score\n Aneuploid - Euploid")+
  ylab("MeanBeta")
# 5x4
# plot.WGL.BetaDelta.pValue_QN_TopDi1qHits.pdf
subset(WGL.ALL.QN.Ts13, Gene%in% HighlightTheseGenes)



# Pathways in Ts13 WGL 
ggplot(WGL.ALL.QN.Ts13, aes(x=meanDiff))+
  geom_density(data = subset(WGL.ALL.QN.Ts13, ! Gene %in% c(Proteosome_Genes$Approved.symbol, RNA_Processing_Spliceosome, MitoElectronChainnAssembly, MitoTranslationTranscription, Ribosomal_Genes_noMT$Approved.symbol, RNA_Processing_rRNA)), color = "black" ) +
  geom_density(data = subset(WGL.ALL.QN.Ts13, Gene %in% Proteosome_Genes$Approved.symbol), color = "chartreuse1") + 
  geom_density(data = subset(WGL.ALL.QN.Ts13, Gene %in% RNA_Processing_Spliceosome), color = "deepskyblue1")+
  geom_density(data = subset(WGL.ALL.QN.Ts13, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_density(data = subset(WGL.ALL.QN.Ts13, Gene %in% MitoTranslationTranscription), color = "gold2") + 
  geom_density(data = subset(WGL.ALL.QN.Ts13, Gene %in% Ribosomal_Genes_noMT$Approved.symbol), color = "#C77CFF")+
  geom_density(data = subset(WGL.ALL.QN.Ts13, Gene %in% RNA_Processing_rRNA), color = "deepskyblue4")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  theme_classic()+
  xlab("Difference in Beta score\n Ts13 - Di13")+
  ylab("Density")
# 4x4
# plot.WGL.Ts13.Density_Pathways.pdf
# Sup Figure XX

t.test(subset(WGL.ALL.QN.Ts13, !Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.Ts13, Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff) # 2.137e-07
t.test(subset(WGL.ALL.QN.Ts13, !Gene %in% RNA_Processing_Spliceosome)$meanDiff, 
       subset(WGL.ALL.QN.Ts13, Gene %in% RNA_Processing_Spliceosome)$meanDiff) # 9.114e-11
t.test(subset(WGL.ALL.QN.Ts13, !Gene %in% MitoElectronChainnAssembly)$meanDiff, 
       subset(WGL.ALL.QN.Ts13, Gene %in% MitoElectronChainnAssembly)$meanDiff) # 0.0002499643
t.test(subset(WGL.ALL.QN.Ts13, !Gene %in% MitoTranslationTranscription)$meanDiff, 
       subset(WGL.ALL.QN.Ts13, Gene %in% MitoTranslationTranscription)$meanDiff) # 3.652714e-08
t.test(subset(WGL.ALL.QN.Ts13, !Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.Ts13, Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff) # 0.0001354
t.test(subset(WGL.ALL.QN.Ts13, !Gene %in% RNA_Processing_rRNA)$meanDiff, 
       subset(WGL.ALL.QN.Ts13, Gene %in% RNA_Processing_rRNA)$meanDiff) # 2.544e-06 




subset(WGL.ALL.QN.Ts13, Gene %in% GeneList13 & meanDiff< -0.4)
# "INTS6" "PCID2" "USPL1" "CUL4A" "IPO5" 
# PCID2 - also significant in depmap, more Ts13 toxic in APC mutants
#           RNA processing: component of TREX2 transcription and export complex 2. 
#           regulates mRNA export. Regulates MAD2 (checkpoint and CIN)! Binds and stabilizes BRCA2 (DNA damage) 


ggplot(subset(WGL.ALL.QN.Ts13, ! Gene%in% GeneList13), aes(x=meanDiff))+
  geom_density()+
  geom_density(data = subset(WGL.ALL.QN.Ts13, Gene%in% GeneList13), color = "Red")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  scale_color_manual(values=c("Black", "Red"), name="Highlight")+ 
  theme_classic()
# 5x4
# plot.Density.Ts13.Screen.pdf



# "NXT1"   "HMGCS1" "DDX10"  "QRSL1"  "ATL2"   "GPN2"  PCID2
GeneOfInterest<- "PCID2"
df<- data.frame(Control=as.numeric(subset(WGL.ALL.QN.Ts13, Gene==GeneOfInterest)[,c(1,3)]), 
                Aneuploid=as.numeric(subset(WGL.ALL.QN.Ts13, Gene==GeneOfInterest)[,c(2,4)]) )
ggpaired(df, cond1 = "Control", cond2 = "Aneuploid", 
         title=GeneOfInterest, fill="condition")+ 
  ylab("Quantile Norm. Beta")+
  scale_fill_manual(values=c("darkolivegreen1", "skyblue"))
# 4x4
# plot.paired.WGL_PCID2.pdf




WGL_TopHits<- c("APC", "IRS2", "CTNNB1", "PCID2") 
WGL.ALL.QN.Ts13.melt<- melt(subset(WGL.ALL.QN.Ts13, Gene %in% WGL_TopHits)[,1:5]) 
ggplot(WGL.ALL.QN.Ts13.melt, aes(x= Gene, y=variable, fill= value))+ 
  geom_tile()+
  ylab("Cell line")+
  scale_fill_gradient2(low="Blue", mid="white", high="Red", midpoint=0, #limits=c(-2, 2),
                       name="Beta Score") +
  theme_classic()
# plot.Tri13_APC.IRS2.pdf


write.csv(WGL.ALL.QN.Ts13, "CRISPRScreens_WGL_13.csv")


    ##### 7p Disomy Specific #### 


WGL.ALL.QN.7pDi<- WGL.ALL.QN[, c(8,10,15)]


# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell (An - EU)
# Where A2780 and A2058 1qDi are "more euploid" than the WT
# negative means more toxic to aneuploid/trisomy cells: 
WGL.ALL.QN.7pDi$meanDiff<- (  (WGL.ALL.QN.7pDi$A2058.Cont.Beta - WGL.ALL.QN.7pDi$A2058.7pDi.Beta) 
)/1
WGL.ALL.QN.7pDi$meanBeta<- (  (WGL.ALL.QN.7pDi$A2058.Cont.Beta + WGL.ALL.QN.7pDi$A2058.7pDi.Beta) 
)/2

WGL.ALL.QN.7pDi<- WGL.ALL.QN.7pDi[order(WGL.ALL.QN.7pDi$meanDiff),]


## Plot difference v mean beta (scatter plot)
WGL_TopHits<- subset(WGL.ALL.QN.7pDi, meanBeta< -0.2 & meanDiff < -0.4)$Gene
ggplot(WGL.ALL.QN.7pDi, aes(y=meanBeta, x=meanDiff, 
                            color=meanBeta < -0.2 & meanDiff< -0.4))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values=c("Black", "red"), name= "Trisomy7p Toxic")+ 
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)\n <- 7p Trisomy toxic | 7p Disomy toxic ->")+
  ylab("Mean Beta Score")
# 5x4
# plot.BetaDelta.MeanBeta.WGL_7pDi.pdf


# Pathways in Ts7p WGL 
ggplot(WGL.ALL.QN.7pDi, aes(x=meanDiff))+
  geom_density(data = subset(WGL.ALL.QN.7pDi, ! Gene %in% c(Proteosome_Genes$Approved.symbol, RNA_Processing_Spliceosome, MitoElectronChainnAssembly, MitoTranslationTranscription, Ribosomal_Genes_noMT$Approved.symbol, RNA_Processing_rRNA)), color = "black" ) +
  geom_density(data = subset(WGL.ALL.QN.7pDi, Gene %in% Proteosome_Genes$Approved.symbol), color = "chartreuse1") + 
  geom_density(data = subset(WGL.ALL.QN.7pDi, Gene %in% RNA_Processing_Spliceosome), color = "deepskyblue1")+
  geom_density(data = subset(WGL.ALL.QN.7pDi, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_density(data = subset(WGL.ALL.QN.7pDi, Gene %in% MitoTranslationTranscription), color = "gold2") + 
  geom_density(data = subset(WGL.ALL.QN.7pDi, Gene %in% Ribosomal_Genes_noMT$Approved.symbol), color = "#C77CFF")+
  geom_density(data = subset(WGL.ALL.QN.7pDi, Gene %in% RNA_Processing_rRNA), color = "deepskyblue4")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  theme_classic()+
  xlab("Difference in Beta score\n Ts7p - Di7p")+
  ylab("Density")
# 4x4
# plot.WGL.Ts7p.Density_Pathways.pdf

t.test(subset(WGL.ALL.QN.7pDi, !Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.7pDi, Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff) # 0.00571 LESS DEPENDENT
t.test(subset(WGL.ALL.QN.7pDi, !Gene %in% RNA_Processing_Spliceosome)$meanDiff, 
       subset(WGL.ALL.QN.7pDi, Gene %in% RNA_Processing_Spliceosome)$meanDiff) # 0.01807 LESS DEPENDENT

x<- t.test(subset(WGL.ALL.QN.7pDi, !Gene %in% MitoElectronChainnAssembly)$meanDiff, 
           subset(WGL.ALL.QN.7pDi, Gene %in% MitoElectronChainnAssembly)$meanDiff) # 6.16245e-17 ***
x<- t.test(subset(WGL.ALL.QN.7pDi, !Gene %in% MitoTranslationTranscription)$meanDiff, 
           subset(WGL.ALL.QN.7pDi, Gene %in% MitoTranslationTranscription)$meanDiff) # 1.795465e-18 ***

t.test(subset(WGL.ALL.QN.7pDi, !Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.7pDi, Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff) # 0.8745 NS
t.test(subset(WGL.ALL.QN.7pDi, !Gene %in% RNA_Processing_rRNA)$meanDiff, 
       subset(WGL.ALL.QN.7pDi, Gene %in% RNA_Processing_rRNA)$meanDiff) # 0.4104 NS



# More toxic to trisomy 
subset(WGL.ALL.QN.7pDi, meanBeta< -0.2 & meanDiff < -1)$Gene


ggplot(WGL.ALL.QN.7pDi, aes(y=meanBeta, x=meanDiff, 
                            color=Gene %in% GeneList7p))+
  geom_point()+
  geom_point(data = subset(WGL.ALL.QN.7pDi, Gene %in% GeneList7p))+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values=c("Black", "red"), name= "On 7p")+ 
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)\n<- 7p Trisomy toxic | 7p Disomy toxic ->")+
  ylab("Mean Beta Score")
# 5x4
# plot.BetaDelta.MeanBeta.WGL_7pDi_Chrm7p.pdf 

subset(WGL.ALL.QN.7pDi, meanBeta< 0 & meanDiff < -.4 & Gene %in% GeneList7p)$Gene
# "TRA2A"  "CYCS"   "MALSU1" "JAZF1" 


write.csv(WGL.ALL.QN.7pDi, "CRISPRScreen_WGL_7p.csv")


    ##### 1q Disomy Specific #### 

setwd(Dependency)
Olfactory_Genes<- read.csv("group-141.csv", header= FALSE)
Olfactory_Genes<- Olfactory_Genes[-1,]
colnames(Olfactory_Genes)<-Olfactory_Genes[1,]
Olfactory_Genes<- Olfactory_Genes[-1,]
Olfactory_Genes$ChrmArm<- gsub("p.*", "p", Olfactory_Genes$Chromosome)
Olfactory_Genes$ChrmArm<- gsub("q.*", "q", Olfactory_Genes$ChrmArm)
Olfactory_Genes_on1q<-subset(Olfactory_Genes, ChrmArm == "1q")


WGL.ALL.QN.1qDi<- WGL.ALL.QN[, c(6,7,8,9,11,12,13,14,15)]
WGL.ALL.QN.1qDi<-  WGL.ALL.QN.1qDi[-grep("neg ", WGL.ALL.QN.1qDi$Gene), ]

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell (An - EU)
# Where A2780 and A2058 1qDi are "more euploid" than the WT
# negative means more toxic to aneuploid/trisomy cells: 
WGL.ALL.QN.1qDi$meanDiff<- (  (WGL.ALL.QN.1qDi$A2780.Cont.Beta - WGL.ALL.QN.1qDi$A2780.Di1q.Beta) +
                                (WGL.ALL.QN.1qDi$A2058.Cont.Beta - WGL.ALL.QN.1qDi$A2058.Di1q.Beta) +
                                (WGL.ALL.QN.1qDi$AGS.Cont.Beta - WGL.ALL.QN.1qDi$AGS.Di1q.Beta)     +
                                (WGL.ALL.QN.1qDi$MCF10A.Cont.Beta - WGL.ALL.QN.1qDi$MCF10A.Di1q.Beta)         # Gain 1q - neutral/loss 1q
)/4
WGL.ALL.QN.1qDi$Diff_TP53WT<- (  (WGL.ALL.QN.1qDi$A2780.Cont.Beta - WGL.ALL.QN.1qDi$A2780.Di1q.Beta) +
                                   (WGL.ALL.QN.1qDi$AGS.Cont.Beta - WGL.ALL.QN.1qDi$AGS.Di1q.Beta)     +
                                   (WGL.ALL.QN.1qDi$MCF10A.Cont.Beta - WGL.ALL.QN.1qDi$MCF10A.Di1q.Beta)         # Gain 1q - neutral/loss 1q
)/3
WGL.ALL.QN.1qDi$meanBeta<- ( WGL.ALL.QN.1qDi$A2780.Cont.Beta + WGL.ALL.QN.1qDi$A2780.Di1q.Beta +
                               WGL.ALL.QN.1qDi$A2058.Cont.Beta + WGL.ALL.QN.1qDi$A2058.Di1q.Beta +
                               WGL.ALL.QN.1qDi$AGS.Cont.Beta + WGL.ALL.QN.1qDi$AGS.Di1q.Beta +
                               WGL.ALL.QN.1qDi$MCF10A.Cont.Beta + WGL.ALL.QN.1qDi$MCF10A.Di1q.Beta
)/8

WGL.ALL.QN.1qDi<- WGL.ALL.QN.1qDi[order(WGL.ALL.QN.1qDi$meanDiff),]


# Now add p-values to difference between all aneuploid and all euploid
WGL.ALL.QN.1qDi$An.minus.Eu.P<- NA
for (i in 1:length(WGL.ALL.QN.1qDi$Gene)){
  test<- t.test(as.numeric(WGL.ALL.QN.1qDi[i,c(1,3,5,7)]), 
                as.numeric(WGL.ALL.QN.1qDi[i,c(2,4,6,8)]), paired=TRUE) # t.test (aneuploid, euploid)
  WGL.ALL.QN.1qDi$An.minus.Eu.P[i]<-test$p.value
} #Paired t-test




## Plot difference v p-value (scatter plot)
WGL_TopHits<- subset(WGL.ALL.QN.1qDi, meanDiff< -0.25 & An.minus.Eu.P < 0.05)$Gene
ggplot(WGL.ALL.QN.1qDi, aes(y=-log2(An.minus.Eu.P), x=meanDiff, 
                            color= An.minus.Eu.P< 0.05 & meanDiff< -0.25))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values=c("Black", "red"), name= "Trisomy1q Toxic")+ 
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)\n <- 1q Trisomy (WT) toxic | 1q Disomy toxic ->")+
  ylab("-log2(p-value)")
# 5x4
# plot.BetaDelta.Pvalue.WGL_1qDi.pdf



ggplot(WGL.ALL.QN.1qDi, aes(y=A2058.Di1q.Beta- A2058.Cont.Beta, x=Diff_TP53WT, 
                            color= (A2058.Di1q.Beta- A2058.Cont.Beta)> 0.1 & Diff_TP53WT< -0.5))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values=c("Black", "red"), name= "Trisomy1q Toxic")+ 
  theme_classic()+
  xlab("Beta delta TP53-WT (aneuploid - euploid)\n <- 1q Trisomy (WT) toxic | 1q Disomy toxic ->")+
  ylab("Beta delta TP53-Mut (aneuploid - euploid)\n")
# 5x4
# plot.BetaDelta.Pvalue.WGL_1qDi_TP53WT.pdf

subset(WGL.ALL.QN.1qDi, (A2058.Di1q.Beta- A2058.Cont.Beta)> 0.1 & Diff_TP53WT< -0.5)$Gene
subset(WGL.ALL.QN.1qDi, (A2058.Di1q.Beta- A2058.Cont.Beta)< -0.1 & Diff_TP53WT< -0.5)$Gene



subset(WGL.ALL.QN.1qDi, Gene == "MDM2")
subset(WGL.ALL.QN.1qDi, meanBeta< -0.1 & meanDiff < -0.25 & An.minus.Eu.P < 0.05)$Gene


# write.csv(WGL.ALL.QN.1qDi, "WGL.QN.1qDi.csv")

# More toxic to trisomy 
subset(WGL.ALL.QN.1qDi, meanBeta< -0.1 & meanDiff < -0.5 & An.minus.Eu.P < 0.05)$Gene


subset(WGL.ALL.QN.1qDi, meanBeta< -0.1 & meanDiff < -0.25 & An.minus.Eu.P < 0.05 & Gene %in% GeneList1q)$Gene
# "XPR1"    # on chromosome 1q! Inorganic ion transporter that mediates phosphate ion export across plasma membrane. phospate regulator
# "IQGAP3"  # on 1q. involved with cytoskeletal. Enables calmodulin binding activity and myosin VI light chain binding activity. Predicted to be involved in regulation of actin cytoskeleton organization
# "DARS2"   # mitochondrial enzyme that specifically aminoacylates aspartyl-tRNA. RNA binding



ggplot(subset(WGL.ALL.QN.1qDi, !Gene %in% GeneList1q), aes(x=meanDiff)) +
  geom_density()+
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% GeneList1q), color = "red")+
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)\n<- 1q Trisomy toxic | 1q Disomy toxic ->")+
  ylab("Density")
# 5x2
# plot.BetaDelta.WGL_1qDi_Chrm1q.pdf 
# Genes on 1q: -0.045
# Genes not on 1q: +0.0016
# Difference = 0.051
# pvalue = 1.159773e-32
# t.test(subset(WGL.ALL.QN.1qDi, Gene %in% GeneList1q)$meanDiff, subset(WGL.ALL.QN.1qDi, !Gene %in% GeneList1q)$meanDiff)


Olfactory_Genes_on1q$`Approved symbol`
GeneList1q_NotOlfactory<- GeneList1q[! GeneList1q %in% Olfactory_Genes_on1q$`Approved symbol`] # 90 olfactory genes on 1q

ggplot(subset(WGL.ALL.QN.1qDi, !Gene %in% GeneList1q), aes(x=meanDiff)) +
  geom_density()+
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% GeneList1q_NotOlfactory), color = "Yellow")+
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% Olfactory_Genes_on1q$`Approved symbol`), color = "Red")+
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)\n<- 1q Trisomy toxic | 1q Disomy toxic ->")+
  ylab("Density")+
  ggtitle("Red= olfactory on 1q \nYellow= other 1q\nBlack=genes not on 1q")

t.test(subset(WGL.ALL.QN.1qDi, !Gene %in% GeneList1q)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% Olfactory_Genes_on1q$`Approved symbol`)$meanDiff )
# p = 5.568e-05, olfactory genes on 1q have negative beta delta score. 
# mean beta delta of genes on olfactory on 1q is -0.060644994 
# genes not on 1q: 0.001460659

t.test(subset(WGL.ALL.QN.1qDi, Gene %in% GeneList1q_NotOlfactory)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% Olfactory_Genes_on1q$`Approved symbol`)$meanDiff )
# p = 0.4177, NS
# non-olfactory beta delta is -0.048
# olfactory beta delta is -0.06064499





# Pathways in Ts1q WGL 
ggplot(WGL.ALL.QN.1qDi, aes(x=meanDiff))+
  geom_density(data = subset(WGL.ALL.QN.1qDi, ! Gene %in% c(Proteosome_Genes$Approved.symbol, RNA_Processing_Spliceosome, MitoElectronChainnAssembly, MitoTranslationTranscription, Ribosomal_Genes_noMT$Approved.symbol, RNA_Processing_rRNA)), color = "black" ) +
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% Proteosome_Genes$Approved.symbol), color = "chartreuse1") + 
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% RNA_Processing_Spliceosome), color = "deepskyblue1")+
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% MitoTranslationTranscription), color = "gold2") + 
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% Ribosomal_Genes_noMT$Approved.symbol), color = "#C77CFF")+
  geom_density(data = subset(WGL.ALL.QN.1qDi, Gene %in% RNA_Processing_rRNA), color = "deepskyblue4")+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  theme_classic()+
  xlab("Difference in Beta score\n Ts1q - Di1q")+
  ylab("Density")
# 4x4
# plot.WGL.Ts1q.Density_Pathways.pdf

t.test(subset(WGL.ALL.QN.1qDi, !Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% Proteosome_Genes$Approved.symbol)$meanDiff) # 0.04679
t.test(subset(WGL.ALL.QN.1qDi, !Gene %in% RNA_Processing_Spliceosome)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% RNA_Processing_Spliceosome)$meanDiff) # 0.2411  NS
t.test(subset(WGL.ALL.QN.1qDi, !Gene %in% MitoElectronChainnAssembly)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% MitoElectronChainnAssembly)$meanDiff) # 7.625e-09
t.test(subset(WGL.ALL.QN.1qDi, !Gene %in% MitoTranslationTranscription)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% MitoTranslationTranscription)$meanDiff) # 2.889e-05
t.test(subset(WGL.ALL.QN.1qDi, !Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$meanDiff) # 0.04321
t.test(subset(WGL.ALL.QN.1qDi, !Gene %in% RNA_Processing_rRNA)$meanDiff, 
       subset(WGL.ALL.QN.1qDi, Gene %in% RNA_Processing_rRNA)$meanDiff) # 0.0003232




# Investigate hit genes
WGL_TopHits<- c("UBE2H", "FOSL1") #MCL1-KO more toxic to Disomy, TP53-KO more beneficial to Disomy, MDM4 NS
ggplot(WGL.ALL.QN.1qDi, aes(y=-log2(An.minus.Eu.P), x=meanDiff, 
                            color= Gene %in% WGL_TopHits))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  geom_point(data = subset(WGL.ALL.QN.1qDi, Gene %in% WGL_TopHits), color = "Red")+
  scale_color_manual(values=c("Black", "Red"), name= "Genes")+ 
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)\n <- 1q Trisomy (WT) toxic | 1q Disomy toxic ->")+
  ylab("-log2(p-value)")
# 5x4
# plot.BetaDelta.MeanBeta.WGL_1qDi.pdf



WGL_TopHits<- c("MCL1", "MDM4","TP53", "MDM2", "TADA3") 
WGL.ALL.QN.1qDi.melt<- melt(subset(WGL.ALL.QN.1qDi, Gene %in% WGL_TopHits)[,1:9]) 
ggplot(WGL.ALL.QN.1qDi.melt, aes(x= Gene, y=variable, fill= value))+ 
  geom_tile()+
  ylab("Cell line")+
  scale_fill_gradient2(low="Red", mid="white", high="Blue", midpoint=0, limits=c(-2, 2),
                       name="Beta Score") +
  theme_classic()
# plot.1qDi_TP53.MCL1.MDM4.pdf
# Conclusion: A2780 1qDi benefits from TP53 and 1qDi suffers from MDM4 and MCL1 KO, but A2059 and AGS are not as sensitive
# MCL1-KO is more detrimental to Disomy, even in TP53-Mut cells (like A2058) 

setwd(ResultsFile)
# write.csv(WGL.ALL.QN.1qDi, "Schukken.WGL.1q.Normalized.Diff.P.csv")
# WGL.ALL.QN.1qDi <- read.csv("Schukken.WGL.1q.Normalized.Diff.P.csv")




    ##### Compare WGL 1q Disomy to DepMap 1q-Gain Specific hits #####

CN1q<- Arm_CN[, c("X", "X1q")]
colnames(CN1q)<- c("Cell_Line", "1q")
Cell_Neutral1q<- subset(CN1q, `1q`=="0")$Cell_Line
Cell_Gain1q<- subset(CN1q, `1q`=="1")$Cell_Line
Cell_Loss1q<- subset(CN1q, `1q`=="-1")$Cell_Line
Cell_Loss.Neutral1q<- subset(CN1q, `1q`!="1")$Cell_Line


CRISPR_DepMap_Ts1q<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch1qGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain1q) #get drug result for all cells with gain for chrm 8
  Drug_Ch1qNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral1q) #get drug result for all cells neutral for chrm 8
  Drug_Ch1qLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral1q) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch1qGain[,2])) >5 & 
      length(na.omit(Drug_Ch1qNeutral[,2])) >5 &
      length(na.omit(Drug_Ch1qLoss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch1qGain[,2], Drug_Ch1qNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch1qGain[,2], Drug_Ch1qLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts1q<- rbind(CRISPR_DepMap_Ts1q, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 

#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts1q<- CRISPR_DepMap_Ts1q[order(CRISPR_DepMap_Ts1q$Gain.Neutral.P),]
CRISPR_DepMap_Ts1q[1:20,]


CRISPR_DepMap_Ts1q$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts1q$GeneID)){
  if (CRISPR_DepMap_Ts1q$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts1q$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts1q$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts1q, SigDiff=="TRUE")
CRISPR_DepMap_Ts1q$Gene<- sub("[..].*", "", CRISPR_DepMap_Ts1q$GeneID)



ggplot(CRISPR_DepMap_Ts1q, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts1q, SigDiff=="TRUE"))+
  theme_classic()+
  geom_density_2d(color = "white")+
  xlab("Difference in Gene Beta Score (DepMap): \nGain 1q minus Neutral/Loss 1q")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "red"), name="P<0.05, Diff < -0.1")
# 5x4 
# plot.CRISPR.Gain1q.LossNeutral.Pvalue.Diff.pdf 
# Minimum 5 samples per condition (gain/loss)


subset(CRISPR_DepMap_Ts1q, Gain.LossNeutral.P < 0.05 & Gain.min.LossNeutral.Diff < -0.1)$Gene
# [1] "OR2T27"   "TDP2"     "NUP50"    "H3"       "MRGBP"    "SGF29"    "PPP1R15B" "HNRNPH1"  "SUMO2"    "MCL1"     "PPIL1"    "TINF2"    "NXT1"     "LCE2C"   "EIF4A1"  


# Top 1q hits less essential to Ts1q cells: 
subset(CRISPR_DepMap_Ts1q, Gain.LossNeutral.P < 0.05 & Gain.min.LossNeutral.Diff > 0.1)$Gene
# [1] "PMVK"     "COPE"     "SREBF1"   "FEN1"     "GGPS1"    "RPL12"    "PCBP1"    "DPAGT1"   "VPS25"    "KRAS"     "ATP6V1B2" "PPP1R12A" "SCAP"    


## Notes about PIK3CA: PIK3CA is more toxic to Ts1q cells if the cell line has PIK3CA- Gain of Function mutation
subset(CRISPR_DepMap_Ts1q, Gene == "PIK3CA")




## Now merge my WGL. CRISPR screen data with Ts1q DepMap data ###
DepMap.Sheltzer.WGL.CRISPR_Ts1q<- merge(WGL.ALL.QN.1qDi,  CRISPR_DepMap_Ts1q, by.x="Gene", by.y="Gene") 





ggplot(DepMap.Sheltzer.WGL.CRISPR_Ts1q, aes(x=meanDiff, y=Gain.min.LossNeutral.Diff))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts1q, 
                           meanDiff < -0.2 & Gain.min.LossNeutral.Diff < -0.1), color="red")+
  #geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts1q, Gene == "EIF4A1"), color="Blue")+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(color = "white")+
  ylab("Difference in Gene Beta Score (DepMap): \nGain 1q minus Neutral/Loss 1q")+
  xlab("Difference in Gene Beta Score (Sheltzer): \nGain 1q minus Loss 1q")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.WGL.CRISPR_Ts1q.pdf 
# Figure XA


subset(DepMap.Sheltzer.WGL.CRISPR_Ts1q, meanDiff < -0.2 & Gain.min.LossNeutral.Diff < -0.1)#$Gene
# "EIF4A1" "NXT1"   "PPIL1" 



setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.WGL.CRISPR_Ts1q, "Chrm1q_Specific_Toxicity_WGL.Sheltzer.DepMap.csv")




    ##### Compare WGL 7p Disomy to DepMap 7p-Gain Specific hits #####

CN7p<- Arm_CN[, c("X", "X7p")]
colnames(CN7p)<- c("Cell_Line", "7p")
Cell_Neutral7p<- subset(CN7p, `7p`=="0")$Cell_Line
Cell_Gain7p<- subset(CN7p, `7p`=="1")$Cell_Line
Cell_Loss7p<- subset(CN7p, `7p`=="-1")$Cell_Line
Cell_Loss.Neutral7p<- subset(CN7p, `7p`!="1")$Cell_Line


CRISPR_DepMap_Ts7p<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch7pGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain7p) #get drug result for all cells with gain for chrm 8
  Drug_Ch7pNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral7p) #get drug result for all cells neutral for chrm 8
  Drug_Ch7pLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral7p) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch7pGain[,2])) >5 & 
      length(na.omit(Drug_Ch7pNeutral[,2])) >5 &
      length(na.omit(Drug_Ch7pLoss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch7pGain[,2], Drug_Ch7pNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch7pGain[,2], Drug_Ch7pLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts7p<- rbind(CRISPR_DepMap_Ts7p, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 

#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts7p<- CRISPR_DepMap_Ts7p[order(CRISPR_DepMap_Ts7p$Gain.Neutral.P),]
CRISPR_DepMap_Ts7p[1:20,]


CRISPR_DepMap_Ts7p$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts7p$GeneID)){
  if (CRISPR_DepMap_Ts7p$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts7p$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts7p$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts7p, SigDiff=="TRUE")
CRISPR_DepMap_Ts7p$Gene<- sub("[..].*", "", CRISPR_DepMap_Ts7p$GeneID)



ggplot(CRISPR_DepMap_Ts7p, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts7p, SigDiff=="TRUE"))+
  theme_classic()+
  geom_density_2d(color = "white")+
  xlab("Difference in Gene Beta Score (DepMap): \nGain 7p minus Neutral/Loss 7p")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.1")
# 5x4 
# plot.CRISPR.Gain7p.LossNeutral.Pvalue.Diff.pdf 


subset(CRISPR_DepMap_Ts7p, Gain.LossNeutral.P < 0.05 & Gain.min.LossNeutral.Diff < -0.09)$Gene
# Depmap top 7p dependencies: 
#[1] "UBR4"     "KCMF1"    "ACTB"     "PPP1R12A" "BIRC6"    "WDR1"     "RIC1"     "HAUS1"    "PSMD14"   "RAC1"     "GRB2"     "SNRPB2"   "NHLRC2"   "INTS6"   
#[15] "PPP1R15B" "WDR36"    "SBDS"     "KIF18A"   "THAP1"    "H2BC15"  

subset(CRISPR_DepMap_Ts7p, Gene == "UBR4")


## Now merge my WGL. CRISPR screen data with Ts7p DepMap data ###
DepMap.Sheltzer.WGL.CRISPR_Ts7p<- merge(WGL.ALL.QN.7pDi,  CRISPR_DepMap_Ts7p, by.x="Gene", by.y="Gene") 

ggplot(DepMap.Sheltzer.WGL.CRISPR_Ts7p, aes(x=(A2058.Cont.Beta - A2058.7pDi.Beta), y=Gain.min.LossNeutral.Diff, 
                                            color= (A2058.Cont.Beta - A2058.7pDi.Beta) < -0.2 & Gain.min.LossNeutral.Diff < -0.1))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts7p, (A2058.Cont.Beta - A2058.7pDi.Beta) < -0.2 & Gain.min.LossNeutral.Diff < -0.1))+
  #geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts7p, Gene == "GRB2"), color ="Blue" )+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(color = "white")+
  ylab("Difference in Gene Beta Score (DepMap): \nGain 7p minus Neutral/Loss 7p")+
  xlab("Difference in Gene Beta Score (Sheltzer): \nGain 7p minus Loss 7p")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.WGL.CRISPR_Ts7p.pdf 

subset(DepMap.Sheltzer.WGL.CRISPR_Ts7p, Gene == "TAF5L")

# potential hits more essential to 7p trisomy cells: 
subset(DepMap.Sheltzer.WGL.CRISPR_Ts7p, (A2058.Cont.Beta - A2058.7pDi.Beta) < -0.2 & Gain.min.LossNeutral.Diff < -0.1) #$Gene
# "ACTB"   "GRB2"   "INTS6"  "SBDS"   "SNRPB2"

# of these, only ACTB is on chrm 7p. Note: EGFR is also on Chrm 7p
# SBDS is on 7q


### Look at genes located ON chromosome 7: 
Chrm7pGenes<- unique(subset(Protein_Info3, ChrmArm=="7p")$`Approved Symbol`)
subset(DepMap.Sheltzer.WGL.CRISPR_Ts7p, Gene %in% Chrm7pGenes)



HighlightTheseGenes<- Chrm7pGenes

ggplot(DepMap.Sheltzer.WGL.CRISPR_Ts7p, aes(y=Gain.min.LossNeutral.Diff, 
                                            x=(A2058.Cont.Beta - A2058.7pDi.Beta),  
                                            color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts7p, Gene %in% HighlightTheseGenes))+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ylab("Difference in Gene Beta Score (DepMap): \nGain 7p minus Neutral/Loss 7p")+
  xlab("Difference in Gene Beta Score (Sheltzer): \nGain 7p minus Neutral 7p, WGL")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.WGL.CRISPR_Ts7p_May2024.pdf 

# Note: ACTB is located on chromosome 7p 

subset(DepMap.Sheltzer.WGL.CRISPR_Ts7p, Gene %in% HighlightTheseGenes &  (A2058.Cont.Beta - A2058.7pDi.Beta) < -0.2 & Gain.min.LossNeutral.Diff < -0.1)


# Top hits for essential to Ts7p cells in Depmap and in A2058 
# "ACTB"    "C3orf38" "CAPZB"   "GRB2"    "INTS6"   "LIN52"   "NUF2"    
# "PSMB7"   "RPL18"   "RPRD1B"  "RUVBL2"  "SNRPB2"  "SOX10"   "SRP68"  


setwd(ResultsFile)
write.csv(DepMap.Sheltzer.WGL.CRISPR_Ts7p, "Chrm7p_Specific_Toxicity_WGL.Sheltzer.DepMap.csv")

# Depmap only chromosome 7 hits: ACTB (7p)  and SNRPB2 (7p)




    ##### Compare WGL 13 Disomy to DepMap 13-Gain Specific hits #####

CN13q<- Arm_CN[, c("X", "X13q")]
colnames(CN13q)<- c("Cell_Line", "13q")
Cell_Neutral13q<- subset(CN13q, `13q`=="0")$Cell_Line
Cell_Gain13q<- subset(CN13q, `13q`=="1")$Cell_Line
Cell_Loss13q<- subset(CN13q, `13q`=="-1")$Cell_Line
Cell_Loss.Neutral13q<- subset(CN13q, `13q`!="1")$Cell_Line


CRISPR_DepMap_Ts13<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch13Gain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain13q) #get drug result for all cells with gain for chrm 8
  Drug_Ch13Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral13q) #get drug result for all cells neutral for chrm 8
  Drug_Ch13Loss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral13q) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch13Gain[,2])) >5 & 
      length(na.omit(Drug_Ch13Neutral[,2])) >5 &
      length(na.omit(Drug_Ch13Loss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch13Gain[,2], Drug_Ch13Neutral[,2], "two.sided")
    test2<-t.test(Drug_Ch13Gain[,2], Drug_Ch13Loss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts13<- rbind(CRISPR_DepMap_Ts13, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 

#length of Genes analyzed:  Genes
CRISPR_DepMap_Ts13<- CRISPR_DepMap_Ts13[order(CRISPR_DepMap_Ts13$Gain.Neutral.P),]
CRISPR_DepMap_Ts13[1:20,]


CRISPR_DepMap_Ts13$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts13$GeneID)){
  if (CRISPR_DepMap_Ts13$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts13$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts13$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts13, SigDiff=="TRUE")
CRISPR_DepMap_Ts13$Gene<- sub("[..].*", "", CRISPR_DepMap_Ts13$GeneID)



ggplot(CRISPR_DepMap_Ts13, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts13, SigDiff == TRUE), color = "Red")+
  theme_classic()+
  geom_density_2d(color = "white")+
  xlab("Difference in Gene Beta Score (DepMap): \nGain 13 minus Neutral/Loss 13")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.1")
# 5x4 
# plot.CRISPR.Gain13.LossNeutral.Pvalue.Diff.pdf 


subset(CRISPR_DepMap_Ts13, Gain.LossNeutral.P < 0.05 & Gain.min.LossNeutral.Diff < -0.1)$Gene
# [1] "IRS2"    "MYBL2"   "FURIN"   "AP2M1"   "BCL2L1"  "DBF4"    "FBXW11"  "SAMD4B"  "CHD4"    "RAB10"   "IGF1R"   "FDPS"    "TFRC"   
#[14] "KRTAP21" "TCF7L2"  "RBM8A"   "AP2S1"   "THOC5"   "RABIF"   "ACLY"    "CYCS"    "IFITM3"  "EFL1"    "TNPO3"   "UPF3A"   "YRDC"   
#[27] "PFAS"    "SOX9"    "DTYMK"   "MYB"     "HMGB1"   "FABP5"   "CTNNB1"  "CNOT3"   "SETDB1"  "ATP6AP1" "DUT"     "MED23"   "ADSL"   
#[40] "ADSS2"   "NXF1"    "SDHB"    "CDC37"   "KLF5"    "DNM2"    "MASTL"   "MMS22L"  "NMT1"    "COX7C"   "H2AZ1"   "GART"    "CTPS1"  
#[53] "SGF29"   "SDHD"    "UROD"    "GMPS"    "SOD2"    "PAICS"   "ATP6AP2" "CHORDC1" "PPAT"    "UMPS"    "SDHC"    "MED12"   "NAMPT"  
#[66] "TINF2"  


subset(CRISPR_DepMap_Ts13, Gene == "IRS2")

# Things less toxic to Ts13 cells: (DepMap) 
subset(CRISPR_DepMap_Ts13, Gain.LossNeutral.P < 0.001 & Gain.min.LossNeutral.Diff > 0.08)$Gene
#      "GOLT1B"   "SLC35B2"  "RPP25L"   "DIAPH3"   "EXT2"     "MPHOSPH8" "SOCS3"    "TRMT6"    "TIMM17A"  "DCTN5"    "DCTN6"    "UBIAD1"   "TBCA"    
#[14] "GTF3A"    "POLR1C"   "MZT1"     "ESCO2"    "NOL12"    "WDR3"     "POMP"     "EXT1"     "EXOSC5"   "PTK2"     "PPP2R1A"  "CENPJ"    "ALG5"    
#[27] "CCDC6"    "TRAPPC13" "SMG6"     "PKN2"     "FERMT2"   "WWTR1"    "NUP35"    "NUP58"    "ZC3H13"  




## Now merge my WGL. CRISPR screen data with Ts13q DepMap data ###
DepMap.Sheltzer.WGL.CRISPR_Ts13<- merge(WGL.ALL.QN.Ts13,  CRISPR_DepMap_Ts13, by.x="Gene", by.y="Gene") 
DepMap.Sheltzer.WGL.CRISPR_Ts13$Diff13<- meanDiff<- ((DepMap.Sheltzer.WGL.CRISPR_Ts13$DLD1.Ts13.Beta - DepMap.Sheltzer.WGL.CRISPR_Ts13$DLD1.Cont.Beta) +
  (DepMap.Sheltzer.WGL.CRISPR_Ts13$Vaco432.Ts13.Beta -DepMap.Sheltzer.WGL.CRISPR_Ts13$Vaco432.Cont.Beta) )/2

ggplot(DepMap.Sheltzer.WGL.CRISPR_Ts13, aes(x=Diff13, y=Gain.min.LossNeutral.Diff, 
                                            color= Diff13 < -0.2 & Gain.min.LossNeutral.Diff < -0.1))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts13, Diff13 < -0.2 & Gain.min.LossNeutral.Diff < -0.1), color = "Red")+
  #geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts13, Gene == "KLF5"), color = "Blue")+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(color = "white")+
  ylab("Difference in Gene Beta Score (DepMap): \nGain 13 minus Neutral/Loss 13")+
  xlab("Difference in Gene Beta Score (Sheltzer): \nGain 13 minus Loss 13")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.WGL.CRISPR_Ts13.pdf 

subset(DepMap.Sheltzer.WGL.CRISPR_Ts13, Diff13 < -0.2 & Gain.min.LossNeutral.Diff < -0.1)$Gene
#"ACLY"     "ADSL"     "AP2S1"    "ATP6AP1"  "ATP6V0C"  "ATP6V0D1" "CFL1"     "CTNNB1"   "CYCS"     "DNM2"     "DUT"      "EFR3A"    "GINS2"    "HSPA8"   
#[15] "MASTL"    "MYBL2"    "NMT1"     "NRBP1"    "NSF"      "RAB10"    "SAMD4B"   "SOX9"     "SUMO2"    "TCF7L2"   "TTC7A"    "TUBG1" 


subset(DepMap.Sheltzer.WGL.CRISPR_Ts13, Gene == "IRS2")




### Look at genes located ON chromosome 13: 
Chrm13Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 13)$`HGNC symbol`)
subset(DepMap.Sheltzer.WGL.CRISPR_Ts13, Gene %in% Chrm13Genes)

Chrm20Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 20)$`HGNC symbol`)

HighlightTheseGenes<- Chrm13Genes


ggplot(DepMap.Sheltzer.WGL.CRISPR_Ts13, aes(y=Gain.min.LossNeutral.Diff, x= Diff13,  
                                            color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL.CRISPR_Ts13, Gene %in% HighlightTheseGenes), color="red")+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Difference in Gene Beta Score (DepMap): \nGain 13 minus Neutral/Loss 13")+
  ylab("Difference in Gene Beta Score (Sheltzer): \nGain 13 minus disomy 13, WGL")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.WGL.CRISPR_Ts13.pdf 
# Figure XX

# Top hits on chromosome 13: IRS2 (mostly DepMap) and INTS6 
# genes on chrm 20 more essential to Ts13 depmap cells: BCL2L1 and MYBL2


subset(DepMap.Sheltzer.WGL.CRISPR_Ts13, Gene %in% HighlightTheseGenes & Gain.min.LossNeutral.Diff < -0.1 )

setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.WGL.CRISPR_Ts13, "Chrm13_Specific_Toxicity_WGL.Sheltzer.DepMap.csv")



### Correlate with Aneuploidy score (WGL only) ####
# Pearson correlation between gene dropout score and aneuploidy score per cell line

colnames(WGL.ALL.QN) 
Merge.WGL.QN.Corr<- WGL.ALL.QN[,1:15]

# Aneuploidy score as defined by number of chromosome arms aneuploid divided by mean ploidy
# Aneuploidy_Scores.Gene is defined by number of genes on aneuploid chromosome, divided by mean ploidy
# I do not count X monosomy as aneuploid

Aneuploidy_ScoresWGL <- c(0, 1/2,  
                        3/2,
                        0, 1/2, 
                        3/2, 2/2, 
                        8/2, 7/2, 7/2, 
                        6/2, 5/2, 
                        4/2, 3/2, 
                        "AneupPloidyScore")
# Background aneuploidy per cell line: (Not including X and Y copy numbers) 
# HCT116 control = 3/2 (8q, 17q, 10q)
# DLD1 control = 0
# A2780 control = 3/2
# A2058 control = 8/2
# SNU1 control = 2/4 # tetraploid
# Vaco432 control = 0 
# AGS control = 6/2
# MCF10A control = 4/2

colnames(Merge.WGL.QN.Corr)

#Aneuploidy_ScoresWGL <- c(0, 1, #number of aneuploid chromosome arms, not including underlying aneuploid arms
#                       0,
#                      0, 1, 
#                      1, 0, 
#                      1, 0, 0, 
#                      1,0,
#                      "AneupPloidyScore")


Merge.WGL.QN.Corr$MeanBeta<- NA
Merge.WGL.QN.Corr$AneuScore.WGL.Corr<- NA
Merge.WGL.QN.Corr$AneuScore.Corr.P<-NA
NumGene<- length(Merge.WGL.QN.Corr$Gene)

for (i in 1:(NumGene)){
  Corr<- cor.test(as.numeric(Aneuploidy_ScoresWGL[1:14]), as.numeric(Merge.WGL.QN.Corr[i,1:14]), method = 'pearson')
  Merge.WGL.QN.Corr$AneuScore.WGL.Corr[i] <-Corr$estimate
  Merge.WGL.QN.Corr$AneuScore.Corr.P[i] <- Corr$p.value
  Merge.WGL.QN.Corr$MeanBeta[i] <- mean(as.numeric(Merge.WGL.QN.Corr[i,1:14]))
}
Merge.WGL.QN.Corr<- Merge.WGL.QN.Corr[order(Merge.WGL.QN.Corr$AneuScore.WGL.Corr),]


Merge.WGL.QN.Corr$Gene[1:30]

subset(Merge.WGL.QN.Corr, MeanBeta < -0.2 & AneuScore.WGL.Corr < -0.2)$Gene

subset(Merge.WGL.QN.Corr, Gene == "ILK")
subset(Merge.WGL.QN.Corr, Gene == "KIF18A")

ggplot(Merge.WGL.QN.Corr, aes(AneuScore.WGL.Corr, MeanBeta)) + 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()+
  geom_point(data = subset(Merge.WGL.QN.Corr, AneuScore.WGL.Corr < -.2 & AneuScore.Corr.P < 0.05 & MeanBeta < -0.25), color = "red") + 
  #geom_point(data = subset(Merge.WGL.QN.Corr, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  #geom_point(data = subset(Merge.WGL.QN.Corr, Gene %in% MitoTranslationTranscription), color = "gold2") + 
  geom_density_2d(color = "white", bins = 10)+
  ylab("Mean Beta Score")+
  xlab("Sheltzer lab: \nAneuploid Score Correlation, WGL")+
  scale_color_manual(values = c("Black", "red"), name = "Top Hits")+
  theme_classic()
# plot.AneuScore.Corr.WGL.pdf
# 5x5 



subset(Merge.WGL.QN.Corr, AneuScore.WGL.Corr < -.5 & AneuScore.Corr.P < 0.005 & MeanBeta < -0.1)$Gene
#  [1] "CDC42"     "G3BP1"     "TTC14"     "PGM3"      "CAPRIN1"   "USP39"     "COPG1"     "LMNA"      "MITF"      "PMF1"      "ZC3H10"   
#[12] "WDHD1"     "CSHL1"     "SOX10"     "USP37"     "AIFM1"     "C20orf196" "DICER1"    "GFPT1"     "DHX8"      "DIEXF"     "PPP1CB"   
#[23] "NSF"       "MICALL1"   "ZC3HC1"    "NDC80"     "TFAP2A"    "XPO5"      "SNRPA"     "TKT"       "NBPF4"     "E2F3"      "SLC25A3"  
#[34] "PPP1R12A"  "ATP7A"     "ARCN1"     "SNRPF"     "AP3D1"     "MPI"       "SKP2"      "SPTSSA"    "SRSF3"     "SLC26A10"  "UGP2"     
#[45] "BRAF"      "UBE2Q1"    "ITPK1"     "POLDIP3"   "RMI1"      "APEX2"     "PIH1D1"    "BRPF1"     "ASUN"      "TRIT1"     "LUC7L3"   
#[56] "PEA15"     "TCEB2"     "HSP90AA1"  "DHX15"     "ORC6"      "YBX1"      "FOXM1"   



subset(Merge.WGL.QN.Corr, AneuScore.WGL.Corr > .5  & MeanBeta < -0.5 & Gene %in% MitoElectronChainnAssembly)$Gene


# top hits in DepMap and Sheltzer if setting ctrl cells to zero:  "EGFR"  "UBE2H".  "VDR"* 
subset(Merge.WGL.QN.Corr, Gene == "UBE2H")


GeneOfInterest<- c("UBC", "AneupPloidyScore")
x<- rbind(Merge.WGL.QN.Corr[,1:15], Aneuploidy_ScoresWGL)
x<- subset(x, Gene %in% GeneOfInterest)
x2<- t(x[,1:14])
x2<- data.frame(x2)

ggplot(x2, aes(y=as.numeric(as.character(x2[,1])), x=as.numeric(as.character(x2[,2]))  )) + 
  geom_point()+ 
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+
  ylab("KIF18A CRISPR score")
# 5x4
# plot.aneuScore.BetaScore.UBC_WGLAneuPloidy.pdf
# CtrlCellZero



### Now compare correlation score with categorical score: 

WGL.categoricalDiff.Correlation<- merge(Merge.WGL.QN.Corr, WGL.ALL.QN.pvalue, by.x="Gene", by.y="Gene")

ggplot(WGL.categoricalDiff.Correlation, aes(x=AneuScore.WGL.Corr, y=meanDiff))+
  geom_point()+
  geom_point(data = subset(WGL.categoricalDiff.Correlation, AneuScore.WGL.Corr < -0.1 & meanDiff < -0.2), color = "red")+
  #geom_point(data = subset(WGL.categoricalDiff.Correlation, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  #geom_point(data = subset(WGL.categoricalDiff.Correlation, Gene %in% MitoTranslationTranscription), color = "gold2") + 
  #geom_point(data = subset(WGL.categoricalDiff.Correlation, Gene %in% c("RAF1", "MAP2K1", "MAPK1")), color = "Red") + 
  theme_classic()+
  geom_density_2d(color = "white")+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("Aneuploidy Burden Correlation: \n Aneuploidy / Ploidy  & Beta Score")+
  ylab("Difference in Beta Score (Schukken et al): \n Aneuploid - Euploid")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.Sheltzer.WGL_AneuScore_betaDelta_RAF.MEK.ERK.pdf

subset(WGL.categoricalDiff.Correlation, Gene %in% c("ARAF", "BRAF", "KSR1", "KSR2", "RAF1", "MAP2K1", "MAPK1") )

subset(WGL.categoricalDiff.Correlation, AneuScore.WGL.Corr < -0.3  & Aneu.CRISPR.Corr.Coef < -0.15)$Gene





   #####  Merge my CRISPR AneuScore Correlation with DepMap aneuscore correlation #####

# Note: BROAD dataset has 20,111 genes (including neg guides, open reading frames ORF)
# DepMap dataset has 17,300 genes after filtering for minimum of 200 cell lines screened
# Joining them leads to 16,559 genes
# anti-join shows  827 genes in DepMap not in my dataset. Identify gene name changes and try to link missing genes: 
DepMap.Sheltzer.WGL_AneuScore1<- merge(CRISPR_DepMap_AneuScore, Merge.WGL.QN.Corr, by.x="Gene", by.y="Gene")   # Merge most genes

length(anti_join(CRISPR_DepMap_AneuScore, Merge.WGL.QN.Corr)$Gene) # 896 genes in Depmap not joined
length(anti_join(Merge.WGL.QN.Corr, CRISPR_DepMap_AneuScore)$Gene) # 3076 genes in BROAD dataset not joined, including neg genes 

AntiJoinDepMapCRISPR<- anti_join(Merge.WGL.QN.Corr, CRISPR_DepMap_AneuScore) # Find genes in my dataset that didn't merge
#AntiJoinDepMapCRISPR<- AntiJoinDepMapCRISPR[-grep("neg ", AntiJoinDepMapCRISPR$Gene), ] # Remove negative controls. genes = 2553

#Change old gene names (BROAD) to newer gene names (used in DepMap) using limma library: 

AntiJoinDepMapCRISPR$Gene2<- AntiJoinDepMapCRISPR$Gene
for (i in 1:length(AntiJoinDepMapCRISPR$Gene)){
  if (length(alias2Symbol(AntiJoinDepMapCRISPR$Gene[i], species = "Hs", expand.symbols = FALSE) )>0){ # If gene can be converted. 
    AntiJoinDepMapCRISPR$Gene2[i]<-alias2Symbol(AntiJoinDepMapCRISPR$Gene[i], species = "Hs", expand.symbols = FALSE) #label to new gene name
  }
} #Is this the best way to run this code? no. but it works. It takes a while to run though... 
# 2512 genes have new names, of the 2553 genes in anti-join list. 

subset(AntiJoinDepMapCRISPR, Gene2 == "SGF29")  # Test! CCDC101 should be converted to SGF29
length(AntiJoinDepMapCRISPR$Gene2) #  3076

#Join DepMap with my data, subset genes with updated gene names: 
Merge2<- merge(CRISPR_DepMap_AneuScore, AntiJoinDepMapCRISPR, by.x="Gene", by.y="Gene2")   # 27
# "Gene" is now the updated gene name! Gene.y is old name

# Change Gene.y to Gene2 to match previous dataset
# Gene2 should be updated (DepMap) gene names. not old (BROAD institute screen) CRISPR names
colnames(Merge2)[24]<-c("Gene2")

Merge2<- Merge2[,c(1:23,25:27)]

# Now merge merged DepMap-MyScreen datasets to get final dataset: 
DepMap.Sheltzer.WGL_AneuScore<- rbind(DepMap.Sheltzer.WGL_AneuScore1, Merge2)


length(DepMap.Sheltzer.WGL_AneuScore$Gene) #17,880
# "Gene" is now updated gene name
# "Gene2 is old name


missingGenes<- anti_join(CRISPR_DepMap_AneuScore, DepMap.Sheltzer.WGL_AneuScore, join_by("Gene"=="Gene"))$Gene
length(missingGenes)# 89
alias2Symbol(missingGenes, species = "Hs", expand.symbols = FALSE)



setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.WGL_AneuScore, "DepMap.Sheltzer.WGL_AneuScoreCorr.csv")
# DepMap.Sheltzer.WGL_AneuScore<- read.csv("DepMap.Sheltzer.WGL_AneuScoreCorr.csv")




ggplot(DepMap.Sheltzer.WGL_AneuScore, aes(x=AneuScore.WGL.Corr, y=Aneu.divPloidy.CRISPR.Corr.Coef))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, # MeanBeta> -1.5 & MeanBeta< -0.1 &
                           AneuScore.WGL.Corr < -0.1  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & AneuScore.Corr.P <0.1), color = "red4") + 
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, # MeanBeta> -1.5 & MeanBeta< -0.1 &
                           AneuScore.WGL.Corr < -0.1  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & AneuScore.Corr.P <0.05), color = "red") + 
  #geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% c("UBE2H")), color = "Blue") + 
  theme_classic()+
  geom_density_2d(color = "white", bins=10)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("Correlation (WGL): \n Aneuploidy Score correlation with Beta Score")+
  ylab("Correlation (DepMap): \n Aneuploidy Score correlation with Effect Score")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.Sheltzer.WGL_AneuScore_AllAneudivPloidy_3.pdf
# SupFig 4B

# top genes correlate with total aneuploid burden: 
subset(DepMap.Sheltzer.WGL_AneuScore, # MeanBeta> -1.5 & MeanBeta< -0.1 &
       AneuScore.WGL.Corr < -0.1  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
         Aneu.divPloidy.CRISPR.Corr.P< 0.005 & AneuScore.Corr.P <0.05)$Gene
#  [1] "ARHGEF7"  "BNIPL"    "CAPN6"    "CCDC105"  "CHEK1"    "CLEC1B"   "COPG1"    "CRTC2"    "CSHL1"    "CYP4A22"  "CYP4F3"   "DDR1"     "DHX15"    "DYNLRB1" 
#[15] "EIF1AX"   "EIF4A3"   "EXOC6B"   "F2RL3"    "FERMT2"   "FOXM1"    "GBF1"     "GH2"      "HAUS3"    "HCAR3"    "HOXC10"   "KIF20A"   "LLGL1"    "LMNA"    
#[29] "LRRC70"   "MARK2"    "MED8"     "MSN"      "NAT8"     "NDC80"    "NEDD1"    "NUDT21"   "NUP160"   "OLIG2"    "OR1I1"    "ORM2"     "PHAX"     "PLEKHF2" 
#[43] "PMPCA"    "POC1A"    "PPP1R12A" "PSMA3"    "PSMD3"    "PYROXD1"  "RBAK"     "RELA"     "RPAP3"    "S100A7"   "SF3A3"    "SKA1"     "SKA3"     "SMU1"    
#[57] "SNRPF"    "SRSF3"    "TEX37"    "TNRC6B"   "TNS3"     "TP53BP2"  "TPR"      "UFL1"     "UFM1"     "URI1"     "VPS37A"   "WDR26"    "WWTR1"    "YTHDC1"  
#[71] "ZMAT3"    "ZNF660"   "ZNF721"   "BCLAF3"   "CFAP298"  "MEIOC" 





# plot.ASdivPloidy.Corr.Sheltzer.DepMap.WGL.pdf

x<- subset(DepMap.Sheltzer.WGL_AneuScore,
           AneuScore.WGL.Corr < -0.1  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & AneuScore.Corr.P <0.05)$Gene
#  [1] "ARHGEF7"  "BNIPL"    "CAPN6"    "CCDC105"  "CHEK1"    "CLEC1B"   "COPG1"    "CRTC2"    "CSHL1"    "CYP4A22"  "CYP4F3"   "DDR1"     "DHX15"   
#[14] "DYNLRB1"  "EIF1AX"   "EIF4A3"   "EXOC6B"   "F2RL3"    "FERMT2"   "FOXM1"    "GBF1"     "GH2"      "HAUS3"    "HCAR3"    "HOXC10"   "KIF20A"  
#[27] "LLGL1"    "LMNA"     "LRRC70"   "MARK2"    "MED8"     "MSN"      "NAT8"     "NDC80"    "NEDD1"    "NUDT21"   "NUP160"   "OLIG2"    "OR1I1"   
#[40] "ORM2"     "PHAX"     "PLEKHF2"  "PMPCA"    "POC1A"    "PPP1R12A" "PSMA3"    "PSMD3"    "PYROXD1"  "RBAK"     "RELA"     "RPAP3"    "S100A7"  
#[53] "SF3A3"    "SKA1"     "SKA3"     "SMU1"     "SNRPF"    "SRSF3"    "TEX37"    "TNRC6B"   "TNS3"     "TP53BP2"  "TPR"      "UFL1"     "UFM1"    
#[66] "URI1"     "VPS37A"   "WDR26"    "WWTR1"    "YTHDC1"   "ZMAT3"    "ZNF660"   "ZNF721"   "BCLAF3"   "CFAP298"  "MEIOC" 

subset(DepMap.Sheltzer.WGL_AneuScore,
       AneuScore.WGL.Corr < -0.1  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
         Aneu.divPloidy.CRISPR.Corr.P< 0.005 & AneuScore.Corr.P <0.1 & ! Gene %in% x)$Gene

#[1] "ACTR1B"   "ADPGK"    "ARHGEF7"  "BAG6"     "BNIPL"    "C11orf24" "CAPN6"    "CCDC105"  "CDK12"    "CHEK1"    "CHMP6"    "CHMP7"    "CLEC1B"  
#[14] "COPB1"    "COPG1"    "CRTC2"    "CSHL1"    "CXADR"    "CYP4A22"  "CYP4F3"   "DDR1"     "DHX15"    "DNAH6"    "DYNLRB1"  "EHD1"     "EIF1AX"  
#[27] "EIF3G"    "EIF4A3"   "EXOC6B"   "F2RL3"    "FERMT2"   "FOXM1"    "GBF1"     "GGNBP2"   "GH2"      "HAUS3"    "HAUS8"    "HCAR3"    "HMGCS1"  
#[40] "HOXC10"   "KIF18A"   "KIF20A"   "LLGL1"    "LMNA"     "LRRC70"   "MARK2"    "MCM3AP"   "MED6"     "MED8"     "MSN"      "NAT8"     "NDC80"   
#[53] "NEDD1"    "NUDT21"   "NUP160"   "NUP98"    "OLIG2"    "OR1I1"    "ORM2"     "P3H2"     "PHAX"     "PLEKHF2"  "PLEKHH1"  "PMPCA"    "POC1A"   
#[66] "POLR2G"   "PPP1R12A" "PPP1R2"   "PRSS38"   "PSG5"     "PSMA3"    "PSMA6"    "PSMD3"    "PYROXD1"  "RACGAP1"  "RBAK"     "RELA"     "RHPN1"   
#[79] "RPAP3"    "S100A7"   "SF3A3"    "SKA1"     "SKA3"     "SMU1"     "SNRPF"    "SRSF3"    "TEX37"    "TIMM8A"   "TNRC6B"   "TNS3"     "TOPBP1"  
#[92] "TP53BP2"  "TPR"      "TROAP"    "TSG101"   "UBC"      "UBE2H"    "UFL1"     "UFM1"     "URI1"     "VPS37A"   "WASL"     "WDR26"    "WWTR1"   
#[105] "YTHDC1"   "ZC3H13"   "ZMAT3"    "ZNF419"   "ZNF45"    "ZNF557"   "ZNF660"   "ZNF682"   "ZNF721"   "BCLAF3"   "CFAP298"  "MEIOC"    "YJU2"    




# Plot pathways in aneuploidy score analysis: 
ggplot(DepMap.Sheltzer.WGL_AneuScore, aes(x=AneuScore.WGL.Corr, y=Aneu.divPloidy.CRISPR.Corr.Coef))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% Proteosome_Genes$Approved.symbol), color = "chartreuse1") + 
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% RNA_Processing_Spliceosome), color = "deepskyblue1")+
  theme_classic()+
  geom_density_2d(color = "white", bins=10)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("Correlation (Sheltzer Lab): \n Aneuploidy Score correlation with Beta Score")+
  ylab("Correlation (DepMap): \n Aneuploidy Score correlation with Effect Score")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.Sheltzer.WGL_AneuScore_AllAneudivPloidy_Proteome.Splice.pdf


# Plot pathways in aneuploidy score analysis: 
ggplot(DepMap.Sheltzer.WGL_AneuScore, aes(x=AneuScore.WGL.Corr, y=Aneu.divPloidy.CRISPR.Corr.Coef))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% MitoTranslationTranscription), color = "gold2") + 
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% Ribosomal_Genes_noMT$Approved.symbol), color = "#C77CFF")+
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% RNA_Processing_rRNA), color = "deepskyblue4")+
  theme_classic()+
  geom_density_2d(color = "white", bins=10)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("Correlation (Sheltzer Lab): \n Aneuploidy Score correlation with Beta Score")+
  ylab("Correlation (DepMap): \n Aneuploidy Score correlation with Effect Score")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.Sheltzer.WGL_AneuScore_AllAneudivPloidy_MitoRibo.pdf



t.test(subset(DepMap.Sheltzer.WGL_AneuScore, !Gene %in% Proteosome_Genes$Approved.symbol)$AneuScore.WGL.Corr, 
       subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% Proteosome_Genes$Approved.symbol)$AneuScore.WGL.Corr) # negative corr. 0.019
t.test(subset(DepMap.Sheltzer.WGL_AneuScore, !Gene %in% RNA_Processing_Spliceosome)$AneuScore.WGL.Corr, 
       subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% RNA_Processing_Spliceosome)$AneuScore.WGL.Corr) # negative corr. 1.052e-13
t.test(subset(DepMap.Sheltzer.WGL_AneuScore, !Gene %in% MitoElectronChainnAssembly)$AneuScore.WGL.Corr, 
       subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% MitoElectronChainnAssembly)$AneuScore.WGL.Corr) # positive corr. 2.272e-12
t.test(subset(DepMap.Sheltzer.WGL_AneuScore, !Gene %in% MitoTranslationTranscription)$AneuScore.WGL.Corr, 
       subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% MitoTranslationTranscription)$AneuScore.WGL.Corr) # positive corr. 1.03e-13
t.test(subset(DepMap.Sheltzer.WGL_AneuScore, !Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$AneuScore.WGL.Corr, 
       subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$AneuScore.WGL.Corr) # positive corr. 0.0003659
t.test(subset(DepMap.Sheltzer.WGL_AneuScore, !Gene %in% RNA_Processing_rRNA)$AneuScore.WGL.Corr, 
       subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% RNA_Processing_rRNA)$AneuScore.WGL.Corr) # negative corr. 0.02107

# plot.ASdivPloidy.Corr.Sheltzer.DepMap.WGL.pdf



# GeneOfInterest Depmap plot
x<- CRISPR_Gene_AneuScore
colnames(x)<-sub("[..].*", "", colnames(x))
x2<- x[, c("Aneuploidy", "UBE2H")]

ggplot(x2, aes(x=as.numeric(as.character(x2[,1])), 
               y=as.numeric(as.character(x2[,2]))  )) + 
  geom_point()+ 
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score (DepMap)")+
  ylab("CRISPR Beta score: UBE2H")
# 5x4
# plot.DepMap.aneuScore.Gene_ILK.pdf

subset(CRISPR_DepMap_AneuScore, Gene == "FOSL1")


HighlightTheseGenes<- c("UBE2H") 
ggplot(DepMap.Sheltzer.WGL_AneuScore, aes(x=AneuScore.WGL.Corr, 
                                          y=Aneu.divPloidy.CRISPR.Corr.Coef, 
                                          color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% HighlightTheseGenes))+
  theme_classic()+
  geom_density_2d(color = "white", bins=10)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("Correlation (Sheltzer Lab): \n Aneuploidy Score correlation with Beta Score")+
  ylab("Correlation (DepMap): \n Aneuploidy Score correlation with Effect Score")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4
# plot.DepMap.Sheltzer.MegaAWGL.Correlation.CRISPR_AneuScore.pdf
subset(DepMap.Sheltzer.WGL_AneuScore, Gene %in% HighlightTheseGenes)

x<- cor.test(DepMap.Sheltzer.WGL_AneuScore$AneuScore.WGL.Corr, DepMap.Sheltzer.WGL_AneuScore$Aneu.divPloidy.CRISPR.Corr.Coef)
x



#write.csv(DepMap.Sheltzer.WGL_AneuScore, "Sheltzer.CRISPR.WGL.correlate.AneuScore.AllAneudivPloidy.DepMapCorr_Update.csv")
#DepMap.Sheltzer.WGL_AneuScore<- read.csv("Sheltzer.CRISPR.WGL.correlate.AneuScore.AllAneudivPloidy.DepMapCorr_Update.csv")



###  Linear regression model DepMap CRISPR with aneuploidy score  ####

## Depmap Aneuploidy score. score is based on number of chromosome arms that are aneuploid. whole number

## We want to do a linear regression model to account for cell OncotreeLineage. 
# So we use an lm() model which is a linear regression model 
# test against t-distribution
# t- estimate is estimate of how extreme of a value it is. 

Aneuploidy_Scores$DepMap_ID
Cell_Info$OncotreeLineage

Aneu_CRISPR_data<- merge(x=Aneuploidy_Scores, y=CRISPR_Gene, by.x="DepMap_ID", by.y="ModelID")
Aneu_CRISPR_data$Score.divPloidy<- Aneu_CRISPR_data$Aneuploidy.score/ round(Aneu_CRISPR_data$Ploidy, digits=0)

Aneu_CRISPR_data_lm<- merge(x=Cell_Info[, c(1,6,18,30) ], y=Aneu_CRISPR_data, by.x="ModelID", by.y= "DepMap_ID")



## Linear regression for aneuploidy score , correcting for lineage: 
CRISPR.DepMap_AneuScore_lm_Lineage<-data.frame(GeneID= character(), 
                                               Gene = character(), 
                                               Count = numeric(), 
                                               Mean.CEG = numeric(),
                                               Aneu.CRISPR.lm.P= numeric(),
                                               Aneu.CRISPR.lm.Coef= numeric()
)

# Run a linear model correcting for lineage: 
for (i in 9:length(Aneu_CRISPR_data_lm)){ #Depmap CRISPR data
  model<- lm(Aneu_CRISPR_data_lm[,i] ~ Score.divPloidy + OncotreeLineage
             , data = Aneu_CRISPR_data_lm ) 
  
  CRISPR.DepMap_AneuScore_lm_Lineage<- rbind(CRISPR.DepMap_AneuScore_lm_Lineage, 
                                             data.frame(GeneID= colnames(Aneu_CRISPR_data_lm)[i],
                                                        Gene = sub("[..].*", "", colnames(Aneu_CRISPR_data_lm)[i] ),
                                                        Count = length(Aneu_CRISPR_data_lm[,i][!is.na(Aneu_CRISPR_data_lm[,i])]), 
                                                        Mean.CEF = mean(na.omit(Aneu_CRISPR_data_lm[,i])),
                                                        AneuDivPloidy.CRISPR.lm.P= summary(model)$coefficients["Score.divPloidy", "Pr(>|t|)"], 
                                                        AneuDivPloidy.CRISPR.lm.Coef= summary(model)$coefficients["Score.divPloidy", "Estimate"]) )
} 



length(CRISPR.DepMap_AneuScore_lm_Lineage$Gene)
# length of Genes analyzed: 17,931 Genes

CRISPR.DepMap_AneuScore_lm_Lineage<- CRISPR.DepMap_AneuScore_lm_Lineage[order(CRISPR.DepMap_AneuScore_lm_Lineage$AneuDivPloidy.CRISPR.lm.P),]
CRISPR.DepMap_AneuScore_lm_Lineage[1:50,]$Gene
subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "ILK")
subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "UBE2H")
subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "FOSL1")
subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "UBC")
subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "KIF18A")

subset(CRISPR.DepMap_AneuScore_lm_Lineage, 
       AneuDivPloidy.CRISPR.lm.P< 0.00015 &
         AneuDivPloidy.CRISPR.lm.Coef < -0.005)$Gene
# Top hits: TP53 and UBC
# [1] "TP53"     "ZNF207"   "PRPF31"   "SMU1"     "POLR2D"   "PSMD8"    "URI1"     
# "EIF4A3"   "PPP1R12A" "PAX8"     "NAPA"     "RBM8A"    "MIS18A"  
# [14] "BUB1B"    "ADAR"     "CCNE1"    "NUP160"   "SON"      "WBP11"    "SF3B2"    
# "POLR2I"   "CDC5L"    "PSMA3"    "UBC"      "LSM2"     "CHERP"   
# [27] "UFM1"     "KIF18A"   "EIF3A"    "IGF1R"   



## Now mark significance: 
# for comparison, TP53 has a lm coefficient of -0.011. and that's the most significant (by p-value) dependency
CRISPR.DepMap_AneuScore_lm_Lineage$SigDiff<- "NS"
for (i in 1:length(CRISPR.DepMap_AneuScore_lm_Lineage$Gene)){
  if (CRISPR.DepMap_AneuScore_lm_Lineage$AneuDivPloidy.CRISPR.lm.P[i]< 0.05/18442 && 
      CRISPR.DepMap_AneuScore_lm_Lineage$AneuDivPloidy.CRISPR.lm.Coef[i] < -0.01) {
    CRISPR.DepMap_AneuScore_lm_Lineage$SigDiff[i]<- "Significant"
  }  else if (CRISPR.DepMap_AneuScore_lm_Lineage$AneuDivPloidy.CRISPR.lm.P[i]< 0.0005 && 
              CRISPR.DepMap_AneuScore_lm_Lineage$AneuDivPloidy.CRISPR.lm.Coef[i] < -0.01) {
    CRISPR.DepMap_AneuScore_lm_Lineage$SigDiff[i]<- "Trending"
  } 
}
subset(CRISPR.DepMap_AneuScore_lm_Lineage, SigDiff=="Significant")#$Gene
subset(CRISPR.DepMap_AneuScore_lm_Lineage, SigDiff=="Trending")#$Gene

CRISPR.DepMap_AneuScore_lm_Lineage[1:20,]$Gene

setwd(ResultsFile)
#write.csv(CRISPR.DepMap_AneuScore_lm_Lineage, "CRISPR.DepMap_AneuScoreDivPloidy_lm_Lineage.csv")
#CRISPR.DepMap_AneuScore_lm_Lineage<-read.csv("CRISPR.DepMap_AneuScoreDivPloidy_lm_Lineage.csv")

#CRISPR.DepMap_AneuScore_lm_Lineage<- CRISPR.DepMap_AneuScore_lm_Lineage[-1,]

ggplot(CRISPR.DepMap_AneuScore_lm_Lineage, aes(x=AneuDivPloidy.CRISPR.lm.Coef, 
                                               y= -log2(AneuDivPloidy.CRISPR.lm.P), 
                                               color = SigDiff))+
  geom_point()+ 
  geom_density2d(color="white")+
  scale_color_manual(values=c("Black", "Red", "Red4"))+ 
  theme_classic()+
  xlab(paste("Linear model coefficient"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("DepMap Linear model: tissue type & aneuploidy score"))
# plot.LM.DepMap_AndivPloidy.pdf
# Figure XXA

# Top hits (red)
#"TP53"   "CHEK2"  "ZNF207" "PRPF31" "SMU1"   "POLR2D" "PSMD8"  "URI1"   "EIF4A3"

subset(CRISPR.DepMap_AneuScore_lm_Lineage, AneuDivPloidy.CRISPR.lm.Coef < -0.01 & AneuDivPloidy.CRISPR.lm.P< 0.0005 )

subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "FOSL1")
subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "UBE2H")
subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene == "UBC")






ggplot(CRISPR.DepMap_AneuScore_lm_Lineage, aes(x=AneuDivPloidy.CRISPR.lm.Coef, 
                                               y= -log2(AneuDivPloidy.CRISPR.lm.P)))+
  geom_point()+ 
  theme_classic()+
  geom_point(data = subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% MitoTranslationTranscription ), color="gold2") + 
  geom_point(data = subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  geom_density2d(color="white")+
  xlab(paste("Linear model coefficient"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("DepMap Linear model: tissue type & aneuploidy score"))
# plot.LM.DepMap_AndivPloidy_pathway2.pdf
# 4x4
# Figure XXC



ggplot(CRISPR.DepMap_AneuScore_lm_Lineage, aes(x=AneuDivPloidy.CRISPR.lm.Coef, 
                                               y= -log2(AneuDivPloidy.CRISPR.lm.P)))+
  geom_point()+ 
  theme_classic()+
  geom_point(data = subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  geom_point(data = subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_density2d(color="white")+
  xlab(paste("Linear model coefficient"))+
  ylab(paste("-log2(p-value)"))+
  ggtitle(paste("DepMap Linear model: tissue type & aneuploidy score"))
# plot.LM.DepMap_AndivPloidy_pathway2.pdf
# 4x4
# Figure XXB


x<- t.test(subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% MitoElectronChainnAssembly)$AneuDivPloidy.CRISPR.lm.Coef, 
           subset(CRISPR.DepMap_AneuScore_lm_Lineage, !Gene %in% MitoElectronChainnAssembly)$AneuDivPloidy.CRISPR.lm.Coef)
x$p.value
x<- t.test(subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% MitoTranslationTranscription)$AneuDivPloidy.CRISPR.lm.Coef, 
           subset(CRISPR.DepMap_AneuScore_lm_Lineage, !Gene %in% MitoTranslationTranscription)$AneuDivPloidy.CRISPR.lm.Coef)
x$p.value
x<- t.test(subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% RNA_Processing_rRNA)$AneuDivPloidy.CRISPR.lm.Coef, 
           subset(CRISPR.DepMap_AneuScore_lm_Lineage, !Gene %in% RNA_Processing_rRNA)$AneuDivPloidy.CRISPR.lm.Coef)
x$p.value
x<- t.test(subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$AneuDivPloidy.CRISPR.lm.Coef, 
           subset(CRISPR.DepMap_AneuScore_lm_Lineage, !Gene %in% Ribosomal_Genes_noMT$Approved.symbol)$AneuDivPloidy.CRISPR.lm.Coef)
x$p.value
x<- t.test(subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% Proteosome_Genes$Approved.symbol)$AneuDivPloidy.CRISPR.lm.Coef, 
           subset(CRISPR.DepMap_AneuScore_lm_Lineage, !Gene %in% Proteosome_Genes$Approved.symbol)$AneuDivPloidy.CRISPR.lm.Coef)
x$p.value
x<- t.test(subset(CRISPR.DepMap_AneuScore_lm_Lineage, Gene %in% RNA_Processing_Spliceosome)$AneuDivPloidy.CRISPR.lm.Coef, 
           subset(CRISPR.DepMap_AneuScore_lm_Lineage, !Gene %in% RNA_Processing_Spliceosome)$AneuDivPloidy.CRISPR.lm.Coef)
x$p.value


# write.csv(CRISPR.DepMap_AneuScore_lm_Lineage, "CRISPR.DepMap_AneuScore_lm_Lineage_Min200.csv")



