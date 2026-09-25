#### CRISPR Screen data analysis of MegaA data for Aneuploid cell line####
## MegaA Screen 1 and 2
##       Goal: identify aneuploid specific genetic dependencies by comparing aneuploid and euploid cells. 
##       CRISPR screen gene dropout results. 
## Klaske Schukken
## October 11, 2022



### INTRO AND LIBRARIES ####
## November 11, 2022

### What this code if for: 
## Import mle files for each Druggable Genome (DG) library screen
## Merge and normalize data 
## visualizing aneuploidy dependency for DG and DG+ WGL
## Arm specific hits. cell line specific hits. tissue specific hits
## G profiler and GSEA

## You will need to open the WGL_all_analysis file and load the WGL mle files to run parts of this script!!! ##


## Folders 
## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening and proteomics datasets: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"


### Libraries
library(ggplot2) 
library(reshape2)
library(tidyverse)
library(readr)
library(xlsx)
library(readxl)
library(ggpubr)
library("cowplot")
library(plyr)
library("viridis")  # color packet
library("Hmisc")
library("PerformanceAnalytics")
library(corrr)
library(corrplot)
library('gprofiler2')


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
### MegaA/ Druggable Genome (DG) Screen data  ##

# MegaA Screen 1 : Gene info
## Update: Get mle data where CDK12 was removed due to CDK12 contamination. 
setwd(SchukkenData)

HCT116.Cont1.MegaA <- read.delim2("mle_HCT116-Cont1.gene_summary.txt", 
                                dec=".", header = TRUE)
#HCT116.Ts13.MegaA <- read.delim2("mle_HCT116-Ts13.gene_summary.txt",  # HCT116 Ts13 is shattered and partially lost
#                              dec=".", header = TRUE)
HCT116.Ts8.MegaA <- read.delim2("mle_HCT116-Ts8.gene_summary.txt", 
                                 dec=".", header = TRUE)
DLD1.Cont1.MegaA <- read.delim2("mle_DLD1-Cont1.gene_summary.txt",
                                  dec=".", header = TRUE)
DLD1.Ts13.MegaA <- read.delim2("mle_DLD1-Ts13.gene_summary.txt",  
                                 dec=".", header = TRUE)
DLD1.Ts8.10.MegaA <- read.delim2("mle_DLD1-Ts8.gene_summary.txt", 
                                dec=".", header = TRUE)

##  CDK12 sgRNA contaminated the sample. removed CDK12. CDK12 sgRNA had >40 million reads instead of ~1000. 

### MegaA Screen 1: guide RNA info
HCT116.Cont1.MegaA.guide <- read.csv("count-HCT116-Cont1-i.f.csv")
HCT116.Ts8.MegaA.guide <- read.csv("count-HCT116-Ts8-i.f.csv") 
DLD1.Cont1.MegaA.guide <- read.csv("count-DLD1-Cont1-i.f.csv")
DLD1.Ts13.MegaA.guide <- read.csv("count-DLD1-Ts13-i.f.csv")
DLD1.Ts8.10.MegaA.guide <- read.csv("count-DLD1-Ts8-i.f.csv")
median(HCT116.Cont1.MegaA.guide$HCT116.Cont1.initial) # 247
median(HCT116.Cont1.MegaA.guide$HCT116.Cont1.final)   # 1422
median(HCT116.Ts8.MegaA.guide$HCT116.Ts8.initial)     # 396
median(HCT116.Ts8.MegaA.guide$HCT116.Ts8.final)       # 1127
median(DLD1.Cont1.MegaA.guide$DLD1.Cont1.initial)     # 1753
median(DLD1.Cont1.MegaA.guide$DLD1.Cont1.final)       # 1253
median(DLD1.Ts13.MegaA.guide$DLD1.Ts13.initial)       # 1823
median(DLD1.Ts13.MegaA.guide$DLD1.Ts13.final)         # 1601
median(DLD1.Ts8.10.MegaA.guide$DLD1.Ts8.initial)      # 1598
median(DLD1.Ts8.10.MegaA.guide$DLD1.Ts8.final)        # 1529


## MegaA Screen 2 : Gene info
# Update  January 30, 2023
#  CDK12 contamination. CDK12 is removed
setwd(SchukkenData)

HCT116.Cont2.MegaA <- read.delim2("mle_HCT116-Cont2.gene_summary.txt", 
                                  dec=".", header = TRUE)
HCT116.Tet5p.Ts5q.MegaA <- read.delim2("mle_HCT116-Ts3.gene_summary.txt", 
                                 dec=".", header = TRUE)
HCT116.Tet5p.MegaA <- read.delim2("mle_HCT116-Ts5.gene_summary.txt", 
                                dec=".", header = TRUE)
### MegaA Screen 2: guide RNA info
HCT116.Cont2.MegaA.guide <- read.csv("count-HCT116-Cont2-i.f.csv") # beautiful curve, ~700 reads/guide 
HCT116.Tet5p.Ts5q.guide <- read.csv("count-HCT116-Ts3-i.f.csv") # beautiful curve, ~700 reads/guide 
HCT116.Tet5p.guide <- read.csv("count-HCT116-Ts5-i.f.csv")  # beautiful curve, ~700 reads/guide 

median(HCT116.Cont2.MegaA.guide$HCT116.Cont2.initial)      # 842
median(HCT116.Cont2.MegaA.guide$HCT116.Cont2.final)        # 618
median(HCT116.Tet5p.Ts5q.guide$HCT116.Ts3.initial)         # 857
median(HCT116.Tet5p.Ts5q.guide$HCT116.Ts3.final)           # 776
median(HCT116.Tet5p.guide$HCT116.Ts5.initial)              # 755
median(HCT116.Tet5p.guide$HCT116.Ts5.final)                # 825



# update March 5, 2023
# very high CDK12. remove CDK12
setwd(SchukkenData)
DLD1.Cont2.MegaA <- read.delim2("mle_DLD1-Cont2_P10.13-stratified.gene_summary.txt", 
                                  dec=".", header = TRUE) # Beta scores do not correlate with any other gene dropout
# Removed due to low quality read counts: 
#DLD1.Cont2.p10.MegaA <- read.delim2("mle_DLD1-Cont2_P10.gene_summary.txt", 
#                                dec=".", header = TRUE) # Beta scores do not correlate with any other gene dropout
#DLD1.Cont2.p13.MegaA <- read.delim2("mle_DLD1-Cont2_P13.gene_summary.txt", 
#                                dec=".", header = TRUE) # Beta scores do not correlate with any other gene dropout

DLD1.Di4.MegaA <- read.delim2("mle_DLD1-Ts4.gene_summary.txt", 
                                dec=".", header = TRUE)  # Beta scores do correlate with HCT116 Cont 2 beta
DLD1.Ts2.18.MegaA <- read.delim2("mle_DLD1-Ts2.18.gene_summary.txt", 
                                dec=".", header = TRUE)

### MegaA Screen 2: guide RNA info
DLD1.Cont2.MegaA.guide1 <- read.csv("count-DLD1-Cont2_P10-i.f.csv") 
DLD1.Cont2.MegaA.guide2 <- read.csv("count-DLD1-Cont2_P13-i.f.csv") 
DLD1.Cont2.MegaA.guide3 <- read.csv("count-DLD1-Cont2_P10.13-i.f-wo-cdk12.csv") 
DLD1.Di4.MegaA.guide <- read.csv("count-DLD1-Ts4-i.f.csv")
DLD1.Ts2.18.MegaA.guide <- read.csv("count-DLD1-Ts2.18-i.f.csv")  
median(DLD1.Di4.MegaA.guide$`DLD1.Ts4.initial`) #931
median(DLD1.Di4.MegaA.guide$`DLD1.Ts4.final`)   #1670
median(DLD1.Ts2.18.MegaA.guide$`DLD1.Ts2.18.initial`)   #1020
median(DLD1.Ts2.18.MegaA.guide$`DLD1.Ts2.18.final`)     #2448


## MegaA Screen 3 : Gene info
# Update  January 30, 2023
setwd(SchukkenData)

DLD1.Cont3.MegaA <- read.delim2("mle_DLD1-Cont3.gene_summary.txt", 
                                  dec=".", header = TRUE)
DLD1.Ts10.21.MegaA <- read.delim2("mle_DLD1-Ts10.21.gene_summary.txt", 
                                dec=".", header = TRUE)
DLD1.Ts5.15.MegaA <- read.delim2("mle_DLD1-Ts5.15.gene_summary.txt", 
                                dec=".", header = TRUE)

### MegaA Screen 3: guide RNA info
DLD1.Cont3.MegaA.guide <- read.csv("count-DLD1-Cont3-i.f.csv") # beautiful curve, ~700 reads/guide 
DLD1.Ts10.21.MegaA.guide <- read.csv("count-DLD1-Ts10.21-i.f.csv") # beautiful curve, ~700 reads/guide 
DLD1.Ts5.15.MegaA.guide <- read.csv("count-DLD1-Ts5.15-i.f.csv")  # beautiful curve, ~700 reads/guide 

median(DLD1.Cont3.MegaA.guide$`DLD1.Cont3.initial`) #1859
median(DLD1.Cont3.MegaA.guide$`DLD1.Cont3.final`)   #1657
median(DLD1.Ts10.21.MegaA.guide$`DLD1.Ts10.21.initial`) #1885
median(DLD1.Ts10.21.MegaA.guide$`DLD1.Ts10.21.final`) #1439
median(DLD1.Ts5.15.MegaA.guide$`DLD1.Ts5.15.initial`) #1668
median(DLD1.Ts5.15.MegaA.guide$`DLD1.Ts5.15.final`)   #1437



## MegaA Screen 4 : A2780 Gene info
# Update  May 24, 2023

setwd(SchukkenData)
# look at A2780 data where you remove sgRNAs that have 0 reads removed. to clean up the data
A2780.Cont1.MegaA <- read.delim2("mle_count-A2780-Cont1.gene_summary.txt", 
                                dec=".", header = TRUE)
A2780.Ts2.MegaA <- read.delim2("mle_count-A2780-Ts2.gene_summary.txt", 
                                  dec=".", header = TRUE)
A2780.Ts10.MegaA <- read.delim2("mle_count-A2780-Ts10.gene_summary.txt", 
                                 dec=".", header = TRUE)
# large Dropouts (Beta scores -60!!!) for some of these. 
# set minimum beta score to -5. 
# MYC and THAP1 are missing from Cont1 and Ts2
# add them back in (set at 0) then after quantile normalization (See below), set them to NA. 
A2780.Cont1.MegaA_forQN<- A2780.Cont1.MegaA
A2780.Cont1.MegaA_forQN$Cont1.beta<- as.numeric(A2780.Cont1.MegaA_forQN$Cont1.beta)
A2780.Cont1.MegaA_forQN[A2780.Cont1.MegaA_forQN$Cont1.beta < -3.5,]$Cont1.beta <- -3.5

A2780.Ts2.MegaA_forQN <- A2780.Ts2.MegaA
A2780.Ts2.MegaA_forQN$Ts2.beta<- as.numeric(A2780.Ts2.MegaA_forQN$Ts2.beta)
A2780.Ts2.MegaA_forQN[A2780.Ts2.MegaA_forQN$Ts2.beta < -3.5,]$Ts2.beta <- -3.5

A2780.Ts10.MegaA_forQN <- A2780.Ts10.MegaA
A2780.Ts10.MegaA_forQN$Ts10.beta<- as.numeric(A2780.Ts10.MegaA_forQN$Ts10.beta)
A2780.Ts10.MegaA_forQN[A2780.Ts10.MegaA_forQN$Ts10.beta < -3.5,]$Ts10.beta <- -3.5

### MegaA Screen 4: guide RNA info

A2780.Cont1.MegaA.guide <- read.csv("count-A2780-Cont1-i.f.csv") # 
A2780.Ts2.MegaA.guide <- read.csv("count-A2780-Ts2-i.f.csv") # 
A2780.Ts10.MegaA.guide <- read.csv("count-A2780-Ts10-i.f.csv")  # 
# There is dropout, but coverage at D=10 is low but acceptable.  

median(A2780.Cont1.MegaA.guide$A2780_Cont1_D0) #1343
median(A2780.Cont1.MegaA.guide$A2780_Cont1_D10) #458
median(A2780.Ts2.MegaA.guide$A2780_Ts2_D0) # 1238
median(A2780.Ts2.MegaA.guide$A2780_Ts2_D10) # 102  *** Very low!  mean = 1409!  Max = 300,000 (Lots of complete losses. low coverage during screen?) 
median(A2780.Ts10.MegaA.guide$A2780_Ts10_D0) # 1461
median(A2780.Ts10.MegaA.guide$A2780_Ts10_D10) # 635



## MegaA Screen 5 : A2780 
# Update  July 5, 2023
setwd(SchukkenData)
A2780.Cont2.MegaA <- read.delim2("mle_count-A2780-Cont2.gene_summary.txt", 
                                 dec=".", header = TRUE)
A2780.Ts8.MegaA <- read.delim2("mle_count-A2780-Ts8.gene_summary.txt",  ### This one has extremely low beta scores as well < -50!
                               dec=".", header = TRUE)
A2780.Ts18.MegaA <- read.delim2("mle_count-A2780-Ts18.gene_summary.txt", 
                                dec=".", header = TRUE)

# large Dropouts (Beta scores -50!!!) for Ts8. likely due to low coverage of D=10 
# set minimum beta score to -3.5. 
A2780.Ts8.MegaA_forQN <- A2780.Ts8.MegaA
A2780.Ts8.MegaA_forQN$Ts8.beta<- as.numeric(A2780.Ts8.MegaA_forQN$Ts8.beta)
A2780.Ts8.MegaA_forQN[A2780.Ts8.MegaA_forQN$Ts8.beta < -3.5,]$Ts8.beta <- -3.5


### MegaA Screen 5: guide RNA info
A2780.Cont2.MegaA.guide <- read.csv("count-A2780-Cont2.csv") # 
A2780.Ts8.MegaA.guide <- read.csv("count-A2780-Ts8.csv") # 
A2780.Ts18.MegaA.guide <- read.csv("count-A2780-Ts18.csv")  # 

median(A2780.Cont2.MegaA.guide$A2780_Cont2_D0) # 1911
median(A2780.Cont2.MegaA.guide$A2780_Cont2_D10) # 1188
median(A2780.Ts8.MegaA.guide$A2780_Ts8_D0) # 1147
median(A2780.Ts8.MegaA.guide$A2780_Ts8_D10) # 573
median(A2780.Ts18.MegaA.guide$A2780_Ts18_D0) # 1678
median(A2780.Ts18.MegaA.guide$A2780_Ts18_D10) # 1096



## MegaA Screen 6 : SNU1 
# Update  September 22, 2023
setwd(SchukkenData)
# C13 = Control/WT SNU1: Gain 20, monosomy X
# C12 =	Ts10.17.5 Ps2 
# C23 = Ts10 Ps6 (D=10 lost it's aneuploidy, do not include)
# C24 = Ts10.1.8.14 Ps5 (D=10 actually labeled C23) 
# C111 = Ts10.6q Ps8.12 
# C114 = Ts9.19 

Snu1.Cont1.C13.MegaA <- read.delim2("mle_count-Snu1C13.gene_summary.txt", 
                                 dec=".", header = TRUE)
Snu1.C12.MegaA <- read.delim2("mle_count-Snu1C12.gene_summary.txt", 
                                dec=".", header = TRUE)
Snu1.C24.MegaA <- read.delim2("mle_count-Snu1C24D0_v_Snu1C23D10.gene_summary.txt",  # Re-do with "C23 D=10" as karyotyping shows that is C24 D10. 
                              dec=".", header = TRUE)
Snu1.C111.MegaA <- read.delim2("mle_count-Snu1C111.gene_summary.txt", 
                              dec=".", header = TRUE)
Snu1.C114.MegaA <- read.delim2("mle_count-Snu1C114.gene_summary.txt", 
                              dec=".", header = TRUE)

### MegaA Screen 6: guide RNA info
Snu1.Cont1.C13.MegaA.guide <- read.csv("count-Snu-1_WT_C13.csv") # 
Snu1.C12.MegaA.guide <- read.csv("count-Snu-1_clone_12.csv")   # 
Snu1.C24.MegaA.guide <- read.csv("count-Snu1C24D0_v_Snu1C23D10.csv")
Snu1.C111.MegaA.guide <- read.csv("count-Snu-1_clone_111.csv") # 
Snu1.C114.MegaA.guide <- read.csv("count-Snu-1_clone_114.csv") # 

median(Snu1.Cont1.C13.MegaA.guide$Snu.1_WT_C13_D.0) #987
median(Snu1.Cont1.C13.MegaA.guide$Snu.1_WT_C13_D.10) #955
median(Snu1.C12.MegaA.guide$Snu.1_clone_12_D.0) #989
median(Snu1.C12.MegaA.guide$Snu.1_clone_12_D.10) #894
median(Snu1.C24.MegaA.guide$Snu.1_clone_24_D.0) #658
median(Snu1.C24.MegaA.guide$Snu.1_clone_23_D.10) #867
median(Snu1.C111.MegaA.guide$Snu.1_clone_111_D.0) #1217
median(Snu1.C111.MegaA.guide$Snu.1_clone_111_D.10) #705
median(Snu1.C114.MegaA.guide$Snu.1_clone_114_D.0) #744
median(Snu1.C114.MegaA.guide$Snu.1_clone_114_D.10) #786



## MegaA Screen 7 : Gene info
# Update  February 14, 2024
setwd(SchukkenData)
DLD1.Cont4.MegaA <- read.delim2("mle_count-DLD1__Ctrl_3.gene_summary.txt", 
                                dec=".", header = TRUE)
DLD1.Ts2.18_try2.MegaA <- read.delim2("mle_count-DLD1__Ts2.18.gene_summary.txt", 
                                  dec=".", header = TRUE)
DLD1.Ts12.17.MegaA <- read.delim2("mle_count-DLD1_Ts12.17.gene_summary.txt", 
                                 dec=".", header = TRUE)

### MegaA Screen 3: guide RNA info
DLD1.Cont4.MegaA.guide <- read.csv("count-DLD1_Ctrl_3.csv") # good coverage (700+ reads/guide), good curve
DLD1.Ts2.18_try2.MegaA.guide <- read.csv("count-DLD1__Ts2.18.csv") # good coverage (700+ reads/guide), good curve
DLD1.Ts12.17.MegaA.guide <- read.csv("count-DLD1_Ts12.17.csv")  # good coverage (500-700 reads/guide), good curve

median(DLD1.Cont4.MegaA.guide$DLD1_Ctrl_3_D0) #933
median(DLD1.Cont4.MegaA.guide$DLD1_Ctrl_3_D10) #855
median(DLD1.Ts2.18_try2.MegaA.guide$DLD1__Ts2.8_D0) #879
median(DLD1.Ts2.18_try2.MegaA.guide$DLD1__Ts2.8_D10) #947
median(DLD1.Ts12.17.MegaA.guide$DLD1_Ts3.12.17_D0) #700
median(DLD1.Ts12.17.MegaA.guide$DLD1_Ts3.12.17_D10) #986





#### Rename all Excel "Date" gene names

#Make function to rename "Sept-1" genes and "SEPT1" to their updated names "SEPTIN1" etc. for all 
#   genes who's names are misinterpretted as dates in Excel. 
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
  
  return(GeneMLEFile)
} 

HCT116.Cont1.MegaA<- RenameExcelGenes(HCT116.Cont1.MegaA) 
HCT116.Ts8.MegaA<- RenameExcelGenes(HCT116.Ts8.MegaA)    
DLD1.Cont1.MegaA<- RenameExcelGenes(DLD1.Cont1.MegaA)     
DLD1.Ts13.MegaA<- RenameExcelGenes(DLD1.Ts13.MegaA)      
DLD1.Ts8.10.MegaA<- RenameExcelGenes(DLD1.Ts8.10.MegaA)   

HCT116.Cont2.MegaA<- RenameExcelGenes(HCT116.Cont2.MegaA)
HCT116.Tet5p.Ts5q.MegaA<- RenameExcelGenes(HCT116.Tet5p.Ts5q.MegaA)
HCT116.Tet5p.MegaA<- RenameExcelGenes(HCT116.Tet5p.MegaA)
DLD1.Ts2.18.MegaA<- RenameExcelGenes(DLD1.Ts2.18.MegaA) #no CDK12

DLD1.Cont3.MegaA<- RenameExcelGenes(DLD1.Cont3.MegaA)
DLD1.Ts10.21.MegaA<- RenameExcelGenes(DLD1.Ts10.21.MegaA)
DLD1.Ts5.15.MegaA<- RenameExcelGenes(DLD1.Ts5.15.MegaA)

A2780.Cont1.MegaA<- RenameExcelGenes(A2780.Cont1.MegaA) # no MYC, THAP1
A2780.Ts2.MegaA<- RenameExcelGenes(A2780.Ts2.MegaA) # no MYC. 
A2780.Ts10.MegaA<- RenameExcelGenes(A2780.Ts10.MegaA) 

A2780.Cont1.MegaA_forQN<-RenameExcelGenes(A2780.Cont1.MegaA_forQN)
A2780.Ts2.MegaA_forQN<-RenameExcelGenes(A2780.Ts2.MegaA_forQN)
A2780.Ts10.MegaA_forQN<-RenameExcelGenes(A2780.Ts10.MegaA_forQN)

A2780.Cont2.MegaA_forQN<-RenameExcelGenes(A2780.Cont2.MegaA)
A2780.Ts8.MegaA_forQN<-RenameExcelGenes(A2780.Ts8.MegaA_forQN)
A2780.Ts18.MegaA_forQN<-RenameExcelGenes(A2780.Ts18.MegaA)


Snu1.Cont1.C13.MegaA<-RenameExcelGenes(Snu1.Cont1.C13.MegaA)
Snu1.C12.MegaA<-RenameExcelGenes(Snu1.C12.MegaA)
Snu1.C24.MegaA<-RenameExcelGenes(Snu1.C24.MegaA)
Snu1.C111.MegaA<-RenameExcelGenes(Snu1.C111.MegaA)
Snu1.C114.MegaA<-RenameExcelGenes(Snu1.C114.MegaA)

DLD1.Cont4.MegaA<-RenameExcelGenes(DLD1.Cont4.MegaA)
DLD1.Ts2.18_try2.MegaA<-RenameExcelGenes(DLD1.Ts2.18_try2.MegaA)
DLD1.Ts12.17.MegaA<-RenameExcelGenes(DLD1.Ts12.17.MegaA)





anti_join(DLD1.Cont3.MegaA, DLD1.Ts8.10.MegaA, by = "Gene")
setdiff(DLD1.Cont1.MegaA$Gene, A2780.Ts10.MegaA$Gene)
setdiff(HCT116.Cont1.MegaA$Gene, A2780.Ts10.MegaA$Gene)
setdiff(HCT116.Cont2.MegaA$Gene, A2780.Ts10.MegaA$Gene)
setdiff(DLD1.Cont3.MegaA$Gene, A2780.Ts10.MegaA$Gene)

setdiff(DLD1.Cont1.MegaA$Gene, DLD1.Ts8.10.MegaA$Gene)
setdiff(HCT116.Cont1.MegaA$Gene, DLD1.Ts8.10.MegaA$Gene)
setdiff(HCT116.Cont2.MegaA$Gene, DLD1.Ts8.10.MegaA$Gene)
setdiff(DLD1.Cont3.MegaA$Gene, DLD1.Ts8.10.MegaA$Gene)
setdiff(A2780.Cont1.MegaA$Gene, DLD1.Ts8.10.MegaA$Gene)





## Depmap: gene info & Chrm number
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
# DepMap, Broad (2022): DepMap 22Q2 Public. figshare. Dataset. https://doi.org/10.6084/m9.figshare.19700056.v2
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







### PLOT RAW CONTROLS ######## 
# Make heatmap for control gene dropout. heatmap beta scores
# then make boxplot of guides per sample read (initial and final per cell line)
# DEFINED FUNCTION
# get count file (from Mageck) as data frames. 

## Define function and positive/negative control RNAs.
PosCont<- c("PCNA", "POLR2D", "POLR2A", "RPL9", "RPL23A", "CDK9", "CDK1", "RPA3")
NegCont<- c("neg") #if you use grepl, you can find all genes containing "neg"

##Define function: for plotting control guide dropout
###Make heatmap of positive and negative controls 
### three negative and three positive control guides? or the genes? 
PlotHeatmap.Cont<-function(DropoutData = list("HCT116 Control 1" =HCT116.Cont1.MegaA, 
                                              "HCT116 Control 2"=HCT116.Cont2.MegaA, 
                                              "HCT116 Ts8"=HCT116.Ts8.MegaA), 
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
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  
  plot.controls
}


##Run function
PlotHeatmap.Cont(DropoutData= list("HCT116 Control 1"= HCT116.Cont1.MegaA, 
                                   "HCT116 Control 2"= HCT116.Cont2.MegaA, 
                                   "HCT116 Ts8"= HCT116.Ts8.MegaA, 
                                   "HCT116 Tet5p_1"= HCT116.Tet5p.Ts5q.MegaA, #labeled Ts3
                                   "HCT116 Tet5p_2"= HCT116.Tet5p.MegaA, 
                                   
                                   "DLD1 Control 1"= DLD1.Cont1.MegaA, 
                                   "DLD1 Control 3"= DLD1.Cont3.MegaA, 
                                   "DLD1 Control 4"= DLD1.Cont4.MegaA, 
                                   "DLD1 Ts2.18"= DLD1.Ts2.18.MegaA,
                                   "DLD1 Ts13"= DLD1.Ts13.MegaA, 
                                   "DLD1 Ts8.10"= DLD1.Ts8.10.MegaA,
                                   "DLD1 Ts10.21"= DLD1.Ts10.21.MegaA, 
                                   "DLD1 Ts5.15"= DLD1.Ts5.15.MegaA, 
                                   "DLD1 Ts2.18_try2"= DLD1.Ts2.18_try2.MegaA, 
                                   "DLD1 Ts12.17"= DLD1.Ts12.17.MegaA, 
                                   
                                   "A2780 Control 1" = A2780.Cont1.MegaA, 
                                   "A2780 Ts2" = A2780.Ts2.MegaA,
                                   "A2780 Ts10" = A2780.Ts10.MegaA, 
                                   "A2780 Control 2" = A2780.Cont2.MegaA, 
                                   "A2780 Ts8" = A2780.Ts8.MegaA,
                                   "A2780 Ts18" = A2780.Ts18.MegaA, 
                                   
                                   "Snu1 Control" = Snu1.Cont1.C13.MegaA, 
                                   "Snu1 C12" = Snu1.C12.MegaA,
                                   "Snu1 C24" = Snu1.C24.MegaA, 
                                   "Snu1 C111" = Snu1.C111.MegaA, 
                                   "Snu1 C114" = Snu1.C114.MegaA),
                 Title="Control guide Beta-scores")
#plot.MegaA1.2.3.4.5.heatmap.controls.pdf
# 8x4

##Run function




#plot boxplot reads per gene: 

df_list<- c(HCT116.Cont1.MegaA.guide, HCT116.Ts8.MegaA.guide,
            DLD1.Cont1.MegaA.guide, DLD1.Ts13.MegaA.guide, DLD1.Ts8.10.MegaA.guide
            )

Reads.MegaA<- merge(HCT116.Cont1.MegaA.guide, HCT116.Ts8.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, DLD1.Cont1.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, DLD1.Ts13.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, DLD1.Ts8.10.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))

Reads.MegaA<- merge(Reads.MegaA, DLD1.Cont2.MegaA.guide1[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, DLD1.Cont2.MegaA.guide2[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))

Reads.MegaA<- merge(Reads.MegaA, DLD1.Cont2.MegaA.guide3[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, DLD1.Di4.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, DLD1.Ts2.18.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))


Reads.MegaA<- merge(Reads.MegaA, A2780.Cont1.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, A2780.Ts2.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, A2780.Ts10.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))

Reads.MegaA<- merge(Reads.MegaA, A2780.Cont2.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, A2780.Ts8.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
Reads.MegaA<- merge(Reads.MegaA, A2780.Ts18.MegaA.guide[,c(1,3,4)], by= c("sgRNA"= "sgRNA"))
#Reads.MegaA<- Reads.MegaA[c(1,3,4,6,7,9,10,12,13,15,16,18,19)]

colnames(Reads.MegaA)<- c( "sgRNA",  "Gene",
                        "HCT116_Cont_initial", "HCT116_Cont_final", 
                        "HCT116_Ts8_initial", "HCT116_Ts8_final", 
                        "DLD1_Cont_initial", "DLD1_Cont_final", 
                        "DLD1_Ts13_initial", "DLD1_Ts13_final", 
                        "DLD1_Ts8.10_initial", "DLD1_Ts8.10_final", 
                        
                        "DLD1_Cont2.10_initial", "DLD1_Cont2.10_final", 
                        "DLD1_Cont2.13_initial", "DLD1_Cont2.13_final", 
                        "DLD1_Cont2.10.13_initial", "DLD1_Cont2.10.13_final", 
                        "DLD1_Ti4_initial", "DLD1_Ti4_final", 
                        "DLD1_Ts2.18_initial", "DLD1_Ts2.18_final", 
                        
                        "A2780_Control_1_initial", "A2780_Control_1_final",
                        "A2780_Ts2_initial",  "A2780_Ts2_final", 
                        "A2780_Ts10_initial", "A2780_Ts10__final", 
                        
                        "A2780_Control_2_initial", "A2780_Control_2_final", 
                        "A2780_Ts8_initial", "A2780_Ts8_final", 
                        "A2780_Ts18_initial", "A2780_Ts18_final")

df<- melt(Reads.MegaA)

ggplot(df, aes(x=variable, y=value))+
  geom_boxplot(outlier.shape=NA)+  #outlier.shape=NA
  #geom_boxplot()+ #with outliers
  xlab("Sample Name")+
  ylab("Reads per sgRNA") +
  ggtitle("Reads per sgRNA per sample")+
  theme(axis.text.x = element_text(angle = 90))+
  ylim(c(0,4500))+
  theme_classic()+
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
## plot.MegaA1.MeanReadssgRNA_noOutliers.pdf
# 8x5

ggplot(Reads.MegaA, aes(x=Reads.MegaA$DLD1_Cont_final))+
  geom_histogram(binwidth = 100)+
  xlim(-100,5000)




### Quantile normalization DG samples #####

## Quantile normalize the dropouts between all cell lines.  
## essentially rank order them, then take the average dropout per rank, give that number per gene. 
## then re-run this data as above, with the new rank-order-averaged dropout scores. 
## to quantile normalize, merge Cont3 and Ts13 data together, set gene name as row name, and look only at BETA SCORES

# Merge MegaA 1: 
MegaA.Merge.ALL<- merge(HCT116.Cont1.MegaA[,c(1,3)], HCT116.Ts8.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, HCT116.Cont2.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, HCT116.Tet5p.Ts5q.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, HCT116.Tet5p.MegaA[,c(1,3)], by=c("Gene", "Gene"))

MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Cont1.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Ts13.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Ts8.10.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Ts2.18.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Cont3.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Ts10.21.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Ts5.15.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Cont4.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Ts2.18_try2.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, DLD1.Ts12.17.MegaA[,c(1,3)], by=c("Gene", "Gene"))

MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, A2780.Cont1.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene")) ### This one is wierd, low coverage leads to extreme negative beta scores
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, A2780.Ts2.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene")) ### This one is wierd, low coverage leads to extreme negative beta scores
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, A2780.Ts10.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, A2780.Cont2.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, A2780.Ts8.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, A2780.Ts18.MegaA[,c(1,3)], by=c("Gene", "Gene"))

MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, Snu1.Cont1.C13.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, Snu1.C12.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, Snu1.C24.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, Snu1.C111.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.ALL<- merge(MegaA.Merge.ALL, Snu1.C114.MegaA[,c(1,3)], by=c("Gene", "Gene"))



colnames(MegaA.Merge.ALL)<- c("Gene", "H.Cont1.beta",  "H.Ts8.beta",  #"H.Ts13.beta", 
                              "H.Cont2.beta", "H.Tet5pTs5q.beta", "H.Tet5p.beta", 
                              "D.Cont1.beta", "D.Ts13.beta", "D.Ts8.10.beta", 
                              #"D.Cont2.beta", "D.Ti4_Loss.beta", 
                              "D.Ts2.18.beta", 
                              "D.Cont3.beta", "D.Ts10.21.beta", "D.Ts5.15.beta",
                              "D.Cont4.beta", "D.Ts2.18_try2.beta", "D.Ts12.17.beta",
                              "A.Cont1.beta", "A.Ts2.beta", "A.Ts10.beta", 
                              "A.Cont2.beta", "A.Ts8.beta", "A.Ts18.beta", 
                              "S.Ctrl.beta", "S.C12.beta", #"S.C23.beta", 
                              "S.C24.beta", "S.C111.beta", "S.C114.beta") 

# Remove all negative controls, except for three negative controls for heatmap plotting purposes. to show positive and negative controls. 
MegaA.Merge.ALL$Gene<- as.character(MegaA.Merge.ALL$Gene)
negControls<- MegaA.Merge.ALL[which(grepl("neg", MegaA.Merge.ALL$Gene)),] #list all negative controls
negControls<- negControls[-c(1:5,13),] #remove 3 negative controls
MegaA.Merge.ALL <- subset(MegaA.Merge.ALL, ! Gene %in% negControls$Gene)   #remove rows with negative controls, except the three above


#MegaA.Merge.ALL[1545, "A.Cont1.beta"]<- 0 # Set MYC A2780 Control 1 beta to 0 instead of NA
#MegaA.Merge.ALL[1545, "A.Ts2.beta"]<- 0 # Set MYC A2780 Ts2 beta to 0 instead of NA
#MegaA.Merge.ALL[2583, "A.Cont1.beta"]<- 0 # Set THAP1 A2780 ontrol1 beta to 0 instead of NA

# write.csv(MegaA.Merge.ALL, "MegaA.Merge.ALL.csv")


### below function will quantile normalize your data
### !!! BUT!!! you need merged control and aneuploid data, and ONLY give function the BETA SCORE COLLUMNS!!!! 
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
  return(as.data.frame(df_final))
}

## Quantile normalize by batch and cell line: 
MegaA.Merge.ALL2<- MegaA.Merge.ALL

# Quantile normalize! 
MegaA.All.Merge.QN<- quantile_normalisation(MegaA.Merge.ALL2[,c(2:length(MegaA.Merge.ALL2))])
MegaA.All.Merge.QN$Gene<- (MegaA.Merge.ALL$Gene)



## Check quantile normalization, they should have the same cutoffs for the 0, 25, 50, 75 and 100 places. 
quantile(MegaA.All.Merge.QN[,1])
quantile(MegaA.All.Merge.QN[,2])
quantile(MegaA.All.Merge.QN[,3])
quantile(MegaA.All.Merge.QN[,10])




### Plot quantile normalized (QN) positive controls
# negative controls have been removed from QN data
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
GeneColumn<- 38


PlotHeatmap.Cont3(DropoutData= list("HCT116 Control 1"= MegaA.All.Merge.QN[,c(1,27)], 
                                    "HCT116 Ts8"= MegaA.All.Merge.QN[,c(2,27)], 
                                    "HCT116 Control 2"= MegaA.All.Merge.QN[,c(3,27)], 
                                    "HCT116 Tet5p"= MegaA.All.Merge.QN[,c(4,27)], 
                                    "HCT116 Tet5p_Try2"= MegaA.All.Merge.QN[,c(5,27)], 
                                    
                                    "DLD1 Control 1"= MegaA.All.Merge.QN[,c(6,27)], 
                                    "DLD1 Ts13"= MegaA.All.Merge.QN[,c(7,27)], 
                                    "DLD1 Ts8.10"= MegaA.All.Merge.QN[,c(8,27)],
                                    "DLD1 Control 3"= MegaA.All.Merge.QN[,c(10,27)], 
                                    "DLD1 Ts2.18"= MegaA.All.Merge.QN[,c(9,27)], 
                                    "DLD1 Ts10.21"= MegaA.All.Merge.QN[,c(11,27)], 
                                    "DLD1 Ts5.15"= MegaA.All.Merge.QN[,c(12,27)], 
                                    
                                    "DLD1 Control 4"= MegaA.All.Merge.QN[,c(13,27)],
                                    "DLD1 Ts2.18_Try2"= MegaA.All.Merge.QN[,c(14,27)],
                                    "DLD1 Ts12.17"= MegaA.All.Merge.QN[,c(15,27)],
                                    
                                    "A2780 Control"= MegaA.All.Merge.QN[,c(16,27)],
                                    "A2780 Ts2"= MegaA.All.Merge.QN[,c(17,27)],
                                    "A2780 Ts10"= MegaA.All.Merge.QN[,c(18,27)],
                                    "A2780 Control 2"= MegaA.All.Merge.QN[,c(19,27)],
                                    "A2780 Ts8"= MegaA.All.Merge.QN[,c(20,27)],
                                    "A2780 Ts18"= MegaA.All.Merge.QN[,c(21,27)],
                                    
                                    "SNU1 Control (Tetraploid)"= MegaA.All.Merge.QN[,c(22,27)],
                                    "SNU1 Ts10.17.5_Penta2"= MegaA.All.Merge.QN[,c(23,27)],
                                    "SNU1 Ts10.1.8.14_Penta5"= MegaA.All.Merge.QN[,c(24,27)],
                                    "SNU1 Ts10.6q_Penta8.12"= MegaA.All.Merge.QN[,c(25,27)],
                                    "SNU1 Ts9.19"= MegaA.All.Merge.QN[,c(26,27)]), 
                                    
                  Title="Control guide Beta-scores")
# plot.MegaA1.2.3.heatmap.controls.QN.pdf






### USE THIS DATASETS: 
MegaA.All.Merge.QN

setwd(ResultsFile)
# write.csv(MegaA.All.Merge.QN, "CRISPRScreens_BetaScores_MegaA_QuantileNorm_April2025.csv")
#MegaA.All.Merge.QN<- read.csv("CRISPRScreens_BetaScores_MegaA_QuantileNorm_April2025.csv")
#MegaA.All.Merge.QN<- MegaA.All.Merge.QN[,-1]



### Plot QN DG Data #####

MegaA.All.Merge.QN_meanDiff<- MegaA.All.Merge.QN

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell 
MegaA.All.Merge.QN_meanDiff$Paired.Diff<- ( 
                             (MegaA.All.Merge.QN_meanDiff$H.Ts8.beta-MegaA.All.Merge.QN_meanDiff$H.Cont1.beta) + 
                             (MegaA.All.Merge.QN_meanDiff$H.Tet5pTs5q.beta-MegaA.All.Merge.QN_meanDiff$H.Cont2.beta) + 
                             (MegaA.All.Merge.QN_meanDiff$H.Tet5p.beta-MegaA.All.Merge.QN_meanDiff$H.Cont2.beta) + 
                              
                             (MegaA.All.Merge.QN_meanDiff$D.Ts13.beta-MegaA.All.Merge.QN_meanDiff$D.Cont1.beta) + 
                             (MegaA.All.Merge.QN_meanDiff$D.Ts8.10.beta-MegaA.All.Merge.QN_meanDiff$D.Cont1.beta) +
                             (MegaA.All.Merge.QN_meanDiff$D.Ts2.18.beta -MegaA.All.Merge.QN_meanDiff$D.Cont1.beta) +
                             (MegaA.All.Merge.QN_meanDiff$D.Ts10.21.beta-MegaA.All.Merge.QN_meanDiff$D.Cont3.beta) +
                             (MegaA.All.Merge.QN_meanDiff$D.Ts5.15.beta -MegaA.All.Merge.QN_meanDiff$D.Cont3.beta) + 
                              
                             (MegaA.All.Merge.QN_meanDiff$A.Ts10.beta -MegaA.All.Merge.QN_meanDiff$A.Cont1.beta) +
                             (MegaA.All.Merge.QN_meanDiff$A.Ts2.beta -MegaA.All.Merge.QN_meanDiff$A.Cont1.beta) +
                              (MegaA.All.Merge.QN_meanDiff$A.Ts8.beta -MegaA.All.Merge.QN_meanDiff$A.Cont2.beta) +
                              (MegaA.All.Merge.QN_meanDiff$A.Ts18.beta -MegaA.All.Merge.QN_meanDiff$A.Cont2.beta) +
                            
                              (MegaA.All.Merge.QN_meanDiff$S.C12.beta -MegaA.All.Merge.QN_meanDiff$S.Ctrl.beta)+ 
                              (MegaA.All.Merge.QN_meanDiff$S.C24.beta -MegaA.All.Merge.QN_meanDiff$S.Ctrl.beta)+ 
                              (MegaA.All.Merge.QN_meanDiff$S.C111.beta -MegaA.All.Merge.QN_meanDiff$S.Ctrl.beta)+ 
                              (MegaA.All.Merge.QN_meanDiff$S.C114.beta -MegaA.All.Merge.QN_meanDiff$S.Ctrl.beta)
                             )/17

# Now add p-values to difference between all aneuploid and all euploid
MegaA.All.Merge.QN_meanDiff$p.paired<- NA
for (i in 1:length(MegaA.All.Merge.QN_meanDiff$Gene)){
  test<- t.test(as.numeric(MegaA.All.Merge.QN_meanDiff[i,c(2, 4,5, 7,8,9, 11,12, 14,15, 17,18, 20,21, 23,24,26)]), 
                as.numeric(MegaA.All.Merge.QN_meanDiff[i,c(1, 3,3, 6,6,6, 10,10, 13,13, 16,16, 19,19, 22,22,22)]), paired=TRUE) # t.test (aneuploid, euploid)
  MegaA.All.Merge.QN_meanDiff$p.paired[i]<-test$p.value
} #paired t-test

MegaA.All.Merge.QN_meanDiff$MeanBeta<- rowMeans(MegaA.All.Merge.QN_meanDiff[,c(1:26)])
MegaA.All.Merge.QN_meanDiff<- MegaA.All.Merge.QN_meanDiff[order(MegaA.All.Merge.QN_meanDiff$p.paired), ]


## Plot difference v pvalue (scatter plot)
ggplot(MegaA.All.Merge.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff, color=p.paired<0.05 & Paired.Diff< -0.2))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values=c("Black", "red"))+ 
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")+
  theme(legend.position="none")
# 5x4
# plot.BetaDelta.pValue.AllMegaA_QN_July2025.pdf

subset(MegaA.All.Merge.QN_meanDiff, p.paired<0.05 & Paired.Diff< -0.2)$Gene
# [1] "BTBD2"   "RAX2"    "EPHA3"   "HOXD10"  "FBXO2"   "TEAD3"   "MYT1"    "MMP21"   "VPS18"   "TRIO"    "ZNF784"  "GLI4"    "KBTBD7"  "ASCL2"   "ZNF768"  "UBE2N"  
# [17] "GATA5"   "FOXS1"   "EBF4"    "TRIM8"   "PRSS33"  "SPI1"    "ZNF771"  "HSF1"    "PHOX2B"  "CDK4"    "CBX7"    "CSNK1E"  "SIX6"    "JUNB"    "PCNA"    "EGR3"   
# [33] "TWIST1"  "ATF6B"   "KLF10"   "NACC1"   "HDAC3"   "UBE2M"   "UBE2H"   "KLF16"   "FOXD4L5"


## Plot Pathways? 
# Very few of the top pathway genes in Druggable genome. 
ggplot(MegaA.All.Merge.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff))+
  geom_point()+
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% Proteosome_Genes$Approved.symbol), color = "chartreuse1") + 
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% RNA_Processing_Spliceosome), color = "deepskyblue1")+
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% MitoElectronChainnAssembly), color = "#F8766D") + 
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% MitoTranslationTranscription), color = "gold2") + 
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% Ribosomal_Genes_noMT$Approved.symbol), color = "#C77CFF")+
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% RNA_Processing_rRNA), color = "deepskyblue4")+
  geom_density_2d(color = "white")+
  theme_classic()+
  xlab("Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")+
  theme(legend.position="none")
# 5x4
# plot.BetaDelta.pValue.AllMegaA_QN_July2025.pdf
# PSMD14, USP49, PA2G4




## highlight top hit gene
HighlightTheseGenes<- c("FOSL1", "UBE2H")

ggplot(MegaA.All.Merge.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff, color=Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  geom_density_2d(color = "white")+
  xlab("Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")+
  theme(legend.position="none")
# 5x4
# plot.MegaA_Only.pairedP.Diff_KinaseHits.pdf
subset(MegaA.All.Merge.QN_meanDiff, Gene %in% HighlightTheseGenes)


GeneOfInterest<- "UBE2H"
df<- data.frame(Control=as.numeric(subset(MegaA.All.Merge.QN_meanDiff, Gene==GeneOfInterest)[,c(1, 3,3, 6,6,6, 10,10, 13,13, 16,16, 19,19, 22,22,22,22)]), 
                Aneuploid=as.numeric(subset(MegaA.All.Merge.QN_meanDiff, Gene==GeneOfInterest)[,c(2, 4,5, 7,8,9, 11,12, 14,15, 17,18, 20,21, 23,24,25,26)]) )
ggpaired(df, cond1 = "Control", cond2 = "Aneuploid", 
         title=GeneOfInterest, fill="condition")+ 
  ylab("Quantile Norm. Beta")+
  scale_fill_manual(values=c("skyblue", "yellow2" ))
# 3x5
# plot.paired.MegaA_GLI4.pdf




# Plot RANK
MegaA.All.Merge.QN_meanDiff$AneuploidRank<- rank(MegaA.All.Merge.QN_meanDiff$Paired.Diff)
length(MegaA.All.Merge.QN_meanDiff$AneuploidRank)*0.1

ggplot(MegaA.All.Merge.QN_meanDiff, aes(x=AneuploidRank, y=Paired.Diff, color= AneuploidRank < 342))+
  geom_point()+
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, AneuploidRank < 342))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  ylab("Difference in Beta Score")+
  xlab("Druggable Genome: Rank")
# 5x3
# plot.MegaA.Rank_Top10.pdf

HighlightTheseGenes<- c("FOSL1")
ggplot(MegaA.All.Merge.QN_meanDiff, aes(x=AneuploidRank, y=Paired.Diff, color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(MegaA.All.Merge.QN_meanDiff, Gene %in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  ylab("Difference in Beta Score")+
  xlab("Druggable Genome: Rank")
# 5x3
# plot.MegaA.Rank_FOSL1.pdf

# write rank data 
# write.csv(subset(MegaA.All.Merge.QN_meanDiff, AneuploidRank < 342)$Gene, "Sheltzer.Druggable.RANK.top10.csv")
#  write.csv(MegaA.All.Merge.QN_meanDiff, "MegaA.All.Merge.QN_meanDiff.csv")



setwd(ResultsFile)
#write.csv(MegaA.All.Merge.QN_meanDiff, "Sheltzer.MegaA.QN_meanDiff.csv")

#MegaA.All.Merge.QN_meanDiff<- read.csv("Sheltzer.MegaA.QN_meanDiff.csv")
#MegaA.All.Merge.QN_meanDiff<- MegaA.All.Merge.QN_meanDiff[,-1]


### Merge DG and WGL results, QN #####
# Import data from WGL_all_Analysis_2.R

MegaA.All.Merge.QN
MegaA.Merge.ALL #from "Quantile normalization ALL" section


MegaA.Merge.ALL

# First merge the WGL data. 
MLE.Merge.ALL3<- merge(DLD1.Cont1.mle[,c(1,3)], DLD1.Ts13.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, Vaco432.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, Vaco432.Ts13.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, A2780.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, A2780.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, A2058.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, A2058.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, A2058.7pDi.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, AGS.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, AGS.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, MCF10A.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3<- merge(MLE.Merge.ALL3, MCF10A.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MLE.Merge.ALL3$Gene<- gsub("neg ", "neg", MLE.Merge.ALL3$Gene)

MLE.Merge.ALL4<- merge(MegaA.Merge.ALL, MLE.Merge.ALL3[,1:14], by.x="Gene", by.y="Gene") 


## 
Non_Merge<- anti_join(MegaA.Merge.ALL, MLE.Merge.ALL3[,1:14], join_by("Gene"=="Gene"))$Gene # Genes in MegaA library not in WGL = 53

MegaA.Merge.ALL_missing<- subset(MegaA.Merge.ALL, Gene %in% Non_Merge)
MegaA.Merge.ALL_missing$Gene2<- MegaA.Merge.ALL_missing$Gene
for (i in 1:length(MegaA.Merge.ALL_missing$Gene)) {
  if (length(alias2Symbol(MegaA.Merge.ALL_missing$Gene[i], species = "Hs", expand.symbols = FALSE) )>0){ # If gene can be converted. 
    MegaA.Merge.ALL_missing$Gene2[i]<-alias2Symbol(MegaA.Merge.ALL_missing$Gene[i], species = "Hs", expand.symbols = FALSE)[1] #label to new gene name
  }
} 

length(MegaA.Merge.ALL_missing$Gene2) # 47

MLE.Merge.ALL3$Gene2<- MLE.Merge.ALL3$Gene
for (i in 1:length(MLE.Merge.ALL3$Gene)){
  Alias<- alias2Symbol(MLE.Merge.ALL3$Gene[i], species = "Hs", expand.symbols = FALSE) # find if gene has alias
  if (length(Alias) > 0){ # If gene can be converted. 
    MLE.Merge.ALL3$Gene2[i]<-Alias[1] #label alias to new gene name
  }
} #Is this the best way to run this code? no. but it works. It takes a while to run though... so try to do it only once. 


Intermediate<- merge(MegaA.Merge.ALL_missing, MLE.Merge.ALL3, by.x="Gene2", by.y="Gene2") # 
length(Intermediate$Gene2) # How many genes? 22

length(unique(Intermediate$Gene2)) # How many unique genes? 21. so one replicate. CRAMP1 as both HN1L and CRAMP1L. I chose CRAMP1L

Intermediate<- Intermediate[-2,] # remove duplicate Gene.y = HN1L 
Intermediate<- Intermediate[, -which(names(Intermediate) %in% c("Gene.x", "Gene.y"))] #remove bad gene names

colnames(Intermediate)[1]<- c("Gene") #rename "Gene2" back to "Gene"
Intermediate<- subset(Intermediate, ! Gene %in% MLE.Merge.ALL4$Gene) # make sure no gene overlaps. KATXXX genes already in MegaA dataset


MLE.Merge.ALL4<- rbind(MLE.Merge.ALL4, Intermediate) # genes merged
length(MLE.Merge.ALL4$Gene) #3377
length(unique(MLE.Merge.ALL4$Gene)) #3377. no replicate gene names. 

# Finding gene symbol for all genes helped merge 13 extra genes 







###



colnames(MLE.Merge.ALL4)<- c("Gene", "HCT116.Cont1.beta", "HCT116.Ts8.beta", # H.Ts13
                                                     "HCT116.Cont2.beta", "HCT116.Tet5p.Ts5q.beta", "HCT116.Tet5.beta", 
                                                     "DLD1.Cont1.beta", "DLD1.Ts13.beta", "DLD1.Ts8.10.beta", 
                                                     "DLD1.Ts2.18.beta", 
                                                     "DLD1.Cont3.beta", "DLD1.Ts10.21.beta", "DLD1.Ts5.15.beta", 
                                                      "DLD1.Cont4.beta", "DLD1.Ts2.18.beta_Try2", "DLD1.Ts12.17.beta",
                             
                                                      "A2780.Cont1.beta", "A2780.Ts2.beta", "A2780.Ts10.beta",
                                                      "A2780.Cont2.beta", "A2780.Ts8.beta", "A2780.Ts18.beta",
                                                      "Snu1.Cont.beta", "Snu1.Ts10.17.5_Ps2.beta", # "Snu1.Ts10_Ps6.beta",
                                                      "Snu1.Ts10.1.8.14_Ps5.beta", "Snu1.Ts10.6q_Ps8.12.beta", "Snu1.Ts9.19.beta",
                             
                            "WGL.DLD1.Cont.beta", "WGL.DLD1.Ts13.beta", 
                            "WGL.Vaco432.Cont.beta", "WGL.Vaco432.Ts13.beta", 
                            # WGL HCT116.control, WGL.HCT116.Ts13
                            "WGL.A2780.Cont.beta", "WGL.A2780.Di1q.beta", 
                            "WGL.A2058.Cont.beta", "WGL.A2058.Di1q.beta", "WGL.A2058.Di7p.beta", 
                            "WGL.AGS.Cont.beta", "WGL.AGS.Di1q.beta", 
                            "WGL.MCF10A.Cont.beta", "WGL.MCF10A.1qDi.beta")


## Quantile normalize by batch and cell line: 
MLE.Merge.ALL4 # quantile normalized across all samples from MegaA and WGL
MLE.Merge.ALL4.2<- MLE.Merge.ALL4

# Since A2780 has several unreliable outliers, turn outliers into "0" 
#    Quantile normalization does not accept NA, so turn NA into 0, then back to NA
MLE.Merge.ALL4.2$A2780.Cont1.beta[is.na(MLE.Merge.ALL4.2$A2780.Cont1.beta)] <- 0
MLE.Merge.ALL4.2$A2780.Ts2.beta[is.na(MLE.Merge.ALL4.2$A2780.Ts2.beta)] <- 0
MLE.Merge.ALL4.2$A2780.Ts10.beta[is.na(MLE.Merge.ALL4.2$A2780.Ts10.beta)] <- 0
MLE.Merge.ALL4.2$A2780.Ts8.beta[is.na(MLE.Merge.ALL4.2$A2780.Ts8.beta)] <- 0

# Quantile normalize! 
Merge.ALL.QN<- quantile_normalisation(MLE.Merge.ALL4.2[,c(2:length(MLE.Merge.ALL4.2))])
Merge.ALL.QN<- as.data.frame(Merge.ALL.QN)
Merge.ALL.QN$Gene<- (MLE.Merge.ALL4.2$Gene)

# Since A2780 has several unreliable outliers, turn outliers into "0" 
#    Quantile normalization does not accept NA, so turn NA into 0, then back to NA
Merge.ALL.QN$A2780.Cont1.beta[is.na(MLE.Merge.ALL3$A2780.Cont1.beta)] <- NA
Merge.ALL.QN$A2780.Ts2.beta[is.na(MLE.Merge.ALL3$A2780.Ts2.beta)] <- NA
Merge.ALL.QN$A2780.Ts10.beta[is.na(MLE.Merge.ALL3$A2780.Ts10.beta)] <- NA
Merge.ALL.QN$A2780.Ts8.beta[is.na(MLE.Merge.ALL3$A2780.Ts8.beta)] <- NA



setwd(ResultsFile)
# write.csv(MLE.Merge.ALL4, "CRISPRScreens_BetaScores_WGLandMEGAA_April2025.csv")
# write.csv(Merge.ALL.QN, "CRISPRScreens_BetaScores_WGLandMEGAA_QuantileNorm_April2025.csv")
# Merge.ALL.QN<- read.csv("CRISPRScreens_BetaScores_WGLandMEGAA_QuantileNorm_April2025.csv")
# Merge.ALL.QN<- Merge.ALL.QN[,-1]


## find difference between aneuploid & control
Merge.ALL.QN2<- Merge.ALL.QN
rownames(Merge.ALL.QN2)<- Merge.ALL.QN2$Gene




### Plot quantile normalized (QN) positive controls
# Negative controls have been removed from QN data
# plot Heatmap quantile normalized data, DG and WGL
GeneColumn<- 40

PlotHeatmap.Cont3(DropoutData= list("HCT116 Control 1"= Merge.ALL.QN2[,c(1,GeneColumn)], 
                                   "HCT116 Ts8"= Merge.ALL.QN2[,c(2,GeneColumn)], 
                                   
                                   "HCT116 Control 2"= Merge.ALL.QN2[,c(3,GeneColumn)], 
                                   "HCT116 Tet5p"= Merge.ALL.QN2[,c(4,GeneColumn)], 
                                   "HCT116 Tet5p_2"= Merge.ALL.QN2[,c(5,GeneColumn)], 
                                   
                                   "DLD1 Control 1"= Merge.ALL.QN2[,c(6,GeneColumn)], 
                                   "DLD1 Ts13"= Merge.ALL.QN2[,c(7,GeneColumn)], 
                                   "DLD1 Ts8.10"= Merge.ALL.QN2[,c(8,GeneColumn)],
                                   "DLD1 Ts2.18"= Merge.ALL.QN2[,c(9,GeneColumn)], 
                                   
                                   "DLD1 Control 3"= Merge.ALL.QN2[,c(10,GeneColumn)], 
                                   "DLD1 Ts10.21"= Merge.ALL.QN2[,c(11,GeneColumn)], 
                                   "DLD1 Ts5.15"= Merge.ALL.QN2[,c(12,GeneColumn)], 
                                   
                                   "DLD1 Control 4"= Merge.ALL.QN2[,c(13,GeneColumn)],
                                   "DLD1 Ts2.18_Try2"= Merge.ALL.QN2[,c(14,GeneColumn)],
                                   "DLD1 Ts12.17"= Merge.ALL.QN2[,c(15,GeneColumn)],
                                   
                                   "A2780 Control"= Merge.ALL.QN2[,c(16,GeneColumn)],
                                   "A2780 Ts2"= Merge.ALL.QN2[,c(17,GeneColumn)],
                                   "A2780 Ts10"= Merge.ALL.QN2[,c(18,GeneColumn)],
                                   
                                   "A2780 Control 2"= Merge.ALL.QN2[,c(19,GeneColumn)],
                                   "A2780 Ts8"= Merge.ALL.QN2[,c(20,GeneColumn)],
                                   "A2780 Ts18"= Merge.ALL.QN2[,c(21,GeneColumn)],
                                   
                                   "SNU1 Control"= Merge.ALL.QN2[,c(22,GeneColumn)],
                                   "SNU1 Ts10.17.5_Ps2"= Merge.ALL.QN2[,c(23,GeneColumn)],
                                   "SNU1 Ts10.1.8.14_Ps5"= Merge.ALL.QN2[,c(24,GeneColumn)],
                                   "SNU1 Ts10.6q_Ps8.12"= Merge.ALL.QN2[,c(25,GeneColumn)],
                                   "SNU1 Ts9.19"= Merge.ALL.QN2[,c(26,GeneColumn)],
                                   
                                   "DLD1 Control WGL"= Merge.ALL.QN2[,c(27,GeneColumn)], 
                                   "DLD1 Ts13 WGL"= Merge.ALL.QN2[,c(28,GeneColumn)],
                                   
                                   "Vaco432 Control WGL"= Merge.ALL.QN2[,c(29,GeneColumn)], 
                                   "Vaco432 Ts13 WGL"= Merge.ALL.QN2[,c(30,GeneColumn)], 
                                   
                                   "A2780 WT WGL"= Merge.ALL.QN2[,c(31,GeneColumn)], 
                                   "A2780 Di1q WGL"= Merge.ALL.QN2[,c(32,GeneColumn)],
                  
                                   "A2058 WT WGL"= Merge.ALL.QN2[,c(33,GeneColumn)], 
                                   "A2058 Di1q WGL"= Merge.ALL.QN2[,c(34,GeneColumn)], 
                                   "A2058 Di7p WGL"= Merge.ALL.QN2[,c(35,GeneColumn)], 
                                   
                                   "AGS WT WGL"= Merge.ALL.QN2[,c(36,GeneColumn)], 
                                   "AGS Di1q WGL"= Merge.ALL.QN2[,c(37,GeneColumn)],   
                  
                                  "MCF10A WT WGL"= Merge.ALL.QN2[,c(38,GeneColumn)], 
                                  "MCF10A Di1q WGL"= Merge.ALL.QN2[,c(39,GeneColumn)]),
                                   Title="Control guide Beta-scores")
# plot.MegaA1.2.3.WGL.heatmap.controls.QN_neg.pdf
# Figure 1E 
# 7x5



## plot single gene. plot gene of interest. categorical camparison
GeneOfInterest<- "UBE2H"
df<- data.frame(Control=as.numeric(subset(Merge.ALL.QN2, Gene==GeneOfInterest)[,c(1, 3,3, 6,6,6, 10,10, 13,13, 16,16, 19,19, 22,22,22,22,   27, 29, 32, 34,35, 37)]), 
                Aneuploid=as.numeric(subset(Merge.ALL.QN2, Gene==GeneOfInterest)[,c(2, 4,5, 7,8,9, 11,12, 14,15, 17,18, 20,21, 23,24,25,26, 28, 30, 31, 33,33, 36)]) )
ggpaired(df, cond1 = "Control", cond2 = "Aneuploid", 
         title=GeneOfInterest, fill="condition")+ 
  ylab("Quantile Norm. Beta")+
  scale_fill_manual(values=c("skyblue", "yellow2"))
# 3x5
# plot.paired.MegaA.WGL_UBE2H.pdf
t.test(as.numeric(subset(Merge.ALL.QN2, Gene==GeneOfInterest)[,c(1, 3,3, 6,6,6, 10,10, 13,13, 16,16, 19,19, 22,22,22,22,   27, 29, 32, 34,35, 37, 39)]),
       as.numeric(subset(Merge.ALL.QN2, Gene==GeneOfInterest)[,c(2, 4,5, 7,8,9, 11,12, 14,15, 17,18, 20,21, 23,24,25,26, 28, 30, 31, 33,33, 36, 38)])  ,
       paired=TRUE) # p= 0.04433, difference = 0.2686855





### Principal Component Analysis (PCA) ####
## https://www.datacamp.com/tutorial/pca-analysis-r

## Step 1: merge CRISPR dropout scores and normalize 
# get Merge.ALL.QN2 from above COMPARE MegaA TO WGL RESULTS, with data from "WGL_all_Analysis.R"
Merge.ALL.QN2


## Step 2:  heatmap clustering WGL + MegaA

# Dendogram sample clustering 
col<- colorRampPalette(c("blue", "white", "red"))(20)       
heatmap(x = as.matrix(Merge.ALL.QN2[,c(1:39)]), col = col)  
# heatmap(x = as.matrix(MLE.Merge.ALL3[,c(2:28)]), col = col) 
# plot all beta scores and correlate gene dropouts for all cell lines. 
# we see DLD1 and HCT116 cluster among themselves
# We see batch specific clustering. 
# plot.QNBeta.Branch.heatmap_WGL.MegaA.pdf
# 12x12
# Supplementary Figure 1A



## Step 3: Principal component analysis (PCA): make correlation matrix
# corr_matrix <- cor(data_normalized)

Corr_matrix_WGL.MegaA<-cor(Merge.ALL.QN2[,1:39])



## PRINCIPAL COMPONENT ANALYSIS (PCA)

#data.pca <- princomp(corr_matrix)
#summary(data.pca)
# from https://www.datacamp.com/tutorial/pca-analysis-r 
# "R PCA summary 
# From the previous screenshot, we notice that nine principal components have been generated (Comp.1 to Comp.9), which also correspond to the number of variables in the data.
# Each component explains a percentage of the total variance in the data set. In the Cumulative Proportion section, the first principal component explains almost 77% of the total variance. This implies that almost two-thirds of the data in the set of 9 variables can be represented by just the first principal component. The second one explains 12.08% of the total variance. 
# The cumulative proportion of Comp.1 and Comp.2 explains nearly 89% of the total variance. This means that the first two principal components can accurately represent the data. 
# It’s great to have the first two components, but what do they really mean? 
# This can be answered by exploring how they relate to each column using the loadings of each principal component." 
WGL.MegaA.pca<-princomp(Corr_matrix_WGL.MegaA)
summary(WGL.MegaA.pca)

#Importance of components:


# Importance of components:
#Comp.1    Comp.2     Comp.3     Comp.4     Comp.5     Comp.6     Comp.7     Comp.8     Comp.9     Comp.10
#Standard deviation     0.6990111 0.3455890 0.23548690 0.14963255 0.12137156 0.11945623 0.10072610 0.09718982 0.08910532 0.082740688
#Proportion of Variance 0.6193690 0.1513914 0.07029344 0.02838138 0.01867305 0.01808835 0.01286072 0.01197355 0.01006442 0.008677998
#Cumulative Proportion  0.6193690 0.7707604 0.84105387 0.86943525 0.88810829 0.90619664 0.91905736 0.93103091 0.94109532 0.949773322



## Step 4: plot PCA

# Percent of explained variance: 
# fviz_eig(data.pca, addlabels = TRUE)
WGL.MegaA.pca2<-as.data.frame(WGL.MegaA.pca$scores)
WGL.MegaA.pca2$Aneuploid<-c(FALSE, TRUE, FALSE, TRUE, TRUE,             # is aneuploid? 
                            FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, 
                            FALSE, TRUE, TRUE, 
                            FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, 
                            FALSE, TRUE, TRUE, TRUE, TRUE, 
                            FALSE, TRUE, 
                            FALSE, TRUE, 
                            TRUE, FALSE, 
                            TRUE, FALSE, FALSE, 
                            TRUE, FALSE, 
                            TRUE, FALSE)
WGL.MegaA.pca2$CellLine<-c("HCT116", "HCT116","HCT116","HCT116","HCT116",
                           "DLD1","DLD1","DLD1","DLD1","DLD1","DLD1","DLD1", "DLD1","DLD1","DLD1",
                           "A2780", "A2780", "A2780", "A2780", "A2780", "A2780", 
                           "SNU1", "SNU1", "SNU1", "SNU1", "SNU1", 
                           "DLD1","DLD1",
                           "Vaco432", "Vaco432", 
                           "A2780", "A2780",
                           "A2058", "A2058","A2058", 
                           "AGS", "AGS", 
                           "MCF10A", "MCF10A")
WGL.MegaA.pca2$Library<-c("MegaA", "MegaA","MegaA","MegaA","MegaA",
                           "MegaA","MegaA","MegaA","MegaA","MegaA","MegaA","MegaA","MegaA","MegaA","MegaA",
                           "MegaA", "MegaA", "MegaA", "MegaA", "MegaA", "MegaA", 
                           "MegaA", "MegaA", "MegaA", "MegaA", "MegaA", 
                           "WGL","WGL",
                           "WGL", "WGL", 
                           "WGL", "WGL",
                           "WGL", "WGL","WGL", 
                          "WGL","WGL", 
                          "WGL","WGL")

ggplot(WGL.MegaA.pca2, aes(Comp.1, Comp.2, color = CellLine, shape = Library))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values=c("black", "Red", "Blue", "Green", "Grey50", "Purple", "Orange", "Cyan"))
# Plot.PCA.Comp1.Comp2_CellLine.pdf
# 1,3,4, etc. lots of cell line specific stuff



### PRINCIPAL COMPONENT ANALYSIS (PCA):  MegaA only 

## Step 3:  make correlation matrix

Corr_matrix_MegaA<-cor(MegaA.All.Merge.QN[,1:26])


MegaA.pca<-princomp(Corr_matrix_MegaA)
summary(MegaA.pca)

#Importance of components:
#Comp.1    Comp.2     Comp.3     Comp.4     Comp.5     Comp.6     Comp.7     Comp.8      Comp.9     Comp.10
#Standard deviation     0.6871232 0.2815905 0.16465418 0.14358090 0.12329104 0.11159080 0.09080035 0.08581259 0.065476018 0.060119367
#Proportion of Variance 0.7059712 0.1185643 0.04053809 0.03082557 0.02272902 0.01861978 0.01232800 0.01101082 0.006410357 0.005404388
#Cumulative Proportion  0.7059712 0.8245354 0.86507353 0.89589910 0.91862812 0.93724790 0.94957590 0.96058673 0.966997084 0.972401472


## Step 4: plot PCA

# Percent of explained variance: 
MegaA.pca2<-as.data.frame(MegaA.pca$scores)
MegaA.pca2$Aneuploid<-c(FALSE, TRUE, FALSE, TRUE, TRUE,             # is aneuploid? 
                            FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, 
                            FALSE, TRUE, TRUE, 
                            FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, 
                            FALSE, TRUE, TRUE, TRUE, TRUE)
MegaA.pca2$CellLine<-c("HCT116", "HCT116","HCT116","HCT116","HCT116",
                           "DLD1","DLD1","DLD1","DLD1","DLD1","DLD1","DLD1", 
                            "DLD1","DLD1","DLD1",
                           "A2780", "A2780", "A2780", "A2780", "A2780", "A2780", 
                           "SNU1", "SNU1", "SNU1", "SNU1", "SNU1")

ggplot(MegaA.pca2, aes(Comp.1, Comp.2, color = CellLine, shape = Aneuploid))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values=c("black", "Red", "Blue", "Green", "Grey50", "Purple", "Orange"))
# Plot.PCA.MegaA.Comp1.Comp2_CellLine.pdf
# 1,3,4, etc. lots of cell line specific stuff
# sup Fig 1C




### PLOT GENE DROPOUT WGL and DG #####

Merge.Beta.WGL.Mega1<- Merge.ALL.QN2

Merge.Beta.WGL.Mega1$p.value<-NA
Merge.Beta.WGL.Mega1$MeanDiff<-NA
Merge.Beta.WGL.Mega1$MeanBeta<-NA
for (i in 1: length(Merge.Beta.WGL.Mega1$p.value)){
             t<- t.test(as.numeric(Merge.Beta.WGL.Mega1[i,c(1, 3,3, 6,6,6, 10,10, 13,13, 16,16, 19,19, 22,22,22,22,   27, 29, 32, 34,35, 37, 39)]), #control
                        as.numeric(Merge.Beta.WGL.Mega1[i,c(2, 4,5, 7,8,9, 11,12, 14,15, 17,18, 20,21, 23,24,25,26, 28, 30, 31, 33,33, 36, 38)]), #aneuploid
                        paired=TRUE)
  Merge.Beta.WGL.Mega1$p.value[i]<-t$p.value
  Merge.Beta.WGL.Mega1$MeanDiff[i]<-t$estimate*-1
  Merge.Beta.WGL.Mega1$MeanBeta[i]<- mean(as.numeric(Merge.Beta.WGL.Mega1[i,1:39]))
}
subset(Merge.Beta.WGL.Mega1, p.value<0.05 & MeanDiff< -0.15)[c(40,41,42,43)]
#
subset(Merge.Beta.WGL.Mega1, MeanDiff< -0.4)
subset(Merge.Beta.WGL.Mega1, p.value< 0.05 &MeanDiff< -0.2)


HighlightTheseGenes<- subset(Merge.Beta.WGL.Mega1, p.value< 0.05 & MeanDiff< -0.2)$Gene
ggplot(Merge.Beta.WGL.Mega1, aes(MeanDiff, -log2(p.value)) ) +
  geom_point()+
  geom_point(data = subset(Merge.Beta.WGL.Mega1, Gene %in% HighlightTheseGenes), color = "red")+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("MegaA + WGL data (Sheltzer) \nBeta Delta (Aneuploid - Euploid)")+
  ylab("-log2 (p-value)")
# plot.MegaA.WGL.BetaDelta.pvalue_TopHits.pdf
#Figure 2B

subset(Merge.Beta.WGL.Mega1, Gene=="MEF2A")

HighlightTheseGenes<- c("FOSL1", "UBE2H")

ggplot(Merge.Beta.WGL.Mega1, aes(MeanDiff, -log2(p.value)) ) +
  geom_point()+
  geom_point(data = subset(Merge.Beta.WGL.Mega1, Gene %in% HighlightTheseGenes), color = "red")+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("MegaA + WGL data (Sheltzer) \nBeta Delta (Aneuploid - Euploid)")+
  ylab("-log2 (p-value)")
# plot.MegaA.WGL.BetaDelta.pvalue_MitoElectronChain_MitoTranscriptionTranslation.pdf

subset(Merge.Beta.WGL.Mega1, Gene %in% HighlightTheseGenes)

# write.csv(Merge.Beta.WGL.Mega1, "WGL.MegaA_Pvalue.diff.csv", row.names=TRUE)

subset(Merge.Beta.WGL.Mega1, p.value<0.05 & MeanDiff < -0.15)$Gene
# [1] "ANAPC11" "ASCL2"   "ATF6B"   "BHLHE40" "BTBD2"   "CARM1"   "CHD8"    "CSNK1E"  "DLX3"    "EGR3"    "FBXO42"  "FOSL1"   "GATA5"   "GLI4"    "GSK3B"  
#[16] "HM13"    "HOXA1"   "HOXD10"  "HSF1"    "MYT1"    "NACC1"   "NR5A1"   "PHOX2B"  "RAX2"    "SIX6"    "SMAD6"   "SPI1"    "SPOP"    "TBL1XR1" "TWIST1" 
#[31] "UBE2H"   "UBE2M"   "UBE2N"   "UQCRC1"  "VPS18"   "ZFP92"   "ZNF771"  "ZNF784"  "ZSCAN31"



### DG AND WGL ANEUPLOIDY ANALYSIS! #####
##  Analyze both MegaA and WGL data 
##  merge all data. Then compare aneuploid vs euploid
MLE.Merge.ALL

Merge.ALL.QN2



# now we want to add a p-value and an average Beta score difference between aneuploid cells and euploid controls
Merge.ALL.QN3<- Merge.ALL.QN2

Merge.ALL.QN3$ANvEU.P<- NA # t-test. two sided, not paired
Merge.ALL.QN3$ANvEU.P.Paired<- NA # t-test. two sided, not paired
Merge.ALL.QN3$ANvEU.Beta<- NA  #difference btween average an and average euploid beta score (aneuploid beta minus euploid beta)
Merge.ALL.QN3$ANvEU.Beta.Paired<- NA 
Merge.ALL.QN3$Mean.An<- NA 
Merge.ALL.QN3$Mean.Eu<- NA 
Merge.ALL.QN3<- Merge.ALL.QN3 %>% relocate(Gene, .after = last_col()) #make sure Gene is at end of columns


for (i in 1:length(Merge.ALL.QN3$Gene)){
  test<-t.test(Merge.ALL.QN3[i,c(1,3,6,10,13,16,19,22,27,29,32,34,35,37,39)], # Euploid
               Merge.ALL.QN3[i,c(2,4,5,7,8,9,11,12,14,15,17,18,20,21,23,24,25,26,28,30,31,33,36,38)], # Aneuploid
               alternative=c("two.sided")) 
  Merge.ALL.QN3$ANvEU.P[i]<- test$p.value
  Merge.ALL.QN3$ANvEU.Beta[i]<- test$estimate[2] - test$estimate[1]
  Merge.ALL.QN3$Mean.An[i]<- test$estimate[2]
  Merge.ALL.QN3$Mean.Eu[i]<- test$estimate[1]
}

for (i in 1:length(Merge.ALL.QN3$Gene)){
  test<-t.test(as.numeric(Merge.ALL.QN3[i, c(1, 3,3, 6,6,6, 10,10, 13,13, 16,16, 19,19, 22,22,22,22, 27,29,32,34,35,37,39)]), 
              as.numeric(Merge.ALL.QN3[i, c(2,  4,5, 7,8,9, 11,12, 14,15, 17,18, 20,21, 23,24,25,26, 28,30,31,33,33,36,38)]), 
               alternative=c("two.sided"), paired=TRUE) 
  Merge.ALL.QN3$ANvEU.P.Paired[i]<- test$p.value
  Merge.ALL.QN3$ANvEU.Beta.Paired[i]<- -1*test$estimate #make negative values more toxic to aneuploid
}
Merge.ALL.QN3<- Merge.ALL.QN3[order(Merge.ALL.QN3$ANvEU.Beta.Paired),]
Merge.ALL.QN3<- Merge.ALL.QN3 %>% relocate("Gene") #more "Gene" gene names back to begining of columns


SigAllAneu<- subset(Merge.ALL.QN3, ANvEU.Beta.Paired< -0.2 & ANvEU.P.Paired< 0.05)$Gene
ggplot(Merge.ALL.QN3, aes(x= ANvEU.Beta.Paired, 
                                 y= -log2(ANvEU.P.Paired), color= Gene %in% SigAllAneu))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  geom_point(data = subset(Merge.ALL.QN3, Gene%in% SigAllAneu))+
  ylab("-log2(P-value)")+
  xlab("Difference in Beta score\nAneuploid - euploid")+
  scale_color_manual(values = c("Black", "red"), name = "P<0.05 \nDiff < -.2")+
  # geom_point(data = subset(Merge.ALL.QN3, Gene%in% c("UBE2H")), color="blue")+
  theme_classic()
# 5x4
# plot.BetaDelta.ANvEU_WGL.MegaA_June2025.pdf

subset(Merge.ALL.QN3, ANvEU.Beta.Paired< -0.2 & ANvEU.P.Paired< 0.05)$Gene
#  "VPS18"!   "SPI1"    "ANAPC11" "ASCL2"   "UBE2H"!   "UBE2M"!   "RAX2"    "GATA5"   "SMAD6"   "SPOP"    "FOSL1"!   "GLI4"!    "HSF1"   "CHD8"    "SIX6" 



# plot Gene of interest: 
subset(Merge.ALL.QN3, Gene == "FOSL1")


HighlightTheseGenes<- c("FOSL1", "JUNB") 
PanAneuploidSubset <- c("PSMD14", "UBE2H", "UBE2N", "PSMB6", "UBC", "PSMB3") # ubiquitin 2 ligase and proteosome


ggplot(Merge.ALL.QN3, aes(x= ANvEU.Beta.Paired, 
                          y= -log2(ANvEU.P.Paired), color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(Merge.ALL.QN3, Gene%in% HighlightTheseGenes))+
  ylab("-log2(P-value)")+
  xlab("WGL + DG data (Sheltzer lab) \nDifference in Beta score\n (Aneuploid - euploid)")+
  #geom_hline(yintercept = -log2(0.05))+
  #geom_vline(xintercept = -0.15)+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values = c("Black", "Red"), name = "Top Hits")+
  theme_classic()
# 4x4
# plot.BetaDelta.ANvEU.WGL.MegaA.pdf
subset(Merge.ALL.QN3, Gene %in% PanAneuploidSubset)




setwd(ResultsFile)
# write.csv(Merge.ALL.QN3, "CRISPRScreens_BetaScores_MegaA.WGL_QN_BetaDelta_2025.csv")
# Merge.ALL.QN3<- read.csv("CRISPRScreens_BetaScores_MegaA.WGL_QN_BetaDelta_2025.csv")
# Merge.ALL.QN3<- Merge.ALL.QN3[-1,]

### Correlate total Aneuploidy burden (DG and WGL) ####
# Pearson correlation between gene dropout score and aneuploidy score per cell line

colnames(Merge.ALL.QN2) 
Merge.ALL.QN.Corr<- Merge.ALL.QN2[,1:40]

# Aneuploidy score as defined by number of chromosome arms aneuploid divided by mean ploidy
# AneuploidyScore.Gene is defined by number of genes on aneuploid chromosome, divided by mean ploidy
# I do not count X monosomy as aneuploid

AneuploidyScore1 <- c(3/2, 5/2, #number of aneuploid chromosome arms divided by ploidy
                      3/2, 5/2, 5/2, 
                      0, 1/2, 4/2, 4/2, 
                      0, 3/2, 3/2, 
                      0, 3/2, 4/2, 
                      3/2, 5/2, 5/2, 
                      3/2, 5/2, 4/2, 
                      2/4, (2+8)/4, (2+9)/4, (2+6)/4, (2+3)/4, 
                      0, 1/2, # WGL DLD1
                      0, 1/2, # WGL Vaco432
                      3/2, 2/2, 
                      8/2, 7/2, 7/2, 
                      6/2, 5/2, 
                      3/2, 2/2,
                      "AneupPloidyScore")
# Background aneuploidy per cell line: (Not including X and Y copy numbers) 
# HCT116 control = 3/2 (8q, 17q, 10q)
# DLD1 control = 0
# A2780 control = 3/2
# A2058 control = 8/2
# SNU1 control = 2/4
# Vaco432 control = 0 
# AGS control = 6/2
# MCF10A control = 4/2

colnames(Merge.ALL.QN.Corr)


Merge.ALL.QN.Corr$MeanBeta<- NA
Merge.ALL.QN.Corr$AneuScore.Corr<- NA
Merge.ALL.QN.Corr$AneuScore.Corr.P<-NA
NumGene<- length(Merge.ALL.QN.Corr$Gene)

for (i in 1:(NumGene)){
  Corr<- cor.test(as.numeric(AneuploidyScore1[1:39]), as.numeric(Merge.ALL.QN.Corr[i,1:39]), method = 'pearson')
  Merge.ALL.QN.Corr$AneuScore.Corr[i] <-Corr$estimate
  Merge.ALL.QN.Corr$AneuScore.Corr.P[i] <- Corr$p.value
  Merge.ALL.QN.Corr$MeanBeta[i] <- mean(as.numeric(Merge.ALL.QN.Corr[i,1:39]))
}
Merge.ALL.QN.Corr<- Merge.ALL.QN.Corr[order(Merge.ALL.QN.Corr$AneuScore.Corr),]


Merge.ALL.QN.Corr$Gene[1:30]

subset(Merge.ALL.QN.Corr, MeanBeta < -0.2 & AneuScore.Corr < -0.2)$Gene
subset(Merge.ALL.QN.Corr, Gene == "UBC")



GeneOfInterest<- c("UBC", "AneupPloidyScore")
x<- rbind(Merge.ALL.QN.Corr[,1:40], AneuploidyScore1)
x<- subset(x, Gene %in% GeneOfInterest)
x2<- t(x[,1:38])
x2<- data.frame(x2)

ggplot(x2, aes(y=as.numeric(as.character(x2[,1])), x=as.numeric(as.character(x2[,2]))  )) + 
  geom_point()+ 
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+
  ylab("UBC CRISPR score")
# 5x4
# plot.aneuScore.BetaScore.UBC_AllAneuPloidy.pdf
# Sup Fig 4A

cor.test(subset(Merge.ALL.QN.Corr, Gene == "UBC")[1:39],  numeric(AneuploidyScore1[1:39]) )



setwd(ResultsFile)
# write.csv(Merge.ALL.QN.Corr, "Merge.ALL.QN.Corr.csv")



### Compare to MegaA categorical analysis: 
MegaA.WGL.QN.Corr.Diff<- merge(Merge.ALL.QN3, Merge.ALL.QN.Corr, by.x="Gene", by.y="Gene") 
ggplot(MegaA.WGL.QN.Corr.Diff, aes(AneuScore.Corr, ANvEU.Beta.Paired)) + 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()+
  #geom_point(data = subset(MegaA.WGL.QN.Corr.Diff, AneuScore.Corr < -.2 & AneuScore.Corr.P < 0.05 & MeanBeta < -0.2), color = "red") + 
  geom_point(data = subset(MegaA.WGL.QN.Corr.Diff, Gene %in% MitoElectronChainnAssembly), color = "red") + 
  geom_point(data = subset(MegaA.WGL.QN.Corr.Diff, Gene %in% MitoTranslationTranscription), color="gold2") + 
  geom_density_2d(color = "white", bins = 10)+
  ylab("Beta Delta: Sheltzer lab MegaA + WGL \n Aneuploid - Euploid")+
  xlab("Correlation Sheltzer lab MegaA + WGL: \nAneuploidy score Correlation with Beta Score")+
  scale_color_manual(values = c("Black", "red"), name = "Top Hits")+
  theme_classic()
# plot.MegaA.WGLAneuScore.Corr.BetaDelta_RedelectronChain_CyanMitoTranslation.pdf
# 5x5 





##  Now compare to DepMap aneuploidy Score
##  Find genes that are more toxic to aneuploid cells in both datasets

CRISPR_Gene_AneuScore<- merge(Aneuploidy_Scores, CRISPR_Gene, by.x="DepMap_ID", by.y="ModelID")


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



## Now merge my CRISPR screen data with  DepMap data ###
DepMap.Sheltzer.CRISPR_AneuScore<- merge(Merge.ALL.QN.Corr, CRISPR_DepMap_AneuScore, by.x="Gene", by.y="Gene") 


ggplot(DepMap.Sheltzer.CRISPR_AneuScore, aes(x=AneuScore.Corr, y=Aneu.divPloidy.CRISPR.Corr.Coef ))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_AneuScore, 
                           AneuScore.Corr < -0.1  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & AneuScore.Corr.P <0.1), color = "red4") + 
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_AneuScore, 
                           AneuScore.Corr < -0.1  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & AneuScore.Corr.P <0.05), color = "red") + 
  #geom_point(data = subset(DepMap.Sheltzer.CRISPR_AneuScore, Gene %in%  c( "UBE2H", "UBC")), color = "Blue") + 
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+ 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("Sheltzer Lab Correlation: \n Aneuploidy Score correlation with Beta Score")+
  ylab("DepMap Correlation: \n Aneuploidy Score correlation with Beta Score")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.Sheltzer.CRISPR_MegaAWGLAneuScore_DepMapAneudivPloidy.pdf
# SupFig 4C

subset(DepMap.Sheltzer.CRISPR_AneuScore, AneuScore.Corr < -0.2  & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1)$Gene






## Plot specific genes: 
HighlightTheseGenes<- c( "UBE2H", "UBC") 

ggplot(DepMap.Sheltzer.CRISPR_AneuScore, aes(y=AneuScore.Corr, x=Aneu.divPloidy.CRISPR.Corr.Coef, 
                                             color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_AneuScore, Gene %in% HighlightTheseGenes))+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_density_2d(color = "white", bin = 10)+ 
  xlab("Correlation between Gene Beta Score (DepMap): \n Aneuploidy Score")+
  ylab("Correlation between Gene Beta Score (Sheltzer): \n Aneuploidy Score")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4
# plot.DepMap.Sheltzer.MegaAWGL.Correlation.CRISPR_AneuScore.pdf
subset(DepMap.Sheltzer.CRISPR_AneuScore, Gene %in% HighlightTheseGenes)


setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.CRISPR_AneuScore, "Sheltzer.CRISPR.correlate.AneuScore.AllAneudivPloidy.DepMapCorr.csv")
# DepMap.Sheltzer.CRISPR_AneuScore<- read.csv("Sheltzer.CRISPR.correlate.AneuScore.AllAneudivPloidy.DepMapCorr.csv")



    ##### Linear regression aneuploidy score & Depmap data ####
# Linear regression analysis 
# linear model, controlling for cell line

colnames(Merge.ALL.QN2) 
Merge.ALL.QN.lm<- Merge.ALL.QN2[,1:40]

# Aneuploidy score as defined by number of chromosome arms aneuploid divided by mean ploidy
# AneuploidyScore.Gene is defined by number of genes on aneuploid chromosome, divided by mean ploidy
# I do not count X monosomy as aneuploid

AneuploidyScore1 <- c(3/2, 5/2, #number of aneuploid chromosome arms divided by ploidy
                      3/2, 5/2, 5/2, 
                      0, 1/2, 4/2, 4/2, 
                      0, 3/2, 3/2, 
                      0, 3/2, 4/2, 
                      3/2, 5/2, 5/2, 
                      3/2, 5/2, 4/2, 
                      2/4, (2+8)/4, (2+9)/4, (2+6)/4, (2+3)/4, 
                      #0, 1/2, # HCT116 Ts13
                      0, 1/2, 
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
# SNU1 control = 2/4
# Vaco432 control = 0 
# AGS control = 6/2
# MCF10A control = 4/2

CellType <- c("HCT116", "HCT116", #number of aneuploid chromosome arms divided by ploidy
              "HCT116", "HCT116", "HCT116", 
              "DLD1", "DLD1", "DLD1", "DLD1", 
              "DLD1", "DLD1", "DLD1", 
              "DLD1", "DLD1", "DLD1", 
              "A2780", "A2780","A2780",
              "A2780","A2780","A2780",
              "SNU1", "SNU1", "SNU1", "SNU1", "SNU1", 
              "DLD1", "DLD1", 
              "Vaco432", "Vaco432",
              "A2780", "A2780", 
              "A2058", "A2058", "A2058", 
              "AGS", "AGS",
              "MCF10A", "MCF10A",
              "CellType")

Merge.ALL.QN.lm<- rbind(CellType, Merge.ALL.QN.lm)
Merge.ALL.QN.lm<- rbind(AneuploidyScore1, Merge.ALL.QN.lm)


Merge.ALL.QN.lm$MeanBeta<- NA
Merge.ALL.QN.lm$AneuScore.lm<- NA
Merge.ALL.QN.lm$AneuScore.lm.R2<-NA
Merge.ALL.QN.lm$AneuScore.lm.P<-NA
NumGene<- length(Merge.ALL.QN.lm$Gene)

for (i in 3:(NumGene)){
  lm_CellLine<- lm(as.numeric(Merge.ALL.QN.lm[i,1:39]) ~ as.numeric(AneuploidyScore1[1:39]) + CellType[1:39]) # linear model (gene Beta score ~ aneuploidy score + Cell type)
  Merge.ALL.QN.lm$AneuScore.lm[i] <- lm_CellLine$coefficients[2]
  Merge.ALL.QN.lm$AneuScore.lm.R2[i] <- summary(lm_CellLine)$r.squared # Multiple R squared. How much does cell types + Aneuploidy score predict protein dropout?
  Merge.ALL.QN.lm$AneuScore.lm.P[i] <- summary(lm_CellLine)$coefficients[2,4] # p-value (column 4) for second attribute (aneuploidy score)
  Merge.ALL.QN.lm$MeanBeta[i] <- mean(as.numeric(Merge.ALL.QN.lm[i,1:39]))
}
Merge.ALL.QN.lm<- Merge.ALL.QN.lm[order(Merge.ALL.QN.lm$AneuScore.lm),]


Merge.ALL.QN.lm$Gene[1:30]
# [1] "ATF4"     "FOXD4L5"  "FXR1"     "PKMYT1"   "KLF10"    "SPI1"     "TRIM28"   "VPS18"    "ADAM10"   "GLI4"     "FIZ1"     "RNF166"  
#[13] "STK11"    "HLX"      "HES1"     "ASCL2"    "PRSS33"   "TRIM61"   "ATF6B"    "MYT1"     "GATA5"    "TMPRSS12" "KLF16"    "MAP4K4"  
#[25] "FOXD4"    "MZF1"     "RNF168"   "LITAF"    "TRIO"     "HSF1"


# Linear model negative coefficient
subset(Merge.ALL.QN.lm, MeanBeta < -0.1 & AneuScore.lm < -0.15 & AneuScore.lm.P < 0.1 & AneuScore.lm.R2 > 0.4)$Gene
#  "FOXD4L5"  "KLF10"    "VPS18"    "GLI4"     "TRIM28"   "PRSS33"   "ATF6B"    "GATA5"    "RNF168"   "TRIM49D1" "NACC1"


#Linear model positive coefficient: 
subset(Merge.ALL.QN.lm, MeanBeta < -0.1 & AneuScore.lm > 0.15 & AneuScore.lm.P < 0.1 & AneuScore.lm.R2 > 0.4)$Gene
#[1] "ANAPC1"  "CDK2"    "ELOC"    "HIPK1"   "COPS5"   "INO80"   "BRD2"    "HUWE1"   "BRAP"    "DCAF13"  "UBE2V1"  "PSMA3"   "DCAF1"   "NAT10"   "ZNF236" 
#[16] "CHD2"    "USPL1"   "ZNF574"  "SETDB1"  "UHRF1"   "DNAJC2"  "SUPT16H" "GTF3C4"  "ZPR1"    "SKP1"    "MIPEP"   "CHD4"    "HINFP"   "MTOR"    "PLK1"   
#[31] "MARCHF5" "POLR2D"  "CDK9"    "SSRP1"   "UBE2I"   "UBA52"   "UBE2D3"  "BUB1B"   "SMG1"    "BPTF"    "METAP2"  "CEBPZ"   "CDK1"    "CUL3"    "PSMA5"  
#[46] "ELP3"    "PSMA6"   "POLR2A"  "ELOB"    "ALG13"   "ZNF626"  "PRMT5"  


subset(Merge.ALL.QN.lm, MeanBeta < -0.2 & AneuScore.lm < -0.2 & AneuScore.lm.P < 0.01)$Gene
# "VPS18" "GLI4" 


ggplot(Merge.ALL.QN.lm, aes(AneuScore.lm, -log2(AneuScore.lm.P))) + 
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()+
  geom_point(data = subset(Merge.ALL.QN.lm, Gene %in% HighlightTheseGenes), color = "Blue") + 
  geom_point(data = subset(Merge.ALL.QN.lm, Gene %in% MitoElectronChainnAssembly), color = "red") + 
  geom_point(data = subset(Merge.ALL.QN.lm, Gene %in% MitoTranslationTranscription), color="gold2") + 
  geom_density_2d(color = "white")+
  ylab("-log2(p-value)")+
  xlab("Linear regression (DG + WGL)\nAneuploidy Score correcting for tissue type")+
  scale_color_manual(values = c("Black", "red"), name = "Top Hits")+
  theme_classic()
# plot.AneuScore.pvalue.lm_Feb2024_CtrlCellType.pdf
# 5x5 


GeneOfInterest<- c("UBC", "AneupPloidyScore")
x<- rbind(Merge.ALL.QN.lm[,1:40], AneuploidyScore1)
x<- rbind(x, CellType)
x<- subset(x, Gene %in% GeneOfInterest)
x2<- t(x[,1:39])
x2<- data.frame(x2)

ggplot(x2, aes(y=as.numeric(as.character(x2[,1])), x=as.numeric(as.character(x2[,2])) )) + 
  geom_point()+ 
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_smooth(method="lm")+
  xlab("Aneuploidy Score")+
  ylab("CRISPR score")
# 5x4
# plot.aneuScore.BetaScore.UBE2H_lmCellType.pdf
# CtrlCellZero 




##  Now compare to DepMap aneuploidy Score
##  Find genes that are more toxic to aneuploid cells in both datasets

setwd(SchukkenData)

CRISPR_DepMap_AneuScore_lm_CancerType<- read.csv("DepMap.CRISPR.AneuScore.lm.Subtype.csv", header=TRUE)
CRISPR_DepMap_AneuScore_lm_CancerType$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_AneuScore_lm_CancerType$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_AneuScore_lm_CancerType$Gene<- sub("[.]", "-", CRISPR_DepMap_AneuScore_lm_CancerType$Gene) # replace single period with -. 



## Now merge my CRISPR screen data with linear regression DepMap data ###
DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine<- merge(Merge.ALL.QN.lm, CRISPR_DepMap_AneuScore_lm_CancerType, by.x="Gene", by.y="Gene") 

ggplot(DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine, aes(x=AneuScore.lm.x, y=AneuScore.lm.y, 
                                                         color= AneuScore.lm.x < -0.05 & AneuScore.lm.y < -0.002))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine, AneuScore.lm.x < -0.05 & AneuScore.lm.y < -0.002))+
  theme_classic()+
  geom_hline(yintercept = 0)+  
  geom_density_2d(color = "white")+
  geom_vline(xintercept = 0)+
  xlab("Linear Regression Beta Score and Aneuploidy Score: \n Sheltzer Lab Screens")+
  ylab("Linear Regression Beta Score and Aneuploidy Score: \n DepMap")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine.pdf


subset(DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine, AneuScore.lm.x < -0.05 & AneuScore.lm.y < -0.002)$Gene
# [1] "AFG3L2"  "ANAPC15" "FBXO42"  "FBXO5"   "FOSL1"   "HOXC10"  "LONP1"   "OLIG3"   "PKMYT1"  "PSMB6"   "PSMD14"  "RELA"    "SYVN1"   "UBE2H"   "USP28"  
# [16] "USP5"    "VDR"     "ZNF763" 


HighlightTheseGenes<- c("FOSL1")

ggplot(DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine, aes(x=AneuScore.lm.x,  y=AneuScore.lm.y, 
                                                         color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine, Gene %in% HighlightTheseGenes))+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_density_2d(color = "white", bins = 9) +
  ylab("Pearson Correlation (DepMap): \n Beta Score and Aneuploidy Score")+
  xlab("Linear regression (Sheltzer): \n Aneuploidy Score contrlling for cell line")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4
# plot.DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine.pdf
subset(DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine, Gene %in% HighlightTheseGenes)

#write.csv(DepMap.Sheltzer.CRISPR_AneuScore_lm_cellLine, "Sheltzer.CRISPR.correlate.AneuScore.lmCellLine.DepMapCorr.csv")




### Paired DG + WGL compared to DepMap correlation ####

### See "correlation with aneuploidy score section to get CRISPR DepMap correlation
setwd(ResultsFile)
# write_csv(CRISPR_DepMap_AneuScore, "DepMap_CRISPR_AneuScore_Correlation.csv")
#  CRISPR_DepMap_AneuScore<- read_csv("DepMap_CRISPR_AneuScore_Correlation.csv")


## Now merge my CRISPR screen data with aneuscore DepMap data ###
DepMap.CRISPR_AneuScore_MegaWGLTtest<- merge(Merge.ALL.QN3, CRISPR_DepMap_AneuScore, by.x="Gene", by.y="Gene") 

ggplot(DepMap.CRISPR_AneuScore_MegaWGLTtest, aes(x=ANvEU.Beta.Paired, y=Aneu.divPloidy.CRISPR.Corr.Coef))+
  geom_point()+
  geom_point(data = subset(DepMap.CRISPR_AneuScore_MegaWGLTtest, 
                           ANvEU.Beta.Paired < -0.2 & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & ANvEU.P.Paired< 0.1 ), color = "red4")+
  geom_point(data = subset(DepMap.CRISPR_AneuScore_MegaWGLTtest, 
                           ANvEU.Beta.Paired < -0.2 & Aneu.divPloidy.CRISPR.Corr.Coef < -0.1 & 
                             Aneu.divPloidy.CRISPR.Corr.P< 0.005 & ANvEU.P.Paired< 0.05 ), color = "red")+
  #geom_point(data = subset(DepMap.CRISPR_AneuScore_MegaWGLTtest, Gene %in% c("UBC")), color = "blue")+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_density_2d(color = "white", bins = 10)+
  ylab("DepMap Correlation:\n Aneuploidy Score correlation with Beta Score")+
  xlab("Beta Delta (Sheltzer lab): \n Aneuploid - Euploid")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.CRISPR_AneuScore_MegaWGLTtest.pdf
# Figure XXB

subset(DepMap.CRISPR_AneuScore_MegaWGLTtest, 
       ANvEU.Beta.Paired < -0.15 & Aneu.divPloidy.CRISPR.Corr.Coef < -0.15 & 
         Aneu.divPloidy.CRISPR.Corr.P< 0.005 & ANvEU.P.Paired< 0.1 )$Gene
# Top Top hits: "FOSL1"  "PCNA"   "UBC"    "UBE2H"  "ZNF207"


#Highlight gene of interest: 
ggplot(DepMap.CRISPR_AneuScore_MegaWGLTtest, aes(x=ANvEU.Beta.Paired, y=Aneu.divPloidy.CRISPR.Corr.Coef))+
  geom_point()+
  #geom_point(data = subset(DepMap.CRISPR_AneuScore_MegaWGLTtest, ANvEU.Beta.Paired < -0.15 & Aneu.divPloidy.CRISPR.Corr.Coef < -0.15), color = "red")+
  geom_point(data = subset(DepMap.CRISPR_AneuScore_MegaWGLTtest, Gene %in% c("UBE2H")), color = "red")+
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_density_2d(color = "white", bins = 10)+
  ylab("DepMap Correlation:\n Aneuploidy Score/ Ploidy correlation with Beta Score")+
  xlab("Beta Delta (Sheltzer lab): \n Aneuploid - Euploid")+
  scale_color_manual(values= c("Black", "Red"), name="Top hits")
# 5x4
# plot.DepMap.CRISPR_AneuScore_MegaWGLTtest.pdf

subset(DepMap.CRISPR_AneuScore_MegaWGLTtest, ANvEU.Beta.Paired < -0.2 & 
         Aneu.divPloidy.CRISPR.Corr.Coef < -0.15 )$Gene
# "FOSL1"  "PCNA"   "UBC"    "UBE2H"  "ZNF207"


# write.csv(DepMap.CRISPR_AneuScore_MegaWGLTtest, "CRISPRScreens_BetaScores_MegaA.WGL_QN_BetaDelta_DepMapCorr_Jun2024.csv")


#### Chromosome specific analysis ####

# CRISPR_Gene 
Chrm10Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 10)$`HGNC symbol`)
Chrm8Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 8)$`HGNC symbol`)
Chrm5Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 5)$`HGNC symbol`)
Chrm2Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 2)$`HGNC symbol`)



    ##### Ts10 GAIN ONLY (3 pairs) (QN) #####
# DLD1 Ts10.21, DLD1 Ts8.10, A2780 Ts10 

MegaA.Merge.ALL #from "Quantile normalization ALL" section


colnames(MegaA.Merge.ALL)
MLE.Merge.Ts10GAIN<- MegaA.Merge.ALL[, c(1,7,9,11,12,17,19)]
colnames(MLE.Merge.Ts10GAIN)<- c("Gene", "MegaA.D.Cont1.beta", "MegaA.D.Ts8.10.beta", 
                             "MegaA.D.Cont3.beta", "MegaA.D.10.21.beta",  
                             "MegaA.A2780.Cont.beta", "MegaA.A2780.Ts10.beta")



## Quantile normalize by batch and cell line: 
Merge.Ts10GAIN.QN # quantile normalized across all samples from MegaA and WGL
Merge.Ts10GAIN.QN<- quantile_normalisation(MLE.Merge.Ts10GAIN[,c(2:length(MLE.Merge.Ts10GAIN))])
Merge.Ts10GAIN.QN<- as.data.frame(Merge.Ts10GAIN.QN)
Merge.Ts10GAIN.QN$Gene<- (MLE.Merge.Ts10GAIN$Gene)
head(Merge.Ts10GAIN.QN)



###
Merge.Ts10.QN3<- Merge.Ts10GAIN.QN
Merge.Ts10.QN3$ANvEU.P<- NA # t-test. two sided, not paired
Merge.Ts10.QN3$ANvEU.Beta<- NA  #difference btween average an and average euploid beta score (aneuploid beta minus euploid beta)
Merge.Ts10.QN3$Mean.An<- NA   
Merge.Ts10.QN3$Mean.Eu<- NA 

for (i in 1:length(Merge.Ts10.QN3$Gene)){
  test<-t.test(as.numeric(Merge.Ts10.QN3[i,c(1,3,5)]), #controls and/or loss 10
               as.numeric(Merge.Ts10.QN3[i,c(2,4,6)]), # Gain10 or neutral (relative to loss)
               paired=TRUE) 
  Merge.Ts10.QN3$ANvEU.P[i]<- test$p.value
  Contbeta<- mean(as.numeric(Merge.Ts10.QN3[i,c(1,3,5)]))
  Ts10beta<- mean(as.numeric(Merge.Ts10.QN3[i,c(2,4,6)]))
  Merge.Ts10.QN3$ANvEU.Beta[i]<- Ts10beta - Contbeta
  Merge.Ts10.QN3$Mean.An[i]<- Ts10beta
  Merge.Ts10.QN3$Mean.Eu[i]<-  Contbeta
}
Merge.Ts10.QN3<-Merge.Ts10.QN3[order(Merge.Ts10.QN3$ANvEU.Beta),]


Mega.WGL.Ts10.Only.TopHits<- subset(Merge.Ts10.QN3, ANvEU.P<0.05 & ANvEU.Beta< -0.2)$Gene

ggplot(Merge.Ts10.QN3, aes(x= ANvEU.Beta, 
                           y= -log2(ANvEU.P), color= Gene %in% Mega.WGL.Ts10.Only.TopHits))+
  geom_point()+
  ylab("-log2(P-value)")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Beta score\n Gain 10 - Neutral 10")+
  scale_color_manual(values = c("Black", "red"), name = "P<0.05 \nDiff < -.2")+
  theme_classic()
# 5x4
# plot.BetaDelta.GAIN10vEU.pdf

subset(Merge.Ts10.QN3, ANvEU.P<0.05 & ANvEU.Beta< -0.2)$Gene


### Look at genes located ON chromosome 10: 
Chrm10Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 10)$`HGNC symbol`)
TopHits<- Chrm10Genes
ggplot(Merge.Ts10.QN3, aes(x= ANvEU.Beta, 
                           y= -log2(ANvEU.P), color= Gene %in% TopHits))+
  geom_point()+
  geom_point(data = subset(Merge.Ts10.QN3, Gene%in% TopHits))+
  ylab("-log2(P-value)")+
  xlab("Difference in Beta score\nTs10 - euploid")+
  scale_color_manual(values = c("Black", "Red"), name = "Top Hits")+
  theme_classic()
# 5x4
# plot.BetaDelta.Ts10vEU.pdf
subset(DepMap.Sheltzer.CRISPR_Ts13, Gene %in% HighlightTheseGenes)


         ###### Ts10 GAIN & DepMap 10q-gain data #####

CN10q<- Arm_CN[, c("X", "X10q")]
colnames(CN10q)<- c("Cell_Line", "10q")
Cell_Neutral10q<- subset(CN10q, `10q`=="0")$Cell_Line
Cell_Gain10q<- subset(CN10q, `10q`=="1")$Cell_Line
Cell_Loss10q<- subset(CN10q, `10q`=="-1")$Cell_Line
Cell_Loss.Neutral10q<- subset(CN10q, `10q`!="1")$Cell_Line


CRISPR_DepMap_Ts10q<-data.frame(GeneID= character(), 
                                Gain.LossNeutral.P = numeric(), 
                                Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch10qGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain10q) #get drug result for all cells with gain for chrm 13
  Drug_Ch10qNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral10q) #get drug result for all cells neutral for chrm 13
  Drug_Ch10qLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral10q) #get drug result for all cells neutral or loss for chrm 13
  
  if (length(na.omit(Drug_Ch10qGain[,2])) >5 & 
      length(na.omit(Drug_Ch10qNeutral[,2])) >5 &
      length(na.omit(Drug_Ch10qLoss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch10qGain[,2], Drug_Ch10qNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch10qGain[,2], Drug_Ch10qLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts10q<- rbind(CRISPR_DepMap_Ts10q, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                                Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                                Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                                Gain.LossNeutral.P = test2$p.value, 
                                                                Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 
#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts10q<- CRISPR_DepMap_Ts10q[order(CRISPR_DepMap_Ts10q$Gain.Neutral.P),]
CRISPR_DepMap_Ts10q[1:20,]


CRISPR_DepMap_Ts10q$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts10q$GeneID)){
  if (CRISPR_DepMap_Ts10q$Gain.LossNeutral.P[i]<0.01 && CRISPR_DepMap_Ts10q$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts10q$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts10q, SigDiff=="TRUE")
CRISPR_DepMap_Ts10q$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_Ts10q$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_Ts10q$Gene<- sub("[.]", "-", CRISPR_DepMap_Ts10q$Gene) # replace single period with -. 

head(CRISPR_DepMap_Ts10q)


## Now merge my CRISPR screen data with Ts10q DepMap data ###
DepMap.SheltzerGAIN.CRISPR_Ts10q<- merge(Merge.Ts10.QN3, CRISPR_DepMap_Ts10q, by.x="Gene", by.y="Gene") 



         ###### Ts10 GAIN & DepMap 10p-gain data #####

CN10p<- Arm_CN[, c("X", "X10p")]
colnames(CN10p)<- c("Cell_Line", "10p")
Cell_Neutral10p<- subset(CN10p, `10p`=="0")$Cell_Line
Cell_Gain10p<- subset(CN10p, `10p`=="1")$Cell_Line
Cell_Loss10p<- subset(CN10p, `10p`=="-1")$Cell_Line
Cell_Loss.Neutral10p<- subset(CN10p, `10p`!="1")$Cell_Line


CRISPR_DepMap_Ts10p<-data.frame(GeneID= character(), 
                                Gain.LossNeutral.P = numeric(), 
                                Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch10pGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain10p) #get drug result for all cells with gain for chrm 13
  Drug_Ch10pNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral10p) #get drug result for all cells neutral for chrm 13
  Drug_Ch10pLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral10p) #get drug result for all cells neutral or loss for chrm 13
  
  if (length(na.omit(Drug_Ch10pGain[,2])) >5 & 
      length(na.omit(Drug_Ch10pNeutral[,2])) >5 &
      length(na.omit(Drug_Ch10pLoss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch10pGain[,2], Drug_Ch10pNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch10pGain[,2], Drug_Ch10pLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts10p<- rbind(CRISPR_DepMap_Ts10p, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                                Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                                Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                                Gain.LossNeutral.P = test2$p.value, 
                                                                Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 
#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts10p<- CRISPR_DepMap_Ts10p[order(CRISPR_DepMap_Ts10p$Gain.Neutral.P),]
CRISPR_DepMap_Ts10p[1:20,]


CRISPR_DepMap_Ts10p$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts10p$GeneID)){
  if (CRISPR_DepMap_Ts10p$Gain.LossNeutral.P[i]<0.01 && CRISPR_DepMap_Ts10p$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts10p$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts10p, SigDiff=="TRUE")
CRISPR_DepMap_Ts10p$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_Ts10p$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_Ts10p$Gene<- sub("[.]", "-", CRISPR_DepMap_Ts10p$Gene) # replace single period with -. 


head(CRISPR_DepMap_Ts10p)



colnames(CRISPR_DepMap_Ts10p)<- c("GeneID", "Gain10p.Neutral.P",  "Gain10p.min.Neutral.Diff",   "Gain10p.LossNeutral.P", "Gain10p.min.LossNeutral.Diff", "SigDiff.DepMap10p", "Gene")
colnames(CRISPR_DepMap_Ts10q)<- c("GeneID", "Gain10q.Neutral.P",  "Gain10q.min.Neutral.Diff",   "Gain10q.LossNeutral.P", "Gain10q.min.LossNeutral.Diff", "SigDiff.DepMap10q", "Gene")

## Now merge my CRISPR screen data with Ts10 DepMap data ###
DepMap.SheltzerGAIN.CRISPR_Ts10pq<- merge(Merge.Ts10.QN3, CRISPR_DepMap_Ts10p, by.x="Gene", by.y="Gene") 
DepMap.SheltzerGAIN.CRISPR_Ts10pq<- merge(DepMap.SheltzerGAIN.CRISPR_Ts10pq, CRISPR_DepMap_Ts10q, by.x="Gene", by.y="Gene") 




setwd(ResultsFile)
# write.csv(DepMap.SheltzerGAIN.CRISPR_Ts10pq, "Chrm10pq_Specific_Toxicity_Sheltzer.DepMap.csv")
# DepMap.SheltzerGAIN.CRISPR_Ts10pq<- read.csv("Chrm10pq_Specific_Toxicity_Sheltzer.DepMap.csv")




### Plot Difference Ts10 with CRISPR values highlighted

DepMap10pHits <- subset(CRISPR_DepMap_Ts10p, Gain10p.min.LossNeutral.Diff< -0.1 & Gain10p.LossNeutral.P< 0.05)$Gene
DepMap10qHits <- subset(CRISPR_DepMap_Ts10q, Gain10q.min.LossNeutral.Diff< -0.1 & Gain10q.LossNeutral.P< 0.05)$Gene

ggplot(Merge.Ts10.QN3, aes(x= ANvEU.Beta, 
                           y= -log2(ANvEU.P)))+
  geom_point()+
  geom_point(data = subset(Merge.Ts10.QN3, Gene %in% DepMap10pHits), color = "orchid1")+ 
  geom_point(data = subset(Merge.Ts10.QN3, Gene %in% DepMap10qHits), color = "orange1")+ 
  #geom_point(data = subset(Merge.Ts10.QN3, Gene %in% DepMap10pHits & Gene %in% DepMap10qHits), color = "red")+ 
  #geom_point(data = subset(Merge.Ts10.QN3, Gene %in% c("PRMT1")), color = "blue")+ 
  ylab("-log2(P-value)")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Beta score\nTs10 - Di10")+
  theme_classic()
# 5x4
# plot.BetaDelta.Ts10vEU.colorDepMap_SigDif.pdf

#Figure 4B



     ##### Ts8 only (4 pairs) (QN) #####
# DLD1 Ts8.10, HCT115 Ts8, A2780 Ts8, SNU1 Ps8.12 Ts10.6p (C111)
### Look at genes located ON chromosome 8: 


MLE.Merge.Ts8<- merge(HCT116.Cont1.MegaA[,c(1,3)], HCT116.Ts8.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MLE.Merge.Ts8<- merge(MLE.Merge.Ts8, DLD1.Cont1.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MLE.Merge.Ts8<- merge(MLE.Merge.Ts8, DLD1.Ts8.10.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MLE.Merge.Ts8<- merge(MLE.Merge.Ts8, A2780.Cont2.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MLE.Merge.Ts8<- merge(MLE.Merge.Ts8, A2780.Ts8.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene"))
MLE.Merge.Ts8<- merge(MLE.Merge.Ts8, Snu1.Cont1.C13.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MLE.Merge.Ts8<- merge(MLE.Merge.Ts8, Snu1.C111.MegaA[,c(1,3)], by=c("Gene", "Gene"))

colnames(MLE.Merge.Ts8)<- c("Gene", "HCT116 Cont1", "HCT116 Ts8", 
                             "DLD1 Cont1", "DLD1 Ts8.10", 
                             "A2780 Cont2", "A2780 Ts8", 
                             "Snu1 Cont2", "Snu1 Ts10.6q Ps8.12")

#write.csv(MLE.Merge.Ts8, "rawData.BetaScore.Gain8Cells.csv")

## Quantile normalize by batch and cell line: 
Merge.Ts8.QN # quantile normalized across all samples from MegaA and WGL
Merge.Ts8.QN<- quantile_normalisation(MLE.Merge.Ts8[,c(2:length(MLE.Merge.Ts8))])
Merge.Ts8.QN<- as.data.frame(Merge.Ts8.QN)
Merge.Ts8.QN$Gene<- (MLE.Merge.Ts8$Gene)
head(Merge.Ts8.QN)

### Plot pos controls of Ts10 only
PlotHeatmap.Cont(DropoutData= list("HCT116 Ctrl MegaA"= Merge.Ts8.QN[,c(9,1,1)],
                                   "HCT116 Ts8 MegaA"= Merge.Ts8.QN[,c(9,1,2)], 
                                   "DLD1 Ctrl MegaA"= Merge.Ts8.QN[,c(9,1,3)], 
                                   "DLD1 Ts8.10 MegaA"= Merge.Ts8.QN[,c(9,1,4)], 
                                   "A2780 Ctrl MegaA"= Merge.Ts8.QN[,c(9,1,5)], 
                                   "A2780 Ts8 MegaA"= Merge.Ts8.QN[,c(9,1,6)], 
                                   "Snu1 Ctrl MegaA"= Merge.Ts8.QN[,c(9,1,7)], 
                                   "Snu1 Ts10.6q Ps8.12 MegaA"= Merge.Ts8.QN[,c(9,1,8)]
),
Title="Control guide Beta-scores")
# 7x4
# plot.Ts8Only.QN.MegaA.WGL.pdf

negControls<- Merge.Ts8.QN[which(grepl("neg", Merge.Ts8.QN$Gene)),]$Gene #list all negative controls

###
Merge.Ts8.QN2<- Merge.Ts8.QN
Merge.Ts8.QN2<- subset(Merge.Ts8.QN2, ! Gene %in% negControls)
Merge.Ts8.QN2$ANvEU.P<- NA # t-test. two sided, not paired
Merge.Ts8.QN2$ANvEU.Beta<- NA  #difference btween average an and average euploid beta score (aneuploid beta minus euploid beta)
Merge.Ts8.QN2$Mean.An<- NA   
Merge.Ts8.QN2$Mean.Eu<- NA 

for (i in 1:length(Merge.Ts8.QN2$Gene)){
  test<-t.test(as.numeric(Merge.Ts8.QN2[i,c(1,3,5,7)]), #controls
               as.numeric(Merge.Ts8.QN2[i,c(2,4,6,8)]), # Ts13
               paired=TRUE) 
  Merge.Ts8.QN2$ANvEU.P[i]<- test$p.value
  Contbeta<- mean(as.numeric(Merge.Ts8.QN2[i,c(1,3,5,7)]))
  Ts8beta<- mean(as.numeric(Merge.Ts8.QN2[i,c(2,4,6,8)]))
  Merge.Ts8.QN2$ANvEU.Beta[i]<- Ts8beta - Contbeta
  Merge.Ts8.QN2$Mean.An[i]<-    Ts8beta
  Merge.Ts8.QN2$Mean.Eu[i]<-    Contbeta
}
Merge.Ts8.QN2<-Merge.Ts8.QN2[order(Merge.Ts8.QN2$ANvEU.Beta),]


Mega.WGL.Ts8.Only.TopHits<- subset(Merge.Ts8.QN2, ANvEU.P<0.05 & ANvEU.Beta< -0.2)$Gene

ggplot(Merge.Ts8.QN2, aes(x= ANvEU.Beta, 
                           y= -log2(ANvEU.P), color= Gene %in% Mega.WGL.Ts8.Only.TopHits))+
  geom_point()+
  ylab("-log2(P-value)")+
  xlab("Difference in Beta score\nTs8 - euploid")+
  geom_density_2d(color = "white", bins = 10)+
  scale_color_manual(values = c("Black", "red"), name = "P<0.05 \nDiff < -.2")+
  theme_classic()
# 5x4
# plot.BetaDelta.Ts8vEU.pdf


subset(Merge.Ts8.QN2, ANvEU.P<0.05 & ANvEU.Beta< -0.2)$Gene



ggplot(subset(Merge.Ts8.QN2, !Gene %in% Chrm8Genes), 
       aes(x=ANvEU.Beta, y= -log2(ANvEU.P),
           color=Gene %in% Chrm8Genes ))+
  geom_point()+ 
  geom_point(data = subset(Merge.Ts8.QN2, Gene %in% Chrm8Genes))+
  scale_color_manual(values=c("Black", "red"), name= "On Chrm 8")+ 
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Beta delta\n <- more toxic gain 8   | more toxic neutral 8 ->")+
  ylab("Density")
# 5x2
# plot.BetaDelta.Density_On8.pdf

x<- t.test(subset(Merge.Ts8.QN2, !Gene %in% Chrm8Genes)$ANvEU.Beta, 
           subset(Merge.Ts8.QN2, Gene %in% Chrm8Genes)$ANvEU.Beta)
# A2780 & AGS: P = NS 




# Image Ts8 MegaA  gene dropout

test<-subset(Merge.Ts8.QN2, Gene=="SH3RF3")
df<- data.frame(Control=as.numeric(test[,c(1,3,5,7)]), Ts8=as.numeric(test[,c(2,4,6,8)]), Name=colnames(test[,c(2,4,6,8)]))
ggpaired(df, cond1 = "Control", cond2 = "Ts8", 
         title="SH3RF3", fill="condition")+ 
  ylab("Quantile Norm. Beta")+
  scale_fill_manual(values=c("darkolivegreen1", "skyblue"))
# 3x5
# plot.pairedLine.Ts8_ADCK5.pdf



setwd(ResultsFile)
# write.csv(Merge.Ts8.QN2, "QN.Beta.Gain8.Difference.pvalue.csv")
# Merge.Ts8.QN2<- read.csv( "QN.Beta.Gain8.Difference.pvalue.csv")
# Merge.Ts8.QN2<- Merge.Ts8.QN2[,-1]

# write.csv(Merge.Ts8.QN2$Gene, "GeneList_MegaA.csv")



         ###### Ts8 correlate with DepMap 8q-gain data #####

CN8q<- Arm_CN[, c("X", "X8q")]
colnames(CN8q)<- c("Cell_Line", "8q")
Cell_Neutral8q<- subset(CN8q, `8q`=="0")$Cell_Line
Cell_Gain8q<- subset(CN8q, `8q`=="1")$Cell_Line
Cell_Loss8q<- subset(CN8q, `8q`=="-1")$Cell_Line
Cell_Loss.Neutral8q<- subset(CN8q, `8q`!="1")$Cell_Line


CRISPR_DepMap_Ts8q<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch8qGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain8q) #get drug result for all cells with gain for chrm 8
  Drug_Ch8qNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral8q) #get drug result for all cells neutral for chrm 8
  Drug_Ch8qLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral8q) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch8qGain[,2])) >5 & 
      length(na.omit(Drug_Ch8qNeutral[,2])) >5 &
      length(na.omit(Drug_Ch8qLoss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch8qGain[,2], Drug_Ch8qNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch8qGain[,2], Drug_Ch8qLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts8q<- rbind(CRISPR_DepMap_Ts8q, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 
#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts8q<- CRISPR_DepMap_Ts8q[order(CRISPR_DepMap_Ts8q$Gain.Neutral.P),]
CRISPR_DepMap_Ts8q[1:20,]


CRISPR_DepMap_Ts8q$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts8q$GeneID)){
  if (CRISPR_DepMap_Ts8q$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts8q$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts8q$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts8q, SigDiff=="TRUE")
CRISPR_DepMap_Ts8q$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_Ts8q$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_Ts8q$Gene<- sub("[.]", "-", CRISPR_DepMap_Ts8q$Gene) # replace single period with -. 



ggplot(CRISPR_DepMap_Ts8q, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, 
                               color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts8q, SigDiff=="TRUE"))+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Gene Beta Score (DepMap): \nGain 8q minus Neutral/Loss 8q")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.1")
# 5x4 
# plot.DepMap.Gain8q.LossNeutral.Pvalue.Diff.pdf 
subset(CRISPR_DepMap_Ts8q, SigDiff=="TRUE")$Gene



## Now merge my CRISPR screen data with Ts8q DepMap data ###
DepMap.Sheltzer.CRISPR_Ts8q<- merge(Merge.Ts8.QN2, CRISPR_DepMap_Ts8q, by.x="Gene", by.y="Gene") 

ggplot(DepMap.Sheltzer.CRISPR_Ts8q, aes(x=ANvEU.Beta, y=Gain.min.LossNeutral.Diff, 
                                        color= ANvEU.Beta < -0.2 & Gain.min.LossNeutral.Diff < -0.05))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_Ts8q, ANvEU.Beta < -0.2 & Gain.min.LossNeutral.Diff < -0.05))+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(color = "white", bins = 10)+
  ylab("Difference in CRISPR Score (DepMap): \nGain 8q minus Neutral/Loss 8q")+
  xlab("Difference in Beta Score (Sheltzer): \nGain 8 minus Neutral/Loss 8")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.CRISPR_Ts8q.pdf 

subset(DepMap.Sheltzer.CRISPR_Ts8q, ANvEU.Beta < -0.2 & Gain.min.LossNeutral.Diff < -0.1)#$Gene
#  "ANAPC15" "PSMB6"   "RFX7"    "SH3RF3"  "UBC"     "UBE2N"  
# top RFX7



HighlightTheseGenes<- Chrm8Genes

ggplot(DepMap.Sheltzer.CRISPR_Ts8q, aes(y=ANvEU.Beta, x=Gain.min.LossNeutral.Diff, 
                                        color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_Ts8q, Gene %in% HighlightTheseGenes))+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Difference in Gene Beta Score (DepMap): \nGain 8q minus Neutral/Loss 8q")+
  ylab("Difference in Gene Beta Score (Sheltzer): \nGain 8 minus Neutral/Loss 8")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.CRISPR_Ts8q.pdf 
subset(DepMap.Sheltzer.CRISPR_Ts8q, Gene %in% HighlightTheseGenes & Gain.min.LossNeutral.Diff< -0.045)


setwd(ResultsFile)
# write.csv(CRISPR_DepMap_Ts8q, "DepMap.Gain8q.NeutralLoss8q_Difference_pvalue.csv")
# CRISPR_DepMap_Ts8q<- read.csv("DepMap.Gain8q.NeutralLoss8q_Difference_pvalue.csv")

setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.CRISPR_Ts8q, "Chrm8q_Specific_Toxicity_Sheltzer.DepMap.csv")


         ###### Ts8 paired & DepMap 8p-gain data #####

CN8p<- Arm_CN[, c("X", "X8p")]
colnames(CN8p)<- c("Cell_Line", "8p")
Cell_Neutral8p<- subset(CN8p, `8p`=="0")$Cell_Line
Cell_Gain8p<- subset(CN8p, `8p`=="1")$Cell_Line
Cell_Loss8p<- subset(CN8p, `8p`=="-1")$Cell_Line
Cell_Loss.Neutral8p<- subset(CN8p, `8p`!="1")$Cell_Line


CRISPR_DepMap_Ts8p<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch8pGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain8p) #get drug result for all cells with gain for chrm 8
  Drug_Ch8pNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral8p) #get drug result for all cells neutral for chrm 8
  Drug_Ch8pLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral8p) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch8pGain[,2])) >5 & 
      length(na.omit(Drug_Ch8pNeutral[,2])) >5 &
      length(na.omit(Drug_Ch8pLoss.Neutral[,2])) >5 )  {
    
    test1<-t.test(Drug_Ch8pGain[,2], Drug_Ch8pNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch8pGain[,2], Drug_Ch8pLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts8p<- rbind(CRISPR_DepMap_Ts8p, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 
#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts8p<- CRISPR_DepMap_Ts8p[order(CRISPR_DepMap_Ts8p$Gain.Neutral.P),]
CRISPR_DepMap_Ts8p[1:20,]


CRISPR_DepMap_Ts8p$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts8p$GeneID)){
  if (CRISPR_DepMap_Ts8p$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts8p$Gain.min.LossNeutral.Diff[i]< -0.2){
    CRISPR_DepMap_Ts8p$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts8p, SigDiff=="TRUE")
CRISPR_DepMap_Ts8p$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_Ts8p$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_Ts8p$Gene<- sub("[.]", "-", CRISPR_DepMap_Ts8p$Gene) # replace single period with -. 




ggplot(CRISPR_DepMap_Ts8p, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts8p, SigDiff=="TRUE"))+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Gene Beta Score (DepMap): \nGain 8p minus Neutral/Loss 8p")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.2")
# 5x4 
# plot.CRISPR.Gain8p.LossNeutral.Pvalue.Diff.pdf 

colnames(CRISPR_DepMap_Ts8p)
colnames(CRISPR_DepMap_Ts8p)<- c("GeneID",  "Gain8p.Neutral.P", "Gain8p.min.Neutral.Diff", "Gain8p.LossNeutral.P", "Gain8p.min.LossNeutral.Diff",
                                 "SigDiff8p" , "Gene")
colnames(CRISPR_DepMap_Ts8q)
colnames(CRISPR_DepMap_Ts8q)<- c("GeneID",  "Gain8q.Neutral.P", "Gain8q.min.Neutral.Diff", "Gain8q.LossNeutral.P", "Gain8q.min.LossNeutral.Diff",
                                 "SigDiff8q" , "Gene")

## Now merge my CRISPR screen data with Ts8p DepMap data ###
DepMap.Sheltzer.CRISPR_Ts8pq<- merge(Merge.Ts8.QN2, CRISPR_DepMap_Ts8p, by.x="Gene", by.y="Gene") 
DepMap.Sheltzer.CRISPR_Ts8pq<- merge(DepMap.Sheltzer.CRISPR_Ts8pq, CRISPR_DepMap_Ts8q, by.x="Gene", by.y="Gene") 

### Plot Difference Ts8 with CRISPR values highlighted

DepMap8pHits <- subset(DepMap.Sheltzer.CRISPR_Ts8pq, Gain8p.min.LossNeutral.Diff< -0.1 & Gain8p.LossNeutral.P< 0.05)$Gene
DepMap8qHits <- subset(CRISPR_DepMap_Ts8q, Gain8q.min.LossNeutral.Diff< -0.1 & Gain8q.LossNeutral.P <0.05)$Gene

ggplot(Merge.Ts8.QN2, aes(x= ANvEU.Beta, 
                          y= -log2(ANvEU.P)))+
  geom_point()+
  geom_point(data = subset(Merge.Ts8.QN2, Gene %in% DepMap8pHits), color = "orchid1")+ 
  geom_point(data = subset(Merge.Ts8.QN2, Gene %in% DepMap8qHits), color = "orange1")+ 
  #geom_point(data = subset(Merge.Ts8.QN2, Gene %in% DepMap8pHits & Gene %in% DepMap8qHits), color = "red")+ 
  #geom_point(data = subset(Merge.Ts8.QN2, Gene %in% c("STK11")), color = "Blue")+ 
  ylab("-log2(P-value)")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Beta score\nTs8 - euploid")+
  theme_classic()
# 5x4
# plot.BetaDelta.Ts8vEU.colorDepMap_sigdif.pdf
# SPI1 8p
# ANAPC15,  UBE2N,  UBC  8q

setwd(ResultsFile)
#write.csv(DepMap.Sheltzer.CRISPR_Ts8pq, "Chrm8pq_Specific_Toxicity_Sheltzer.DepMap.csv")
# DepMap.Sheltzer.CRISPR_Ts8pq<- read.csv("Chrm8pq_Specific_Toxicity_Sheltzer.DepMap.csv")



     ##### Ts5p only (4 pairs) (QN) #####
## HCT116 Tet 5 (2x), DLD1 Ts5.15, SNU1 Ps5
# C13 = Control/WT SNU1: Gain 20, monosomy X 
# C24 = Ts10.1.8.14 Ps5

MegaA.Merge.ALL #from "Quantile normalization ALL" section
colnames(MegaA.Merge.ALL)

MLE.Merge.Ts5<- MegaA.Merge.ALL[, c(1,4,5,6,11,13,23,25)]
colnames(MLE.Merge.Ts5)<- c("Gene", "MegaA.H.Cont2.beta", "MegaA.H.Tet5p.Tri5q.beta", "MegaA.H.Tet5p.beta", 
                            "MegaA.D.Cont2.beta", "MegaA.D.5.15.beta", 
                            "MegaA.Snu1.Cont1.beta", "MegaA.Snu1.Ts10.1.8.14.Ps5.beta") #c24



## Quantile normalize by batch and cell line: 
Merge.Ts5.QN # quantile normalized across all samples from MegaA and WGL
Merge.Ts5.QN<- quantile_normalisation(MLE.Merge.Ts5[,c(2:length(MLE.Merge.Ts5))])
Merge.Ts5.QN<- as.data.frame(Merge.Ts5.QN)
Merge.Ts5.QN$Gene<- (MLE.Merge.Ts5$Gene)
head(Merge.Ts5.QN)

### Plot pos controls of Ts5 only
PlotHeatmap.Cont(DropoutData= list("HCT116 Ctrl MegaA"= Merge.Ts5.QN[,c(8,1,1)],
                                   "HCT116 Tet5p.Tri5q MegaA"= Merge.Ts5.QN[,c(8,1,2)], 
                                   "HCT116 Tet5p MegaA"= Merge.Ts5.QN[,c(8,1,3)], 
                                   "DLD1 Ctrl MegaA"= Merge.Ts5.QN[,c(8,1,4)], 
                                   "DLD1 Ts5.15 MegaA"= Merge.Ts5.QN[,c(8,1,5)], 
                                   "Snu1 Ctrl MegaA"= Merge.Ts5.QN[,c(8,1,6)], 
                                   "Snu1 Snu1.Ts10.1.8.14 Ps5 MegaA"= Merge.Ts5.QN[,c(8,1,7)]
),
Title="Control guide Beta-scores")
# 7x4
# plot.Ts5Only.QN.MegaA.WGL.pdf


###
Merge.Ts5.QN2<- Merge.Ts5.QN
Merge.Ts5.QN2$ANvEU.P<- NA # t-test. two sided, not paired
Merge.Ts5.QN2$ANvEU.Beta<- NA  #difference btween average an and average euploid beta score (aneuploid beta minus euploid beta)
Merge.Ts5.QN2$Mean.An<- NA   
Merge.Ts5.QN2$Mean.Eu<- NA 

for (i in 1:length(Merge.Ts5.QN2$Gene)){
  test<-t.test(as.numeric(Merge.Ts5.QN2[i,c(1,1,4,6)]), #controls
               as.numeric(Merge.Ts5.QN2[i,c(2,3,5,7)]), # Ts5
               paired=TRUE) 
  Merge.Ts5.QN2$ANvEU.P[i]<- test$p.value
  Contbeta<- mean(as.numeric(Merge.Ts5.QN2[i,c(1,1,4,6)]))
  Ts5beta<- mean(as.numeric(Merge.Ts5.QN2[i,c(2,3,5,7)]))
  Merge.Ts5.QN2$ANvEU.Beta[i]<- Ts5beta - Contbeta
  Merge.Ts5.QN2$Mean.An[i]<-    Ts5beta
  Merge.Ts5.QN2$Mean.Eu[i]<-    Contbeta
}
Merge.Ts5.QN2<-Merge.Ts5.QN2[order(Merge.Ts5.QN2$ANvEU.Beta),]


Mega.WGL.Ts5.Only.TopHits<- subset(Merge.Ts5.QN2, ANvEU.P<0.1 & ANvEU.Beta< -0.2)$Gene
ggplot(Merge.Ts5.QN2, aes(x= ANvEU.Beta, 
                          y= -log2(ANvEU.P), color= Gene %in% Mega.WGL.Ts5.Only.TopHits))+
  geom_point()+
  ylab("-log2(P-value)")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Beta score\nTs5 - euploid")+
  scale_color_manual(values = c("Black", "red"), name = "P<0.05 \nDiff < -.2")+
  theme_classic()
# 5x4
# plot.BetaDelta.Ts5vEU.pdf

subset(Merge.Ts5.QN2, ANvEU.P< 0.05 & ANvEU.Beta< -0.2)$Gene
# [1] "CDK4"    "TGIF2LX" "TCF7L2"  "UBE2H"   "KCTD2"   "MYLK2"   "ELMSAN1" "TRIM8"   "PCNA"    "ILK"     "KALRN"   "RPS6KB1" "TRIM61" 
# [14] "SPI1"    "ZNF99"   "ENDOU"   "MSX1"    "PRKDC"   "PRKCZ"   "ZNF784"  "CPZ"     "GFI1"    "MAPK12" 


TopHits<- c("PPP1R12A") 
TopHits<- Chrm5Genes

ggplot(Merge.Ts5.QN2, aes(x= ANvEU.Beta, 
                          y= -log2(ANvEU.P), color= Gene %in% TopHits))+
  geom_point()+
  geom_point(data = subset(Merge.Ts5.QN2, Gene%in% TopHits))+
  ylab("-log2(P-value)")+
  xlab("Difference in Beta score\nTs5 - euploid")+
  scale_color_manual(values = c("Black", "Red"), name = "Top Hits")+
  theme_classic()
# 5x4
# plot.MegaA.BetaDelta.pValue_QN_TopTs5qHits.pdf 

subset(Merge.Ts5.QN2, Gene %in% TopHits & ANvEU.P<0.1)


# Image Ts5 MegaA  gene dropout

test<-subset(Merge.Ts5.QN2, Gene=="TRIO")
df<- data.frame(Control=as.numeric(test[,c(1,1,4,6)]), Ts5=as.numeric(test[,c(2,3,5,7)]), Name=colnames(test[,c(2,3,5,7)]))
ggpaired(df, cond1 = "Control", cond2 = "Ts5", 
         title="TRIO", fill="condition")+ 
  ylab("Quantile Norm. Beta")+
  scale_fill_manual(values=c("darkolivegreen1", "skyblue"))
# 3x5
# plot.pairedLine.Ts5_CERS2.pdf


### Look at genes located ON chromosome 5: 
Chrm5Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 5)$`HGNC symbol`)
subset(Merge.Ts5.QN2, Gene %in% Chrm5Genes)
subset(Merge.Ts5.QN2, Gene %in% Chrm5Genes & ANvEU.P < 0.1)
# of genes within MegaA and on chromosome 5


         ###### Ts5 correlate with DepMap 5p-gain data #####

CN5p<- Arm_CN[, c("X", "X5p")]
colnames(CN5p)<- c("Cell_Line", "5p")
Cell_Neutral5p<- subset(CN5p, `5p`=="0")$Cell_Line
Cell_Gain5p<- subset(CN5p, `5p`=="1")$Cell_Line
Cell_Loss5p<- subset(CN5p, `5p`=="-1")$Cell_Line
Cell_Loss.Neutral5p<- subset(CN5p, `5p`!="1")$Cell_Line


CRISPR_DepMap_Ts5p<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch5pGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain5p) #get drug result for all cells with gain for chrm 8
  Drug_Ch5pNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral5p) #get drug result for all cells neutral for chrm 8
  Drug_Ch5pLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral5p) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch5pGain[,2])) >5 & 
      length(na.omit(Drug_Ch5pNeutral[,2])) >5 &
      length(na.omit(Drug_Ch5pLoss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch5pGain[,2], Drug_Ch5pNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch5pGain[,2], Drug_Ch5pLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts5p<- rbind(CRISPR_DepMap_Ts5p, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 
#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts5p<- CRISPR_DepMap_Ts5p[order(CRISPR_DepMap_Ts5p$Gain.Neutral.P),]
CRISPR_DepMap_Ts5p[1:20,]


CRISPR_DepMap_Ts5p$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts5p$GeneID)){
  if (CRISPR_DepMap_Ts5p$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts5p$Gain.min.LossNeutral.Diff[i]< -0.2){
    CRISPR_DepMap_Ts5p$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts5p, SigDiff=="TRUE")
CRISPR_DepMap_Ts5p$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_Ts5p$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_Ts5p$Gene<- sub("[.]", "-", CRISPR_DepMap_Ts5p$Gene) # replace single period with -. 




ggplot(CRISPR_DepMap_Ts5p, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts5p, SigDiff=="TRUE"))+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in CRISPR Score (DepMap): \nGain 5p minus Neutral/Loss 5p")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.1")
# 5x4 
# plot.CRISPR.Gain5p.LossNeutral.Pvalue.Diff.pdf 


colnames(CRISPR_DepMap_Ts5p)<- c("GeneID",  "Gain5p.Neutral.P","Gain5p.min.Neutral.Diff","Gain5p.LossNeutral.P", "Gain5p.min.LossNeutral.Diff",
                                 "SigDiff5p","Gene")

## Now merge my CRISPR screen data with Ts5p DepMap data ###
DepMap.Sheltzer.CRISPR_Ts5p<- merge(Merge.Ts5.QN2, CRISPR_DepMap_Ts5p, by.x="Gene", by.y="Gene") 
anti_join(Merge.Ts5.QN2, CRISPR_DepMap_Ts5p)$Gene # 175 did not merge (some are negative controls)

AntiJoinMegaATs5<- anti_join(Merge.Ts5.QN2, CRISPR_DepMap_Ts5p) # Find genes in my dataset that didn't merge
AntiJoinMegaATs5<- AntiJoinMegaATs5[-grep("neg ", AntiJoinMegaATs5$Gene), ] # Remove negative controls. genes = 173


#Change old gene names (BROAD) to newer gene names (used in DepMap) using limma library: 

AntiJoinMegaATs5$Gene2<- AntiJoinMegaATs5$Gene
for (i in 1:length(AntiJoinMegaATs5$Gene)){
  if (length(alias2Symbol(AntiJoinMegaATs5$Gene[i], species = "Hs", expand.symbols = FALSE) )>0){ # If gene can be converted. 
    AntiJoinMegaATs5$Gene2[i]<-alias2Symbol(AntiJoinMegaATs5$Gene[i], species = "Hs", expand.symbols = FALSE)[1] #label to new gene name
  }
} #Is this the best way to run this code? no. but it works. It takes a while to run though... 
# 2512 genes have new names, of the 2553 genes in anti-join list. 

subset(AntiJoinMegaATs5, Gene2 == "IRF4")  # Test! MUM1 should be converted to IRF4
length(AntiJoinMegaATs5$Gene2) # 173

#Join DepMap with my data, subset genes with updated gene names: 
AntiJoinMegaATs5_CRISPR<- merge(AntiJoinMegaATs5, CRISPR_DepMap_Ts5p, by.x="Gene2", by.y="Gene")   # 49 
DepMap.Sheltzer.CRISPR_Ts5p$Gene2<- DepMap.Sheltzer.CRISPR_Ts5p$Gene
DepMap.Sheltzer.CRISPR_Ts5p<- rbind(DepMap.Sheltzer.CRISPR_Ts5p, AntiJoinMegaATs5_CRISPR) #length= 3281
# "Gene" is now the updated gene name! Gene.y is old name




ggplot(DepMap.Sheltzer.CRISPR_Ts5p, aes(x=ANvEU.Beta, y=Gain.min.LossNeutral.Diff, 
                                        color= ANvEU.Beta < -0.2 & Gain.min.LossNeutral.Diff < -0.05))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_Ts5p, ANvEU.Beta < -0.2 & Gain.min.LossNeutral.Diff < -0.05))+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(color = "white", bins = 10)+
  ylab("Difference in CRISPR Score (DepMap): \nGain 5p minus Neutral/Loss 5p")+
  xlab("Difference in Beta Score (Sheltzer): \nGain 5 minus Neutral/Loss 5")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.CRISPR_Ts5p.pdf 
subset(DepMap.Sheltzer.CRISPR_Ts5p, ANvEU.Beta < -0.2 & Gain.min.LossNeutral.Diff < -0.05)#$Gene
# [1] "BIRC6"  "EGFR"   "FOSL1"  "ILK"    "PCNA"   "PKMYT1" "UBE2N"  "ZNF207"



HighlightTheseGenes<-  c("EGFR")
ggplot(DepMap.Sheltzer.CRISPR_Ts5p, aes(x=ANvEU.Beta, y=Gain.min.LossNeutral.Diff, 
                                        color= Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(DepMap.Sheltzer.CRISPR_Ts5p, Gene %in% HighlightTheseGenes))+
  theme_classic()+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(color = "white", bins = 10)+
  ylab("Difference in Gene Beta Score (DepMap): \nGain 5p minus Neutral/Loss 5p")+
  xlab("Difference in Gene Beta Score (Sheltzer): \nGain 5p minus Neutral/Loss 5p")+
  scale_color_manual(values= c("Black", "Red"), name="")
# 5x4 
# plot.DepMap.Sheltzer.CRISPR_Ts5p.pdf 
subset(DepMap.Sheltzer.CRISPR_Ts5p, Gene %in% HighlightTheseGenes & Gain.min.LossNeutral.Diff < -0.025)


setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.CRISPR_Ts5p, "Chrm5p_Specific_Toxicity_Sheltzer.DepMap.csv")



### Plot Difference 5p with CRISPR values highlighted

DepMap5pHits <- subset(DepMap.Sheltzer.CRISPR_Ts5p, 
                       Gain5p.min.LossNeutral.Diff< -0.1 & Gain5p.LossNeutral.P<0.05)$Gene

ggplot(Merge.Ts5.QN2, aes(x= ANvEU.Beta, 
                          y= -log2(ANvEU.P)))+ 
  geom_point()+
  geom_point(data = subset(Merge.Ts5.QN2, Gene %in% DepMap5pHits), color = "orchid1")+ 
  #geom_point(data = subset(Merge.Ts5.QN2, Gene %in% c("UBE2H")), color = "blue")+ 
  ylab("-log2(P-value)")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Beta score\nTs5p - euploid")+
  theme_classic()
# 4x4
# plot.BetaDelta.Ts5qvEU.colorDepMap_SigDif.pdf


subset(Merge.Ts5.QN2, Gene %in% DepMap5pHits & ANvEU.Beta < 0.2)



     ##### Ts2 only (4 pairs) (QN) #####
## A2780 Ts2, DLD1 Ts2.18, S.Ts10.17.5_Ps2, DLD1 Ts2.18_try2

# C13 = Control/WT SNU1: Gain 20, monosomy X
# C12 =	Ts10.17.5 Ps2 
# C23 = Ts10 Ps6 (D=10 lost it's aneuploidy, do not include)
# C24 = Ts10.1.8.14 Ps5 (D=10 actually labeled C23) 
# C111 = Ts10.6q Ps8.12 
# C114 = Ts9.19 

MegaA.Merge.ALL #from "Quantile normalization ALL" section
colnames(MegaA.Merge.ALL) 

Merge.Ts2.QN<- MegaA.Merge.ALL[, c(1,7,10,14,15,17,18,23,24)]
colnames(Merge.Ts2.QN)<- c("Gene", "DLD1 WT", "DLD1 Ts2.18", "DLD1 WT", "DLD1 Ts2.18", 
                           "A2780 WT", "A2780 Ts2", "SNU1 WT", "SNU1 Ts10.17.5 Ps2")



## Quantile normalize by batch and cell line: 
Merge.Ts2.QN # quantile normalized across all samples from MegaA and WGL
Merge.Ts2.QN<- quantile_normalisation(Merge.Ts2.QN[,c(2:length(Merge.Ts2.QN))])
Merge.Ts2.QN<- as.data.frame(Merge.Ts2.QN)
Merge.Ts2.QN$Gene<- (MegaA.Merge.ALL$Gene)
head(Merge.Ts2.QN)

### Plot pos controls of Ts2 only
PlotHeatmap.Cont(DropoutData= list("DLD1 WT MegaA"= Merge.Ts2.QN[,c(9,1,1)],
                                   "DLD1 Ts2.18 MegaA"= Merge.Ts2.QN[,c(9,1,2)], 
                                   "DLD1 WT MegaA_try2"= Merge.Ts2.QN[,c(9,1,3)],
                                   "DLD1 Ts2.18 MegaA_try2"= Merge.Ts2.QN[,c(9,1,4)], 
                                   "A2780 WT MegaA"= Merge.Ts2.QN[,c(9,1,5)], 
                                   "A2780 Ts2 MegaA"= Merge.Ts2.QN[,c(9,1,6)], 
                                   "SNU1 WT MegaA"= Merge.Ts2.QN[,c(9,1,7)], 
                                   "Snu1 Ts10.17.5 Ps2 MegaA"= Merge.Ts2.QN[,c(9,1,8)]
),
Title="Control guide Beta-scores")
# 7x4
# plot.Ts2.QN.MegaA.WGL.pdf


###
Merge.Ts2.QN2<- Merge.Ts2.QN
Merge.Ts2.QN2$ANvEU.P <-   NA # t-test. two sided, not paired
Merge.Ts2.QN2$ANvEU.Beta<- NA  #difference btween average an and average euploid beta score (aneuploid beta minus euploid beta)
Merge.Ts2.QN2$Mean.An<- NA   
Merge.Ts2.QN2$Mean.Eu<- NA 

for (i in 1:length(Merge.Ts2.QN2$Gene)){
  test<-t.test(as.numeric(Merge.Ts2.QN2[i,c(1,3,5,7)]), #controls
               as.numeric(Merge.Ts2.QN2[i,c(2,4,6,8)]), # Ts5
               paired=TRUE) 
  Merge.Ts2.QN2$ANvEU.P[i]<- test$p.value
  Contbeta<- mean(as.numeric(Merge.Ts2.QN2[i,c(1,3,5,7)]))
  Ts2beta<- mean(as.numeric(Merge.Ts2.QN2[i,c(2,4,6,8)]))
  Merge.Ts2.QN2$ANvEU.Beta[i]<- Ts2beta - Contbeta
  Merge.Ts2.QN2$Mean.An[i]<-    Ts2beta
  Merge.Ts2.QN2$Mean.Eu[i]<-    Contbeta
}
Merge.Ts2.QN2<-Merge.Ts2.QN2[order(Merge.Ts2.QN2$ANvEU.Beta),]


Mega.WGL.Ts2.TopHits<- subset(Merge.Ts2.QN2, ANvEU.P<0.05 & ANvEU.Beta< -0.2)$Gene
ggplot(Merge.Ts2.QN2, aes(x= ANvEU.Beta, 
                          y= -log2(ANvEU.P), color= Gene %in% Mega.WGL.Ts2.TopHits))+
  geom_point()+
  ylab("-log2(P-value)")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Beta score\nTs2 - euploid")+
  scale_color_manual(values = c("Black", "red"), name = "P<0.05 \nDiff < -.2")+
  theme_classic()
# 5x4
# plot.BetaDelta.Ts2vEU.pdf

subset(Merge.Ts2.QN2, ANvEU.P<0.05 & ANvEU.Beta< -0.2 )$Gene

# [1] "BTBD2"     "RPS6KL1"   "SMAD6"     "AURKC"     "ZNF135"    "HOXD10"    "GSK3B"     "RAX2"      "USP38"     "UQCRC1"    "NEUROD2"  
#[12] "CBLL1"     "NR1D2"     "MAK"       "ZNF528"    "INSM1"     "ZNF786"    "RNF207"    "IRAK4"     "PLAGL1"    "POU2F3"    "TDRD3"    
#[23] "ZNF711"    "ZFP92"     "GLI4"      "PRKCQ"     "ZNF841"    "EPHA3"     "PBX4"      "ZNF789"    "ZNF146"    "MARCHF11"  "KCNG1"    
#[34] "KLHL36"    "TCF23"     "HNF4A"     "TMPRSS11D" "MYF5"      "ZNF324B"   "LCK"       "REPIN1"    "ESR2"      "NHLRC1"    "RNF130"  




TopHits<- c("PKN3", "HOXC9") 
TopHits<- Chrm2Genes
ggplot(Merge.Ts2.QN2, aes(x= ANvEU.Beta, 
                          y= -log2(ANvEU.P), color= Gene %in% TopHits))+
  geom_point()+
  geom_point(data = subset(Merge.Ts2.QN2, Gene%in% TopHits))+
  ylab("-log2(P-value)")+
  xlab("Difference in Beta score\nTrisomy 2 - Disomy 2")+
  scale_color_manual(values = c("Black", "Red"), name = "Top Hits")+
  theme_classic()
# 5x4
# plot.MegaA.BetaDelta.pValue_QN_TopTs2Hits.pdf 

subset(Merge.Ts2.QN2, Gene %in% TopHits & ANvEU.P<0.05 & ANvEU.Beta < -0.1)


test<-subset(Merge.Ts2.QN2, Gene=="TCF23")
df<- data.frame(Control=as.numeric(test[,c(1,3,5,7)]), Aneuploid=as.numeric(test[,c(2,4,6,8)]), Name=colnames(test[,c(2,4,6,8)]))
ggpaired(df, cond1 = "Control", cond2 = "Aneuploid", 
         title="TCF23", fill="condition")+ 
  ylab("Quantile Norm. Beta")+
  scale_fill_manual(values=c("darkolivegreen1", "skyblue"))
# 3x5
# plot.paired.Ts2_TCF23.pdf



### Look at genes located ON chromosome 2: 
Chrm2Genes<- unique(subset(Gene_Info, `Chromosome/scaffold name` == 2)$`HGNC symbol`)
subset(Merge.Ts2.QN2, Gene %in% Chrm2Genes)
subset(Merge.Ts2.QN2, Gene %in% Chrm2Genes & ANvEU.P < 0.05 & ANvEU.Beta < -0.2)
#  HOXD10, TCF23



         ###### Ts2p correlate with DepMap 2p-gain data #####

CN2p<- Arm_CN[, c("X", "X2p")]
colnames(CN2p)<- c("Cell_Line", "2p")
Cell_Neutral2p<- subset(CN2p, `2p`=="0")$Cell_Line
Cell_Gain2p<- subset(CN2p, `2p`=="1")$Cell_Line
Cell_Loss2p<- subset(CN2p, `2p`=="-1")$Cell_Line
Cell_Loss.Neutral2p<- subset(CN2p, `2p`!="1")$Cell_Line


CRISPR_DepMap_Ts2p<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch2pGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain2p) #get drug result for all cells with gain for chrm 8
  Drug_Ch2pNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral2p) #get drug result for all cells neutral for chrm 8
  Drug_Ch2pLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral2p) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch2pGain[,2])) >5 & 
      length(na.omit(Drug_Ch2pNeutral[,2])) >5 &
      length(na.omit(Drug_Ch2pLoss.Neutral[,2])) >5)  {
    
    test1<-t.test(Drug_Ch2pGain[,2], Drug_Ch2pNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch2pGain[,2], Drug_Ch2pLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts2p<- rbind(CRISPR_DepMap_Ts2p, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 
#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts2p<- CRISPR_DepMap_Ts2p[order(CRISPR_DepMap_Ts2p$Gain.Neutral.P),]
CRISPR_DepMap_Ts2p[1:20,]


CRISPR_DepMap_Ts2p$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts2p$GeneID)){
  if (CRISPR_DepMap_Ts2p$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts2p$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts2p$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts2p, SigDiff=="TRUE")
CRISPR_DepMap_Ts2p$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_Ts2p$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_Ts2p$Gene<- sub("[.]", "-", CRISPR_DepMap_Ts2p$Gene) # replace single period with -. 



ggplot(CRISPR_DepMap_Ts2p, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts2p, SigDiff=="TRUE"))+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Gene Beta Score (DepMap): \nGain 2p minus Neutral/Loss 2p")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.1")
# 5x4 
# plot.CRISPR.Gain2p.LossNeutral.Pvalue.Diff.pdf 



setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.CRISPR_Ts2p, "Chrm2p_Specific_Toxicity_Sheltzer.DepMap.csv")





         ###### Ts2q correlate with DepMap 2q-gain data #####

CN2q<- Arm_CN[, c("X", "X2q")]
colnames(CN2q)<- c("Cell_Line", "2q")
Cell_Neutral2q<- subset(CN2q, `2q`=="0")$Cell_Line
Cell_Gain2q<- subset(CN2q, `2q`=="1")$Cell_Line
Cell_Loss2q<- subset(CN2q, `2q`=="-1")$Cell_Line
Cell_Loss.Neutral2q<- subset(CN2q, `2q`!="1")$Cell_Line


CRISPR_DepMap_Ts2q<-data.frame(GeneID= character(), 
                               Gain.LossNeutral.P = numeric(), 
                               Gain.min.LossNeutral.Diff= numeric()
)

for (i in 2:length(CRISPR_Gene)){
  Drug_Ch2qGain<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Gain2q) #get drug result for all cells with gain for chrm 8
  Drug_Ch2qNeutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Neutral2q) #get drug result for all cells neutral for chrm 8
  Drug_Ch2qLoss.Neutral<- subset(CRISPR_Gene[,c(1,i)], ModelID %in% Cell_Loss.Neutral2q) #get drug result for all cells neutral or loss for chrm 8
  
  if (length(na.omit(Drug_Ch2qGain[,2])) >5 & 
      length(na.omit(Drug_Ch2qNeutral[,2])) >5 &
      length(na.omit(Drug_Ch2qLoss.Neutral[,2])) >5 )  {
    
    test1<-t.test(Drug_Ch2qGain[,2], Drug_Ch2qNeutral[,2], "two.sided")
    test2<-t.test(Drug_Ch2qGain[,2], Drug_Ch2qLoss.Neutral[,2], "two.sided")
    
    CRISPR_DepMap_Ts2q<- rbind(CRISPR_DepMap_Ts2q, data.frame(GeneID= colnames(CRISPR_Gene)[i],
                                                              Gain.Neutral.P= test1$p.value, #p-value between drug scores with gain of 13 and drug score in cells neutral for 13
                                                              Gain.min.Neutral.Diff= test1$estimate[1]- test1$estimate[2], # Mean Drug score Gain13  minus  mean drugscore Neutral 13
                                                              Gain.LossNeutral.P = test2$p.value, 
                                                              Gain.min.LossNeutral.Diff= test2$estimate[1]- test2$estimate[2]) ) # Mean Drug score Gain13  minus  mean drugscore Neutral or loss 13
  }
} 
#length of Genes analyzed: 17386 Genes
CRISPR_DepMap_Ts2q<- CRISPR_DepMap_Ts2q[order(CRISPR_DepMap_Ts2q$Gain.Neutral.P),]
CRISPR_DepMap_Ts2q[1:20,]


CRISPR_DepMap_Ts2q$SigDiff<- "FALSE"
for (i in 1:length(CRISPR_DepMap_Ts2q$GeneID)){
  if (CRISPR_DepMap_Ts2q$Gain.LossNeutral.P[i]<0.05 && CRISPR_DepMap_Ts2q$Gain.min.LossNeutral.Diff[i]< -0.1){
    CRISPR_DepMap_Ts2q$SigDiff[i]<- "TRUE"
  }
}
subset(CRISPR_DepMap_Ts2q, SigDiff=="TRUE")
CRISPR_DepMap_Ts2q$Gene<- sub("\\.\\..*", "", CRISPR_DepMap_Ts2q$GeneID) #Two periods \\.\\. followed by any character .* replace with nothing
CRISPR_DepMap_Ts2q$Gene<- sub("[.]", "-", CRISPR_DepMap_Ts2q$Gene) # replace single period with -. 
# ITGAV, KIF2C, FERMT2


ggplot(CRISPR_DepMap_Ts2q, aes(y=-log(Gain.LossNeutral.P,2), x=Gain.min.LossNeutral.Diff, color= SigDiff))+
  geom_point()+
  geom_point(data = subset(CRISPR_DepMap_Ts2p, SigDiff=="TRUE"))+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Gene Beta Score (DepMap): \nGain 2q minus Neutral/Loss 2q")+
  ylab("-log2(P-value)")+
  scale_color_manual(values= c("Black", "Red"), name="P<0.05, Diff < -0.1")
# 5x4 
# plot.CRISPR.Gain2p.LossNeutral.Pvalue.Diff.pdf 



setwd(ResultsFile)
# write.csv(DepMap.Sheltzer.CRISPR_Ts2q, "Chrm2q_Specific_Toxicity_Sheltzer.DepMap.csv")




### Plot Difference 2p with CRISPR values highlighted
colnames(CRISPR_DepMap_Ts2p) <- c("GeneID",                    "Gain2p.Neutral.P",            "Gain2p.min.Neutral.Diff",     "Gain2p.LossNeutral.P",        "Gain2p.min.LossNeutral.Diff",
                                  "SigDiff2p",                   "Gene"  )
colnames(CRISPR_DepMap_Ts2q)  <- c("GeneID",                    "Gain2q.Neutral.P",            "Gain2q.min.Neutral.Diff",     "Gain2q.LossNeutral.P",        "Gain2q.min.LossNeutral.Diff",
                                   "SigDiff2q",                   "Gene"  )
Merge.Ts2.QN2.DepMappq<- merge(Merge.Ts2.QN2, CRISPR_DepMap_Ts2p)
Merge.Ts2.QN2.DepMappq<- merge(Merge.Ts2.QN2.DepMappq, CRISPR_DepMap_Ts2q)
write.csv(Merge.Ts2.QN2.DepMappq, "Merge.Ts2.QN2.DepMappq.csv")


DepMap2pHits<- subset(CRISPR_DepMap_Ts2p,  Gain2p.min.LossNeutral.Diff< -0.1 & Gain2p.LossNeutral.P< 0.05)$Gene
DepMap2qHits <- subset(CRISPR_DepMap_Ts2q, Gain2q.min.LossNeutral.Diff< -0.1 & Gain2q.LossNeutral.P< 0.05)$Gene

ggplot(Merge.Ts2.QN2, aes(x= ANvEU.Beta, 
                          y= -log2(ANvEU.P)))+
  geom_point()+
  geom_point(data = subset(Merge.Ts2.QN2, Gene %in% DepMap2pHits), color = "orchid1")+ 
  geom_point(data = subset(Merge.Ts2.QN2, Gene %in% DepMap2qHits), color = "orange1")+ 
  geom_point(data = subset(Merge.Ts2.QN2, Gene %in% DepMap2qHits & Gene %in% DepMap2pHits), color = "red")+ 
  #geom_point(data = subset(Merge.Ts2.QN2, Gene == "NRBP1"), color = "blue")+ 
  ylab("-log2(P-value)")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Difference in Beta score\nTs2 - euploid")+
  theme_classic()
# 5x4
# plot.BetaDelta.Ts2vEU.colorDepMap_SigDif.pdf
# 2p

#PSMB3



### Top hits per cell type: DLD1, HCT116, A2780 and SNU1 ####
### below function will quantile normalize your data
### !!! BUT!!! you need merged control and aneuploid data, and ONLY give function the BETA SCORE COLLUMNS!!!! 
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
  return(as.data.frame(df_final))
}

     ##### HCT116 (3 pairs) #####
## Quantile normalize the dropouts per cell type (HCT116)
## essentially rank order them, then take the average dropout per rank, give that number per gene. 
## then re-run this data as above, with the new rank-order-averaged dropout scores. 
## to quantile normalize, merge Cont3 and Ts13 data together, set gene name as row name, and look only at BETA SCORES
setwd(ResultsFile)

# Merge MegaA HCT116: 
MegaA.Merge.HCT116<- merge(HCT116.Cont1.MegaA[,c(1,3)], HCT116.Ts8.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.HCT116<- merge(MegaA.Merge.HCT116, HCT116.Cont2.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.HCT116<- merge(MegaA.Merge.HCT116, HCT116.Tet5p.Ts5q.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.HCT116<- merge(MegaA.Merge.HCT116, HCT116.Tet5p.MegaA[,c(1,3)], by=c("Gene", "Gene"))


colnames(MegaA.Merge.HCT116)<- c("Gene", "H.Cont1.beta",  "H.Ts8.beta", 
                              "H.Cont2.beta", "H.Tet5pTs5q.beta", "H.Tet5p.beta"
                              ) 

# Remove all negative controls, except for three negative controls for heatmap plotting purposes. to show positive and negative controls. 
MegaA.Merge.HCT116$Gene<- as.character(MegaA.Merge.HCT116$Gene)
negControls<- MegaA.Merge.HCT116[which(grepl("neg", MegaA.Merge.HCT116$Gene)),] #list all negative controls
negControls<- negControls[-c(1,3,4),] #remove 3 negative controls
MegaA.Merge.HCT116 <- subset(MegaA.Merge.HCT116, ! Gene %in% negControls$Gene)   #remove rows with negative controls, except the three above

## Quantile normalize by batch and cell line: 
MegaA.Merge.HCT116.QN<- MegaA.Merge.HCT116

# Quantile normalize! 
MegaA.Merge.HCT116.QN<- quantile_normalisation(MegaA.Merge.HCT116[,c(2:length(MegaA.Merge.HCT116))])
MegaA.Merge.HCT116.QN<- as.data.frame(MegaA.Merge.HCT116.QN)
MegaA.Merge.HCT116.QN$Gene <- as.character(MegaA.Merge.HCT116$Gene)


## Check quantile normalization, they should have the same cutoffs for the 0, 25, 50, 75 and 100 places. 
quantile(MegaA.Merge.HCT116.QN[,1])
quantile(MegaA.Merge.HCT116.QN[,2])


# Plot Quantile normalized positive genes
dfm <- melt(MegaA.Merge.HCT116.QN)
dmf2<- subset(dfm, Gene %in% c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10", "neg11", "neg100"))
dmf2$Gene <- factor(dmf2$Gene, levels= c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10","neg11",  "neg100"))
ggplot(dmf2, aes(x = Gene, y = variable, fill = value)) + 
  geom_tile()+
  xlab("Gene")+
  ylab("Cell line") +
  scale_fill_gradient2(low="Red", mid="white", high="Blue", midpoint=0, 
                       name="Beta-Score") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.MegaA_HCT116.heatmap.controls.QuantileNorm.pdf
# 6x5


MegaA.Merge.HCT116.QN
setwd(ResultsFile)
# write.csv(MegaA.Merge.HCT116.QN, "CRISPRScreens_BetaScores_MegaA_HCT116.csv")




## Now make p-value and average difference for HCT116

MegaA.Merge.HCT116.QN_meanDiff<- MegaA.Merge.HCT116.QN

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell 
MegaA.Merge.HCT116.QN_meanDiff$Paired.Diff<- ((MegaA.Merge.HCT116.QN_meanDiff$H.Ts8.beta-MegaA.Merge.HCT116.QN_meanDiff$H.Cont1.beta) + 
                                              (MegaA.Merge.HCT116.QN_meanDiff$H.Tet5pTs5q.beta-MegaA.Merge.HCT116.QN_meanDiff$H.Cont2.beta) + 
                                              (MegaA.Merge.HCT116.QN_meanDiff$H.Tet5p.beta-MegaA.Merge.HCT116.QN_meanDiff$H.Cont2.beta)  
)/3

# Now add p-values to difference between all aneuploid and all euploid
MegaA.Merge.HCT116.QN_meanDiff$p.paired<- NA
for (i in 1:length(MegaA.Merge.HCT116.QN_meanDiff$Gene)){
  test<- t.test(as.numeric(MegaA.Merge.HCT116.QN_meanDiff[i,c(2,4,5)]), 
                as.numeric(MegaA.Merge.HCT116.QN_meanDiff[i,c(1,3,3)]), paired=TRUE) # t.test (aneuploid, euploid)
  MegaA.Merge.HCT116.QN_meanDiff$p.paired[i]<-test$p.value
} #paired t-test

MegaA.Merge.HCT116.QN_meanDiff$MeanBeta<- rowMeans(MegaA.Merge.HCT116.QN_meanDiff[,c(1:5)])


## Plot difference v pvalue (scatter plot)
ggplot(MegaA.Merge.HCT116.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff, 
                                           color=p.paired<0.05 & Paired.Diff< -0.15))+
  geom_point()+
  scale_color_manual(values=c("Black", "red"))+ 
  theme_classic()+
  xlab("HCT116: Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.BetaDelta.pValueHCT116_QN.pdf
subset(MegaA.Merge.HCT116.QN_meanDiff, p.paired<0.061 & Paired.Diff< -0.149 & MeanBeta < -0.2)$Gene
subset(MegaA.Merge.HCT116.QN_meanDiff, p.paired<0.01 & Paired.Diff< -0.3 & MeanBeta < -0.2)


#  "ANAPC15" "ATF6B"   "CARM1"   "EIF3H"   "GTF3C4"  "JMJD8"   "MAPK1"   "MBD2"    "NSMCE2"  "OVOL1"   "PKMYT1"  "SMARCA4"
#  "STK11"   "TCF7L2" 

#  Top hit?   NSMCE2 



## highlight top hits

HighlightTheseGenes<- subset(MegaA.Merge.HCT116.QN_meanDiff, p.paired<0.05 & Paired.Diff< -0.2)$Gene
ggplot(MegaA.Merge.HCT116.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff, color=Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(MegaA.Merge.HCT116.QN_meanDiff, Gene %in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"))+ 
  #geom_point(data = subset(MegaA.Merge.HCT116.QN_meanDiff, Gene %in% c("CARM1")), color="Blue")+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  theme(legend.position="none")+
  xlab("HCT116: Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.MegaA_HCT116.pairedP.Diff_Highlight.pdf
# Figure 4C
subset(MegaA.Merge.HCT116.QN_meanDiff, p.paired<0.005 & Paired.Diff< -0.25)


setwd(ResultsFile)
# write.csv(MegaA.Merge.HCT116.QN_meanDiff, "CRISPRScreens_BetaScores_MegaA_HCT116.csv")



     ##### DLD1 (8 pairs) ####


setwd(SchukkenData)
DLD1.Cont1.mle <- read.delim2("mle_DLD1_Cont.gene_summary.txt", 
                              dec=".", header = TRUE)
DLD1.Ts13.mle<- read.delim2("mle_DLD1_Ts13.gene_summary.txt", 
                            dec=".", header = TRUE)

MegaA.Merge.DLD1<- merge(DLD1.Cont1.MegaA[,c(1,3)],  DLD1.Ts13.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Ts8.10.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Ts2.18.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Cont3.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Ts10.21.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Ts5.15.MegaA[,c(1,3)], by=c("Gene", "Gene"))

MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Cont4.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Ts2.18_try2.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Ts12.17.MegaA[,c(1,3)], by=c("Gene", "Gene"))


MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MegaA.Merge.DLD1<- merge(MegaA.Merge.DLD1, DLD1.Ts13.mle[,c(1,3)], by.x="Gene", by.y="Gene")


colnames(MegaA.Merge.DLD1)<- c("Gene", 
                              "D.Cont1.beta", "D.Ts13.beta", "D.Ts8.10.beta", 
                              "D.Ts2.18.beta", 
                              "D.Cont3.beta", "D.Ts10.21.beta", "D.Ts5.15.beta", 
                              "D.Cont4.beta", "D.Ts2.18..beta", "D.Ts12.17.beta", 
                              "D.Cont1_WGL.beta", "D.Ts13_WGL.beta"
                              ) 

# Remove all negative controls, except for three negative controls for heatmap plotting purposes. to show positive and negative controls. 
MegaA.Merge.DLD1$Gene<- as.character(MegaA.Merge.DLD1$Gene)
negControls<- MegaA.Merge.DLD1[which(grepl("neg", MegaA.Merge.DLD1$Gene)),] #list all negative controls
negControls<- negControls[-c(1,3,4),] #remove 3 negative controls
MegaA.Merge.DLD1 <- subset(MegaA.Merge.DLD1, ! Gene %in% negControls$Gene)   #remove rows with negative controls, except the three above


## Quantile normalize by batch and cell line: 
MegaA.Merge.DLD1.QN<- MegaA.Merge.DLD1

# Quantile normalize! 
MegaA.Merge.DLD1.QN<- quantile_normalisation(MegaA.Merge.DLD1.QN[,c(2:length(MegaA.Merge.DLD1.QN))])
MegaA.Merge.DLD1.QN$Gene<- (MegaA.Merge.DLD1$Gene)


## Check quantile normalization, they should have the same cutoffs for the 0, 25, 50, 75 and 100 places. 
quantile(MegaA.Merge.DLD1.QN[,1])
quantile(MegaA.Merge.DLD1.QN[,2])


# Plot Quantile normalized positive genes
dfm <- melt(MegaA.Merge.DLD1.QN)
dmf2<- subset(dfm, Gene %in% c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10", "neg11", "neg100"))
dmf2$Gene <- factor(dmf2$Gene, levels= c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10", "neg11", "neg100"))
ggplot(dmf2, aes(x = Gene, y = variable, fill = value)) + 
  geom_tile()+
  xlab("Gene")+
  ylab("Cell line") +
  scale_fill_gradient2(low="Red", mid="white", high="Blue", midpoint=0, 
                       name="Beta-Score") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.MegaA_DLD1.heatmap.controls.QuantileNorm.pdf
# 6x5

MegaA.Merge.DLD1.QN


setwd(ResultsFile)
#write.csv(MegaA.Merge.DLD1.QN, "CRISPRScreens_BetaScores_MegaA_DLD1.csv")





## Now make p-valeu and average difference for DLD1

MegaA.Merge.DLD1.QN_meanDiff<- MegaA.Merge.DLD1.QN

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell 
MegaA.Merge.DLD1.QN_meanDiff$Paired.Diff<- ( (MegaA.Merge.DLD1.QN_meanDiff$D.Ts13.beta- MegaA.Merge.DLD1.QN_meanDiff$D.Cont1.beta) + 
                                                 (MegaA.Merge.DLD1.QN_meanDiff$D.Ts8.10.beta-MegaA.Merge.DLD1.QN_meanDiff$D.Cont1.beta) + 
                                                 (MegaA.Merge.DLD1.QN_meanDiff$D.Ts2.18.beta-MegaA.Merge.DLD1.QN_meanDiff$D.Cont1.beta) + 
                                                 (MegaA.Merge.DLD1.QN_meanDiff$D.Ts10.21.beta-MegaA.Merge.DLD1.QN_meanDiff$D.Cont3.beta) +
                                               (MegaA.Merge.DLD1.QN_meanDiff$D.Ts2.18..beta-MegaA.Merge.DLD1.QN_meanDiff$D.Cont4.beta) +
                                               (MegaA.Merge.DLD1.QN_meanDiff$D.Ts12.17.beta-MegaA.Merge.DLD1.QN_meanDiff$D.Cont4.beta) +
                                               (MegaA.Merge.DLD1.QN_meanDiff$D.Ts5.15.beta-MegaA.Merge.DLD1.QN_meanDiff$D.Cont3.beta) +
                                               (MegaA.Merge.DLD1.QN_meanDiff$D.Ts13_WGL.beta-MegaA.Merge.DLD1.QN_meanDiff$D.Cont1_WGL.beta) 
)/8

# Now add p-values to difference between all aneuploid and all euploid
MegaA.Merge.DLD1.QN_meanDiff$p.paired<- NA
for (i in 1:length(MegaA.Merge.DLD1.QN_meanDiff$Gene)){
  test<- t.test(as.numeric(MegaA.Merge.DLD1.QN_meanDiff[i,c(2,3,4,6,7,9,10,12)]), #aneuploid
                as.numeric(MegaA.Merge.DLD1.QN_meanDiff[i,c(1,1,1,5,5,8,8,11)]), paired=TRUE) # t.test (aneuploid, euploid)
  MegaA.Merge.DLD1.QN_meanDiff$p.paired[i]<-test$p.value
} #paired t-test

MegaA.Merge.DLD1.QN_meanDiff$MeanBeta<- rowMeans(MegaA.Merge.DLD1.QN_meanDiff[,c(1:12)])



## highlight hits

HighlightTheseGenes<- subset(MegaA.Merge.DLD1.QN_meanDiff, p.paired<0.05 & Paired.Diff< -0.2)$Gene
ggplot(MegaA.Merge.DLD1.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff, color=Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(MegaA.Merge.DLD1.QN_meanDiff, Gene %in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"))+ 
  #geom_point(data = subset(MegaA.Merge.DLD1.QN_meanDiff, Gene %in% c("PRKDC")), color="Blue")+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  theme(legend.position="none")+
  xlab("DLD1: Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.MegaA_DLD1.pairedP.Diff_Highlight.pdf
# Figure 4C

subset(MegaA.Merge.DLD1.QN_meanDiff, p.paired<0.01 & Paired.Diff< -0.4)



setwd(ResultsFile)
# write.csv(MegaA.Merge.DLD1.QN_meanDiff, "CRISPRScreens_BetaScores_MegaA_DLD1.csv")



     ##### A2780 (5 pairs) ####

# A2780  cells WGL data: 
setwd(SchukkenData)

A2780.Cont1.mle <- read.delim2("mle_count-A2780-Cont1-d0-d10.gene_summary.txt", 
                               dec=".", header = TRUE)
A2780.Di1q.mle<- read.delim2("mle_count-A2780-1qDi.gene_summary.txt", 
                             dec=".", header = TRUE)


MegaA.Merge.A2780<- merge(A2780.Cont1.MegaA_forQN[,c(1,3)], A2780.Ts2.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene")) ### This one is wierd, low coverage leads to extreme negative beta scores
MegaA.Merge.A2780<- merge(MegaA.Merge.A2780, A2780.Ts10.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.A2780<- merge(MegaA.Merge.A2780, A2780.Cont2.MegaA[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.A2780<- merge(MegaA.Merge.A2780, A2780.Ts8.MegaA_forQN[,c(1,3)], by=c("Gene", "Gene"))
MegaA.Merge.A2780<- merge(MegaA.Merge.A2780, A2780.Ts18.MegaA[,c(1,3)], by=c("Gene", "Gene"))

MegaA.Merge.A2780<- merge(MegaA.Merge.A2780, A2780.Cont1.mle[,c(1,3)], by.x="Gene", by.y="Gene")
MegaA.Merge.A2780<- merge(MegaA.Merge.A2780, A2780.Di1q.mle[,c(1,3)], by.x="Gene", by.y="Gene")

colnames(MegaA.Merge.A2780)<- c("Gene", 
                              "A.Cont1.beta", "A.Ts2.beta", "A.Ts10.beta", 
                              "A.Cont2.beta", "A.Ts8.beta", "A.Ts18.beta", 
                              "A.Cont1_WGL.beta", "A.Di1q_WGL.beta") 


# Remove all negative controls, except for three negative controls for heatmap plotting purposes. to show positive and negative controls. 
MegaA.Merge.A2780$Gene<- as.character(MegaA.Merge.A2780$Gene)
negControls<- MegaA.Merge.A2780[which(grepl("neg", MegaA.Merge.A2780$Gene)),] #list all negative controls
negControls<- negControls[-c(1,3,4),] #remove 3 negative controls
MegaA.Merge.A2780 <- subset(MegaA.Merge.A2780, ! Gene %in% negControls$Gene)   #remove rows with negative controls, except the three above


## Quantile normalize by batch and cell line: 
MegaA.Merge.A2780.QN<- MegaA.Merge.A2780


# Quantile normalize! 
MegaA.Merge.A2780.QN<- quantile_normalisation(MegaA.Merge.A2780.QN[,c(2:length(MegaA.Merge.A2780.QN))])
MegaA.Merge.A2780.QN$Gene<- (MegaA.Merge.A2780$Gene)


## Check quantile normalization, they should have the same cutoffs for the 0, 25, 50, 75 and 100 places. 
quantile(MegaA.Merge.A2780.QN[,1])
quantile(MegaA.Merge.A2780.QN[,2])


# Plot Quantile normalized positive genes
dfm <- melt(MegaA.Merge.A2780.QN)
dmf2<- subset(dfm, Gene %in% c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10", "neg11", "neg100"))
dmf2$Gene <- factor(dmf2$Gene, levels= c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10", "neg11", "neg100"))
ggplot(dmf2, aes(x = Gene, y = variable, fill = value)) + 
  geom_tile()+
  xlab("Gene")+
  ylab("Cell line") +
  scale_fill_gradient2(low="Red", mid="white", high="Blue", midpoint=0, 
                       name="Beta-Score") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.MegaA_DLD1.heatmap.controls.QuantileNorm.pdf
# 6x5

MegaA.Merge.A2780.QN

setwd(ResultsFile)
#write.csv(MegaA.Merge.A2780.QN, "CRISPRScreens_BetaScores_MegaA_A2780.csv")




## Now make p-valeu and average difference for A2780

MegaA.Merge.A2780.QN_meanDiff<- MegaA.Merge.A2780.QN

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell 
MegaA.Merge.A2780.QN_meanDiff$Paired.Diff<- ( (MegaA.Merge.A2780.QN_meanDiff$A.Ts2.beta - MegaA.Merge.A2780.QN_meanDiff$A.Cont1.beta) + 
                                               (MegaA.Merge.A2780.QN_meanDiff$A.Ts10.beta - MegaA.Merge.A2780.QN_meanDiff$A.Cont1.beta) + 
                                               (MegaA.Merge.A2780.QN_meanDiff$A.Ts8.beta - MegaA.Merge.A2780.QN_meanDiff$A.Cont2.beta) + 
                                               (MegaA.Merge.A2780.QN_meanDiff$A.Ts18.beta - MegaA.Merge.A2780.QN_meanDiff$A.Cont2.beta) +
                                                
                                                (MegaA.Merge.A2780.QN_meanDiff$A.Cont1_WGL.beta - MegaA.Merge.A2780.QN_meanDiff$A.Di1q_WGL.beta) 
)/5

# Now add p-values to difference between all aneuploid and all euploid
MegaA.Merge.A2780.QN_meanDiff$p.paired<- NA
for (i in 1:length(MegaA.Merge.A2780.QN_meanDiff$Gene)){
  test<- t.test(as.numeric(MegaA.Merge.A2780.QN_meanDiff[i,c(2,3,5,6,7)]), 
                as.numeric(MegaA.Merge.A2780.QN_meanDiff[i,c(1,1,4,4,8)]), paired=TRUE) # t.test (aneuploid, euploid)
  MegaA.Merge.A2780.QN_meanDiff$p.paired[i]<-test$p.value
} #paired t-test

MegaA.Merge.A2780.QN_meanDiff$MeanBeta<- rowMeans(MegaA.Merge.A2780.QN_meanDiff[,c(1:8)])
MegaA.Merge.A2780.QN_meanDiff<- MegaA.Merge.A2780.QN_meanDiff[order(MegaA.Merge.A2780.QN_meanDiff$MeanBeta),]


## Highlight hits
HighlightTheseGenes <- subset(MegaA.Merge.A2780.QN_meanDiff, p.paired<0.05 & Paired.Diff< -0.2)$Gene

ggplot(MegaA.Merge.A2780.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff, color=Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(MegaA.Merge.A2780.QN_meanDiff, Gene %in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"))+ 
  #geom_point(data = subset(MegaA.Merge.A2780.QN_meanDiff, Gene %in% c("PRSS42")), color="Blue")+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  theme(legend.position="none")+
  xlab("A2780: Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.MegaA_A2780.pairedP.Diff_Highlight.pdf

subset(MegaA.Merge.A2780.QN_meanDiff, Gene %in% HighlightTheseGenes)


subset(MegaA.Merge.A2780.QN_meanDiff, p.paired<0.005 & Paired.Diff< -0.4)$Gene


setwd(ResultsFile)
#write.csv(MegaA.Merge.A2780.QN_meanDiff, "CRISPRScreens_BetaScores_MegaA_A2780.csv")




     ##### SNU1 (4 pairs) ####
# SNU1  cells data: #

MegaA.Merge.SNU1<- MegaA.Merge.ALL[,c(1, 20:24)]

colnames(MegaA.Merge.SNU1)<- c("Gene", "SNU1 WT", "SNU1 C12", "SNU1 C24", "SNU1 C111", "SNU1 C114") 


# Remove all negative controls, except for three negative controls for heatmap plotting purposes. to show positive and negative controls. 
MegaA.Merge.SNU1$Gene<- as.character(MegaA.Merge.SNU1$Gene)
negControls<- MegaA.Merge.SNU1[which(grepl("neg", MegaA.Merge.SNU1$Gene)),] #list all negative controls. Not needed: neg controls already removed above
negControls<- negControls[-c(1,3,4),] #remove 3 negative controls for plot. 
MegaA.Merge.SNU1 <- subset(MegaA.Merge.SNU1, ! Gene %in% negControls$Gene)   #remove rows with negative controls, except the three above


## Quantile normalize by batch and cell line: 
MegaA.Merge.SNU1.QN<- MegaA.Merge.SNU1

# Quantile normalize! 
MegaA.Merge.SNU1.QN<- quantile_normalisation(MegaA.Merge.SNU1.QN[,c(2:length(MegaA.Merge.SNU1.QN))])
MegaA.Merge.SNU1.QN<- as.data.frame(MegaA.Merge.SNU1.QN)
MegaA.Merge.SNU1.QN$Gene<- MegaA.Merge.SNU1$Gene


## Check quantile normalization, they should have the same cutoffs for the 0, 25, 50, 75 and 100 places. 
quantile(MegaA.Merge.SNU1.QN[,1])
quantile(MegaA.Merge.SNU1.QN[,2])


# Plot Quantile normalized positive genes
dfm <- melt(MegaA.Merge.SNU1.QN)
dmf2<- subset(dfm, Gene %in% c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10", "neg11", "neg100"))
dmf2$Gene <- factor(dmf2$Gene, levels= c("CDK1", "CDK9", "PCNA", "POLR2A", "POLR2D", "RPA3", "RPL23A", "neg1", "neg10", "neg11", "neg100"))
ggplot(dmf2, aes(x = Gene, y = variable, fill = value)) + 
  geom_tile()+
  xlab("Gene")+
  ylab("Cell line") +
  scale_fill_gradient2(low="Red", mid="white", high="Blue", midpoint=0, 
                       name="Beta-Score") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# plot.MegaA_SNU1.heatmap.controls.QuantileNorm.pdf
# 6x5

MegaA.Merge.SNU1.QN

setwd(ResultsFile)
# write.csv(MegaA.Merge.SNU1.QN, "CRISPRScreens_BetaScores_MegaA_SNU1.csv")




## Now make p-value and average difference for SNU1

MegaA.Merge.SNU1.QN_meanDiff<- MegaA.Merge.SNU1.QN

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell 
MegaA.Merge.SNU1.QN_meanDiff$Paired.Diff<- ( (MegaA.Merge.SNU1.QN_meanDiff$`SNU1 C12` - MegaA.Merge.SNU1.QN_meanDiff$`SNU1 WT`) + 
                                                #(MegaA.Merge.SNU1.QN_meanDiff$`SNU1 C23` - MegaA.Merge.SNU1.QN_meanDiff$`SNU1 WT`) + 
                                                (MegaA.Merge.SNU1.QN_meanDiff$`SNU1 C24` - MegaA.Merge.SNU1.QN_meanDiff$`SNU1 WT`) + 
                                                (MegaA.Merge.SNU1.QN_meanDiff$`SNU1 C111` - MegaA.Merge.SNU1.QN_meanDiff$`SNU1 WT`) + 
                                               (MegaA.Merge.SNU1.QN_meanDiff$`SNU1 C114` - MegaA.Merge.SNU1.QN_meanDiff$`SNU1 WT`)
)/4

# Now add p-values to difference between all aneuploid and all euploid
MegaA.Merge.SNU1.QN_meanDiff$p.paired<- NA
for (i in 1:length(MegaA.Merge.SNU1.QN_meanDiff$Gene)){
  test<- t.test(as.numeric(MegaA.Merge.SNU1.QN_meanDiff[i,c(2,3,4,5)]), 
                as.numeric(MegaA.Merge.SNU1.QN_meanDiff[i,c(1,1,1,1)]), paired=TRUE) # t.test (aneuploid, euploid)
  MegaA.Merge.SNU1.QN_meanDiff$p.paired[i]<-test$p.value
} #paired t-test

MegaA.Merge.SNU1.QN_meanDiff$MeanBeta<- rowMeans(MegaA.Merge.SNU1.QN_meanDiff[,c(1:5)])
MegaA.Merge.SNU1.QN_meanDiff<- MegaA.Merge.SNU1.QN_meanDiff[order(MegaA.Merge.SNU1.QN_meanDiff$MeanBeta),]


## Highlight hits
HighlightTheseGenes<-  c("FOSL1", "UBE2H")
HighlightTheseGenes<-  subset(MegaA.Merge.SNU1.QN_meanDiff, Paired.Diff< -0.2 & p.paired<0.05)$Gene
ggplot(MegaA.Merge.SNU1.QN_meanDiff, aes(y=-log(p.paired,2), x=Paired.Diff, color=Gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(MegaA.Merge.SNU1.QN_meanDiff, Gene %in% HighlightTheseGenes))+
  scale_color_manual(values=c("Black", "Red"))+ 
  #geom_point(data = subset(MegaA.Merge.SNU1.QN_meanDiff, Gene %in% c("FURIN")), color = "blue")+
  theme_classic()+
  geom_density_2d(color = "white", bins = 10)+
  theme(legend.position="none")+
  xlab("A2780: Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.MegaA_A2780.pairedP.Diff_Highlight.pdf
# Figure XXC
subset(MegaA.Merge.SNU1.QN_meanDiff, Gene %in% HighlightTheseGenes)


setwd(ResultsFile)
# write.csv(MegaA.Merge.SNU1.QN_meanDiff, "CRISPRScreens_BetaScores_MegaA_SNU1.csv")


     ##### Top 10% and 20% in all four cell lines groups #### 

# find overall in all four cell line specific analysis. 
HCT116.Top10<-MegaA.Merge.HCT116.QN_meanDiff$Gene[MegaA.Merge.HCT116.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.HCT116.QN_meanDiff$Paired.Diff, 0.1, na.rm = TRUE)]
DLD1.Top10<-MegaA.Merge.DLD1.QN_meanDiff$Gene[MegaA.Merge.DLD1.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.DLD1.QN_meanDiff$Paired.Diff, 0.1, na.rm = TRUE)]
A2780.Top10<-MegaA.Merge.A2780.QN_meanDiff$Gene[MegaA.Merge.A2780.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.A2780.QN_meanDiff$Paired.Diff, 0.1, na.rm = TRUE)]
SNU1.Top10<-MegaA.Merge.SNU1.QN_meanDiff$Gene[MegaA.Merge.SNU1.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.SNU1.QN_meanDiff$Paired.Diff, 0.1, na.rm = TRUE)]

common_genes_4CellLines <- Reduce(
  intersect,
  list(HCT116.Top10, DLD1.Top10, A2780.Top10, SNU1.Top10)
)
common_genes_4CellLines
# Top 10% in all four cell line groups: 
# "CAPN1"  "FBXO41" "PKMYT1" "POU2F1" "PSMD14" "SPI1" 


# find overlap in top 20% of aneuploid dependencies 
HCT116.Top20<-MegaA.Merge.HCT116.QN_meanDiff$Gene[MegaA.Merge.HCT116.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.HCT116.QN_meanDiff$Paired.Diff, 0.2, na.rm = TRUE)]
DLD1.Top20<-MegaA.Merge.DLD1.QN_meanDiff$Gene[MegaA.Merge.DLD1.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.DLD1.QN_meanDiff$Paired.Diff, 0.2, na.rm = TRUE)]
A2780.Top20<-MegaA.Merge.A2780.QN_meanDiff$Gene[MegaA.Merge.A2780.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.A2780.QN_meanDiff$Paired.Diff, 0.2, na.rm = TRUE)]
SNU1.Top20<-MegaA.Merge.SNU1.QN_meanDiff$Gene[MegaA.Merge.SNU1.QN_meanDiff$Paired.Diff <= quantile(MegaA.Merge.SNU1.QN_meanDiff$Paired.Diff, 0.2, na.rm = TRUE)]

common_genes_4CellLines <- Reduce(
  intersect,
  list(HCT116.Top20, DLD1.Top20, A2780.Top20, SNU1.Top20)
)
common_genes_4CellLines

# Top 20% in all four cell line groups: 
#  [1] "ARID3C" "CAPN1"  "CARM1"  "DLX3"   "FBXO41" "FOSL1"  "HSF1"   "LITAF"  "PKMYT1" "PLK3"   "POU2F1" "PRSS3"  "PSMD14" "RFX1"  
# [15] "SPI1"   "STK32C" "THAP10" "UBE2H"  "ULK4"   "VDR"    "ZFPL1"  "ZNF703" "ZNF784"


# random overlap for these datasets: 
Random10<- (341*(341/3411)^4) # 0.03406003 genes would randomly overlap in top 10%
Random20<- (682*(682/3411)^4) # 1.089921 genes would randomly overlap in top 20%



     ##### Colorectal Cancer: HCT116 DLD1 and Vaco432 ####
## Quantile normalize the dropouts for colorectal cancers (HCT116, DLD1 and Vaco432)


Merge.CRC<- Merge.ALL.QN2

# Remove all negative controls, except for three negative controls for heatmap plotting purposes. to show positive and negative controls. 
Merge.CRC$Gene<- as.character(Merge.CRC$Gene)
negControls<- Merge.CRC[which(grepl("neg", Merge.CRC$Gene)),] #list all negative controls
negControls<- negControls[-c(1,3,4),] #remove 3 negative controls
Merge.CRC <- subset(Merge.CRC, ! Gene %in% negControls$Gene)   #remove rows with negative controls, except the three above




Merge.CRC$CRC.p.paired<-NA
Merge.CRC$CRC.Diff.paired<-NA
Merge.CRC$CRC.MeanBeta<-NA
for (i in 1: length(Merge.CRC$CRC.p.paired)){
  t<- t.test(as.numeric(Merge.CRC[i,c(1, 3,3, 6,6,6, 10,10, 13,13, 27, 29)]), #control
             as.numeric(Merge.CRC[i,c(2, 4,5, 7,8,9, 11,12, 14,15, 28, 30)]), #aneuploid
             paired=TRUE)
  Merge.CRC$CRC.p.paired[i]<-t$p.value
  Merge.CRC$CRC.MeanBeta[i]<- mean(as.numeric(Merge.CRC[i,c(1:15,27:30)]))
}
Merge.CRC$CRC.Diff.paired<-( (Merge.CRC$HCT116.Ts8.beta-Merge.CRC$HCT116.Cont1.beta) + 
                                  (Merge.CRC$HCT116.Tet5p.Ts5q.beta-Merge.CRC$HCT116.Cont2.beta) + 
                                  (Merge.CRC$HCT116.Tet5.beta-Merge.CRC$HCT116.Cont2.beta) +
                                  (Merge.CRC$DLD1.Ts13.beta- Merge.CRC$DLD1.Cont1.beta) +
                                  (Merge.CRC$DLD1.Ts8.10.beta- Merge.CRC$DLD1.Cont1.beta) +
                                  (Merge.CRC$DLD1.Ts2.18.beta- Merge.CRC$DLD1.Cont1.beta) +
                                  (Merge.CRC$DLD1.Ts10.21.beta- Merge.CRC$DLD1.Cont3.beta) +
                                  (Merge.CRC$DLD1.Ts5.15.beta- Merge.CRC$DLD1.Cont3.beta) +
                                  (Merge.CRC$DLD1.Ts2.18.beta_Try2- Merge.CRC$DLD1.Cont4.beta) +
                                  (Merge.CRC$DLD1.Ts12.17.beta- Merge.CRC$DLD1.Cont4.beta) +
                                  (Merge.CRC$WGL.DLD1.Ts13.beta- Merge.CRC$WGL.DLD1.Cont.beta) +
                                  (Merge.CRC$WGL.Vaco432.Ts13.beta- Merge.CRC$WGL.Vaco432.Cont.beta) 
)/12


subset(Merge.CRC, CRC.p.paired<0.05 & CRC.Diff.paired< -0.15)[c(40,41,42,43)]

HighlightTheseGenes<- subset(Merge.CRC, CRC.p.paired< 0.05 &CRC.Diff.paired< -0.2)$Gene


ggplot(Merge.CRC, aes(y=-log(CRC.p.paired,2), x=CRC.Diff.paired, 
                      color= CRC.p.paired<0.05 
                      & CRC.Diff.paired< -0.2))+
  geom_point()+
  scale_color_manual(values=c("Black", "red"))+ 
  # geom_point(data = subset(Merge.CRC, Gene == "PKMYT1"), color="blue")+
  theme_classic()+
  theme(legend.position="none")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Colorectal Cancer Only: Beta delta\n aneuploid - euploid")+
  ylab("-log2 (p-value)")
# plot.MegaA.WGL.CRC_TopHits.pdf 




# Plot RANK CRC

Merge.CRC$AneuploidRank<- rank(Merge.CRC$CRC.Diff.paired)
ggplot(Merge.CRC, aes(x=AneuploidRank, y=CRC.Diff.paired, color= Gene == "ILK"))+
  geom_point()+
  geom_point(data = subset(Merge.CRC, Gene == "ILK"))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  ylab("Beta delta (aneuploid - euploid)")+
  xlab("Mega-A Colorectal: Rank")
# 5x3
# plot.MegaA.Colorectal.Rank_ILK.pdf
subset(Merge.CRC, Gene == "UBE2H")

# 
# top hits: "CARM1"   "CDK11B"  "FOSL1"   "HM13"    "ILK"     "PKMYT1"  "PRKDC"   "PRKRIR"  "RNF166"  "SMARCA4" 
# "STK11"   "TFAM"    "TOX4"    "UBE2N"   "UCHL5"   "ZSCAN10"


## highlight hits
HighlightTheseGenes<- c("UBE2H") #
HighlightTheseGenes<- subset(CRISPR_DepMap_AneuScore, SigDiff==TRUE)$Gene


ggplot(Merge.CRC, aes(y=-log(CRC.p.paired,2), x=CRC.Diff.paired))+
  geom_point()+
  geom_point(data = subset(Merge.CRC, CRC.p.paired< 0.05 & CRC.Diff.paired< -0.15), color = "red")+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  #geom_vline(xintercept = 0)+
  geom_density_2d(color = "white", bins = 10)+
  theme(legend.position="none")+ 
  xlab("CRC: Beta delta (aneuploid - euploid)")+ 
  ylab("-log2(p-value)")
# 5x4
# plot.MegaA_CRC.pairedP.Diff_ILK.pdf
subset(Merge.CRC, Gene %in% HighlightTheseGenes)


setwd(ResultsFile)
# write.csv(Merge.CRC, "CRISPRScreens_BetaScores_MegaA_Colorectal.csv")
# Merge.CRC<- read.csv("CRISPRScreens_BetaScores_MegaA_Colorectal.csv")


## Compare the Colorectal cancer analysis in DepMap.
# get data from DepMap_Drug_AneuPloidyScore.R

CRISPR_DepMap_AneuScore_Colorectal<- read.csv("DepMap_AneuScore.CRISPR.Corr_Colorectal.csv")

CRC.MegaA.WGL.DepMap<- merge( CRISPR_DepMap_AneuScore_Colorectal, Merge.CRC, by.x = "Gene", by.y= "Gene")

HighlightTheseGenes<- c("HPN"  ,  "RNF166", "TOX4" ,  "UCHL5" )
HighlightTheseGenes<- subset(CRC.MegaA.WGL.DepMap, Aneu.CRISPR.Corr.Coef< -0.1 & CRC.Diff.paired< -0.2)$Gene


ggplot(CRC.MegaA.WGL.DepMap, aes(y=Aneu.CRISPR.Corr.Coef, x=CRC.Diff.paired))+
  geom_point()+
  geom_point(data = subset(CRC.MegaA.WGL.DepMap, Gene %in% HighlightTheseGenes), color = "red")+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Colorectal Cancer Only: \n Diff Beta Sheltzer lab (aneuploid - euploid)")+
  ylab("Colorectal Cancer Only: \n DepMap correlation Score")
# 5x4
# plot.CRC.MegaA.WGL_DepMap.pdf


subset(CRC.MegaA.WGL.DepMap, Aneu.CRISPR.Corr.Coef< -0.1 & CRC.Diff.paired< -0.2)$Gene
# "ASCL2"  "FOSL1"  "HDAC3"  "MYT1"   "NR0B1"  "PAX3"   "RNF166" "TRIM61" "ZFP92" 

subset(CRC.MegaA.WGL.DepMap, Aneu.CRISPR.Corr.Coef < -0.25 & CRC.Diff.paired < -0.1 & 
         Mean.Beta < -0.1 & Aneu.CRISPR.Corr.P< 0.05 & CRC.p.paired< 0.05)




         ###### Non Colorectal cancer #### 
# I want to compare CRC specific hits to non-CRC hits to differentiate between pan cancer aneuploid hits and CRC specific hits

Merge.nonCRC<- Merge.ALL.QN2

Merge.nonCRC$ANvEU.P.Paired.nonCRC<- NA  #difference btween average an and average euploid beta score (aneuploid beta minus euploid beta)
Merge.nonCRC$ANvEU.Beta.Paired.nonCRC<- NA 
Merge.nonCRC<- Merge.nonCRC %>% relocate(Gene, .after = last_col()) #make sure Gene is at end of columns


for (i in 1:length(Merge.nonCRC$Gene)){
  test<-t.test(as.numeric(Merge.nonCRC[i, c(16,16, 19,19, 22,22,22,22, 32, 34,35, 37, 39)]), 
               as.numeric(Merge.nonCRC[i, c(17,18, 20,21, 23,24,25,26, 31, 33,33, 36, 38)]), 
               alternative=c("two.sided"), paired=TRUE) 
  Merge.nonCRC$ANvEU.P.Paired.nonCRC[i]<- test$p.value
  Merge.nonCRC$ANvEU.Beta.Paired.nonCRC[i]<- -1*test$estimate #make negative values more toxic to aneuploid
}
Merge.nonCRC<- Merge.nonCRC[order(Merge.nonCRC$ANvEU.Beta.Paired.nonCRC),]
Merge.nonCRC<- Merge.nonCRC %>% relocate("Gene") #more "Gene" gene names back to begining of columns


SigAllAneu<- subset(Merge.nonCRC, ANvEU.Beta.Paired.nonCRC< -0.2 & ANvEU.P.Paired.nonCRC< 0.05)$Gene
ggplot(Merge.nonCRC, aes(x= ANvEU.Beta.Paired.nonCRC, 
                          y= -log2(ANvEU.P.Paired.nonCRC), color= Gene %in% SigAllAneu))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  geom_point(data = subset(Merge.nonCRC, Gene%in% SigAllAneu))+
  ylab("-log2(P-value)")+
  xlab("Difference in Beta score\nAneuploid - euploid")+
  scale_color_manual(values = c("Black", "red"), name = "P<0.05 \nDiff < -.1")+
  theme_classic()
# 5x4
# plot.BetaDelta.ANvEU_WGL.nonCRC.pdf


CRC.vs.nonCRC<- merge(Merge.CRC, Merge.nonCRC[, c(1,41,42)], by.x = "Gene", by.y = "Gene")

SigAllAneu<- subset(CRC.vs.nonCRC, CRC.Diff.paired< -0.2 & 
                      CRC.p.paired< 0.05 & 
                      ANvEU.Beta.Paired.nonCRC > 0)$Gene

SigAllAneu<-subset(CRC.vs.nonCRC, CRC.Diff.paired< -0.2 & 
            CRC.p.paired< 0.05 & 
              ANvEU.Beta.Paired.nonCRC > 0)$Gene


ggplot(CRC.vs.nonCRC, aes(x= CRC.Diff.paired, 
                          y= ANvEU.Beta.Paired.nonCRC, 
                          color= Gene %in% SigAllAneu))+
  geom_point()+
  geom_hline(yintercept =0)+
  geom_vline(xintercept =0)+
  geom_density_2d(color = "white", bins = 10)+
  geom_point(data = subset(CRC.vs.nonCRC, Gene%in% SigAllAneu))+
  #geom_point(data = subset(CRC.vs.nonCRC, Gene == "TOX4"), color = "blue")+
  ylab("Difference in Beta score (Non-CRC)")+
  xlab("Difference in Beta score (CRC)")+
  scale_color_manual(values = c("Black", "red"), 
                     name = "CRC specific")+
  theme_classic()
# 5x4
# plot.BetaDelta.CRC.nonCRC.pdf 

#  "HPN"    "RNF166" "TOX4"   "UCHL5" 

subset(CRC.vs.nonCRC, CRC.Diff.paired< -0.2 & 
         CRC.p.paired< 0.05 & 
         ANvEU.Beta.Paired.nonCRC > 0)$Gene
#  "HPN"    "RNF166" "TOX4"!   "UCHL5" 


subset(CRC.vs.nonCRC, CRC.Diff.paired< -0.15 & 
         CRC.p.paired< 0.05 & 
         ANvEU.Beta.Paired.nonCRC < -0.15 &
         ANvEU.P.Paired.nonCRC< 0.05)$Gene #Both CRC & non CRC
# "CARM1" "FOSL1" "RAX2"  "VPS18"
# "ASCL2"  "CARM1"  "CHD8"   "CSNK1E" "FBXO42" "FOSL1"  "GATA5"  "GLI4"   "HOXA1"  "HSF1"   "MYT1"   "PKMYT1" "RAX2"   "SIX6"   "SPI1"   "SPOP"  "VPS18" 

# write.csv(CRC.vs.nonCRC, "CRC.vs.nonCRC.csv")


     ##### Stomach Cancer: SNU1 and AGS (5 pairs) ####
## Quantile normalize the dropouts for Stomach cancers (SNU1 and AGS)


Merge.St<- Merge.ALL.QN2

# Remove all negative controls, except for three negative controls for heatmap plotting purposes. to show positive and negative controls. 
Merge.St$Gene<- as.character(Merge.St$Gene)
negControls<- Merge.St[which(grepl("neg", Merge.St$Gene)),] #list all negative controls
negControls<- negControls[-c(1,3,4),] #remove 3 negative controls
Merge.St <- subset(Merge.St, ! Gene %in% negControls$Gene)   #remove rows with negative controls, except the three above




Merge.St$Stomach.p.paired<-NA
Merge.St$Stomach.Diff.paired<-NA
Merge.St$Stomach.MeanBeta<-NA
for (i in 1: length(Merge.St$Stomach.p.paired)){
  t<- t.test(as.numeric(Merge.St[i,c(22,22,22,22,37)]), #control
             as.numeric(Merge.St[i,c(23,24,25,26,36)]), #aneuploid
             paired=TRUE)
  Merge.St$Stomach.p.paired[i]<-t$p.value
  Merge.St$Stomach.MeanBeta[i]<- mean(as.numeric(Merge.St[i,c(22:26)]))
}
Merge.St$Stomach.Diff.paired<-( (Merge.St$Snu1.Ts10.17.5_Ps2.beta-Merge.St$Snu1.Cont.beta) + 
                               (Merge.St$Snu1.Ts10.1.8.14_Ps5.beta-Merge.St$Snu1.Cont.beta) + 
                               (Merge.St$Snu1.Ts10.6q_Ps8.12.beta-Merge.St$Snu1.Cont.beta) +
                               (Merge.St$Snu1.Ts9.19.beta- Merge.St$Snu1.Cont.beta) +
                               (Merge.St$WGL.AGS.Cont.beta- Merge.St$WGL.AGS.Di1q.beta)
                               
)/5


subset(Merge.St, Stomach.p.paired<0.05 & Stomach.Diff.paired< -0.15)[c(40,41,42,43)]
#


HighlightTheseGenes<- subset(Merge.St, Stomach.p.paired< 0.05 & Stomach.Diff.paired< -0.2)$Gene


ggplot(Merge.St, aes(y=-log(Stomach.p.paired,2), x=Stomach.Diff.paired, 
                      color= Stomach.p.paired<0.05 
                      & Stomach.Diff.paired< -0.2))+
  geom_point()+
  scale_color_manual(values=c("Black", "red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  geom_density_2d(color = "white", bins = 10)+
  xlab("Stomach Cancer Only: Beta delta\n aneuploid - euploid")+
  ylab("-log2 (p-value)")
# plot.MegaA.WGL.Stomach_TopHits.pdf 





## highlight hits
HighlightTheseGenes<- c("PLK4") #
HighlightTheseGenes<- subset(CRISPR_DepMap_AneuScore, SigDiff==TRUE)$Gene


ggplot(Merge.St, aes(y=-log(Stomach.p.paired,2), x=Stomach.Diff.paired))+
  geom_point()+
  geom_point(data = subset(Merge.St, Stomach.p.paired< 0.05 & Stomach.Diff.paired< -0.15), color = "red")+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  #geom_vline(xintercept = 0)+
  geom_density_2d(color = "white", bins = 10)+
  theme(legend.position="none")+ 
  xlab("Stomach Cancer: Beta delta (aneuploid - euploid)")+ 
  ylab("-log2(p-value)")
# 5x4
# plot.MegaA_Stomach.pairedP.Diff.pdf
subset(Merge.St, Gene %in% HighlightTheseGenes)


setwd(ResultsFile)
# write.csv(Merge.St, "CRISPRScreens_BetaScores_MegaA_Stomach.csv")
# Merge.St<- read.csv("CRISPRScreens_BetaScores_MegaA_Stomach.csv")


         ###### Non Stomach cancer #### 
# I want to compare Stomach specific hits to non-Stomach hits to differentiate between pan cancer aneuploid hits and Stomach specific hits

Merge.nonStomach<- Merge.ALL.QN2

Merge.nonStomach$ANvEU.P.Paired.nonStomach<- NA  #difference btween average an and average euploid beta score (aneuploid beta minus euploid beta)
Merge.nonStomach$ANvEU.Beta.Paired.nonStomach<- NA 
Merge.nonStomach<- Merge.nonStomach %>% relocate(Gene, .after = last_col()) #make sure Gene is at end of columns


for (i in 1:length(Merge.nonStomach$Gene)){
  test<-t.test(as.numeric(Merge.nonStomach[i, c(1, 3,3, 6,6,6, 10,10, 13,13, 16,16, 19,19, 27, 29, 32, 34,35, 39)]), 
               as.numeric(Merge.nonStomach[i, c(2, 4,5, 7,8,9, 11,12, 14,15, 17,18, 20,21, 18, 30, 31, 33,33, 38)]), 
               alternative=c("two.sided"), paired=TRUE) 
  Merge.nonStomach$ANvEU.P.Paired.nonStomach[i]<- test$p.value
  Merge.nonStomach$ANvEU.Beta.Paired.nonStomach[i]<- -1*test$estimate #make negative values more toxic to aneuploid
}
Merge.nonStomach<- Merge.nonStomach[order(Merge.nonStomach$ANvEU.Beta.Paired.nonStomach),]
Merge.nonStomach<- Merge.nonStomach %>% relocate("Gene") #more "Gene" gene names back to begining of columns


SigAllAneu<- subset(Merge.nonStomach, ANvEU.Beta.Paired.nonStomach< -0.2 & ANvEU.P.Paired.nonStomach< 0.05)$Gene

ggplot(Merge.nonStomach, aes(x= ANvEU.Beta.Paired.nonStomach, 
                         y= -log2(ANvEU.P.Paired.nonStomach), color= Gene %in% SigAllAneu))+
  geom_point()+
  geom_density_2d(color = "white", bins = 10)+
  geom_point(data = subset(Merge.nonStomach, Gene%in% SigAllAneu))+
  ylab("-log2(P-value)")+
  xlab("Difference in Beta score\nAneuploid - euploid")+
  scale_color_manual(values = c("Black", "red"), name = "P<0.05 \nDiff < -.1")+
  theme_classic()
# 5x4
# plot.BetaDelta.ANvEU_WGL.nonStomach.pdf


Stomach.vs.nonStomach<- merge(Merge.St, Merge.nonStomach[,c(1,41,42)], by.x = "Gene", by.y = "Gene")

subset(Stomach.vs.nonStomach, Stomach.Diff.paired< -0.5 & 
                      Stomach.p.paired< 0.05 & 
                      ANvEU.Beta.Paired.nonStomach < 0)$Gene


SigAllAneu<-subset(Stomach.vs.nonStomach, Stomach.Diff.paired< -0.2 & 
                     Stomach.p.paired< 0.05 & 
                     ANvEU.Beta.Paired.nonStomach > 0)$Gene

ggplot(Stomach.vs.nonStomach, aes(x= Stomach.Diff.paired, 
                          y= ANvEU.Beta.Paired.nonStomach, 
                          color= Gene %in% SigAllAneu))+
  geom_point()+
  geom_hline(yintercept =0)+
  geom_vline(xintercept =0)+
  geom_point(data = subset(Stomach.vs.nonStomach, Gene%in% SigAllAneu))+
  scale_color_manual(values = c("black", "red"), 
                     name = "Stomach specific")+
  #geom_point(data = subset(Stomach.vs.nonStomach, Gene == "PLK4"), color = "blue")+
  geom_density_2d(color = "white", bins = 10)+
  ylab("Difference in Beta score (Non-Stomach)")+
  xlab("Difference in Beta score (Stomach)")+
  theme_classic()
# 5x4
# plot.BetaDelta.Stomach.nonStomach.pdf 



subset(Stomach.vs.nonStomach, Stomach.Diff.paired< -0.5 & 
         Stomach.p.paired< 0.05 & 
         ANvEU.Beta.Paired.nonStomach > 0)[,c(1,41,42,43,83,84)]
#  


subset(Stomach.vs.nonStomach, Stomach.Diff.paired< -0.2 & 
         Stomach.p.paired< 0.05 & 
         ANvEU.Beta.Paired.nonStomach < 0 &
         ANvEU.P.Paired.nonStomach< 0.05)[,c(1,41,42,43,83,84)] #Both Stomach & non Stomach

# write.csv(Stomach.vs.nonStomach, "Stomach.vs.nonStomach.csv")

