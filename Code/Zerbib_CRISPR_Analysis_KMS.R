#### CRISPR Screen data analysis from Zerbib et al. 2024 ####
##   Ben_David_CRISPR_Analysis_CC.R
##  From Zerbib J, et al. (2024) 
##  Uri Ben-David and Stefano Santaguida labs


# Data from: 
# Zerbib J, Ippolito MR, Eliezer Y, De Feudis G, Reuveni E, Savir Kadmon A, Martin S, 
# Viganò S, Leor G, Berstler J, Muenzner J, Mülleder M, Campagnolo EM, Shulman ED, Chang T, 
# Rubolino C, Laue K, Cohen-Sharir Y, Scorzoni S, Taglietti S, Ratti A, Stossel C, Golan T, 
# Nicassio F, Ruppin E, Ralser M, Vazquez F, Ben-David U, Santaguida S. Human aneuploid cells 
# depend on the RAF/MEK/ERK pathway for overcoming increased DNA damage. Nat Commun. 
# 2024 Sep 9;15(1):7772. doi: 10.1038/s41467-024-52176-x. PMID: 39251587; PMCID: PMC11385192.



### INTRO AND LIBRARIES ####
# Folders
## !! Update these locations to pathway where data was downloaded: 
# Folder with CRISPR screening and proteomics datasets: 
SchukkenData<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Data"
# Folder with dependency datasets: 
Dependency<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Dependency files"
# Folder with results: 
ResultsFile<- "/Volumes/Schukken_SSD/CRISPR SCREEN PAPER/R Code Data Files/Results"



### libraries
# mle analysis data
library(ggplot2) 
library(reshape2)
library(tidyverse)
library(readxl)
library(ggpubr)
library("cowplot")
library(plyr)
library('gprofiler2')
library(xlsx)
library("viridis")  # color packet
library("Hmisc")
library(corrplot)
library("PerformanceAnalytics")



### GET DATA #####
### Ben-David & Santaguida labs Screen data 

# Data from: 
# Zerbib J, Ippolito MR, Eliezer Y, De Feudis G, Reuveni E, Savir Kadmon A, Martin S, 
# Viganò S, Leor G, Berstler J, Muenzner J, Mülleder M, Campagnolo EM, Shulman ED, Chang T, 
# Rubolino C, Laue K, Cohen-Sharir Y, Scorzoni S, Taglietti S, Ratti A, Stossel C, Golan T, 
# Nicassio F, Ruppin E, Ralser M, Vazquez F, Ben-David U, Santaguida S. Human aneuploid cells 
# depend on the RAF/MEK/ERK pathway for overcoming increased DNA damage. Nat Commun. 
# 2024 Sep 9;15(1):7772. doi: 10.1038/s41467-024-52176-x. PMID: 39251587; PMCID: PMC11385192.

# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://pmc.ncbi.nlm.nih.gov/articles/PMC11385192/pdf/41467_2024_Article_52176.pdf


## Screen: screened with Avana library. 
# which contains 73,372 guides with an average of 4 guides per gene
# CRISPR dependency scores (CERES scores) were calculated as previously described 103 
# and were integrated with the data from all the cell lines screened as part of the Cancer Dependency Map, 21Q3 release

# Supplementary files 1 (cell line data) and 8 (CRISPR screen data) 

#Cell line names and aneuploidy status:

# RPE1-SS48	    Near diploid (Gain 10q)
# RPE1-SS77	    Near diploid (Gain 10q) 
# RPE1-SS6	    Ts7
# RPE1-SS119	  Ts8
# RPE1-SS51     Ts7.22


setwd(Dependency)

Dropouts.Santaguida <- read_xlsx("41467_2024_52176_MOESM10_ESM.xlsx")
# CERES Scores

colnames(Dropouts.Santaguida)<- c("gene", "RPE1", "RPE1_try2", "RPE1_Ts7", "RPE1_Ts8", "RPE1_Ts7.22")

#### Rename all Excel "Date" gene names

#Make function to rename "Sept-1" genes and "SEPT1" to their updated names "SEPTIN1" etc. for all 
#   genes who's names are mis-interpretted as dates in Excel. 
RenameExcelGenes2<- function(GeneMLEFile){
  GeneMLEFile$gene<- as.character(GeneMLEFile$gene)
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '1-Dec'] <- "DELEC1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'DEC1'] <- "DELEC1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '37226'] <- "DELEC1" #European Excel date format 1-December "2001"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '45627'] <- "DELEC1" #USA Excel date format, 2024
  
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '1-Sep'] <- "SEPTIN1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '2-Sep'] <- "SEPTIN2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '3-Sep'] <- "SEPTIN3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '4-Sep'] <- "SEPTIN4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '5-Sep'] <- "SEPTIN5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '6-Sep'] <- "SEPTIN6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '7-Sep'] <- "SEPTIN7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '8-Sep'] <- "SEPTIN8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '9-Sep'] <- "SEPTIN9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '10-Sep'] <- "SEPTIN10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '11-Sep'] <- "SEPTIN11"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '12-Sep'] <- "SEPTIN12"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '14-Sep'] <- "SEPTIN14"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '15-Sep'] <- "SELENOF"
  
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT1'] <- "SEPTIN1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT2'] <- "SEPTIN2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT3'] <- "SEPTIN3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT4'] <- "SEPTIN4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT5'] <- "SEPTIN5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT6'] <- "SEPTIN6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT7'] <- "SEPTIN7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT8'] <- "SEPTIN8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT9'] <- "SEPTIN9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT10'] <- "SEPTIN10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT11'] <- "SEPTIN11"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT12'] <- "SEPTIN12"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEPT14'] <- "SEPTIN14"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'SEP15'] <- "SELENOF"
  
  #European date format SEPT9 -> 1-September, 2009, -> numeric: 40057
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '37135'] <- "SEPTIN1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '37500'] <- "SEPTIN2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '37865'] <- "SEPTIN3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '38231'] <- "SEPTIN4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '38596'] <- "SEPTIN5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '38961'] <- "SEPTIN6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '39326'] <- "SEPTIN7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '39692'] <- "SEPTIN8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '40057'] <- "SEPTIN9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '40422'] <- "SEPTIN10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '40787'] <- "SEPTIN11"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '41153'] <- "SEPTIN12"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '41883'] <- "SEPTIN14"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '42248'] <- "SELENOF"
  
  
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '1-Mar'] <- "MARCHF1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '2-Mar'] <- "MARCHF2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '3-Mar'] <- "MARCHF3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '4-Mar'] <- "MARCHF4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '5-Mar'] <- "MARCHF5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '6-Mar'] <- "MARCHF6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '7-Mar'] <- "MARCHF7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '8-Mar'] <- "MARCHF8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '9-Mar'] <- "MARCHF9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '10-Mar'] <- "MARCHF10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '11-Mar'] <- "MARCHF11"
  
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH1'] <- "MARCHF1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH2'] <- "MARCHF2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH3'] <- "MARCHF3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH4'] <- "MARCHF4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH5'] <- "MARCHF5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH6'] <- "MARCHF6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH7'] <- "MARCHF7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH8'] <- "MARCHF8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH9'] <- "MARCHF9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH10'] <- "MARCHF10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == 'MARCH11'] <- "MARCHF11"
  
  # European date format MAR8 => 1-March 2008 ->
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '36951'] <- "MARCHF1"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '37316'] <- "MARCHF2"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '37681'] <- "MARCHF3"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '38047'] <- "MARCHF4"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '38412'] <- "MARCHF5"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '38777'] <- "MARCHF6"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '39142'] <- "MARCHF7"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '39508'] <- "MARCHF8"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '39873'] <- "MARCHF9"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '40238'] <- "MARCHF10"
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '40603'] <- "MARCHF11"
  
  #USA format saved in 2019. excel assumes the date you save is in year you saved it. so MAR1 = MARCH 1, YearOfSaving
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43525'] <- 'MARCHF1' 
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43534'] <- 'MARCHF10'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43535'] <- 'MARCHF11'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43526'] <- 'MARCHF2'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43527'] <- 'MARCHF3'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43528'] <- 'MARCHF4'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43529'] <- 'MARCHF5'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43530'] <- 'MARCHF6'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43531'] <- 'MARCHF7'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43532'] <- 'MARCHF8'
  GeneMLEFile["gene"][GeneMLEFile["gene"]  == '43533'] <- 'MARCHF9'
  
  
  
  return(GeneMLEFile)
} 

Dropouts.Santaguida<- RenameExcelGenes2(Dropouts.Santaguida) 

## Get pathway gene groups: 
setwd(Dependency)
Proteosome_Genes<- read.csv("group-690.csv", header=TRUE)
Ribosomal_Genes<- read.csv("group-1054.csv", skip = 1 , header=TRUE) #HGNC Ribosomal Proteins
Ribosomal_Genes_noMT<- subset(Ribosomal_Genes, Group %in% c("L ribosomal proteins", "S ribosomal proteins"))
Olfactory_Genes<- read.csv("group-141.csv", skip = 1 , header=TRUE)

RNA_Processing<- read_tsv("GO_term_summary_20250720_134600.txt")
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




### Plot Zerbib difference p.value #####
length(Dropouts.Santaguida$gene) # 18,119 genes


Dropouts.Santaguida_DiffP<- Dropouts.Santaguida
Dropouts.Santaguida_DiffP$MeanDi<- rowMeans(Dropouts.Santaguida_DiffP[,c(2,3)])

# $meanDiff is mean difference in gene dropout between an aneuploid cell line and it's corresponding control cell 
Dropouts.Santaguida_DiffP$Paired.Diff<- ( 
  (Dropouts.Santaguida_DiffP$RPE1_Ts7-Dropouts.Santaguida_DiffP$MeanDi) + 
    (Dropouts.Santaguida_DiffP$RPE1_Ts8-Dropouts.Santaguida_DiffP$MeanDi) + 
    (Dropouts.Santaguida_DiffP$RPE1_Ts7.22-Dropouts.Santaguida_DiffP$MeanDi) 
)/3

Dropouts.Santaguida_DiffP$MeanCRONOS<- rowMeans(Dropouts.Santaguida_DiffP[,c(2:6)])


# Now add p-values to difference between all aneuploid and all euploid
Dropouts.Santaguida_DiffP$pvalue<- NA
for (i in 1:length(Dropouts.Santaguida_DiffP$gene)){
  test<- t.test(as.numeric(Dropouts.Santaguida_DiffP[i,c(2,3)]), 
                as.numeric(Dropouts.Santaguida_DiffP[i,c(4,5,6)]), paired=FALSE) # t.test (aneuploid, euploid)
  Dropouts.Santaguida_DiffP$pvalue[i]<-test$p.value
} #paired t-test

Dropouts.Santaguida_DiffP<- Dropouts.Santaguida_DiffP[order(Dropouts.Santaguida_DiffP$pvalue), ]


## Plot difference v pvalue (scatter plot), highlight top hits: 
ggplot(Dropouts.Santaguida_DiffP, aes(y=-log(pvalue,2), x=Paired.Diff, 
                                      color=pvalue<0.05 & Paired.Diff< -0.25))+
  geom_point()+
  geom_density_2d(color = "white")+
  scale_color_manual(values=c("Black", "red"))+ 
  theme_classic()+
  theme(legend.position = "none")+
  xlab("CHRONOS delta (aneuploid - diploid)")+
  ylab("-log2(p-value)") 
# 5x4
# plot.Zerbib.CHRONOSDelta.pValue_RPE1.pdf

subset(Dropouts.Santaguida_DiffP, pvalue<0.05 & Paired.Diff< -0.25)$gene
# "VMA21"  "AMBRA1" "MVD"    "TAF13"  "KRTAP1" "NRG3"   "RB1"    "KMT2D"  "MBTPS2" "RIC1"   "SRP54"  "BICRA" 




## Highlight hits
HighlightTheseGenes<- c("PSMD14")
HighlightTheseGenes<- c("UBE2H")
HighlightTheseGenes<- c("FOSL1")
HighlightTheseGenes<- subset(Dropouts.Santaguida_DiffP, pvalue<0.05 & Paired.Diff< -0.2)$gene

ggplot(Dropouts.Santaguida_DiffP, aes(y=-log(pvalue,2), x=Paired.Diff, 
                                      color=gene %in% HighlightTheseGenes))+
  geom_point()+
  geom_point(data = subset(Dropouts.Santaguida_DiffP, gene %in% HighlightTheseGenes))+
  geom_density_2d(color = "white")+
  scale_color_manual(values=c("Black", "red"))+ 
  theme_classic()+
  theme(legend.position = "none")+
  xlab("Beta delta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.Zerbib.pairedP.Diff_TopHits.pdf


subset(Dropouts.Santaguida_DiffP, gene %in% HighlightTheseGenes)$gene



# Plot aneuploid pathways: 
ggplot(Dropouts.Santaguida_DiffP, aes(y=-log(pvalue,2), x=Paired.Diff))+
  geom_point()+
  geom_point(data = subset(Dropouts.Santaguida_DiffP, gene %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_point(data = subset(Dropouts.Santaguida_DiffP, gene %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  geom_point(data = subset(Dropouts.Santaguida_DiffP, gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_point(data = subset(Dropouts.Santaguida_DiffP, gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  geom_point(data = subset(Dropouts.Santaguida_DiffP, gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_point(data = subset(Dropouts.Santaguida_DiffP, gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  geom_density_2d(color = "white")+
  theme(legend.position="none")+
  theme_classic()+
  xlab("Difference in Beta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 5x4
# plot.Zerbib.pairedP.Diff_AllPathway.pdf
# Figure 4D

# density plot of pathways: 
ggplot(Dropouts.Santaguida_DiffP, aes(x=Paired.Diff))+
  geom_density(size=2, data = subset(Dropouts.Santaguida_DiffP, gene %in% MitoElectronChainnAssembly), color="#F8766D")+ #F8766D
  geom_density(size=2, data = subset(Dropouts.Santaguida_DiffP, gene %in% MitoTranslationTranscription), color="gold2")+   #7CAE00
  geom_density(size=2, data = subset(Dropouts.Santaguida_DiffP, gene %in% RNA_Processing_rRNA ), color="deepskyblue4") + #00BFC4
  geom_density(size=2, data = subset(Dropouts.Santaguida_DiffP, gene %in% Ribosomal_Genes_noMT$Approved.symbol ), color="#C77CFF") + #C77CFF
  geom_density(size=2, data = subset(Dropouts.Santaguida_DiffP, gene %in% Proteosome_Genes$Approved.symbol ), color="chartreuse1") + #7CAE00
  geom_density(size=2, data = subset(Dropouts.Santaguida_DiffP, gene %in% RNA_Processing_Spliceosome ), color="deepskyblue1") + #00BFC4
  geom_density(size=2)+
  theme(legend.position="none")+
  theme_classic()+
  xlab("Difference in Beta (aneuploid - euploid)")+
  ylab("-log2(p-value)")
# 3x4
# plot.Zerbib.DENSITY_AllPathway.pdf

# pathways t-test: 
t.test(subset(Dropouts.Santaguida_DiffP, ! gene %in% MitoElectronChainnAssembly)$Paired.Diff, 
       subset(Dropouts.Santaguida_DiffP, gene %in% MitoElectronChainnAssembly)$Paired.Diff ) # 1.63e-05
t.test(subset(Dropouts.Santaguida_DiffP, ! gene %in% MitoTranslationTranscription)$Paired.Diff, 
       subset(Dropouts.Santaguida_DiffP, gene %in% MitoTranslationTranscription)$Paired.Diff ) # 1.232e-14
t.test(subset(Dropouts.Santaguida_DiffP, ! gene %in% RNA_Processing_rRNA)$Paired.Diff, 
       subset(Dropouts.Santaguida_DiffP, gene %in% RNA_Processing_rRNA)$Paired.Diff ) # 0.01417
t.test(subset(Dropouts.Santaguida_DiffP, ! gene %in% Ribosomal_Genes_noMT$Approved.symbol)$Paired.Diff, 
       subset(Dropouts.Santaguida_DiffP, gene %in% Ribosomal_Genes_noMT$Approved.symbol)$Paired.Diff ) # 1.711743e-19
t.test(subset(Dropouts.Santaguida_DiffP, ! gene %in% Proteosome_Genes$Approved.symbol)$Paired.Diff, 
       subset(Dropouts.Santaguida_DiffP, gene %in% Proteosome_Genes$Approved.symbol)$Paired.Diff ) # 1.025e-05
t.test(subset(Dropouts.Santaguida_DiffP, ! gene %in% RNA_Processing_Spliceosome)$Paired.Diff, 
       subset(Dropouts.Santaguida_DiffP, gene %in% RNA_Processing_Spliceosome)$Paired.Diff ) # 1.885e-07




#### RANK GENES Zerbib ####
Dropouts.Santaguida_DiffP$AneuploidRank<- rank(Dropouts.Santaguida_DiffP$Paired.Diff)

# Highlight top 10% of aneuploid genes: 
length(Dropouts.Santaguida_DiffP$AneuploidRank) # 18,119 genes
HighlightTheseGenes<- subset(Dropouts.Santaguida_DiffP, AneuploidRank < 1812 )$gene

ggplot(Dropouts.Santaguida_DiffP, aes(x=AneuploidRank, y=Paired.Diff, 
                                      color= AneuploidRank < 1812))+
  geom_point()+
  geom_point(data = subset(Dropouts.Santaguida_DiffP, AneuploidRank < 1812))+
  scale_color_manual(values=c("Black", "Red"))+ 
  theme_classic()+
  theme(legend.position="none")+
  ylab("Difference in Beta score \n(aneuploid - euploid)")+
  xlab("Zerbib J et al: Rank")
# 5x3
# plot.Zerbib.Rank_Top10.pdf
# Figure 4B
subset(Dropouts.Santaguida_DiffP, gene %in% HighlightTheseGenes)

mean(subset(Dropouts.Santaguida_DiffP, gene %in% MitoTranslationTranscription)$Paired.Diff ) # -0.06498772
mean(subset(Dropouts.Santaguida_DiffP, gene %in% MitoElectronChainnAssembly)$Paired.Diff ) # -0.02763292


# write data as csv: 
setwd(ResultsFile)
write.csv(Dropouts.Santaguida_DiffP, "BenDavid.CRISPR_Diff.Pvalue.csv")
write.csv(subset(Dropouts.Santaguida_DiffP, AneuploidRank < 1812 )$gene, "Zerbib.BenDavid.RANK.Top10.csv")


### G Profiler: BenDavid #####

Dropouts.Santaguida_DiffP


### terms enriched in toxic to aneuploid cells
g.Antoxic.BD<- gost(
  subset(Dropouts.Santaguida_DiffP, Paired.Diff< -0.2)$gene,
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
  custom_bg = as.character(Dropouts.Santaguida_DiffP$gene),
  numeric_ns = "",
  sources = NULL,
  as_short_link = FALSE
)
g.Antoxic.BD_result<-g.Antoxic.BD$result


# terms enriched in toxic to aneuploid cells
P.AntoxicBD<-g.Antoxic.BD_result[order(g.Antoxic.BD_result$p_value),]
P.AntoxicBD<- subset(P.AntoxicBD, ! source %in% c("HPA"))
P.AntoxicBD$term_name[1:400]
P.AntoxicBD<-P.AntoxicBD[c(1,2,11,  3,4,  9,10,32,  51,56,58),]
P.AntoxicBD$term_name<- factor(P.AntoxicBD$term_name, levels= P.AntoxicBD$term_name)
P.AntoxicBD$Termtype<-c("RNA processing", "RNA processing", "RNA processing", 
                        "Protein Complex", "Protein Complex", 
                        "Ribosome","Ribosome","Ribosome",
                        "Metabolism", "Metabolism","Metabolism") 

ggplot(P.AntoxicBD, 
       aes(x=term_name, y=-log2(p_value), fill=Termtype))+
  geom_bar(stat="identity")+
  ylab("-log2 (P-Value)") +
  xlab("Key RNA terms")+
  ggtitle ("Top biological terms enriched: toxic to Aneuploid, Ben-David")+
  theme_classic()+
  coord_flip()
# plot.GProfiler.WGL.BenDavid.ANvEU.pdf
# 8x5

setwd(ResultsFile)
write.csv(g.Antoxic.BD_result[,1:13], "gprofiler.AneuploidDependency.BenDavid.2024.csv")



